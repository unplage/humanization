#!/usr/bin/env python3
"""Adapter: legacy ProteinMPNN CLI -> modern protein_mpnn_run.py CLI.

The humanization pipeline (scripts/humanize/mpnn.py) invokes ProteinMPNN with
the pre-2023 flags ``--chains_to_design`` / ``--fixed_residues``. Upstream
renamed the entry point to protein_mpnn_run.py and replaced those flags with
``--pdb_path_chains`` + ``--fixed_positions_jsonl``.

This adapter accepts the legacy CLI, builds the required JSONL files, calls the
real protein_mpnn_run.py, and copies the resulting FASTA back where the
pipeline's glob expects it (directly under --out_folder, not --out_folder/seqs/).

Usage (invoked by the pipeline, not by hand):
    python proteinmpnn_adapter.py --pdb_path in.pdb --out_folder out/ \
        --num_seq_per_target 32 --batch_size 8 --chains_to_design A \
        --fixed_residues 1,2,3,5 --seed 0 --suppress_print 1

Point the pipeline at it:
    PROTEINMPNN_DIR=/path/to/ProteinMPNN \
    python scripts/humanize/cli.py run --input ab.fasta --outdir out \
      --mpnn-mode local --mpnn-script tools/mpnn/proteinmpnn_adapter.py \
      --mpnn-python /path/to/proteinMPNN/bin/python

Fixed-residue numbering: --fixed_residues are 1-based positions in the chain
*sequence* (the same space the pipeline's NumberedChain.index uses), which is
exactly what ProteinMPNN's fixed_positions_jsonl expects.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

# Fallback search order when PROTEINMPNN_DIR is unset.
_DEFAULT_MPNN_DIRS = [
    os.path.expanduser("~/ProteinMPNN"),
    os.path.expanduser("~/ProteinMPNN-main"),
    os.path.expanduser("~/ProteinMPNN/vanilla_proteinmpnn"),
    "/home/jiemiaoxing/miniconda3/envs/proteinMPNN/ProteinMPNN-main",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="ProteinMPNN legacy->modern CLI adapter",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    # Legacy flags the pipeline passes.
    p.add_argument("--pdb_path", required=True, help="Input PDB file")
    p.add_argument("--out_folder", default="mpnn_out", help="Output folder")
    p.add_argument("--num_seq_per_target", type=int, default=32)
    p.add_argument("--batch_size", type=int, default=8)
    p.add_argument("--chains_to_design", default="A",
                   help="Space- or comma-separated chain IDs to design")
    p.add_argument("--fixed_residues", default="",
                   help="Comma-separated 1-based chain-sequence positions to fix")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--suppress_print", type=int, default=1)
    # Adapter-only controls.
    p.add_argument("--proteinmpnn-dir", default=os.environ.get("PROTEINMPNN_DIR", ""),
                   help="Path to the ProteinMPNN repo (or set PROTEINMPNN_DIR)")
    p.add_argument("--model-name", default="v_48_020")
    p.add_argument("--use-soluble-model", action="store_true")
    p.add_argument("--sampling-temp", default="0.1")
    p.add_argument("--keep-native", action="store_true",
                   help="Keep the native record in the output FASTA "
                        "(default: drop it so the pipeline does not count WT as a design)")
    p.add_argument("--strict-gaps", action="store_true",
                   help="Abort if the PDB chain has missing residue numbers "
                        "(ProteinMPNN inserts GAP slots, which can shift fixed positions)")
    return p.parse_args()


def find_mpnn_repo(explicit: str) -> str:
    """Locate the ProteinMPNN repo and return the dir containing protein_mpnn_run.py."""
    candidates = []
    if explicit:
        candidates.append(explicit)
    candidates.extend(_DEFAULT_MPNN_DIRS)
    for d in candidates:
        d = os.path.abspath(os.path.expanduser(d))
        if os.path.isfile(os.path.join(d, "protein_mpnn_run.py")):
            return d
        # Accept a path that points straight at the script.
        if os.path.basename(d) == "protein_mpnn_run.py" and os.path.isfile(d):
            return os.path.dirname(d)
    raise FileNotFoundError(
        "protein_mpnn_run.py not found. Set PROTEINMPNN_DIR or pass "
        "--proteinmpnn-dir. Searched: " + ", ".join(candidates)
    )


def _observe_residues(pdb_path: str) -> dict:
    """Return {chain: {(resseq, icode): resname}} from ATOM (+ HETATM MSE) records."""
    obs = {}
    with open(pdb_path, "r", errors="ignore") as fh:
        for line in fh:
            rec = line[:6].strip()
            resname = line[17:20].strip().upper()
            if rec == "HETATM" and resname != "MSE":
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            cid = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26].strip()
            obs.setdefault(cid, {})[(resseq, icode)] = resname
    return obs


def _span_and_observed(residues: dict):
    """Return (span_length, observed_count) using upstream's counting rules.

    Upstream (parse_PDB_biounits) iterates min..max residue number and emits
    one slot per insertion code, filling absent numbers with a GAP slot. So the
    parsed chain length is the numbered *span*, which can exceed the number of
    observed residues when the structure has missing residues.
    """
    if not residues:
        return 0, 0
    by_num = {}
    for (resseq, icode) in residues:
        by_num.setdefault(resseq, []).append(icode)
    span = 0
    for n in range(min(by_num), max(by_num) + 1):
        span += len(by_num.get(n, [None]))
    return span, len(residues)


def pdb_chains_and_lengths(pdb_path: str) -> dict:
    """Return {chain_id: parsed_length} using upstream's counting rules.

    Getting this right matters because ProteinMPNN indexes fixed positions into
    a chain_length array, so an over-count silently lets out-of-range positions
    through and crashes later in tied_featurize.
    """
    return {cid: _span_and_observed(res)[0]
            for cid, res in _observe_residues(pdb_path).items()}


def chain_span_gaps(pdb_path: str) -> dict:
    """Return {chain_id: (span, observed)}; span > observed means gaps exist."""
    return {cid: _span_and_observed(res)
            for cid, res in _observe_residues(pdb_path).items()}


def _is_cif(path: str) -> bool:
    return path.lower().endswith((".cif", ".mmcif"))


def convert_cif_to_pdb(cif_path: str):
    """Convert an mmCIF file to a temporary PDB. Returns (pdb_path, temp_path)."""
    try:
        from Bio.PDB.MMCIFParser import MMCIFParser
        from Bio.PDB.PDBIO import PDBIO
    except ImportError as exc:  # noqa: BLE001
        raise RuntimeError(
            "biopython is required to read mmCIF input; install it in the "
            "ProteinMPNN env or pass a PDB file"
        ) from exc
    structure = MMCIFParser(QUIET=True).get_structure("input", cif_path)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp.close()
    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp.name)
    return tmp.name, tmp.name


def target_name(pdb_path: str) -> str:
    """Mirror upstream: basename with the last 4 chars (extension) stripped."""
    return os.path.basename(pdb_path)[:-4]


def parse_fixed(res_str: str) -> list:
    """'1,2,3' -> [1, 2, 3]; tolerates whitespace and empty string."""
    return sorted({int(x) for x in res_str.replace(" ", "").split(",") if x.strip()})


def _run(args, work_pdb: str) -> int:
    repo = find_mpnn_repo(args.proteinmpnn_dir)
    run_py = os.path.join(repo, "protein_mpnn_run.py")

    chains_lengths = pdb_chains_and_lengths(work_pdb)
    chains_all = list(chains_lengths)
    if not chains_all:
        print(f"ERROR: no chains found in {work_pdb}", file=sys.stderr)
        return 1

    # Empty / "auto" means: design the only chain, or demand an explicit choice.
    raw_designed = args.chains_to_design.strip()
    if raw_designed.lower() in ("", "auto"):
        if len(chains_all) == 1:
            designed = [chains_all[0]]
        else:
            print(
                f"ERROR: PDB has multiple chains {chains_all}; specify "
                f"--chains_to_design explicitly (or --mpnn-chains on the pipeline)",
                file=sys.stderr,
            )
            return 1
    else:
        designed = raw_designed.replace(",", " ").split()
    bad = [c for c in designed if c not in chains_all]
    if bad:
        print(f"ERROR: chains {bad} not in PDB (available: {chains_all})",
              file=sys.stderr)
        return 1
    # Upstream derives the target name from the PDB it is given (work_pdb), so
    # the JSONL key must match that. Keep the output FASTA named after the
    # original input for readability.
    jsonl_name = target_name(work_pdb)
    out_stem = target_name(args.pdb_path)

    # ProteinMPNN treats seed=0 as "pick a random seed"; map it to a fixed
    # value so the pipeline's --seed 0 is reproducible.
    seed = args.seed if args.seed else 37

    # num_seq_per_target must be a multiple of batch_size or NUM_BATCHES=0.
    num_seqs = max(1, args.num_seq_per_target)
    batch = max(1, min(args.batch_size, num_seqs))
    if num_seqs % batch:
        num_seqs = (num_seqs // batch + 1) * batch

    os.makedirs(args.out_folder, exist_ok=True)

    fixed_positions = parse_fixed(args.fixed_residues)
    # Clip fixed positions to each designed chain's length. ProteinMPNN indexes
    # fixed_positions as np.array(pos)-1 into a chain_length array, so an
    # out-of-range position crashes with an opaque IndexError. A mismatch here
    # means the pipeline's chain numbering does not match the PDB (e.g. the
    # design PDB is not the grafted variant structure).
    clipped = {}
    for c in designed:
        n = chains_lengths[c]
        keep = [p for p in fixed_positions if 1 <= p <= n]
        dropped = [p for p in fixed_positions if p not in keep]
        clipped[c] = keep
        if dropped:
            print(
                f"WARNING: dropped {len(dropped)} fixed position(s) outside "
                f"chain {c} length {n}: {dropped[:10]}"
                f"{'...' if len(dropped) > 10 else ''}",
                file=sys.stderr,
            )
    with tempfile.TemporaryDirectory(prefix="mpnn_adapter_") as tmp:
        # fixed_positions_jsonl: {name: {chain: [1-based fixed positions]}}
        fp_path = os.path.join(tmp, "fixed_positions.jsonl")
        fp_dict = {jsonl_name: {}}
        for c in chains_all:
            fp_dict[jsonl_name][c] = list(clipped[c]) if c in designed else []
        with open(fp_path, "w") as fh:
            fh.write(json.dumps(fp_dict) + "\n")

        cmd = [
            sys.executable, run_py,
            "--pdb_path", work_pdb,
            "--pdb_path_chains", " ".join(designed),
            "--out_folder", args.out_folder,
            "--num_seq_per_target", str(num_seqs),
            "--batch_size", str(batch),
            "--sampling_temp", args.sampling_temp,
            "--model_name", args.model_name,
            "--fixed_positions_jsonl", fp_path,
            "--seed", str(seed),
            "--suppress_print", str(args.suppress_print),
        ]
        if args.use_soluble_model:
            cmd.append("--use_soluble_model")

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            sys.stderr.write(proc.stdout[-4000:])
            sys.stderr.write(proc.stderr[-4000:])
            print(f"ERROR: protein_mpnn_run.py exited {proc.returncode}",
                  file=sys.stderr)
            return 1
        if not args.suppress_print and proc.stdout:
            sys.stdout.write(proc.stdout)

        # Upstream writes <out_folder>/seqs/<name>.fa; the pipeline globs
        # <out_folder>/*fa. Copy the FASTA up one level.
        produced = sorted(glob.glob(
            os.path.join(args.out_folder, "seqs", f"{jsonl_name}.fa")))
        if not produced:
            produced = sorted(glob.glob(
                os.path.join(args.out_folder, "seqs", "*.fa")))
        if not produced:
            print("ERROR: protein_mpnn_run.py produced no FASTA", file=sys.stderr)
            return 1

        for fa in produced:
            if len(produced) == 1:
                dest = os.path.join(args.out_folder, out_stem + ".fa")
            else:
                dest = os.path.join(args.out_folder, os.path.basename(fa))
            if args.keep_native:
                shutil.copyfile(fa, dest)
            else:
                _drop_native_record(fa, dest)
    return 0


def main() -> int:
    args = parse_args()

    work_pdb = args.pdb_path
    temp_pdb = None
    try:
        if _is_cif(work_pdb):
            try:
                work_pdb, temp_pdb = convert_cif_to_pdb(work_pdb)
            except Exception as exc:  # noqa: BLE001
                print(f"ERROR: failed to convert {args.pdb_path} to PDB: {exc}",
                      file=sys.stderr)
                return 1

        # Missing residue numbers make ProteinMPNN insert GAP slots, which can
        # shift fixed positions that lie after the gap. Warn (or abort).
        gapped = {c: (span, obs)
                  for c, (span, obs) in chain_span_gaps(work_pdb).items()
                  if span > obs}
        if gapped:
            detail = ", ".join(f"{c}:{obs}/{span}"
                               for c, (span, obs) in gapped.items())
            msg = (f"WARNING: chain(s) have missing residue numbers "
                   f"(observed/span: {detail}); ProteinMPNN inserts GAP slots, "
                   f"so fixed positions after a gap may shift")
            if args.strict_gaps:
                print("ERROR: " + msg, file=sys.stderr)
                return 1
            print(msg, file=sys.stderr)

        return _run(args, work_pdb)
    finally:
        if temp_pdb and os.path.exists(temp_pdb):
            try:
                os.unlink(temp_pdb)
            except OSError:
                pass


def _drop_native_record(src: str, dest: str) -> None:
    """Copy FASTA to dest, dropping the native (WT) record.

    Upstream always writes the native first, tagged with designed_chains=.
    Generated samples are tagged with T=<temp>. Keep only sample records.
    """
    with open(src) as fh, open(dest, "w") as out:
        header = None
        seq_lines = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                _emit(header, seq_lines, out)
                header, seq_lines = line, []
            elif header is not None:
                if line.strip():
                    seq_lines.append(line.strip())
        _emit(header, seq_lines, out)


def _emit(header, seq_lines, out) -> None:
    if header is None:
        return
    if "designed_chains=" in header:   # native / WT record
        return
    out.write(header + "\n")
    for s in seq_lines:
        out.write(s + "\n")


if __name__ == "__main__":
    sys.exit(main())
