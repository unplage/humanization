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


def pdb_chains_and_lengths(pdb_path: str) -> dict:
    """Return {chain_id: n_residues} using upstream's counting rules.

    Mirrors ProteinMPNN's parse_PDB_biounits(): only ATOM records count (plus
    HETATM MSE re-read as MET); residues are keyed by (resSeq, iCode); and
    missing residue numbers *inside* the min..max range still occupy a slot, so
    the length equals the numbered span, not the number of observed residues.
    Getting this right matters because ProteinMPNN indexes fixed positions into
    a chain_length array, so an over-count silently lets out-of-range
        positions through and crashes later in tied_featurize.
    """
    # chain -> {(resseq, icode): resname}
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

    lengths = {}
    for cid, residues in obs.items():
        # Group insertion codes under their residue number (upstream sorts them).
        by_num = {}
        for (resseq, icode) in residues:
            by_num.setdefault(resseq, []).append(icode)
        if not by_num:
            lengths[cid] = 0
            continue
        # Upstream iterates min..max and emits one slot per insertion code,
        # filling absent numbers with a GAP slot.
        total = 0
        for n in range(min(by_num), max(by_num) + 1):
            total += len(by_num.get(n, [None]))
        lengths[cid] = total
    return lengths


def target_name(pdb_path: str) -> str:
    """Mirror upstream: basename with the last 4 chars (extension) stripped."""
    return os.path.basename(pdb_path)[:-4]


def parse_fixed(res_str: str) -> list:
    """'1,2,3' -> [1, 2, 3]; tolerates whitespace and empty string."""
    return sorted({int(x) for x in res_str.replace(" ", "").split(",") if x.strip()})


def main():
    args = parse_args()

    repo = find_mpnn_repo(args.proteinmpnn_dir)
    run_py = os.path.join(repo, "protein_mpnn_run.py")

    chains_lengths = pdb_chains_and_lengths(args.pdb_path)
    chains_all = list(chains_lengths)
    if not chains_all:
        print(f"ERROR: no chains found in {args.pdb_path}", file=sys.stderr)
        return 1

    designed = args.chains_to_design.replace(",", " ").split()
    bad = [c for c in designed if c not in chains_all]
    if bad:
        print(f"ERROR: chains {bad} not in PDB (available: {chains_all})",
              file=sys.stderr)
        return 1
    name = target_name(args.pdb_path)

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
        fp_dict = {name: {}}
        for c in chains_all:
            fp_dict[name][c] = list(clipped[c]) if c in designed else []
        with open(fp_path, "w") as fh:
            fh.write(json.dumps(fp_dict) + "\n")

        cmd = [
            sys.executable, run_py,
            "--pdb_path", args.pdb_path,
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
            os.path.join(args.out_folder, "seqs", f"{name}.fa")))
        if not produced:
            produced = sorted(glob.glob(
                os.path.join(args.out_folder, "seqs", "*.fa")))
        if not produced:
            print("ERROR: protein_mpnn_run.py produced no FASTA", file=sys.stderr)
            return 1

        for fa in produced:
            dest = os.path.join(args.out_folder, os.path.basename(fa))
            if args.keep_native:
                shutil.copyfile(fa, dest)
            else:
                _drop_native_record(fa, dest)
    return 0


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
