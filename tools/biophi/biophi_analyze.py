#!/usr/bin/env python3
"""Standalone BioPhi wrapper for Sapiens and OASis humanness analysis.

This tool is intentionally independent of ``scripts/humanize``. It shells out
to a dedicated BioPhi conda environment so the portable pipeline keeps its
stdlib-only core and does not get pinned to BioPhi's old ANARCI/abnumber stack.

Examples
--------
python3 tools/biophi/biophi_analyze.py check \
    --oasis-db data/biophi/OASis_9mers_v1.db

python3 tools/biophi/biophi_analyze.py sapiens \
    --input data/examples/mouse_4d5_fab.fasta --outdir outputs/biophi

python3 tools/biophi/biophi_analyze.py oasis \
    --input data/examples/mouse_4d5_fab.fasta \
    --oasis-db data/biophi/OASis_9mers_v1.db --outdir outputs/biophi
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

TOOL_VERSION = "0.1.0"
DEFAULT_ENV = "biophi"
DEFAULT_DB = "data/biophi/OASis_9mers_v1.db"
HMMER_BINARIES = ("hmmsearch", "hmmbuild", "hmmscan")
SAPIENS_CACHE_KEYS = (
    "models--prihodad--biophi-sapiens1-vh",
    "models--prihodad--biophi-sapiens1-vl",
    "models--prihodad--biophi-sapiens1-tokenizer",
)


class ToolError(RuntimeError):
    pass


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _md5_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    seq_parts: List[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                if name is None:
                    raise ToolError(f"FASTA sequence before header in {path}")
                seq_parts.append(line.replace(" ", "").upper())
    if name is not None:
        records.append((name, "".join(seq_parts)))
    if not records:
        raise ToolError(f"No FASTA records found in {path}")
    seen: Dict[str, int] = {}
    for record_name, seq in records:
        if not seq:
            raise ToolError(f"Empty sequence for record {record_name!r} in {path}")
        seen[record_name] = seen.get(record_name, 0) + 1
    duplicates = sorted(k for k, v in seen.items() if v > 1)
    if duplicates:
        raise ToolError(
            "Duplicate FASTA IDs are not supported by BioPhi: "
            + ", ".join(duplicates[:5])
        )
    return records


def _conda_run(env: str, cmd: Sequence[str], timeout: int,
               extra_env: Optional[Dict[str, str]] = None) -> subprocess.CompletedProcess:
    conda = shutil.which("conda")
    if not conda:
        raise ToolError("conda not found and biophi is not on PATH")
    full = [conda, "run", "-n", env, *cmd]
    env_vars = os.environ.copy()
    if extra_env:
        env_vars.update(extra_env)
    try:
        return subprocess.run(
            full, capture_output=True, text=True, timeout=timeout, env=env_vars
        )
    except subprocess.TimeoutExpired as exc:
        raise ToolError(
            f"BioPhi timed out after {timeout}s: {' '.join(cmd)}"
        ) from exc


def _run_biophi(args: Sequence[str], env: str, timeout: int,
                hf_endpoint: Optional[str], verbose: bool = False) -> subprocess.CompletedProcess:
    extra_env: Dict[str, str] = {}
    if hf_endpoint:
        extra_env["HF_ENDPOINT"] = hf_endpoint
    else:
        # The server has no direct HuggingFace access. Local cache is populated.
        extra_env["HF_HUB_OFFLINE"] = "1"
        extra_env["TRANSFORMERS_OFFLINE"] = "1"
    binary = shutil.which(args[0])
    if binary:
        proc_env = os.environ.copy()
        proc_env.update(extra_env)
        if verbose:
            # Stream output in real-time
            print(f"  Running: {' '.join(args)}", file=sys.stderr, flush=True)
            proc = subprocess.Popen(
                list(args), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=proc_env
            )
            stdout_lines = []
            assert proc.stdout is not None
            for line in proc.stdout:
                print(f"  {line.rstrip()}", file=sys.stderr, flush=True)
                stdout_lines.append(line)
            proc.wait()
            return subprocess.CompletedProcess(args, proc.returncode, "".join(stdout_lines), "")
        else:
            try:
                return subprocess.run(
                    list(args), capture_output=True, text=True,
                    timeout=timeout, env=proc_env
                )
            except subprocess.TimeoutExpired as exc:
                raise ToolError(f"BioPhi timed out after {timeout}s: {' '.join(args)}") from exc
    return _conda_run(env, args, timeout=timeout, extra_env=extra_env)


def _find_env_python(env: str) -> Optional[str]:
    """Locate the interpreter of a named conda environment."""
    candidates: List[Path] = []
    conda = shutil.which("conda")
    if conda:
        conda_root = Path(conda).resolve().parent.parent
        candidates.append(conda_root / "envs" / env / "bin" / "python")
        try:
            out = subprocess.run(
                [conda, "info", "--envs"], capture_output=True, text=True, timeout=60
            )
            for line in out.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == env:
                    candidates.append(Path(parts[1]) / "bin" / "python")
        except Exception:
            pass
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix and Path(prefix).name == env:
        candidates.append(Path(prefix) / "bin" / "python")
    for cand in candidates:
        if cand.exists() and os.access(str(cand), os.X_OK):
            return str(cand)
    return None


def _run_env_python(env: str, script_path: Path, script_args: Sequence[str],
                    timeout: int, hf_endpoint: Optional[str],
                    verbose: bool = False) -> subprocess.CompletedProcess:
    """Run a helper script with the conda env's Python (never system Python)."""
    extra_env: Dict[str, str] = {}
    if hf_endpoint:
        extra_env["HF_ENDPOINT"] = hf_endpoint
    else:
        extra_env["HF_HUB_OFFLINE"] = "1"
        extra_env["TRANSFORMERS_OFFLINE"] = "1"
    env_python = _find_env_python(env)
    cmd: List[str]
    if env_python:
        cmd = [env_python, str(script_path), *script_args]
        proc_env = os.environ.copy()
        proc_env.update(extra_env)
        if verbose:
            print(f"  Running: {' '.join(cmd)}", file=sys.stderr, flush=True)
            proc = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                text=True, bufsize=1, env=proc_env,
            )
            lines: List[str] = []
            assert proc.stdout is not None
            for line in proc.stdout:
                print(f"  {line.rstrip()}", file=sys.stderr, flush=True)
                lines.append(line)
            proc.wait()
            return subprocess.CompletedProcess(cmd, proc.returncode, "".join(lines), "")
        try:
            return subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout, env=proc_env
            )
        except subprocess.TimeoutExpired as exc:
            raise ToolError(
                f"Helper script timed out after {timeout}s: {script_path.name}"
            ) from exc
    return _conda_run(
        env, ["python", str(script_path), *script_args],
        timeout=timeout, extra_env=extra_env,
    )


def _check_process(proc: subprocess.CompletedProcess, label: str) -> None:
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        if len(detail) > 1500:
            detail = detail[-1500:]
        raise ToolError(f"{label} failed (exit {proc.returncode}):\n{detail}")


def _env_info(env: str) -> Dict[str, object]:
    code = (
        "import json;"
        "import biophi, sapiens, torch, transformers, abnumber;"
        "print(json.dumps({"
        "'biophi': biophi.__version__,"
        "'sapiens': sapiens.__version__,"
        "'torch': torch.__version__,"
        "'transformers': transformers.__version__,"
        "'abnumber': abnumber.__version__"
        "}))"
    )
    proc = _conda_run(env, ["python", "-c", code], timeout=120)
    if proc.returncode != 0:
        return {"error": (proc.stderr or proc.stdout or "").strip()[-1000:]}
    try:
        return json.loads(proc.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return {"error": proc.stdout.strip()[-1000:]}


def _hf_cache_info() -> Dict[str, object]:
    cache = Path(os.environ.get("HF_HOME", Path.home() / ".cache" / "huggingface"))
    hub = cache / "hub"
    present = {key: (hub / key).is_dir() for key in SAPIENS_CACHE_KEYS}
    return {"hub": str(hub), "models": present, "complete": all(present.values())}


def _db_info(path: Path, with_md5: bool = False) -> Dict[str, object]:
    info: Dict[str, object] = {
        "path": str(path.resolve()),
        "exists": path.exists(),
    }
    if path.exists():
        stat = path.stat()
        info.update({"size_bytes": stat.st_size, "mtime": stat.st_mtime})
        if path.suffix == ".gz":
            info["error"] = "OASis database must be unzipped (gunzip the .db.gz first)"
        if with_md5:
            info["md5"] = _md5_file(path)
    return info


def _write_json(path: Path, payload: Dict[str, object]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: Iterable[Dict[str, object]],
               fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def cmd_check(args: argparse.Namespace) -> int:
    problems: List[str] = []
    warnings: List[str] = []
    report: Dict[str, object] = {
        "tool": "tools/biophi",
        "tool_version": TOOL_VERSION,
        "timestamp": _now(),
        "env": args.env,
    }

    if shutil.which("conda") is None:
        problems.append("conda not found (needed to run the isolated BioPhi env)")
    else:
        report["env_packages"] = _env_info(args.env)
        if "error" in report["env_packages"]:
            problems.append("BioPhi environment is not usable: "
                            + str(report["env_packages"]["error"]))

    missing_hmmer = [name for name in HMMER_BINARIES if shutil.which(name) is None]
    report["hmmer"] = {"missing": missing_hmmer}
    if missing_hmmer:
        problems.append("HMMER binaries missing: " + ", ".join(missing_hmmer))

    hf = _hf_cache_info()
    report["sapiens_cache"] = hf
    if not hf["complete"]:
        warnings.append(
            "Sapiens VH/VL/tokenizer are not fully cached; use --hf-endpoint "
            "https://hf-mirror.com once, or pre-download the files listed in README"
        )

    db = Path(args.oasis_db)
    report["oasis_db"] = _db_info(db, with_md5=args.md5)
    if not db.exists():
        warnings.append(f"OASis database not found: {db}")
    elif db.suffix == ".gz":
        problems.append(f"OASis database is still gzipped: {db} (run gunzip)")
    elif args.md5 and report["oasis_db"].get("md5") and db.suffix == ".gz":
        expected = "28fa18135b3ef1710970f72b5953cbc5"
        if report["oasis_db"]["md5"] != expected:
            warnings.append(
                f"OASis DB MD5 differs from the published v1 file ({expected})"
            )

    report["problems"] = problems
    report["warnings"] = warnings
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 1 if problems else 0


def _sapiens_outputs(records: List[Tuple[str, str]], scores: Dict[str, Dict[str, object]],
                     outdir: Path, args: argparse.Namespace,
                     input_path: Path, env_info: Dict[str, object],
                     detail: Optional[Dict[str, Dict[str, Any]]] = None) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    detail = detail or {}
    rows = []
    for name, seq in records:
        score = scores.get(name, {})
        rec_detail = detail.get(name, {})
        rows.append({
            "id": name,
            "length": len(seq),
            "chain": score.get("chain") or rec_detail.get("chain", ""),
            "sapiens_score": score.get("sapiens_score", ""),
            "sapiens_mean_api": rec_detail.get("mean"),
            "n_positions": len(rec_detail.get("per_position") or []),
        })

    payload = {
        "tool": "tools/biophi",
        "tool_version": TOOL_VERSION,
        "mode": "sapiens",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
            "n_records": len(records),
        },
        "parameters": {
            "numbering_scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
        },
        "environment": env_info,
        "results": rows,
        "per_position": {
            rid: d.get("per_position") for rid, d in detail.items()
        },
    }
    _write_json(outdir / f"{stem}_biophi_sapiens.json", payload)
    _write_csv(outdir / f"{stem}_biophi_sapiens.csv", rows,
               ["id", "length", "chain", "sapiens_score", "sapiens_mean_api",
                "n_positions"])

    # Full 20-aa matrix (BioPhi --scores-only equivalent)
    matrix_rows = []
    for rid, d in detail.items():
        for entry in d.get("matrix", []):
            row = {"id": rid, "chain": d.get("chain", "")}
            row.update(entry)
            matrix_rows.append(row)
    if matrix_rows:
        _write_csv(
            outdir / f"{stem}_biophi_scores.csv",
            matrix_rows,
            ["id", "chain", "position", "numbering", "input_aa",
             "sapiens_score", *list("ACDEFGHIKLMNPQRSTVWY")],
        )

    def _cell(value: object) -> str:
        if value is None or value == "":
            return "N/A"
        return str(value).replace("|", "\\|")

    lines = [
        f"# BioPhi Sapiens report: {input_path.name}",
        "",
        f"- Generated: {payload['timestamp']}",
        f"- BioPhi: {env_info.get('biophi', 'unknown')}",
        f"- Sapiens: {env_info.get('sapiens', 'unknown')}",
        f"- Numbering: {args.scheme} / CDR definition: {args.cdr_definition}",
        "",
        "| ID | Length | Chain | Mean Sapiens score |",
        "|----|--------|-------|--------------------|",
    ]
    for row in rows:
        lines.append(
            f"| {_cell(row['id'])} | {row['length']} | {_cell(row['chain'])} | "
            f"{_cell(row['sapiens_score'])} |"
        )
    lines += [
        "",
        "The Sapiens score is the mean per-residue probability under the human "
        "antibody language model. Higher is more human-like; it is not a "
        "germline identity percentage.",
        "",
        "## Per-record details",
        "",
    ]
    for rid, d in detail.items():
        lines.append(f"### {_cell(rid)}")
        lines.append(
            f"- mean={_cell(d.get('mean'))} min={_cell(d.get('min'))} "
            f"max={_cell(d.get('max'))} positions<0.5={d.get('n_below_0_5')}"
        )
        hotspots = d.get("hotspots") or []
        if hotspots:
            hot = ", ".join(
                f"{h.get('numbering', h.get('position'))}({h.get('residue')})"
                f"={h.get('score')}"
                for h in hotspots
            )
            lines.append(f"- hotspots (lowest): {hot}")
        if d.get("best_score_sequence"):
            lines.append(
                f"- best-score sequence: `{d['best_score_sequence']}`"
            )
        if d.get("masked_sequence"):
            lines.append(f"- masked sequence: `{d['masked_sequence']}`")
        lines.append("")
    (outdir / f"{stem}_biophi_sapiens.md").write_text("\n".join(lines), encoding="utf-8")


def cmd_sapiens(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    records = read_fasta(input_path)
    outdir = Path(args.outdir)

    with tempfile.TemporaryDirectory(prefix="biophi_sapiens_") as td:
        scores_csv = Path(td) / "sapiens_scores.csv"
        cmd = [
            "biophi", "sapiens", str(input_path),
            "--mean-score-only",
            "--output", str(scores_csv),
            "--scheme", args.scheme,
            "--cdr-definition", args.cdr_definition,
        ]
        proc = _run_biophi(cmd, args.env, args.timeout, args.hf_endpoint)
        _check_process(proc, "biophi sapiens")
        if not scores_csv.exists():
            raise ToolError("biophi sapiens did not write the expected CSV")
        scores: Dict[str, Dict[str, object]] = {}
        with scores_csv.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                record_id = row.get("id", "")
                raw = row.get("sapiens_score", "")
                try:
                    value: object = round(float(raw), 6)
                except (TypeError, ValueError):
                    value = raw
                scores[record_id] = {
                    "chain": row.get("chain", ""),
                    "sapiens_score": value,
                }

        if args.humanize:
            outdir.mkdir(parents=True, exist_ok=True)
            humanized = outdir / f"{input_path.stem}_biophi_humanized.fa"
            humanize_cmd = [
                "biophi", "sapiens", str(input_path),
                "--fasta-only", "--output", str(humanized),
                "--scheme", args.scheme,
                "--cdr-definition", args.cdr_definition,
                "--iterations", str(args.iterations),
            ]
            if args.humanize_cdrs:
                humanize_cmd.append("--humanize-cdrs")
            humanize_proc = _run_biophi(
                humanize_cmd, args.env, args.timeout, args.hf_endpoint
            )
            _check_process(humanize_proc, "biophi sapiens --fasta-only")

    # Per-position + full 20-aa matrix (default: on, fullest analysis)
    detail: Dict[str, Dict[str, Any]] = {}
    if not args.no_per_position:
        numbering_map: Optional[Dict[str, List[str]]] = None
        metadata_path = getattr(args, "metadata", None)
        if metadata_path and Path(metadata_path).exists():
            with open(metadata_path) as fh:
                meta = json.load(fh)
            variants_map = meta.get("variants", meta) if isinstance(meta, dict) else {}
            numbering_map = {}
            for rid, _seq in records:
                parts = rid.split("|")
                key = parts[1] if len(parts) > 1 else rid
                labels = variants_map.get(key) or variants_map.get(rid)
                if isinstance(labels, list):
                    numbering_map[rid] = labels
        detail = _sapiens_per_position_python_api(
            records, env=args.env, timeout=args.timeout,
            hf_endpoint=args.hf_endpoint,
            scheme=args.scheme, cdr_definition=args.cdr_definition,
            numbering_map=numbering_map,
            verbose=getattr(args, "verbose", False),
        )

    env_info = _env_info(args.env)
    _sapiens_outputs(records, scores, outdir, args, input_path, env_info, detail)

    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.json')}")
    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.md')}")
    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.csv')}")
    if detail:
        print(f"Wrote {outdir / (input_path.stem + '_biophi_scores.csv')}")
    return 0


def _clean_cell(value: object) -> object:
    # BioPhi 1.0.11 sometimes writes a one-element tuple such as "(61,)",
    # which xlsxwriter serialises as the literal string "(61,)".
    if isinstance(value, tuple) and len(value) == 1:
        return value[0]
    if isinstance(value, str):
        match = re.fullmatch(r"\((\d+),\)", value.strip())
        if match:
            return int(match.group(1))
    return value


def _read_oasis_overview(xlsx: Path) -> List[Dict[str, object]]:
    try:
        import openpyxl  # type: ignore
    except ImportError:
        try:
            import pandas as pd  # type: ignore
        except ImportError:
            raise ToolError(
                "Reading BioPhi XLSX requires openpyxl or pandas in the current "
                "Python environment"
            )
        df = pd.read_excel(xlsx, sheet_name="Overview")
        df = df.where(df.notna(), "")
        return [
            {key: _clean_cell(value) for key, value in row.items()}
            for row in df.reset_index().to_dict(orient="records")
        ]

    workbook = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
    if "Overview" not in workbook.sheetnames:
        raise ToolError(f"BioPhi report has no Overview sheet: {xlsx}")
    sheet = workbook["Overview"]
    rows = list(sheet.iter_rows(values_only=True))
    if not rows:
        return []
    header = [str(h) if h is not None else "" for h in rows[0]]
    records: List[Dict[str, object]] = []
    for raw in rows[1:]:
        if all(cell is None for cell in raw):
            continue
        record = {}
        for key, value in zip(header, raw):
            if key == "" and value is not None:
                key = "Antibody"
            record[key] = _clean_cell(value) if value is not None else ""
        records.append(record)
    workbook.close()
    return records


# Helper executed inside the biophi conda env (system Python has no sapiens).
_PER_POSITION_HELPER = r'''
import json, os, sys
os.environ.setdefault("HF_HUB_OFFLINE", "1")
os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")
from sapiens.predict import predict_scores

records_path, out_path, scheme, cdr_definition = sys.argv[1:5]
numbering_path = sys.argv[5] if len(sys.argv) > 5 else ""
with open(records_path) as fh:
    records = json.load(fh)
numbering_map = {}
if numbering_path and os.path.exists(numbering_path):
    with open(numbering_path) as fh:
        numbering_map = json.load(fh)

def chain_type_for(record_id, seq):
    try:
        from abnumber import Chain
        chain = Chain(seq, scheme=scheme, cdr_definition=cdr_definition)
        if chain.chain_type in ("H", "K", "L"):
            return chain.chain_type
    except Exception:
        pass
    rid = record_id.upper()
    for suffix, token in (
        ("_VH", "H"), ("_VL", "L"), ("_HC", "H"), ("_LC", "L"),
        (".VH", "H"), (".VL", "L"), (".HC", "H"), (".LC", "L"),
        ("_H", "H"), ("_L", "L"), (".H", "H"), (".L", "L"),
    ):
        if rid.endswith(suffix):
            return token
    if "|H_" in record_id or rid.startswith("H"):
        return "H"
    return "L"

results = {}
for record_id, seq in records:
    ct = chain_type_for(record_id, seq)
    labels = numbering_map.get(record_id)
    try:
        df = predict_scores(seq, chain_type=ct)
        per_pos = []
        matrix = []
        observed = []
        best_chars = []
        for idx, residue in enumerate(seq):
            row = df.iloc[idx]
            score = float(row[residue]) if residue in df.columns else None
            if score is not None:
                observed.append(score)
            best_aa = str(row.idxmax())
            best_chars.append(best_aa)
            label = labels[idx] if labels and idx < len(labels) else str(idx + 1)
            per_pos.append({
                "position": idx + 1,
                "numbering": label,
                "residue": residue,
                "score": round(score, 6) if score is not None else None,
            })
            entry = {
                "position": idx + 1,
                "numbering": label,
                "input_aa": residue,
                "sapiens_score": round(score, 6) if score is not None else None,
            }
            for aa in "ACDEFGHIKLMNPQRSTVWY":
                entry[aa] = round(float(row[aa]), 6) if aa in df.columns else None
            matrix.append(entry)
        best_seq = "".join(best_chars)
        masked_seq = "".join(
            b if a in "*X" else a for a, b in zip(seq, best_seq)
        )
        hotspots = sorted(
            (p for p in per_pos if p["score"] is not None),
            key=lambda p: p["score"],
        )[:5]
        results[record_id] = {
            "chain": ct,
            "mean": round(sum(observed) / len(observed), 6) if observed else None,
            "min": round(min(observed), 6) if observed else None,
            "max": round(max(observed), 6) if observed else None,
            "n_below_0_5": sum(1 for s in observed if s < 0.5),
            "per_position": per_pos,
            "matrix": matrix,
            "best_score_sequence": best_seq,
            "masked_sequence": masked_seq,
            "hotspots": hotspots,
        }
    except Exception as exc:
        results[record_id] = {
            "chain": ct, "mean": None, "per_position": [], "matrix": [],
            "error": str(exc),
        }

with open(out_path, "w") as fh:
    json.dump(results, fh)
'''


def _sapiens_per_position_python_api(
    records: List[Tuple[str, str]], env: str = DEFAULT_ENV,
    timeout: int = 3600, hf_endpoint: Optional[str] = None,
    scheme: str = "kabat", cdr_definition: str = "kabat",
    numbering_map: Optional[Dict[str, List[str]]] = None,
    verbose: bool = False,
) -> Dict[str, Dict[str, Any]]:
    """Per-position + full 20-aa matrix via the Sapiens Python API.

    Runs inside the biophi conda env (the BioPhi CLI ``--scores-only`` has a
    known bug that leaves position/residue fields empty).

    Returns dict: {record_id: {"chain", "mean", "per_position", "matrix",
    "best_score_sequence", "masked_sequence", "hotspots", ...}}.
    """
    with tempfile.TemporaryDirectory(prefix="biophi_perpos_") as td:
        td_path = Path(td)
        helper = td_path / "per_position_helper.py"
        helper.write_text(_PER_POSITION_HELPER, encoding="utf-8")
        records_json = td_path / "records.json"
        out_json = td_path / "per_position.json"
        records_json.write_text(
            json.dumps([[rid, seq] for rid, seq in records]), encoding="utf-8"
        )
        cmd_args = [str(records_json), str(out_json), scheme, cdr_definition]
        if numbering_map:
            numbering_json = td_path / "numbering.json"
            numbering_json.write_text(
                json.dumps(numbering_map), encoding="utf-8"
            )
            cmd_args.append(str(numbering_json))
        proc = _run_env_python(
            env, helper, cmd_args,
            timeout=timeout, hf_endpoint=hf_endpoint, verbose=verbose,
        )
        _check_process(proc, "Sapiens per-position helper")
        if not out_json.exists():
            raise ToolError("Sapiens per-position helper did not write its output")
        with out_json.open() as fh:
            raw: Dict[str, Dict[str, Any]] = json.load(fh)
    return raw


def _oasis_outputs(rows: List[Dict[str, object]], outdir: Path, args: argparse.Namespace,
                   input_path: Path, xlsx: Path, db: Path,
                   env_info: Dict[str, object]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    kept = outdir / f"{stem}_biophi_oasis.xlsx"
    shutil.copyfile(xlsx, kept)

    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    if not fieldnames:
        fieldnames = ["Antibody"]

    payload = {
        "tool": "tools/biophi",
        "tool_version": TOOL_VERSION,
        "mode": "oasis",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
        },
        "parameters": {
            "min_percent_subjects": args.min_percent_subjects,
            "scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
        },
        "oasis_db": _db_info(db, with_md5=args.md5),
        "environment": env_info,
        "results": rows,
        "raw_report": str(kept.resolve()),
    }
    _write_json(outdir / f"{stem}_biophi_oasis.json", payload)
    _write_csv(outdir / f"{stem}_biophi_oasis.csv", rows, fieldnames)

    def _fmt(key: str) -> str:
        values = []
        for row in rows:
            if key in row and row[key] != "":
                values.append(str(row[key]))
        return values[0] if values else "-"

    lines = [
        f"# BioPhi OASis report: {input_path.name}",
        "",
        f"- Generated: {payload['timestamp']}",
        f"- BioPhi: {env_info.get('biophi', 'unknown')}",
        f"- OASis DB: `{payload['oasis_db'].get('path', db)}`",
        f"- Threshold: >={args.min_percent_subjects}% of OAS subjects",
        f"- Numbering: {args.scheme} / CDR definition: {args.cdr_definition}",
        f"- Raw XLSX: `{kept.resolve()}`",
        "",
        "| Antibody | OASis Identity | OASis Percentile | Germline Content "
        "| Heavy Identity | Light Identity |",
        "|----------|----------------|------------------|------------------"
        "|----------------|----------------|",
    ]
    for row in rows:
        lines.append(
            f"| {row.get('Antibody', '-')} | {row.get('OASis Identity', '-')} | "
            f"{row.get('OASis Percentile', '-')} | {row.get('Germline Content', '-')} | "
            f"{row.get('Heavy OASis Identity', '-')} | "
            f"{row.get('Light OASis Identity', '-')} |"
        )
    lines += [
        "",
        "OASis Identity is the fraction of 9-mer peptides found in the OAS "
        "repertoires of at least the chosen percentage of human subjects. "
        "OASis Percentile compares that identity with background antibodies.",
        "",
    ]
    (outdir / f"{stem}_biophi_oasis.md").write_text("\n".join(lines), encoding="utf-8")


def cmd_oasis(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    read_fasta(input_path)
    db = Path(args.oasis_db)
    if not db.exists():
        raise ToolError(f"OASis database not found: {db}")
    if db.suffix == ".gz":
        raise ToolError(
            f"OASis database must be unzipped first: gunzip {db}"
        )
    outdir = Path(args.outdir)

    with tempfile.TemporaryDirectory(prefix="biophi_oasis_") as td:
        xlsx = Path(td) / "oasis.xlsx"
        cmd = [
            "biophi", "oasis", str(input_path),
            "--oasis-db", str(db.resolve()),
            "--output", str(xlsx),
            "--scheme", args.scheme,
            "--cdr-definition", args.cdr_definition,
            "--min-percent-subjects", str(args.min_percent_subjects),
        ]
        if args.summary:
            cmd.append("--summary")
        proc = _run_biophi(cmd, args.env, args.timeout, args.hf_endpoint)
        _check_process(proc, "biophi oasis")
        if not xlsx.exists():
            raise ToolError("biophi oasis did not write the expected XLSX")
        rows = _read_oasis_overview(xlsx)
        env_info = _env_info(args.env)
        _oasis_outputs(rows, outdir, args, input_path, xlsx, db, env_info)

    stem = input_path.stem
    print(f"Wrote {outdir / (stem + '_biophi_oasis.json')}")
    print(f"Wrote {outdir / (stem + '_biophi_oasis.md')}")
    print(f"Wrote {outdir / (stem + '_biophi_oasis.csv')}")
    print(f"Wrote {outdir / (stem + '_biophi_oasis.xlsx')}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Standalone BioPhi wrapper (Sapiens + OASis) for humanness "
                    "cross-checks. Does not modify the humanization pipeline."
    )
    parser.add_argument("--tool-version", action="version", version=TOOL_VERSION)
    sub = parser.add_subparsers(dest="mode", required=True)

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--env", default=DEFAULT_ENV,
                       help=f"conda env containing biophi (default: {DEFAULT_ENV})")
        p.add_argument("--timeout", type=int, default=3600,
                       help="BioPhi subprocess timeout in seconds (default: 3600)")
        p.add_argument("--hf-endpoint", default=None,
                       help="HuggingFace endpoint/mirror; omit for offline cache mode")
        p.add_argument("--md5", action="store_true",
                       help="Compute MD5 of the OASis DB (slow for ~22 GB)")

    check = sub.add_parser("check", help="Check env, HMMER, Sapiens cache and OASis DB")
    check.add_argument("--oasis-db", default=DEFAULT_DB)
    add_common(check)
    check.set_defaults(func=cmd_check)

    sapiens = sub.add_parser("sapiens", help="Mean Sapiens humanness score")
    sapiens.add_argument("--input", required=True)
    sapiens.add_argument("--outdir", default="outputs/biophi")
    sapiens.add_argument("--scheme", default="kabat",
                         choices=["kabat", "chothia", "imgt", "aho"])
    sapiens.add_argument("--cdr-definition", default="kabat",
                         choices=["kabat", "chothia", "imgt", "north"])
    sapiens.add_argument("--humanize", action="store_true",
                         help="Also write a Sapiens-humanized FASTA")
    sapiens.add_argument("--metadata", default=None,
                         help="variants_numbering.json for Kabat-mapped "
                              "positions in scores/MD")
    sapiens.add_argument("--no-per-position", dest="no_per_position",
                         action="store_true",
                         help="Skip per-position scores + full matrix "
                              "(default: compute them)")
    sapiens.add_argument("--iterations", type=int, default=1,
                         help="Sapiens humanization iterations (default: 1)")
    sapiens.add_argument("--humanize-cdrs", action="store_true",
                         help="Allow humanizing mutations inside CDRs")
    sapiens.add_argument("--verbose", action="store_true",
                         help="Verbose subprocess output")
    add_common(sapiens)
    sapiens.set_defaults(func=cmd_sapiens)

    oasis = sub.add_parser("oasis", help="OASis 9-mer humanness identity")
    oasis.add_argument("--input", required=True)
    oasis.add_argument("--oasis-db", default=DEFAULT_DB)
    oasis.add_argument("--outdir", default="outputs/biophi")
    oasis.add_argument("--min-percent-subjects", type=float, default=10.0,
                       help="Minimum percent of OAS subjects (1-90, default: 10)")
    oasis.add_argument("--scheme", default="kabat",
                       choices=["kabat", "chothia", "imgt", "aho"])
    oasis.add_argument("--cdr-definition", default="kabat",
                       choices=["kabat", "chothia", "imgt", "north"])
    oasis.add_argument("--summary", action="store_true",
                       help="Ask BioPhi for a summary-only XLSX")
    add_common(oasis)
    oasis.set_defaults(func=cmd_oasis)

    evaluate = sub.add_parser("evaluate-variants",
                              help="Comprehensive BioPhi evaluation of pipeline variants")
    evaluate.add_argument("--input", required=True,
                          help="Variants FASTA from pipeline (variants.fasta)")
    evaluate.add_argument("--metadata", required=True,
                          help="Variants numbering JSON (variants_numbering.json)")
    evaluate.add_argument("--pipeline-result",
                          help="Pipeline result JSON (humanization_result.json)")
    evaluate.add_argument("--outdir", default="outputs/biophi_eval")
    evaluate.add_argument("--oasis-db", default=DEFAULT_DB)
    evaluate.add_argument("--scheme", default="kabat",
                          choices=["kabat", "chothia", "imgt", "aho"])
    evaluate.add_argument("--cdr-definition", default="kabat",
                          choices=["kabat", "chothia", "imgt", "north"])
    evaluate.add_argument("--run-sapiens", action="store_true", default=True,
                          help="Run Sapiens scoring (default: True)")
    evaluate.add_argument("--no-sapiens", dest="run_sapiens", action="store_false",
                          help="Skip Sapiens scoring")
    evaluate.add_argument("--run-oasis", action="store_true", default=True,
                          help="Run OASis 9-mer identity (requires --oasis-db, default: on)")
    evaluate.add_argument("--no-oasis", dest="run_oasis", action="store_false",
                          help="Skip OASis 9-mer identity")
    evaluate.add_argument("--run-humanize", action="store_true", default=True,
                          help="Generate Sapiens-humanized FASTA (default: on)")
    evaluate.add_argument("--no-humanize", dest="run_humanize", action="store_false",
                          help="Skip Sapiens-humanized FASTA")
    evaluate.add_argument("--heatmap", action="store_true", default=True,
                          help="Generate Sapiens per-position heatmap JSON (default: on)")
    evaluate.add_argument("--no-heatmap", dest="heatmap", action="store_false",
                          help="Skip heatmap JSON")
    evaluate.add_argument("--min-percent-subjects", type=float, default=10.0,
                          help="OASis minimum percent of subjects (1-90, default: 10)")
    evaluate.add_argument("--iterations", type=int, default=1,
                          help="Sapiens humanization iterations (default: 1)")
    evaluate.add_argument("--humanize-cdrs", action="store_true",
                          help="Allow humanizing mutations inside CDRs")
    evaluate.add_argument("--verbose", action="store_true",
                          help="Verbose output (show BioPhi subprocess output)")
    add_common(evaluate)
    evaluate.set_defaults(func=cmd_evaluate_variants)

    return parser


def cmd_evaluate_variants(args: argparse.Namespace) -> int:
    """Comprehensive BioPhi evaluation of pipeline variants."""
    input_path = Path(args.input)
    metadata_path = Path(args.metadata)
    pipeline_result_path = Path(args.pipeline_result) if args.pipeline_result else None
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read pipeline inputs
    records = read_fasta(input_path)
    with metadata_path.open() as fh:
        metadata = json.load(fh)
    pipeline_result = None
    if pipeline_result_path and pipeline_result_path.exists():
        with pipeline_result_path.open() as fh:
            pipeline_result = json.load(fh)

    stem = input_path.stem
    env_info = _env_info(args.env)
    numbering_variants = metadata.get("variants", {}) if isinstance(metadata, dict) else {}

    # Build a pipeline index: variant name -> chain/variant metrics
    pipeline_index: Dict[str, Dict[str, object]] = {}
    if pipeline_result:
        for chain in pipeline_result.get("chains", []):
            chain_bms = {
                bm.get("position"): bm
                for bm in chain.get("backmutations", [])
                if isinstance(bm, dict)
            }
            for variant in chain.get("variants", []):
                name = variant.get("name")
                if not name:
                    continue
                bm_positions = [
                    p for p in (variant.get("backmutations") or [])
                    if isinstance(p, str)
                ]
                composites = []
                tiers: Dict[str, int] = {}
                for p in bm_positions:
                    bm = chain_bms.get(p) or {}
                    if bm.get("composite") is not None:
                        try:
                            composites.append(float(bm["composite"]))
                        except (TypeError, ValueError):
                            pass
                    tier = bm.get("tier") or "unknown"
                    tiers[tier] = tiers.get(tier, 0) + 1
                pipeline_index[name] = {
                    "chain_name": chain.get("name"),
                    "chain_type": chain.get("chain_type"),
                    "human_likeness": chain.get("human_likeness", {}),
                    "cvi_homology": chain.get("cvi_homology"),
                    "germline_v": (chain.get("germline") or {}).get("v_gene"),
                    "n_backmutations": len(bm_positions),
                    "backmutations": bm_positions,
                    "backmutation_composite_mean": (
                        round(sum(composites) / len(composites), 4)
                        if composites else None
                    ),
                    "tier_distribution": tiers,
                }

    # ---------------------------------------------------------------------
    # 1. Sapiens scoring (mean + per-position + full matrix)
    # ---------------------------------------------------------------------
    sapiens_results: Dict[str, Any] = {}
    sapiens_detail: Dict[str, Dict[str, Any]] = {}
    humanized_fasta_path = None

    if args.run_sapiens:
        print("[1/5] Running Sapiens scoring...", file=sys.stderr, flush=True)
        with tempfile.TemporaryDirectory(prefix="biophi_eval_sapiens_") as td:
            td_path = Path(td)
            scores_csv = td_path / "sapiens_scores.csv"

            # Mean scores (authoritative BioPhi CLI numbers)
            cmd = [
                "biophi", "sapiens", str(input_path),
                "--mean-score-only", "--output", str(scores_csv),
                "--scheme", args.scheme,
                "--cdr-definition", args.cdr_definition,
            ]
            proc = _run_biophi(cmd, args.env, args.timeout, args.hf_endpoint, args.verbose)
            _check_process(proc, "biophi sapiens (mean-score-only)")
            if scores_csv.exists():
                with scores_csv.open(newline="", encoding="utf-8") as fh:
                    for row in csv.DictReader(fh):
                        rid = row.get("id", "")
                        raw = row.get("sapiens_score", "")
                        try:
                            sapiens_results[rid] = round(float(raw), 6)
                        except (TypeError, ValueError):
                            sapiens_results[rid] = raw

            # Per-position + full 20-aa matrix via Python API inside the env
            # (bypasses the broken CLI --scores-only)
            numbering_map: Dict[str, List[str]] = {}
            for rid, _seq in records:
                parts = rid.split("|")
                key = parts[1] if len(parts) > 1 else rid
                labels = numbering_variants.get(key) or numbering_variants.get(rid)
                if isinstance(labels, list):
                    numbering_map[rid] = labels
            sapiens_detail = _sapiens_per_position_python_api(
                records, env=args.env, timeout=args.timeout,
                hf_endpoint=args.hf_endpoint,
                scheme=args.scheme, cdr_definition=args.cdr_definition,
                numbering_map=numbering_map,
                verbose=args.verbose,
            )

            # Persist the full score matrix (BioPhi --scores-only equivalent)
            matrix_rows = []
            for rid, detail in sapiens_detail.items():
                for entry in detail.get("matrix", []):
                    row = {"id": rid, "chain": detail.get("chain", "")}
                    row.update(entry)
                    matrix_rows.append(row)
            if matrix_rows:
                fieldnames = ["id", "chain", "position", "numbering",
                              "input_aa",
                              "sapiens_score", *list("ACDEFGHIKLMNPQRSTVWY")]
                _write_csv(
                    outdir / f"{stem}_biophi_scores.csv", matrix_rows, fieldnames
                )

            # Humanized FASTA
            if args.run_humanize:
                humanized_path = outdir / f"{stem}_biophi_humanized.fa"
                cmd = [
                    "biophi", "sapiens", str(input_path),
                    "--fasta-only", "--output", str(humanized_path),
                    "--scheme", args.scheme,
                    "--cdr-definition", args.cdr_definition,
                    "--iterations", str(args.iterations),
                ]
                if args.humanize_cdrs:
                    cmd.append("--humanize-cdrs")
                proc = _run_biophi(cmd, args.env, args.timeout, args.hf_endpoint, args.verbose)
                _check_process(proc, "biophi sapiens (humanize)")
                if humanized_path.exists():
                    humanized_fasta_path = str(humanized_path.resolve())

    # ---------------------------------------------------------------------
    # 2. OASis 9-mer identity
    # ---------------------------------------------------------------------
    oasis_rows = []
    if args.run_oasis:
        print("[2/5] Running OASis 9-mer identity...", file=sys.stderr, flush=True)
        db = Path(args.oasis_db)
        if not db.exists():
            raise ToolError(f"OASis database not found: {db}")
        if db.suffix == ".gz":
            raise ToolError(f"OASis database must be unzipped first: gunzip {db}")

        with tempfile.TemporaryDirectory(prefix="biophi_eval_oasis_") as td:
            xlsx = Path(td) / "oasis.xlsx"
            cmd = [
                "biophi", "oasis", str(input_path),
                "--oasis-db", str(db.resolve()),
                "--output", str(xlsx),
                "--scheme", args.scheme,
                "--cdr-definition", args.cdr_definition,
                "--min-percent-subjects", str(args.min_percent_subjects),
            ]
            proc = _run_biophi(cmd, args.env, args.timeout, args.hf_endpoint, args.verbose)
            _check_process(proc, "biophi oasis")
            if xlsx.exists():
                oasis_rows = _read_oasis_overview(xlsx)
                kept = outdir / f"{stem}_biophi_oasis.xlsx"
                shutil.copyfile(xlsx, kept)

    # ---------------------------------------------------------------------
    # 3. Assemble unified results
    # ---------------------------------------------------------------------
    print("[3/5] Assembling unified report...", file=sys.stderr, flush=True)

    unified: Dict[str, Any] = {
        "tool": "tools/biophi",
        "tool_version": TOOL_VERSION,
        "mode": "evaluate-variants",
        "timestamp": _now(),
        "input": {
            "variants_fasta": str(input_path.resolve()),
            "variants_sha256": _sha256_file(input_path),
            "metadata_json": str(metadata_path.resolve()),
            "pipeline_result_json": (
                str(pipeline_result_path.resolve()) if pipeline_result_path else None
            ),
        },
        "parameters": {
            "scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
            "run_sapiens": args.run_sapiens,
            "run_oasis": args.run_oasis,
            "run_humanize": args.run_humanize,
            "heatmap": args.heatmap,
            "min_percent_subjects": args.min_percent_subjects,
            "iterations": args.iterations,
            "humanize_cdrs": args.humanize_cdrs,
        },
        "environment": env_info,
        "variants": [],
        "cross_validation": {},
    }

    oasis_by_id = {
        str(row.get("Antibody") or row.get("id")): row for row in oasis_rows
    }

    for rid, seq in records:
        variant_data: Dict[str, object] = {
            "id": rid,
            "sequence": seq,
            "length": len(seq),
        }

        # Extract variant info from ID (format: "4D5_VH|H_V2|...")
        parts = rid.split("|")
        variant_name = parts[1] if len(parts) > 1 else ""
        if len(parts) >= 2:
            variant_data["base_name"] = parts[0]
            variant_data["variant_name"] = variant_name
            variant_data["description"] = (
                "|".join(parts[2:]) if len(parts) > 2 else ""
            )

        if rid in sapiens_results:
            variant_data["sapiens_mean"] = sapiens_results[rid]
        detail = sapiens_detail.get(rid)
        if detail:
            variant_data["sapiens_chain"] = detail.get("chain")
            variant_data["sapiens_per_position"] = detail.get("per_position")
            if detail.get("mean") is not None:
                variant_data["sapiens_mean_api"] = detail.get("mean")
            if detail.get("best_score_sequence"):
                variant_data["sapiens_best_score_sequence"] = (
                    detail["best_score_sequence"]
                )
            if detail.get("masked_sequence"):
                variant_data["sapiens_masked_sequence"] = detail["masked_sequence"]
            if detail.get("hotspots"):
                variant_data["sapiens_hotspots"] = detail["hotspots"]

        oasis_row = oasis_by_id.get(rid)
        if oasis_row:
            variant_data["oasis"] = oasis_row

        pipe = pipeline_index.get(variant_name)
        if pipe:
            variant_data["pipeline"] = pipe

        if variant_name and variant_name in numbering_variants:
            variant_data["numbering"] = numbering_variants[variant_name]
        elif rid in numbering_variants:
            variant_data["numbering"] = numbering_variants[rid]

        unified["variants"].append(variant_data)

    # Cross-validation
    unified["cross_validation"] = _build_cross_validation(
        unified["variants"], pipeline_index, sapiens_results, oasis_by_id
    )

    # ---------------------------------------------------------------------
    # 4. Heatmap (Kabat-mapped per-position scores)
    # ---------------------------------------------------------------------
    if args.heatmap:
        print("[4/5] Writing heatmap...", file=sys.stderr, flush=True)
        heatmap: Dict[str, List[Dict[str, Any]]] = {}
        for v in unified["variants"]:
            rid = str(v["id"])
            labels = v.get("numbering")
            per_position = v.get("sapiens_per_position") or []
            mapped = []
            for p in per_position:
                if p.get("score") is None:
                    continue
                idx = p.get("position", 0) - 1
                label = (
                    labels[idx]
                    if isinstance(labels, list) and 0 <= idx < len(labels)
                    else p.get("position")
                )
                mapped.append({
                    "position": p.get("position"),
                    "numbering": label,
                    "residue": p.get("residue"),
                    "score": p.get("score"),
                })
            if mapped:
                heatmap[rid] = mapped
        _write_json(
            outdir / f"{stem}_biophi_heatmap.json",
            {"heatmap": heatmap, "timestamp": _now()},
        )
        # Long-format CSV companion for easy plotting / spreadsheets
        csv_rows = []
        for rid, points in heatmap.items():
            for p in points:
                csv_rows.append({
                    "id": rid,
                    "position": p.get("position"),
                    "numbering": p.get("numbering"),
                    "residue": p.get("residue"),
                    "score": p.get("score"),
                })
        if csv_rows:
            _write_csv(
                outdir / f"{stem}_biophi_heatmap.csv",
                csv_rows,
                ["id", "position", "numbering", "residue", "score"],
            )

    # ---------------------------------------------------------------------
    # 5. Write outputs
    # ---------------------------------------------------------------------
    print("[5/5] Writing outputs...", file=sys.stderr, flush=True)

    _write_json(outdir / f"{stem}_biophi_eval.json", unified)

    # CSV summary
    csv_rows = []
    for v in unified["variants"]:
        oasis = v.get("oasis") or {}
        pipe = v.get("pipeline") or {}
        row = {
            "id": v["id"],
            "variant": v.get("variant_name", ""),
            "sapiens_mean": v.get("sapiens_mean", ""),
            "oasis_identity": oasis.get("OASis Identity", ""),
            "oasis_percentile": oasis.get("OASis Percentile", ""),
            "pipeline_n_backmutations": pipe.get("n_backmutations", ""),
            "pipeline_cvi_homology": pipe.get("cvi_homology", ""),
            "humanized_fasta": (
                Path(humanized_fasta_path).name if humanized_fasta_path else ""
            ),
        }
        csv_rows.append(row)
    _write_csv(
        outdir / f"{stem}_biophi_eval.csv",
        csv_rows,
        ["id", "variant", "sapiens_mean", "oasis_identity", "oasis_percentile",
         "pipeline_n_backmutations", "pipeline_cvi_homology", "humanized_fasta"],
    )

    # Markdown report
    _write_eval_markdown(unified, outdir / f"{stem}_biophi_eval.md")

    print(f"Wrote {outdir / f'{stem}_biophi_eval.json'}", file=sys.stderr, flush=True)
    print(f"Wrote {outdir / f'{stem}_biophi_eval.md'}", file=sys.stderr, flush=True)
    print(f"Wrote {outdir / f'{stem}_biophi_eval.csv'}", file=sys.stderr, flush=True)
    if (outdir / f"{stem}_biophi_scores.csv").exists():
        print(f"Wrote {outdir / f'{stem}_biophi_scores.csv'}", file=sys.stderr, flush=True)
    if args.heatmap:
        print(f"Wrote {outdir / f'{stem}_biophi_heatmap.json'}", file=sys.stderr, flush=True)
        if (outdir / f"{stem}_biophi_heatmap.csv").exists():
            print(f"Wrote {outdir / f'{stem}_biophi_heatmap.csv'}", file=sys.stderr, flush=True)
    if humanized_fasta_path:
        print(f"Wrote {humanized_fasta_path}", file=sys.stderr, flush=True)

    return 0


def _build_cross_validation(variants, pipeline_index, sapiens_results, oasis_by_id):
    """Build cross-validation section comparing BioPhi vs pipeline metrics.

    ``pipeline_index`` maps variant names (e.g. ``H_V2``) to pipeline metrics.
    ``oasis_by_id`` maps FASTA record ids to BioPhi OASis Overview rows.
    """
    def _spearman(xs, ys):
        pairs = [
            (float(x), float(y))
            for x, y in zip(xs, ys)
            if x is not None and y is not None
        ]
        n = len(pairs)
        if n < 3:
            return None

        def _rank(values):
            order = sorted(range(n), key=lambda i: values[i])
            ranks = [0.0] * n
            i = 0
            while i < n:
                j = i
                while j + 1 < n and values[order[j + 1]] == values[order[i]]:
                    j += 1
                avg = (i + j) / 2.0 + 1.0
                for k in range(i, j + 1):
                    ranks[order[k]] = avg
                i = j + 1
            return ranks

        rx = _rank([p[0] for p in pairs])
        ry = _rank([p[1] for p in pairs])
        d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
        denom = n * (n * n - 1)
        return round(1.0 - 6.0 * d2 / denom, 4) if denom else None

    per_variant = []
    for v in variants:
        vid = v["id"]
        vname = v.get("variant_name", "")
        sapiens_score = sapiens_results.get(vid)
        if sapiens_score is None:
            sapiens_score = v.get("sapiens_mean_api")
        oasis_data = oasis_by_id.get(vid) or {}
        pipe = pipeline_index.get(vname)
        row = {
            "variant": vid,
            "variant_name": vname,
            "biophi_sapiens": sapiens_score,
            "biophi_oasis_identity": oasis_data.get("OASis Identity"),
            "biophi_oasis_percentile": oasis_data.get("OASis Percentile"),
            "pipeline": None,
        }
        if pipe:
            row["pipeline"] = {
                "chain_name": pipe.get("chain_name"),
                "human_likeness": pipe.get("human_likeness"),
                "cvi_homology": pipe.get("cvi_homology"),
                "germline_v": pipe.get("germline_v"),
                "n_backmutations": pipe.get("n_backmutations"),
                "backmutation_composite_mean": pipe.get(
                    "backmutation_composite_mean"
                ),
                "tier_distribution": pipe.get("tier_distribution"),
            }
        per_variant.append(row)

    sapiens_vals = [r["biophi_sapiens"] for r in per_variant]
    oasis_vals = [
        _maybe_float(r["biophi_oasis_identity"]) for r in per_variant
    ]
    bm_vals = [
        (r["pipeline"] or {}).get("n_backmutations") for r in per_variant
    ]
    cvi_vals = [(r["pipeline"] or {}).get("cvi_homology") for r in per_variant]

    return {
        "per_variant": per_variant,
        "spearman": {
            "sapiens_vs_oasis": _spearman(sapiens_vals, oasis_vals),
            "sapiens_vs_pipeline_n_backmutations": _spearman(
                sapiens_vals, bm_vals
            ),
            "oasis_vs_pipeline_n_backmutations": _spearman(oasis_vals, bm_vals),
            "sapiens_vs_pipeline_cvi": _spearman(sapiens_vals, cvi_vals),
        },
        "notes": [
            "Sapiens/OASis are BioPhi humanness metrics; pipeline "
            "human_likeness/CVI are germline-identity metrics (different "
            "estimators, reported side by side).",
            "Spearman rho measures rank agreement across variants.",
        ],
    }


def _maybe_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _md_cell(value):
    if value is None or value == "":
        return "N/A"
    if isinstance(value, float):
        text = f"{value:.4f}"
    else:
        text = str(value)
    return text.replace("|", "\\|")


def _write_eval_markdown(unified, path):
    """Write markdown evaluation report."""
    lines = [
        f"# BioPhi Comprehensive Evaluation: {unified['input']['variants_fasta']}",
        "",
        f"- Generated: {unified['timestamp']}",
        f"- BioPhi: {unified['environment'].get('biophi', 'unknown')}",
        f"- Sapiens: {unified['environment'].get('sapiens', 'unknown')}",
        f"- Numbering: {unified['parameters']['scheme']} / CDR: {unified['parameters']['cdr_definition']}",
        f"- Variants evaluated: {len(unified['variants'])}",
        "",
        "## Variant Summary",
        "",
        "| ID | Variant | Sapiens Mean | OASis Identity | OASis Percentile "
        "| Pipeline back-mutations |",
        "|----|---------|--------------|----------------|------------------"
        "|-------------------------|",
    ]

    for v in unified["variants"]:
        oasis = v.get("oasis") or {}
        pipe = v.get("pipeline") or {}
        lines.append(
            f"| {_md_cell(v['id'])} | {_md_cell(v.get('variant_name', ''))} | "
            f"{_md_cell(v.get('sapiens_mean'))} | "
            f"{_md_cell(oasis.get('OASis Identity'))} | "
            f"{_md_cell(oasis.get('OASis Percentile'))} | "
            f"{_md_cell(pipe.get('n_backmutations'))} |"
        )

    lines += [
        "",
        "## Per-Variant Details",
        "",
    ]

    for v in unified["variants"]:
        lines += [
            f"### {_md_cell(v['id'])}",
            f"- Variant: {v.get('variant_name', 'N/A')}",
            f"- Description: {v.get('description', 'N/A')}",
            f"- Length: {v.get('length', 'N/A')}",
            f"- Sapiens Mean: {v.get('sapiens_mean', 'N/A')}",
        ]

        oasis = v.get("oasis") or {}
        if oasis:
            lines += [
                f"- OASis Identity: {oasis.get('OASis Identity', 'N/A')}",
                f"- OASis Percentile: {oasis.get('OASis Percentile', 'N/A')}",
                f"- Germline Content: {oasis.get('Germline Content', 'N/A')}",
                f"- Heavy Identity: {oasis.get('Heavy OASis Identity', 'N/A')}",
                f"- Light Identity: {oasis.get('Light OASis Identity', 'N/A')}",
            ]

        if v.get("sapiens_per_position"):
            valid_positions = [
                p for p in v["sapiens_per_position"]
                if p.get("position") and p.get("score") is not None
            ]
            lines.append(
                f"- Sapiens per-position scores: {len(valid_positions)} valid "
                f"positions (out of {len(v['sapiens_per_position'])} total)"
            )
            if valid_positions:
                hot = sorted(valid_positions, key=lambda p: p["score"])[:5]
                hot_txt = ", ".join(
                    f"{p.get('numbering', p.get('position'))}"
                    f"({p.get('residue')})={p.get('score')}"
                    for p in hot
                )
                lines.append(f"- Sapiens hotspots (lowest): {hot_txt}")
        if v.get("sapiens_best_score_sequence"):
            lines.append(
                f"- Sapiens best-score sequence: "
                f"`{v['sapiens_best_score_sequence']}`"
            )
        if v.get("sapiens_masked_sequence"):
            lines.append(
                f"- Sapiens masked sequence: `{v['sapiens_masked_sequence']}`"
            )

        pipe = v.get("pipeline")
        if pipe:
            lines.append(
                f"- Pipeline backmutations: "
                f"{', '.join(pipe.get('backmutations') or []) or 'none'}"
            )
            hl = (pipe.get("human_likeness") or {}).get("kabat")
            lines.append(
                f"- Pipeline human-likeness (kabat): {hl} | "
                f"CVI: {pipe.get('cvi_homology')} | "
                f"germline: {pipe.get('germline_v')}"
            )

        lines.append("")

    # Cross-validation
    cv = unified.get("cross_validation", {})
    if cv.get("per_variant"):
        lines += [
            "## Cross-Validation: BioPhi vs Pipeline",
            "",
            "| Variant | BioPhi Sapiens | BioPhi OASis identity "
            "| Pipeline human-likeness (kabat) | Pipeline CVI "
            "| Pipeline back-mutations | Composite mean |",
            "|---------|----------------|----------------------"
            "|-------------------------------|--------------"
            "|-------------------------|----------------|",
        ]
        for row in cv["per_variant"]:
            pipe = row.get("pipeline") or {}
            hl = (pipe.get("human_likeness") or {}).get("kabat")
            lines.append(
                f"| {_md_cell(row['variant'])} | "
                f"{_md_cell(row.get('biophi_sapiens'))} | "
                f"{_md_cell(_maybe_float(row.get('biophi_oasis_identity')))} | "
                f"{_md_cell(hl)} | {_md_cell(pipe.get('cvi_homology'))} | "
                f"{_md_cell(pipe.get('n_backmutations'))} | "
                f"{_md_cell(pipe.get('backmutation_composite_mean'))} |"
            )
        lines += ["", "### Spearman rank correlations", ""]
        for key, value in (cv.get("spearman") or {}).items():
            lines.append(f"- {key}: {_md_cell(value)}")
        lines += ["", "### Notes", ""]
        for note in cv.get("notes", []):
            lines.append(f"- {note}")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "min_percent_subjects", None) is not None:
        if not 1 <= args.min_percent_subjects <= 90:
            parser.error("--min-percent-subjects must be between 1 and 90")
    try:
        return int(args.func(args))
    except ToolError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
