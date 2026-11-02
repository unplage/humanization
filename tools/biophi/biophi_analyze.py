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
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

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
                hf_endpoint: Optional[str]) -> subprocess.CompletedProcess:
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
        try:
            return subprocess.run(
                list(args), capture_output=True, text=True,
                timeout=timeout, env=proc_env
            )
        except subprocess.TimeoutExpired as exc:
            raise ToolError(f"BioPhi timed out after {timeout}s: {' '.join(args)}") from exc
    return _conda_run(env, args, timeout=timeout, extra_env=extra_env)


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
                     input_path: Path, env_info: Dict[str, object]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    rows = []
    for name, seq in records:
        score = scores.get(name, {})
        rows.append({
            "id": name,
            "length": len(seq),
            "chain": score.get("chain", ""),
            "sapiens_score": score.get("sapiens_score", ""),
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
    }
    _write_json(outdir / f"{stem}_biophi_sapiens.json", payload)
    _write_csv(outdir / f"{stem}_biophi_sapiens.csv", rows,
               ["id", "length", "chain", "sapiens_score"])

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
            f"| {row['id']} | {row['length']} | {row['chain']} | "
            f"{row['sapiens_score']} |"
        )
    lines += [
        "",
        "The Sapiens score is the mean per-residue probability under the human "
        "antibody language model. Higher is more human-like; it is not a "
        "germline identity percentage.",
        "",
    ]
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
            ]
            humanize_proc = _run_biophi(
                humanize_cmd, args.env, args.timeout, args.hf_endpoint
            )
            _check_process(humanize_proc, "biophi sapiens --fasta-only")

    env_info = _env_info(args.env)
    _sapiens_outputs(records, scores, outdir, args, input_path, env_info)

    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.json')}")
    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.md')}")
    print(f"Wrote {outdir / (input_path.stem + '_biophi_sapiens.csv')}")
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

    return parser


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
