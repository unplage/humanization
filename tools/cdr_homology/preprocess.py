#!/usr/bin/env python3
"""Preprocess TheraSAbDab: extract CDRs with ANARCI IMGT numbering.

Usage:
    python3 tools/cdr_homology/preprocess.py \
        --database data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx \
        --output data/therasabdab_cdrs.json

Writes two files:
  * ``--output``            CDR database (JSON)
  * ``therasabdab_stats.json`` (next to it) with database statistics and the
    provenance of the source spreadsheet.
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cdr_core as core  # noqa: E402

# Backwards-compatible names: the numbering/CDR helpers moved into cdr_core.
_extract_cdrs_from_posmap = core.extract_cdrs_from_posmap
extract_cdrs_from_posmap = core.extract_cdrs_from_posmap
extract_cdrs = core.extract_cdrs
number_sequence_imgt = core.number_sequence_imgt

# Columns that have moved between TheraSAbDab releases are located by prefix.
COLUMNS = {
    "name": ("Therapeutic",),
    "format": ("Format",),
    "genetics": ("Genetics (Bispecifics delimited with semicolon)",),
    "target": ("Target",),
    "tech": ("Development Tech",),
    "vd_lc": ("VD LC",),
    "est_status": ("Est. Status",),
    "conditions": ("Conditions Approved",),
    "companies": ("Companies",),
}
COLUMN_PREFIXES = {
    "trial": "Highest_Clin_Trial",
}
SEQUENCE_COLUMNS = {
    "H": "HeavySequence",
    "L": "LightSequence",
    "H_bispec": "HeavySequence(ifbispec)",
    "L_bispec": "LightSequence(ifbispec)",
}


def _clean(value) -> str:
    """Return a trimmed string; NaN/None become ``''`` rather than ``'nan'``."""
    try:
        import pandas as pd
        if value is None or pd.isna(value):
            return ""
    except (ImportError, TypeError, ValueError):
        if value is None:
            return ""
    text = str(value).strip()
    return "" if text.lower() in ("na", "none", "nan", "<na>") else text


def _find_column(df, candidates=(), prefix=None):
    """Locate a spreadsheet column by exact name, then by prefix."""
    for name in candidates:
        if name in df.columns:
            return name
    if prefix:
        for name in df.columns:
            if str(name).startswith(prefix):
                return name
    return None


def _extract_chain(sequence: str, chain_type: str):
    """Number a sequence and extract its CDRs, or return ``None``."""
    posmap = core.number_sequence_imgt(sequence, chain_type)
    if not posmap:
        return None
    return core.extract_cdrs_from_posmap(posmap, chain_type)


def preprocess_excel(database_path: str, output_path: str):
    """Read TheraSAbDab Excel, extract CDRs, write the JSON index."""
    try:
        import pandas as pd
    except ImportError:
        print("Error: pandas is required (pip install pandas openpyxl)",
              file=sys.stderr)
        raise SystemExit(1)

    print(f"Reading {database_path} ...")
    df = pd.read_excel(database_path)
    print(f"  {len(df)} rows loaded")

    resolved = {key: _find_column(df, names) for key, names in COLUMNS.items()}
    trial_col = _find_column(df, prefix=COLUMN_PREFIXES["trial"])
    missing = [k for k, v in resolved.items() if v is None]
    if missing:
        print(f"  [WARN] columns not found: {', '.join(missing)}", file=sys.stderr)
    if trial_col is None:
        print("  [WARN] no Highest_Clin_Trial column found; trial phase will be blank",
              file=sys.stderr)
    else:
        print(f"  trial column: {trial_col!r}")

    database = []
    skipped = 0
    numbering_failed = 0
    n_primary = 0
    n_bispecific = 0

    for idx, row in df.iterrows():
        def cell(key, default=""):
            col = resolved.get(key)
            if not col:
                return default
            return _clean(row.get(col, default))

        entry = {
            "name": cell("name") or f"entry_{idx}",
            "format": cell("format"),
            "genetics": cell("genetics"),
            "target": cell("target"),
            "development_tech": cell("tech"),
            "vd_lc": cell("vd_lc"),
            "highest_trial": _clean(row.get(trial_col, "")) if trial_col else "",
            "est_status": cell("est_status"),
            "conditions": cell("conditions"),
            "companies": cell("companies"),
            "chains": {},
        }

        for chain_type, column in (("H", SEQUENCE_COLUMNS["H"]),
                                   ("L", SEQUENCE_COLUMNS["L"])):
            sequence = _clean(row.get(column)) if column in df.columns else ""
            if not sequence:
                continue
            cdrs = _extract_chain(sequence, chain_type)
            if not cdrs:
                numbering_failed += 1
                print(f"  [WARN] row {idx} ({entry['name']}): ANARCI could not "
                      f"number the {'heavy' if chain_type == 'H' else 'light'} chain",
                      file=sys.stderr)
                continue
            chain_key = "VH" if chain_type == "H" else "VL"
            entry["chains"][chain_key] = {"sequence": sequence, **cdrs}
            n_primary += 1

        bispecific = {}
        for chain_type, column in (("H", SEQUENCE_COLUMNS["H_bispec"]),
                                   ("L", SEQUENCE_COLUMNS["L_bispec"])):
            sequence = _clean(row.get(column)) if column in df.columns else ""
            if not sequence:
                continue
            cdrs = _extract_chain(sequence, chain_type)
            if not cdrs:
                numbering_failed += 1
                print(f"  [WARN] row {idx} ({entry['name']}): ANARCI could not "
                      f"number the bispecific chain", file=sys.stderr)
                continue
            chain_key = "VH" if chain_type == "H" else "VL"
            bispecific[chain_key] = {"sequence": sequence, **cdrs}
            n_bispecific += 1
        if bispecific:
            entry["bispecific"] = bispecific

        if entry["chains"] or entry.get("bispecific"):
            database.append(entry)
        else:
            skipped += 1

    payload = database
    core.write_json(output_path, payload)

    print(f"Done: {len(database)} entries written, {skipped} skipped, "
          f"{numbering_failed} chain(s) unnumbered")
    print(f"  chains: {n_primary} primary, {n_bispecific} bispecific")
    return database, {
        "source_file": os.path.basename(database_path),
        "source_path": os.path.abspath(database_path),
        "source_sha256": core.file_sha256(database_path),
        "source_rows": int(len(df)),
        "generated_at": core.timestamp(),
        "tool_version": core.TOOL_VERSION,
        "libraries": core.library_versions(),
        "n_primary_chains": n_primary,
        "n_bispecific_chains": n_bispecific,
        "n_skipped_rows": skipped,
        "n_unnumbered_chains": numbering_failed,
    }


def build_statistics(database: list) -> dict:
    """Compute summary statistics from the preprocessed database.

    CDR length distributions cover **every** chain in the database, primary
    and bispecific alike, so the appendix sample sizes match the real corpus.
    """
    stats = {
        "schema_version": 2,
        "total_entries": len(database),
        "format_counts": {},
        "genetics_counts": {},
        "vd_lc_counts": {},
        "trial_counts": {},
        "target_counts": {},
        "cdr_lengths": {f"{c}_{r}": [] for c in ("VH", "VL")
                        for r in ("CDR1", "CDR2", "CDR3")},
        "cdr_lengths_primary": {f"{c}_{r}": [] for c in ("VH", "VL")
                                for r in ("CDR1", "CDR2", "CDR3")},
    }

    for entry in database:
        stats["format_counts"][entry.get("format") or "Unknown"] = \
            stats["format_counts"].get(entry.get("format") or "Unknown", 0) + 1

        for value in (entry.get("genetics") or "").split(";"):
            value = value.strip()
            if value:
                stats["genetics_counts"][value] = \
                    stats["genetics_counts"].get(value, 0) + 1

        for value in (entry.get("vd_lc") or "").split(";"):
            value = value.strip()
            if value:
                stats["vd_lc_counts"][value] = stats["vd_lc_counts"].get(value, 0) + 1

        trial = entry.get("highest_trial") or "Unknown"
        stats["trial_counts"][trial] = stats["trial_counts"].get(trial, 0) + 1

        target = entry.get("target") or ""
        if target:
            stats["target_counts"][target] = stats["target_counts"].get(target, 0) + 1

        for group, is_primary in ((entry.get("chains", {}), True),
                                  (entry.get("bispecific", {}), False)):
            for chain_key, chain in group.items():
                for cdr in ("CDR1", "CDR2", "CDR3"):
                    key = f"{chain_key}_{cdr}"
                    sequence = chain.get(cdr, "")
                    if not sequence:
                        continue
                    stats["cdr_lengths"][key].append(len(sequence))
                    if is_primary:
                        stats["cdr_lengths_primary"][key].append(len(sequence))

    stats["cdr_length_summary"] = _summarise_lengths(stats["cdr_lengths"])
    stats["cdr_length_summary_primary"] = _summarise_lengths(
        stats["cdr_lengths_primary"])
    del stats["cdr_lengths"], stats["cdr_lengths_primary"]
    return stats


def _summarise_lengths(lengths_by_key: dict) -> dict:
    """min/max/mean/median/count per CDR (median via statistics.median)."""
    summary = {}
    for key, lengths in lengths_by_key.items():
        if not lengths:
            continue
        summary[key] = {
            "min": min(lengths),
            "max": max(lengths),
            "mean": round(sum(lengths) / len(lengths), 1),
            "median": statistics.median(lengths),
            "count": len(lengths),
        }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Preprocess TheraSAbDab: extract CDRs with ANARCI IMGT numbering")
    parser.add_argument(
        "--database",
        default="data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx",
        help="Path to the TheraSAbDab Excel file")
    parser.add_argument(
        "--output",
        default="data/therasabdab_cdrs.json",
        help="Output CDR database JSON")
    args = parser.parse_args()

    if not os.path.exists(args.database):
        print(f"Error: database not found: {args.database}", file=sys.stderr)
        return 1
    if not core.anarci_available():
        print("Error: ANARCI is required (pip install anarci) and needs the "
              "HMMER binaries (hmmsearch/hmmbuild) on PATH", file=sys.stderr)
        return 1

    database, provenance = preprocess_excel(args.database, args.output)

    stats = build_statistics(database)
    stats.update(provenance)
    stats_path = os.path.join(os.path.dirname(os.path.abspath(args.output)),
                              "therasabdab_stats.json")
    core.write_json(stats_path, stats)
    print(f"Statistics saved to {stats_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
