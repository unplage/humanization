#!/usr/bin/env python3
"""Preprocess TheraSAbDab database: extract CDRs using ANARCI IMGT numbering.

Usage:
    python3 tools/cdr_homology/preprocess.py \
        --database data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx \
        --output data/therasabdab_cdrs.json
"""

import argparse
import json
import sys
import os
from pathlib import Path

import pandas as pd


def _position_number(label: str):
    """Numeric IMGT position of a label like 'H105' or 'H111A'."""
    digits = "".join(c for c in label[1:] if c.isdigit())
    return int(digits) if digits else None


def extract_cdrs_from_posmap(posmap: dict, chain_type: str) -> dict:
    """Extract CDR sequences from an IMGT-numbered position map.

    IMGT CDR boundaries (Lefranc 2003; identical to NATIVE_IMGT_CDRS in
    ``scripts/humanize/numbering.py``):
        CDR1: 27-38, CDR2: 56-65, CDR3: 105-117

    The posmap is built in ANARCI sequence order. Iterating it in that order
    (instead of numeric order with ascending insertion letters) is required
    for CDR3: IMGT places insertions as 111A, 111B, ... then 112C, 112B,
    112A and finally the 112 base, so a numeric walk transposes the loop.
    """
    imgt_cdr_ranges = {
        "CDR1": (27, 38),
        "CDR2": (56, 65),
        "CDR3": (105, 117),
    }
    cdrs = {name: "" for name in imgt_cdr_ranges}
    for label, aa in posmap.items():
        if not aa or aa == "-" or not label.startswith(chain_type):
            continue
        pos_num = _position_number(label)
        if pos_num is None:
            continue
        for cdr_name, (start, end) in imgt_cdr_ranges.items():
            if start <= pos_num <= end:
                cdrs[cdr_name] += aa
                break
    return cdrs


def number_sequence_anarci(sequence: str, chain_type: str):
    """Number a sequence using ANARCI IMGT, return posmap dict."""
    try:
        from anarci import anarci
    except ImportError:
        print("  [WARN] ANARCI not available, skipping sequence", file=sys.stderr)
        return None

    seq = sequence.upper().strip()
    try:
        result = anarci(
            [("query", seq)],
            scheme="imgt",
            output=False,
        )
    except Exception as e:
        print(f"  [WARN] ANARCI failed: {e}", file=sys.stderr)
        return None

    try:
        numbering_list = result[0][0][0][0]
    except (IndexError, TypeError):
        return None

    if not numbering_list:
        return None

    posmap = {}
    for (pos_num, ins_flag), aa in numbering_list:
        if aa == "-":
            continue
        if ins_flag and ins_flag != " ":
            label = f"{chain_type}{pos_num}{ins_flag}"
        else:
            label = f"{chain_type}{pos_num}"
        posmap[label] = aa.upper()
    return posmap


def preprocess_excel(database_path: str, output_path: str):
    """Read TheraSAbDab Excel, extract CDRs, write JSON index."""
    print(f"Reading {database_path} ...")
    df = pd.read_excel(database_path)
    print(f"  {len(df)} rows loaded")

    database = []
    skipped = 0

    for idx, row in df.iterrows():
        name = str(row.get("Therapeutic", f"entry_{idx}"))
        fmt = str(row.get("Format", ""))
        genetics = str(row.get("Genetics (Bispecifics delimited with semicolon)", ""))
        target = str(row.get("Target", ""))
        tech = str(row.get("Development Tech", ""))
        vd_lc = str(row.get("VD LC", ""))
        highest_trial = str(row.get("Highest_Clin_Trial (Feb '25)", ""))
        est_status = str(row.get("Est. Status", ""))
        conditions = str(row.get("Conditions Approved", ""))
        companies = str(row.get("Companies", ""))

        entry = {
            "name": name,
            "format": fmt,
            "genetics": genetics,
            "target": target,
            "development_tech": tech,
            "vd_lc": vd_lc,
            "highest_trial": highest_trial,
            "est_status": est_status,
            "conditions": conditions,
            "companies": companies,
            "chains": {},
        }

        # Primary chains
        for chain_type, col_name in [("H", "HeavySequence"), ("L", "LightSequence")]:
            val = row.get(col_name)
            if pd.notna(val) and str(val).strip().lower() not in ("na", "none", ""):
                sequence = str(val).strip()
                posmap = number_sequence_anarci(sequence, chain_type)
                if posmap:
                    cdrs = extract_cdrs_from_posmap(posmap, chain_type)
                    chain_key = "VH" if chain_type == "H" else "VL"
                    entry["chains"][chain_key] = {
                        "sequence": sequence,
                        **cdrs,
                    }

        # Bispecific chains
        bispec = {}
        for chain_type, col_name in [("H", "HeavySequence(ifbispec)"), ("L", "LightSequence(ifbispec)")]:
            val = row.get(col_name)
            if pd.notna(val) and str(val).strip().lower() not in ("na", "none", ""):
                sequence = str(val).strip()
                posmap = number_sequence_anarci(sequence, chain_type)
                if posmap:
                    cdrs = extract_cdrs_from_posmap(posmap, chain_type)
                    chain_key = "VH" if chain_type == "H" else "VL"
                    bispec[chain_key] = {
                        "sequence": sequence,
                        **cdrs,
                    }
        if bispec:
            entry["bispecific"] = bispec

        if entry["chains"]:
            database.append(entry)
        else:
            skipped += 1

    # Write output
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(database, f, indent=2, ensure_ascii=False)

    print(f"Done: {len(database)} entries written, {skipped} skipped")
    return database


def build_statistics(database: list) -> dict:
    """Compute summary statistics from the preprocessed database."""
    stats = {
        "total_entries": len(database),
        "format_counts": {},
        "genetics_counts": {},
        "vd_lc_counts": {},
        "trial_counts": {},
        "target_list": [],
        "cdr_lengths": {
            "VH_CDR1": [],
            "VH_CDR2": [],
            "VH_CDR3": [],
            "VL_CDR1": [],
            "VL_CDR2": [],
            "VL_CDR3": [],
        },
    }

    for entry in database:
        # Format
        fmt = entry.get("format", "Unknown")
        stats["format_counts"][fmt] = stats["format_counts"].get(fmt, 0) + 1

        # Genetics
        genetics = entry.get("genetics", "Unknown")
        for g in genetics.split(";"):
            g = g.strip()
            if g:
                stats["genetics_counts"][g] = stats["genetics_counts"].get(g, 0) + 1

        # VD LC
        vd_lc = entry.get("vd_lc", "Unknown")
        for v in vd_lc.split(";"):
            v = v.strip()
            if v:
                stats["vd_lc_counts"][v] = stats["vd_lc_counts"].get(v, 0) + 1

        # Clinical trial
        trial = entry.get("highest_trial", "Unknown")
        stats["trial_counts"][trial] = stats["trial_counts"].get(trial, 0) + 1

        # Target
        target = entry.get("target", "")
        if target and target not in ("nan", "None", ""):
            stats["target_list"].append(target)

        # CDR lengths
        for chain_key in ["VH", "VL"]:
            chain = entry.get("chains", {}).get(chain_key)
            if not chain:
                continue
            for cdr in ["CDR1", "CDR2", "CDR3"]:
                key = f"{chain_key}_{cdr}"
                seq = chain.get(cdr, "")
                if seq:
                    stats["cdr_lengths"][key].append(len(seq))

    # Deduplicate targets
    stats["unique_targets"] = sorted(set(stats["target_list"]))
    stats["target_counts"] = {}
    for t in stats["target_list"]:
        stats["target_counts"][t] = stats["target_counts"].get(t, 0) + 1
    del stats["target_list"]

    # CDR length summary
    stats["cdr_length_summary"] = {}
    for key, lengths in stats["cdr_lengths"].items():
        if lengths:
            stats["cdr_length_summary"][key] = {
                "min": min(lengths),
                "max": max(lengths),
                "mean": round(sum(lengths) / len(lengths), 1),
                "median": sorted(lengths)[len(lengths) // 2],
                "count": len(lengths),
            }

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Preprocess TheraSAbDab: extract CDRs with ANARCI IMGT numbering"
    )
    parser.add_argument(
        "--database",
        default="data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx",
        help="Path to TheraSAbDab Excel file",
    )
    parser.add_argument(
        "--output",
        default="data/therasabdab_cdrs.json",
        help="Output JSON path",
    )
    args = parser.parse_args()

    database = preprocess_excel(args.database, args.output)

    # Also save statistics
    stats = build_statistics(database)
    stats_path = Path(args.output).with_name("therasabdab_stats.json")
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"Statistics saved to {stats_path}")


if __name__ == "__main__":
    main()
