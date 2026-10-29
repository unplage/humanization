#!/usr/bin/env python3
"""Create paired VH/VL FASTA files for CDR homology / humanness analysis.

Reads a humanization pipeline ``variants.fasta`` and emits one VH/VL pair per
variant-ladder index (V0, V1, V2a, V2b, V2, V3, Vmin, ...), so the standalone
CDR-homology and humanness tools can consume it.

Pipeline ``variants.fasta`` header convention (see ``report.write_all``)::

    >AMG110.H|H_V0|pure graft: human framework + donor CDRs
    >AMG110.H|H_V1|V0 + Tier-1 back-mutations (structural pillars)
    ...
    >AMG110.L|L_V0|pure graft: ...
    >AMG110.L|L_V1|...

Output headers pair the same ladder index across chains::

    >AMG110_V0_VH
    >AMG110_V0_VL
    >AMG110_V1_VH
    ...

A single-chain input (e.g. a VHH) yields only ``_VH`` records; the downstream
tools treat them as single-chain queries.

Grouping key
------------
The grouping key is the *antibody name* plus the ladder index, so two
antibodies in one file (``>AMG110.H|H_V0`` and ``>OTHER.H|H_V0``) stay
separate instead of silently overwriting each other.  A duplicate chain type
for the same group is reported as an error rather than dropped.

Usage:
    python3 tools/cdr_homology/create_paired_fasta.py \
        --input outputs/amg110_humanization_step3/step3/variants.fasta \
        --output outputs/amg110_humanization_step3/step3/analysis_input/variants_paired.fasta \
        --name AMG110
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cdr_core as core  # noqa: E402

_CHAIN_SUFFIXES = ("_VH", "_VL", "_HEAVY", "_LIGHT")


def parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Return ``[(header, sequence), ...]`` (header without the leading '>')."""
    return core.parse_fasta(path)


def _classify(header: str) -> Tuple[str, str, str]:
    """Return ``(chain_type, antibody, index)`` for a FASTA header.

    Pipeline output ``CHAIN|ladder|desc`` (``AMG110.H|H_V2a|...``) yields the
    antibody name from the first field and the ladder index from the second.
    Any other header is split into an optional ``_VH``/``_VL`` suffix (chain
    type) and the remainder (antibody + index).
    """
    parts = [p.strip() for p in header.split("|")]
    if len(parts) >= 2 and re.match(r"^[HL]_", parts[1]):
        variant = parts[1]
        chain_type = variant[0]
        index = variant.split("_", 1)[1] if "_" in variant else variant
        antibody = core.base_name(parts[0])
        # Pipeline first field is "<antibody>.<H|L>"; drop the chain suffix.
        if re.search(r"\.[HL]$", antibody, re.IGNORECASE):
            antibody = antibody[:-2]
        return chain_type, antibody, index

    # Header has no ladder field: "<name>[_VH|_VL][ description]".
    token = parts[0].split()[0] if parts[0].split() else parts[0]
    upper = token.upper()
    if upper.endswith(("_VL", "_LIGHT")):
        chain_type = "L"
    elif upper.endswith(("_VH", "_HEAVY")):
        chain_type = "H"
    elif "VL" in upper or "LIGHT" in upper or "KAPPA" in upper or "LAMBDA" in upper:
        chain_type = "L"
    elif "VH" in upper or "HEAVY" in upper:
        chain_type = "H"
    else:
        chain_type = "H"

    antibody = core.base_name(token)
    return chain_type, antibody, ""


def build_pairs(records: List[Tuple[str, str]],
                errors: List[str] | None = None) -> Dict[str, Dict[str, str]]:
    """Group ``[(header, seq)]`` into ``{group: {chain_type: seq}}``.

    The group key is ``"<antibody>_<index>"`` (or just the antibody when there
    is no ladder index).  Insertion order follows first appearance.
    """
    if errors is None:
        errors = []
    variants: Dict[str, Dict[str, str]] = {}
    seen: Dict[Tuple[str, str], str] = {}
    for header, sequence in records:
        chain_type, antibody, index = _classify(header)
        group = f"{antibody}_{index}" if index else antibody
        slot = variants.setdefault(group, {})
        if chain_type in slot:
            errors.append(
                f"duplicate {chain_type} chain for group '{group}' "
                f"(headers '{seen[(group, chain_type)]}' and '{header}'); "
                f"the second one is ignored")
            continue
        slot[chain_type] = sequence
        seen[(group, chain_type)] = header
    return variants


def write_paired_fasta(variants: Dict[str, Dict[str, str]],
                       output_path: str,
                       name: str = "query") -> str:
    """Write ``{group: {H/L: seq}}`` as a paired FASTA; return the path."""
    outdir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(outdir, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        for group, chains in variants.items():
            base = group
            if name and not base.upper().startswith(name.upper()):
                base = f"{name}_{group}"
            if "H" in chains:
                fh.write(f">{base}_VH\n{chains['H']}\n")
            if "L" in chains:
                fh.write(f">{base}_VL\n{chains['L']}\n")
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a paired VH/VL FASTA from a humanization variants.fasta")
    parser.add_argument("--input", "-i", required=True,
                        help="pipeline variants.fasta (or any multi-sequence FASTA)")
    parser.add_argument("--output", "-o", required=True,
                        help="output paired FASTA path")
    parser.add_argument("--name", default="query",
                        help="base name for output records (default: query)")
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        parser.error(f"input FASTA not found: {args.input}")

    records = parse_fasta(args.input)
    if not records:
        parser.error(f"no sequences found in {args.input}")

    errors: List[str] = []
    variants = build_pairs(records, errors)
    path = write_paired_fasta(variants, args.output, args.name)

    n_pairs = sum(1 for c in variants.values() if "H" in c and "L" in c)
    n_single = sum(1 for c in variants.values() if ("H" in c) != ("L" in c))
    print(f"Variants: {', '.join(variants.keys())}")
    print(f"  paired VH/VL: {n_pairs}, single-chain: {n_single}")
    for message in errors:
        print(f"  [WARN] {message}", file=sys.stderr)
    print(f"Paired FASTA written: {path}")
    return 2 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
