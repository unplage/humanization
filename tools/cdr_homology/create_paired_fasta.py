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

Output headers pair the same ladder index across chains::

    >AMG110_V0_VH
    >AMG110_V0_VL
    >AMG110_V1_VH
    ...

A single-chain input (e.g. a VHH) yields only ``_VH`` records; the downstream
tools treat them as single-chain queries.

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
from typing import Dict, List, Tuple


def parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Return ``[(header, sequence), ...]`` (header without the leading '>')."""
    records: List[Tuple[str, str]] = []
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:], []
            else:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


_CHAIN_SUFFIXES = ("_VH", "_VL", "_HEAVY", "_LIGHT")


def _sanitize(text: str) -> str:
    """Make a header token safe for a FASTA record name."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", text).strip("_") or "query"


def _classify(header: str) -> Tuple[str, str]:
    """Return ``(chain_type, index)`` for a FASTA header.

    Pipeline output uses ``CHAIN|ladder|desc`` (``AMG110.H|H_V2a|...``) and is
    grouped by the ladder index (``V2a``). Any other header is sanitized and
    used as the index with a known ``_VH``/``_VL``/``_Heavy``/``_Light`` suffix
    stripped (so ``4D5_VH`` + ``4D5_VL`` pair as ``4D5``).
    """
    parts = [p.strip() for p in header.split("|")]
    if len(parts) >= 2 and re.match(r"^[HL]_", parts[1]):
        variant = parts[1]
        return variant[0], (variant.split("_", 1)[1] if "_" in variant else variant)

    up = header.upper()
    if "VH" in up or "HEAVY" in up:
        chain_type = "H"
    elif "VL" in up or "LIGHT" in up or "KAPPA" in up or "LAMBDA" in up:
        chain_type = "L"
    else:
        chain_type = "H"

    base = _sanitize(parts[0].split()[0])
    for suffix in _CHAIN_SUFFIXES:
        if base.upper().endswith(suffix):
            base = base[: -len(suffix)]
            break
    return chain_type, base


def build_pairs(records: List[Tuple[str, str]]) -> "Dict[str, Dict[str, str]]":
    """Group ``[(header, seq)]`` into ``{index: {chain_type: seq}}``.

    Insertion order follows the first appearance, which for pipeline output is
    the ladder order (V0, V1, V2a, ...).
    """
    variants: Dict[str, Dict[str, str]] = {}
    for header, sequence in records:
        chain_type, index = _classify(header)
        slot = variants.setdefault(index, {})
        slot[chain_type] = sequence
    return variants


def write_paired_fasta(
    variants: Dict[str, Dict[str, str]],
    output_path: str,
    name: str = "query",
) -> str:
    """Write ``{index: {H/L: seq}}`` as a paired FASTA; return the path."""
    outdir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(outdir, exist_ok=True)
    with open(output_path, "w") as fh:
        for index, chains in variants.items():
            if not name or index.upper().startswith(name.upper()):
                base = index
            else:
                base = f"{name}_{index}"
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

    variants = build_pairs(records)
    path = write_paired_fasta(variants, args.output, args.name)

    n_pairs = sum(1 for c in variants.values() if "H" in c and "L" in c)
    n_single = sum(1 for c in variants.values() if ("H" in c) != ("L" in c))
    print(f"Variants: {', '.join(variants.keys())}")
    print(f"  paired VH/VL: {n_pairs}, single-chain: {n_single}")
    print(f"Paired FASTA written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
