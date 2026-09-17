#!/usr/bin/env python3
"""Create a combined VH/VL FASTA (and per-combo Kabat numbering) for
immunogenicity analysis of combinatorial variants.

Reads the FASTA files produced by
``scripts/humanize/combinatorial_igfold.py`` (``combo_<H>_<L>.fasta``) and emits
one combined file whose records are named ``VH|combo_<H>_<L>`` and
``VL|combo_<H>_<L>``. The immunogenicity analyzer groups these by chain/variant
and the structural-risk tool pairs the variant name with ``combo_<H>_<L>.pdb``.

Optionally emits a chain-aware per-combo Kabat numbering JSON
(``{"combo_V0_V1|H": [...], "combo_V0_V1|L": [...]}``) from the pipeline's
``variants_numbering.json`` so FR4 / insertion-letter positions match exactly.

Usage:
    python3 tools/immunogenicity/create_combined_fasta.py \
        --input-dir af3/amg110_step3_combinatorial \
        --output outputs/.../immunogenicity_input/all_49_variants.fasta \
        --numbering-source outputs/.../step3/variants_numbering.json \
        --numbering-output outputs/.../immunogenicity_input/combo_numbering.json
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
from typing import Dict, List, Optional, Tuple


def parse_fasta(path: str) -> List[Tuple[str, str]]:
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


def read_combo(path: str) -> Tuple[str, Optional[str], Optional[str], Optional[str], Optional[str]]:
    """Return ``(combo_name, h_variant, vh_seq, l_variant, vl_seq)``.

    ``combo_name`` is the file stem without the ``combo_`` prefix; the H/L
    variant names come from the pipeline header (``H|H_V0|...``).
    """
    combo = re.sub(r"^combo_", "", os.path.splitext(os.path.basename(path))[0])
    h_var = l_var = None
    vh = vl = None
    for hdr, seq in parse_fasta(path):
        parts = [p.strip() for p in hdr.split("|")]
        chain = parts[0].upper() if parts else ""
        variant = parts[1] if len(parts) > 1 else None
        if chain.startswith("H"):
            vh, h_var = seq, variant
        elif chain.startswith("L"):
            vl, l_var = seq, variant
    return combo, h_var, vh, l_var, vl


def build_combined(input_dir: str, output: str) -> List[Tuple[str, Optional[str], Optional[str]]]:
    """Write the combined FASTA; return ``[(combo, h_var, l_var), ...]``."""
    files = sorted(glob.glob(os.path.join(input_dir, "combo_*.fasta")))
    if not files:
        raise SystemExit(f"no combo_*.fasta found in {input_dir}")

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    pairs: List[Tuple[str, Optional[str], Optional[str]]] = []
    with open(output, "w") as out:
        for path in files:
            combo, h_var, vh, l_var, vl = read_combo(path)
            if vh:
                out.write(f">VH|combo_{combo}|{combo} heavy chain\n{vh}\n")
            if vl:
                out.write(f">VL|combo_{combo}|{combo} light chain\n{vl}\n")
            pairs.append((combo, h_var, l_var))
    return pairs


def build_combo_numbering(
    pairs: List[Tuple[str, Optional[str], Optional[str]]],
    numbering_source: str, output: str,
) -> int:
    """Map per-chain pipeline labels onto per-combo, chain-aware keys."""
    with open(numbering_source) as fh:
        base = json.load(fh)

    combo_num: Dict[str, List[str]] = {}
    for combo, h_var, l_var in pairs:
        if h_var and h_var in base:
            combo_num[f"combo_{combo}|H"] = base[h_var]
        if l_var and l_var in base:
            combo_num[f"combo_{combo}|L"] = base[l_var]

    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w") as fh:
        json.dump(combo_num, fh, indent=2)
    return len(combo_num)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build combined FASTA + per-combo numbering for immunogenicity")
    parser.add_argument("--input-dir", "-i", required=True,
                        help="directory with combo_*.fasta")
    parser.add_argument("--output", "-o", required=True,
                        help="output combined FASTA path")
    parser.add_argument("--numbering-source",
                        help="pipeline variants_numbering.json (per-chain labels)")
    parser.add_argument("--numbering-output",
                        help="output chain-aware per-combo numbering JSON")
    args = parser.parse_args()

    pairs = build_combined(args.input_dir, args.output)
    n_h = sum(1 for _, h, _l in pairs if h)
    n_l = sum(1 for _, _h, l in pairs if l)
    print(f"Combos: {len(pairs)} (VH: {n_h}, VL: {n_l})")
    print(f"Combined FASTA written: {args.output}")

    if args.numbering_source:
        if not os.path.isfile(args.numbering_source):
            raise SystemExit(f"numbering source not found: {args.numbering_source}")
        out = args.numbering_output or os.path.splitext(args.output)[0] + "_numbering.json"
        n = build_combo_numbering(pairs, args.numbering_source, out)
        print(f"Per-combo numbering ({n} entries) written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
