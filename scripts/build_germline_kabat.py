#!/usr/bin/env python3
"""Build the bundled human germline dataset (Kabat-space maps).

Input : data/germline/abnumber_human_imgt.json  (extracted from abnumber,
        MIT licensed; the same human IMGT germline set used by IgBLAST)
Output: data/germline/human_germline_kabat.json
        {V: {"IGHV1-69*01": {"H1": "Q", ...}}, J: {...}}
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from humanize.numbering import number_heavy, number_light

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data", "germline")
SRC = os.path.join(DATA, "abnumber_human_imgt.json")
DST = os.path.join(DATA, "human_germline_kabat.json")

IMGT_TO_CHAIN = {"H": "H", "K": "L", "L": "L"}
# generic J tails used only to give the numbering engine its FR4 anchor;
# the FR4 residues are dropped from the output map.
J_TAIL = {"H": "WGQGTLVTVSS", "L": "FGQGTKVEIK"}


def ungapped(mm):
    """Reconstruct the amino-acid sequence from a position->aa map."""
    items = sorted(mm.items(), key=lambda kv: (kv[0][0], int("".join(c for c in kv[0] if c.isdigit())), kv[0]))
    return "".join(aa for _, aa in items)


def number_j_gene(gene: str, mm: dict):
    """Map a J germline (partial CDR3 + FR4) to FR4 Kabat positions.
    Heavy: anchor W at H103 + 10 residues; Light: anchor F at L98 + 9."""
    import re
    seq = "".join(aa for _, aa in sorted(mm.items(), key=lambda kv: (kv[0][0], int("".join(c for c in kv[0] if c.isdigit())))) if aa not in ".-")
    if gene.startswith("IGHJ"):
        m = re.search(r"W[GRASLVQFKY][GRASLVQFKT]G", seq)
        if not m:
            return None
        start = m.start()
        anchor = "H103"
    else:
        m = re.search(r"FG[GQSTPE][GTKQE]", seq)
        if not m:
            return None
        start = m.start()
        anchor = "L98"
    out = {}
    base = int(anchor[1:])
    for i, aa in enumerate(seq[start : start + 11 if gene.startswith("IGHJ") else start + 10]):
        out[anchor[0] + str(base + i)] = aa
    return out


def main():
    with open(SRC) as fh:
        blob = json.load(fh)
    out = {"V": {}, "J": {}}
    skipped = []
    for kind in ("V", "J"):
        for gene, mm in blob[kind].items():
            if kind == "J":
                mapped = number_j_gene(gene, mm)
                if mapped is None:
                    skipped.append(gene)
                    continue
                out[kind][gene] = mapped
                continue
            seq = "".join(aa for _, aa in sorted(mm.items(), key=lambda kv: (kv[0][0], int("".join(c for c in kv[0] if c.isdigit())))) if aa not in ".-")
            ctype = gene[2] if gene.startswith("IG") else ""
            chain_type = IMGT_TO_CHAIN.get(ctype, "")
            if chain_type is None:
                skipped.append(gene)
                continue
            try:
                chain = number_heavy(seq + J_TAIL["H"]) if chain_type == "H" else number_light(seq + J_TAIL["L"])
            except ValueError as e:
                skipped.append(f"{gene} ({e})")
                continue
            # keep V-region residues only (drop FR4 / J tail)
            fr4_num = 103 if chain_type == "H" else 98
            out[kind][gene] = {
                r.pos: r.aa for r in chain.residues
                if int("".join(c for c in r.pos if c.isdigit())) < fr4_num and r.region != "FR4"
            }
    with open(DST, "w") as fh:
        json.dump(out, fh, indent=0)
    print(f"wrote {DST}: V={len(out['V'])} J={len(out['J'])}")
    if skipped:
        print(f"skipped {len(skipped)}: {skipped[:10]}")


if __name__ == "__main__":
    main()
