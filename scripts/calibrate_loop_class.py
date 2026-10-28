#!/usr/bin/env python3
"""Calibrate the per-antibody loop-class thresholds (D3) on real structures.

The loop-class correction (scripts/humanize/loopclass.py) uses three
distance cutoffs:

  * ``LOOP_CLASS["shell"]`` (8 A) - a literature canonical claim is credible
    only if the loop it claims to support is within this distance;
  * ``LOOP_CLASS["far"]``  (12 A) - a literature vernier claim is dropped if
    the position is farther than this from every CDR;
  * ``LOOP_CLASS["direct"]`` (4.5 A) - side-chain contact cutoff separating a
    real contact (already reported as ``cdr_contact_sc``) from the second
    shell (the weak ``loop_shell`` evidence).

This script measures those relations on a set of PDB files so the cutoffs can
be justified (and re-checked) with data instead of intuition. It reports, for
every antibody chain it can number:

  * canonical support positions -> distance to the loop they support;
  * vernier positions -> distance to the nearest CDR;
  * how many framework positions would receive the weak ``loop_shell``
    evidence (the precision of the addition rule).

Usage::

    python3 scripts/calibrate_loop_class.py af3/humab25_parents/*.pdb
    python3 scripts/calibrate_loop_class.py --glob 'af3/*/*.pdb'

The distribution the shipped thresholds were set from (86 chains):
canonical median 3.2 A (3% beyond 8 A), vernier median 3.1 A (0.4% beyond
12 A), loop-shell candidates 26% of all framework positions.
"""
from __future__ import annotations

import argparse
import glob as globmod
import os
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))

from humanize import config
from humanize.config import (CANONICAL, CANONICAL_SUPPORT_BY_LOOP,
                             INTERFACE_CORE, INTERFACE_EXTENDED,
                             VERNIER_ZONE)
from humanize.fr_indel import FR_REGIONS
from humanize.graft import is_cdr_loop_position, loop_of_position
from humanize.numbering import number_heavy, number_light
from humanize.sequences import _anarci_chain_type
from humanize.structure import compute_hints, load_model

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def _num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def _in_fr(chain_type: str, n: int) -> bool:
    return any(lo <= n <= hi for lo, hi in FR_REGIONS[chain_type].values())


def _chain_sequence(atoms):
    res = {a.resseq: a.resname for a in atoms}
    return "".join(AA3.get(res[r], "X") for r in sorted(res))


def analyse(pdb_path: str):
    """Return (canonical_dists, vernier_dists, fr_total, shell_candidates)."""
    model = load_model(pdb_path)
    if not model:
        return None
    chains = {}
    for a in model.atoms:
        chains.setdefault(a.chain, []).append(a)

    canon, vernier = [], []
    fr_total = shell_candidates = 0
    for chain_label, atoms in chains.items():
        seq = _chain_sequence(atoms)
        if len(seq) < 80:
            continue
        # chain typing: ANARCI when available (most reliable), otherwise try
        # heavy and light and keep the interpretation with sensible CDRs
        anarci_type = _anarci_chain_type(seq)
        if anarci_type in ("H", "K"):
            order = (("H", number_heavy),)
        elif anarci_type == "L":
            order = (("L", number_light),)
        else:
            order = (("H", number_heavy), ("L", number_light))
        for chain_type, numberer in order:
            try:
                numbered = numberer(seq)
            except Exception:
                continue
            if numbered is None or len(numbered.residues) < 80:
                continue
            pos_to_resseq = {r.pos: r.index + 1 for r in numbered.residues}
            cdrs = {p: n for p, n in pos_to_resseq.items()
                    if is_cdr_loop_position(chain_type, _num(p), "kabat", "kabat")}
            if not cdrs:
                continue
            loops = {}
            for p in cdrs:
                loop = loop_of_position(chain_type, _num(p), "kabat", "kabat")
                if loop:
                    loops[p] = loop
            hints = compute_hints(model, chain_label, pos_to_resseq, cdrs,
                                  cdr_loops=loops)
            if not hints.has_loop_geometry():
                break
            for res in numbered.residues:
                pos, n = res.pos, _num(res.pos)
                for loop, positions in CANONICAL_SUPPORT_BY_LOOP.get(chain_type, {}).items():
                    if n in positions:
                        d = hints.loop_min_dist(pos).get(loop)
                        if d is not None:
                            canon.append(d)
                if n in VERNIER_ZONE.get(chain_type, set()):
                    dd = hints.loop_min_dist(pos)
                    if dd:
                        vernier.append(min(dd.values()))
                if (n in CANONICAL.get(chain_type, set())
                        or n in VERNIER_ZONE.get(chain_type, set())
                        or n in INTERFACE_CORE.get(chain_type, set())
                        or n in INTERFACE_EXTENDED.get(chain_type, set())):
                    continue
                if not _in_fr(chain_type, n):
                    continue
                fr_total += 1
                sc = hints.loop_min_dist_sc(pos)
                if not any(config.LOOP_CLASS["direct"] < d <= config.LOOP_CLASS["shell"]
                           for d in sc.values()):
                    continue
                buried = hints.buried(chain_type, pos)
                rel = hints.rel_sasa(pos)
                if not (buried is True
                        or (rel is not None and rel <= config.LOOP_CLASS["shell_rel_sasa"])):
                    continue
                shell_candidates += 1
            break
    return canon, vernier, fr_total, shell_candidates


def _report(name, values, cutoffs):
    if not values:
        print(f"  {name}: no data")
        return
    xs = sorted(values)
    q = lambda f: xs[min(len(xs) - 1, int(f * len(xs)))]
    print(f"  {name}: n={len(xs)} median={statistics.median(xs):.1f} "
          f"q25={q(.25):.1f} q75={q(.75):.1f} min={xs[0]:.1f} max={xs[-1]:.1f}")
    for c in cutoffs:
        over = sum(1 for x in xs if x > c)
        print(f"      > {c:>5} A: {over}/{len(xs)} = {100.0 * over / len(xs):.1f}%")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdbs", nargs="*", help="PDB files to analyse")
    ap.add_argument("--glob", default=None, help="glob pattern for PDB files")
    args = ap.parse_args()

    paths = list(args.pdbs)
    if args.glob:
        paths += sorted(globmod.glob(args.glob))
    paths = [p for p in dict.fromkeys(paths) if os.path.exists(p)]
    if not paths:
        print("no PDB files given (use positional arguments or --glob)")
        return 1

    canon, vernier = [], []
    fr_total = shell_candidates = 0
    used = 0
    for p in paths:
        res = analyse(p)
        if not res:
            continue
        c, v, ft, sc = res
        canon += c
        vernier += v
        fr_total += ft
        shell_candidates += sc
        used += 1

    print(f"analysed {used} structures ({len(canon)} canonical / "
          f"{len(vernier)} vernier support relations)")
    print(f"thresholds: shell={config.LOOP_CLASS['shell']} A  "
          f"far={config.LOOP_CLASS['far']} A  "
          f"direct={config.LOOP_CLASS['direct']} A")
    _report("canonical -> loop it supports", canon,
            [config.LOOP_CLASS["shell"], config.LOOP_CLASS["far"]])
    _report("vernier -> nearest CDR", vernier,
            [config.LOOP_CLASS["shell"], config.LOOP_CLASS["far"]])
    if fr_total:
        print(f"  loop-shell candidates: {shell_candidates}/{fr_total} = "
              f"{100.0 * shell_candidates / fr_total:.1f}% of framework "
              f"positions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
