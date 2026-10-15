#!/usr/bin/env python3
"""Derive IMGT-space position sets from the paired germline files.

Input : data/germline/human_germline_kabat.json  (Kabat-space maps)
        data/germline/abnumber_human_imgt.json   (IMGT-space maps)
Output: data/germline/kabat_imgt_positions.json

Method: both files describe the same genes residue-for-residue (the Kabat
file was generated from the IMGT file by scripts/build_germline_kabat.py).
For every paired V gene the two posmaps are ordered as sequence (Kabat label
sort; IMGT label sort incl. '-' gap placeholders skipped) and zipped, giving
a per-residue (Kabat label -> IMGT label) vote across 373 genes.

Kabat<->IMGT label correspondence is NOT a global bijection: IMGT CDR gaps
(31-34, 58-60, ...) and Kabat insertions make CDR-boundary positions
context-dependent. Majority votes are therefore used as follows:

  * frac >= 0.95 and a bare-int target  -> that single IMGT int
  * otherwise                           -> union of all observed bare-int
                                           images (letter images collapse to
                                           their base int), i.e. a conservative
                                           over-inclusive membership set
  * J/FR4 positions (no V-gene votes)   -> definition rule
                                           H: imgt = kabat - 103 + 118
                                           L: imgt = kabat -  98 + 118

Only the config position sets are translated (they are the only bare-int
sets the engine probes); the raw votes for every translated label are kept
for audit.
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import humanize.config as config  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "germline")
SRC_KABAT = os.path.join(DATA, "human_germline_kabat.json")
SRC_IMGT = os.path.join(DATA, "abnumber_human_imgt.json")
DST = os.path.join(DATA, "kabat_imgt_positions.json")

THRESHOLD = 0.95
VERSION = 1


def _order(mm):
    return sorted(mm.items(),
                  key=lambda kv: (kv[0][0],
                                  int("".join(c for c in kv[0] if c.isdigit())),
                                  kv[0]))


def build_votes():
    with open(SRC_KABAT) as fh:
        kabat = json.load(fh)
    with open(SRC_IMGT) as fh:
        imgt = json.load(fh)
    votes = collections.defaultdict(collections.Counter)
    n_paired = 0
    skipped = []
    for gene in sorted(set(kabat["V"]) & set(imgt["V"])):
        kk = [(p, a) for p, a in _order(kabat["V"][gene]) if a not in ".-"]
        mm = [(p, a) for p, a in _order(imgt["V"][gene]) if a not in ".-"]
        if len(kk) != len(mm):
            skipped.append(gene)
            continue
        for (kp, _), (ip, _) in zip(kk, mm):
            votes[kp][ip] += 1
        n_paired += 1
    return votes, n_paired, skipped


def translate(chain: str, x: int, votes) -> dict:
    """Return {"targets": [imgt ints], "votes": {...}, "rule": str}."""
    # FR4 / J positions never appear in V genes: definition rule.
    if chain == "H" and x >= 103:
        return {"targets": [118 + (x - 103)], "votes": {}, "rule": "j_region_def"}
    if chain == "L" and x >= 98:
        return {"targets": [118 + (x - 98)], "votes": {}, "rule": "j_region_def"}

    lab = f"{chain}{x}"
    v = votes.get(lab)
    if not v:
        return {"targets": [x], "votes": {}, "rule": "passthrough_no_evidence"}
    tot = sum(v.values())
    top, n = v.most_common(1)[0]
    raw = {k: c for k, c in v.most_common()}
    frac = n / tot
    top_int = int("".join(c for c in top if c.isdigit())) if any(
        c.isdigit() for c in top) else None
    if frac >= THRESHOLD and top[1:].isdigit():
        return {"targets": [int(top[1:])], "votes": raw, "rule": "majority"}
    # union of bare images (letters collapse to their base int)
    union = []
    for a in v:
        b = int("".join(c for c in a if c.isdigit())) if any(
            c.isdigit() for c in a) else None
        if b is not None and b not in union:
            union.append(b)
    if not union and top_int is not None:
        union = [top_int]
    return {"targets": sorted(union), "votes": raw,
            "rule": "union", "top_frac": round(frac, 3)}


def translate_set(name: str, values, votes):
    out = {}
    for chain in ("H", "L"):
        vals = values.get(chain) if isinstance(values, dict) else None
        if vals is None:
            continue
        seen = []
        audit = {}
        for x in sorted(vals):
            r = translate(chain, x, votes)
            audit[str(x)] = r
            for t in r["targets"]:
                if t not in seen:
                    seen.append(t)
        out[chain] = sorted(seen)
        out.setdefault("_audit", {})[chain] = audit
    return out


def main():
    votes, n_paired, skipped = build_votes()
    sets_in = {
        "VERNIER_ZONE": config.VERNIER_ZONE,
        "INTERFACE_CORE": config.INTERFACE_CORE,
        "INTERFACE_EXTENDED": config.INTERFACE_EXTENDED,
        "CANONICAL": config.CANONICAL,
        "EMPIRICAL_NO_EFFECT": config.EMPIRICAL_NO_EFFECT,
    }
    out_sets = {}
    audit = {}
    for name, vals in sets_in.items():
        tr = translate_set(name, vals, votes)
        out_sets[name] = {k: v for k, v in tr.items() if not k.startswith("_")}
        audit[name] = tr.get("_audit", {})
    out_sets["VHH_HALLMARK"] = translate_set(
        "VHH_HALLMARK", {"H": sorted(config.VHH_HALLMARK)}, votes
    )["H"]

    blob = {
        "version": VERSION,
        "method": f"residue-zip majority vote over {n_paired} paired V genes; "
                  f"frac>={THRESHOLD} single target else union of images",
        "threshold": THRESHOLD,
        "paired_genes": n_paired,
        "skipped_genes": skipped,
        "sets": out_sets,
        "audit": audit,
    }
    with open(DST, "w") as fh:
        json.dump(blob, fh, indent=1, sort_keys=True)
    print(f"wrote {DST} (paired={n_paired}, skipped={len(skipped)})")
    # summary of weak translations
    weak = 0
    for name, ch in audit.items():
        for chain, positions in ch.items():
            for x, r in positions.items():
                if r["rule"] != "majority" and r["rule"] != "j_region_def":
                    weak += 1
                    print(f"  {name} {chain}{x}: {r['rule']} -> {r['targets']}"
                          + (f" frac={r.get('top_frac')}" if "top_frac" in r else ""))
    print(f"weak/passthrough translations: {weak}")


if __name__ == "__main__":
    main()
