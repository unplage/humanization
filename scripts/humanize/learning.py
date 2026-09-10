"""Closed-loop learning from experimental humanization data.

After each humanization round you record (parent sequence, variant sequences,
measured affinities). This module converts those records into empirical
per-position effects and recalibrates the scoring:

  1. parse experiments (JSON): parent + variants + KD (optional Tm/expr)
  2. per-position effect: within-parent contrast
     effect(p) = mean(ΔΔG of variants carrying the DONOR residue at p)
               - mean(ΔΔG of variants carrying the HUMAN residue at p)
     ΔΔG = RT*ln(KD_var/KD_parent), RT = 0.593 kcal/mol (25 C)
  3. shrink effects toward 0 by sample size (empirical Bayes-lite)
  4. write calibration.json: per-position effects + feature/tier adjustments
  5. pipeline consumes calibration: composite modifiers + tier overrides +
     report column "empirical ΔΔG"

Design rules for the experimental data (see docs/learning_loop.md):
  * include V0 (pure graft) and V2 (T1+T2) at minimum, plus single-position
    variants for high-value positions;
  * KD must be measured in the same assay/format as the parent;
  * record the format (IgG/Fab/VHH-Fc) - effects are format-relative.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .germline import GermlineDB
from .numbering import number_heavy, number_light

RT_KCAL = 0.593   # RT at 25 C

# ---------------------------------------------------------------------------
# data model
# ---------------------------------------------------------------------------

@dataclass
class ExperimentRecord:
    name: str
    parent_vh: str
    parent_vl: Optional[str]
    parent_kd: float                # nM
    variants: List[dict] = field(default_factory=list)  # {name, seq, kd, tm?, expression?}

    def variant(self, name: str) -> Optional[dict]:
        for v in self.variants:
            if v.get("name") == name:
                return v
        return None


def parse_experiments(path: str) -> List[ExperimentRecord]:
    """Load experiments JSON.

    Schema:
      [{
        "name": "Ab1",
        "parent_vh": "...", "parent_vl": "...",
        "parent_kd": 0.15,
        "variants": [
          {"name": "Ab1_hV2", "seq": "...", "kd": 0.20, "tm": 68.0}
        ]
      }]
    seq = VH sequence for VHH; for Fab provide "vh"+"vl" keys instead.
    """
    with open(path) as fh:
        data = json.load(fh)
    records = []
    for item in data:
        variants = []
        for v in item.get("variants", []):
            if not (v.get("seq") or v.get("vh") or v.get("vl")):
                raise ValueError(f"variant {v.get('name')} missing sequence")
            variants.append(dict(v))
        records.append(ExperimentRecord(
            name=item["name"],
            parent_vh=item.get("parent_vh", ""),
            parent_vl=item.get("parent_vl"),
            parent_kd=float(item["parent_kd"]),
            variants=variants,
        ))
    return records


# ---------------------------------------------------------------------------
# per-position effect estimation
# ---------------------------------------------------------------------------

@dataclass
class PositionEffect:
    position: str
    donor_aa: str
    human_aa: str
    effect: float            # kcal/mol; >0 => donor is better (revert helps)
    n: int
    raw_n: int = 0


def _ddg(kd_var: float, kd_parent: float) -> float:
    if kd_var <= 0 or kd_parent <= 0:
        return 0.0
    return RT_KCAL * math.log(kd_var / kd_parent)


def _solve_linear(A: List[List[float]], b: List[float]) -> List[float]:
    """Gauss-Jordan solve of A x = b (small dense system, pure stdlib)."""
    n = len(b)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(M[r][col]))
        if abs(M[piv][col]) < 1e-12:
            continue
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [x / pv for x in M[col]]
        for r in range(n):
            if r != col and abs(M[r][col]) > 1e-15:
                f = M[r][col]
                M[r] = [M[r][k] - f * M[col][k] for k in range(n + 1)]
    return [M[i][n] for i in range(n)]


def _ridge_solve(rows: List[List[str]], y: List[float],
                 pos_list: List[str], lam: float = 1.0) -> List[float]:
    """Ridge regression of ddG on the 'carries human residue' indicators.

    Deconvolves co-occurring framework substitutions: with a single V2-style
    variant a per-position mean would smear the whole ddG across every changed
    position; least squares attributes it to the positions actually present in
    the variant set.
    """
    p = len(pos_list)
    idx = {pos: i for i, pos in enumerate(pos_list)}
    A = [[0.0] * p for _ in range(p)]
    b = [0.0] * p
    for row, yy in zip(rows, y):
        cols = [idx[pos] for pos in row]
        for i in cols:
            b[i] += yy
            for j in cols:
                A[i][j] += 1.0
    for i in range(p):
        A[i][i] += lam
    return _solve_linear(A, b)


def compute_position_effects(
    records: List[ExperimentRecord],
    db: GermlineDB,
) -> Tuple[Dict[str, PositionEffect], List[str]]:
    """Estimate per-position effects. Returns (effects, warnings).

    For each experiment the ddG of every variant is regressed on which
    framework positions carry the humanized residue (ridge, lam=1). When the
    variant set does not identify the positions (fewer variants than changed
    positions, e.g. a single V0), it falls back to attributing each variant's
    ddG equally to its changed positions. Per-experiment effects are then
    combined across experiments weighted by sample size and shrunk toward 0
    (empirical-Bayes-lite).
    """
    warnings: List[str] = []
    accum: Dict[str, Dict] = {}

    for exp in records:
        try:
            parent = number_heavy(exp.parent_vh) if exp.parent_vh else None
            pl = number_light(exp.parent_vl) if exp.parent_vl else None
        except ValueError as e:
            warnings.append(f"[{exp.name}] parent numbering failed: {e}")
            continue
        pd = {}
        if parent:
            pd.update(parent.posmap())
        if pl:
            pd.update(pl.posmap())

        pos_pairs: Dict[str, Tuple[str, str]] = {}
        rows: List[List[str]] = []
        ys: List[float] = []
        for v in exp.variants:
            vh = vl = None
            if v.get("vh"):
                try:
                    vh = number_heavy(v["vh"])
                except ValueError:
                    vh = None
            if v.get("vl"):
                try:
                    vl = number_light(v["vl"])
                except ValueError:
                    vl = None
            if vh is None and vl is None and v.get("seq"):
                try:
                    vh = number_heavy(v["seq"])
                except ValueError:
                    try:
                        vl = number_light(v["seq"])
                    except ValueError:
                        vl = None
            if vh is None and vl is None:
                warnings.append(f"[{exp.name}] variant {v['name']} not numberable")
                continue
            vd = {}
            if vh:
                vd.update(vh.posmap())
            if vl:
                vd.update(vl.posmap())
            ddg = _ddg(float(v["kd"]), exp.parent_kd)
            row = []
            for pos in sorted(set(pd) & set(vd)):
                num = int("".join(c for c in pos if c.isdigit()))
                if num >= (103 if pos[0] == "H" else 98):
                    continue
                if pd[pos] == vd[pos]:
                    continue
                # resolve the region against the chain that OWNS the position
                num_chain = vh if pos.startswith("H") else vl
                reg = num_chain.region_of(pos) or "" if num_chain else ""
                if not reg.startswith("FR"):
                    continue
                pos_pairs[pos] = (pd[pos], vd[pos])
                row.append(pos)
            rows.append(row)
            ys.append(ddg)

        pos_list = sorted(pos_pairs)
        if not pos_list:
            continue
        if len(rows) >= len(pos_list) and len(pos_list) > 1:
            beta = _ridge_solve(rows, ys, pos_list, lam=1.0)
            contribs = [(pos, b, len(rows)) for pos, b in zip(pos_list, beta)]
        else:
            sums: Dict[str, float] = {}
            counts: Dict[str, int] = {}
            for row, yy in zip(rows, ys):
                for pos in row:
                    sums[pos] = sums.get(pos, 0.0) + yy
                    counts[pos] = counts.get(pos, 0) + 1
            contribs = [(pos, sums[pos] / counts[pos], counts[pos])
                        for pos in sums if counts[pos] > 0]

        for pos, eff, n in contribs:
            a = accum.setdefault(pos, {"num": 0.0, "den": 0.0, "raw_n": 0,
                                       "best_n": 0, "donor": "", "human": ""})
            a["num"] += eff * n
            a["den"] += n
            a["raw_n"] += n
            if n > a["best_n"]:
                a["best_n"] = n
                a["donor"], a["human"] = pos_pairs[pos]

    effects: Dict[str, PositionEffect] = {}
    for pos, a in accum.items():
        if a["den"] <= 0:
            continue
        combined = a["num"] / a["den"]
        # shrinkage toward 0 (empirical-Bayes-lite): weight by total sample size
        shrunk = combined * a["den"] / (a["den"] + 1.0)
        effects[pos] = PositionEffect(
            position=pos,
            donor_aa=a["donor"],
            human_aa=a["human"],
            effect=round(shrunk, 3),
            n=1,
            raw_n=a["raw_n"],
        )
    return effects, warnings


# ---------------------------------------------------------------------------
# calibration artifact
# ---------------------------------------------------------------------------

def write_calibration(path: str, effects: Dict[str, PositionEffect],
                      meta: Optional[dict] = None) -> None:
    payload = {
        "schema": "humanize-calibration-v1",
        "meta": meta or {},
        "position_effects": {
            p: {
                "donor_aa": e.donor_aa,
                "human_aa": e.human_aa,
                "ddG_kcal": e.effect,
                "n_experiments": e.n,
                "n_variants": e.raw_n,
            }
            for p, e in sorted(effects.items())
        },
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)


def load_calibration(path: str) -> Dict[str, dict]:
    with open(path) as fh:
        data = json.load(fh)
    if data.get("schema") != "humanize-calibration-v1":
        raise ValueError(f"not a humanize calibration file: {path}")
    return data["position_effects"]


def effect_thresholds(ddG: float, n_variants: int) -> str:
    """Map an empirical effect to a scoring adjustment class.

    >0: donor residue empirically better (reverting helps affinity).
    <0: human residue empirically fine (reverting is unnecessary).
    """
    if n_variants == 0:
        return "none"
    if ddG >= 0.41:        # KD ~2x better with donor
        return "keep_donor"
    if ddG >= 0.20:
        return "promote"
    if ddG <= -0.20:
        return "demote"
    return "neutral"
