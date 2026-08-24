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


def compute_position_effects(
    records: List[ExperimentRecord],
    db: GermlineDB,
) -> Tuple[Dict[str, PositionEffect], List[str]]:
    """Estimate per-position effects. Returns (effects, warnings)."""
    warnings: List[str] = []
    # donor/human aa at each FR position, per experiment
    state: Dict[str, Dict[str, dict]] = {}   # exp -> pos -> {donor: [ddg], human: [ddg]}
    pos_info: Dict[str, Dict] = {}           # exp -> pos -> (donor_aa, human_aa)

    for exp in records:
        try:
            parent = number_heavy(exp.parent_vh) if exp.parent_vh else None
            pl = number_light(exp.parent_vl) if exp.parent_vl else None
        except ValueError as e:
            warnings.append(f"[{exp.name}] parent numbering failed: {e}")
            continue
        pd = (parent.posmap() if parent else {}) | (pl.posmap() if pl else {})
        pidx = {}
        if parent:
            pidx.update({r.pos: r.index for r in parent.residues})
        if pl:
            pidx.update({r.pos: r.index for r in pl.residues})

        exp_state: Dict[str, Dict[str, list]] = {}
        exp_info: Dict[str, Tuple[str, str]] = {}
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
            vd = (vh.posmap() if vh else {}) | (vl.posmap() if vl else {})
            ddg = _ddg(float(v["kd"]), exp.parent_kd)
            for pos in set(pd) & set(vd):
                num = int("".join(c for c in pos if c.isdigit()))
                if num >= (103 if pos[0] == "H" else 98):
                    continue
                if pd[pos] == vd[pos]:
                    continue
                # only framework positions (structure of the CDR is constant)
                num_chain = vh if vh is not None else vl
                reg = num_chain.region_of(pos) or "" if num_chain else ""
                if not reg.startswith("FR"):
                    continue
                donor_aa, human_aa = pd[pos], vd[pos]
                if donor_aa == human_aa:
                    continue
                exp_state.setdefault(pos, {"donor": [], "human": []})
                exp_info[pos] = (donor_aa, human_aa)
                # is the variant carrying the donor residue at pos?
                key = "donor" if donor_aa == vd[pos] else "human"
                # note: for a humanized variant, vd[pos] is usually the human
                # residue; a back-mutated variant carries the donor residue.
                # The contrast below pools both directions within the parent.
                exp_state[pos][key].append(ddg)
        state[exp.name] = exp_state
        pos_info[exp.name] = exp_info

    effects: Dict[str, PositionEffect] = {}
    for exp_name, exp_state in state.items():
        for pos, sides in exp_state.items():
            donor_ddg = sides.get("donor", [])
            human_ddg = sides.get("human", [])
            n = len(donor_ddg) + len(human_ddg)
            if n == 0:
                continue
            mean_donor = sum(donor_ddg) / len(donor_ddg) if donor_ddg else 0.0
            mean_human = sum(human_ddg) / len(human_ddg) if human_ddg else 0.0
            # effect > 0  <=>  variants carrying the HUMAN residue have worse
            # affinity  <=>  reverting to the donor residue helps
            effect = mean_human - mean_donor
            # shrinkage toward 0 (empirical-Bayes-lite): weight by n
            shrunk = effect * n / (n + 1.0)
            donor_aa, human_aa = pos_info[exp_name][pos]
            if pos not in effects or effects[pos].raw_n < n:
                effects[pos] = PositionEffect(
                    position=pos,
                    donor_aa=donor_aa,
                    human_aa=human_aa,
                    effect=round(shrunk, 3),
                    n=1,
                    raw_n=n,
                )
            else:
                e = effects[pos]
                e.effect = round((e.effect * e.raw_n + shrunk) / (e.raw_n + n), 3)
                e.raw_n += n
                e.n += 1
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
