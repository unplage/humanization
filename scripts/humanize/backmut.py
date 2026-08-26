"""Back-mutation candidate identification and scoring.

For every framework position (FR1-FR3) where the donor differs from the
chosen human germline, we classify the position by structural role and
immunogenicity, then assign:

  * structural_score  - how much the donor residue supports CDR/fold
  * benefit_score     - how much humanizing the position buys (immunogenicity)
  * chemical_score    - developability-motif effects of the substitution
  * composite score   - blended 0-100
  * tier              - T1 must / T2 recommended / T3 optional /
                        KEEP_HUMAN / KEEP_DONOR
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .config import (
    CANONICAL,
    EMPIRICAL_NO_EFFECT,
    EMPIRICAL_NO_EFFECT_NOTE,
    INTERFACE_CORE,
    INTERFACE_EXTENDED,
    VERNIER_ZONE,
    VHH_HALLMARK,
    WEIGHTS,
)
from .learning import effect_thresholds
from .germline import GermlineDB, GermlineGene
from .numbering import NumberedChain

FR_REGIONS = ("FR1", "FR2", "FR3")
J_ANCHOR = {"H": 103, "L": 98}


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


@dataclass
class BackMutationCandidate:
    position: str              # e.g. "H67"
    donor_aa: str
    human_aa: str
    features: List[str]
    structural_score: float
    benefit_score: float
    chemical_score: float
    composite: float
    tier: str                  # T1/T2/T3/KEEP_HUMAN/KEEP_DONOR
    rationale: List[str] = field(default_factory=list)
    buried: Optional[bool] = None
    cdr_contact: Optional[bool] = None
    antigen_contact: Optional[bool] = None
    empirical_ddG: Optional[float] = None   # kcal/mol from experiments
    empirical_n: int = 0
    empirical_note: str = ""


@dataclass
class BackMutationResult:
    chain_type: str
    candidates: List[BackMutationCandidate]
    v_gene: GermlineGene
    germline_conservation: Dict[str, float] = field(default_factory=dict)

    def by_tier(self, tier: str) -> List[BackMutationCandidate]:
        return [c for c in self.candidates if c.tier == tier]

    def revert_positions(self, tiers=("T1", "T2")) -> List[str]:
        return [c.position for c in self.candidates if c.tier in tiers]


class StructureHints:
    """Structural annotations from AF3 (optional). None = unknown.

    data keys:
      buried          {pos: bool}
      cdr_contact     {pos: bool}      framework residue contacts any CDR
      antigen_contact {pos: bool}      framework/CDR residue contacts antigen
      cdr_partners    {fr_pos: [cdr_pos, ...]}  which CDR residues a
                       framework residue contacts (heavy atom < 4.5 A)
    """

    def __init__(self, data: Optional[dict] = None):
        self.data = data or {}

    def buried(self, chain: str, pos: str) -> Optional[bool]:
        d = self.data.get("buried")
        return d.get(pos) if d else None

    def cdr_contact(self, chain: str, pos: str) -> Optional[bool]:
        d = self.data.get("cdr_contact")
        return d.get(pos) if d else None

    def antigen_contact(self, chain: str, pos: str) -> Optional[bool]:
        d = self.data.get("antigen_contact")
        return d.get(pos) if d else None

    def cdr_partners(self, pos: str) -> set:
        d = self.data.get("cdr_partners") or {}
        v = d.get(pos) or []
        return set(v)

    def exposure(self, pos: str) -> float:
        b = self.buried("", pos)
        if b is None:
            return 0.5
        return 0.15 if b else 0.85


def _motif_hit(seq: str, pattern: str) -> bool:
    return bool(re.search(pattern, seq))


def analyze_backmutations(
    donor: NumberedChain,
    v_gene: GermlineGene,
    is_vhh: bool = False,
    structure: Optional[StructureHints] = None,
    top_germlines: Optional[List[Tuple[GermlineGene, dict]]] = None,
    calibration: Optional[Dict[str, dict]] = None,
) -> BackMutationResult:
    """Score all framework positions where donor != chosen germline."""
    chain_type = donor.chain_type
    if v_gene.numbered is None:
        raise ValueError(f"[{chain_type}] germline gene without numbering")
    dmap = donor.posmap()
    gmap = v_gene.numbered.posmap()
    structure = structure or StructureHints()

    # conservation of the donor residue across the top germlines
    conservation: Dict[str, float] = {}
    for p, aa in dmap.items():
        hits = 0
        n = 0
        for g, _s in (top_germlines or []):
            gm = g.numbered.posmap() if g.numbered else {}
            if p in gm and gm[p]:
                n += 1
                if gm[p] == aa:
                    hits += 1
        conservation[p] = hits / n if n else 0.0

    candidates: List[BackMutationCandidate] = []
    for pos in sorted(set(dmap) & set(gmap), key=lambda p: (_pos_num(p), p)):
        num = _pos_num(pos)
        if num >= J_ANCHOR[chain_type]:
            continue
        if donor.region_of(pos) not in FR_REGIONS:
            continue
        # H93/H94 carry the first two CDR3-loop residues (strict-Kabat FR3
        # labels, IMGT CDR3 105-106). They are grafted from the donor as part
        # of the loop and must never become back-mutation candidates.
        if chain_type == "H" and num in (93, 94):
            continue
        donor_aa = dmap[pos].upper()
        human_aa = gmap[pos].upper()
        if donor_aa == human_aa or donor_aa in ("", "X"):
            continue

        features: List[str] = []
        if num in INTERFACE_CORE[chain_type]:
            features.append("interface_core")
        if num in INTERFACE_EXTENDED[chain_type]:
            features.append("interface_extended")
        if num in VERNIER_ZONE[chain_type]:
            features.append("vernier")
        if num in CANONICAL[chain_type]:
            features.append("canonical")
        if is_vhh and chain_type == "H" and num in VHH_HALLMARK:
            features.append("vhh_hallmark")
        if donor_aa == "C":
            features.append("disulfide_cys")

        buried = structure.buried(chain_type, pos)
        cdr_contact = structure.cdr_contact(chain_type, pos)
        ag_contact = structure.antigen_contact(chain_type, pos)
        if buried:
            features.append("buried")
        if cdr_contact:
            features.append("cdr_contact")
        if ag_contact:
            features.append("antigen_contact")

        # ---- structural score ----
        w = WEIGHTS["structural"]
        struct_scores = [w[f] for f in features if f in w]
        # buried + any structural feature boosts confidence
        structural = max(struct_scores) if struct_scores else 0.0
        if features and buried is True:
            structural = max(structural, 0.7)
        if features and buried is False:
            structural = max(structural * 0.85, 0.0)

        # ---- immunogenicity benefit ----
        exposure = structure.exposure(pos)
        conservation_val = conservation.get(pos, 0.5)
        # rare donor residue among germlines -> more human-like to revert
        rare = 1.0 - conservation_val
        benefit = 0.3 + 0.5 * exposure * rare
        if not features:
            benefit = min(benefit, 0.45)   # exposed non-structural: low value

        # ---- chemical score (developability) ----
        # WARNING: 仅基于序列模式检测，未考虑结构暴露状态；
        # 埋藏位点实际风险较低，表面暴露位点风险更高。
        # 已排除保守 Cys（VH 22/92, VL 23/88）和双计风险。
        wc = WEIGHTS["chemical"]
        chem = 0.0
        # effect of reverting donor->human at this position: scan a window
        rpos = donor.residue(pos)
        if rpos is None:
            ridx = max(0, donor.sequence.find(donor_aa))
        else:
            ridx = rpos.index
        win_d = donor.sequence[max(0, ridx - 2): ridx + 3]

        # N-glycan (NxS/T) - 优先检测，排除重叠
        has_nglycan = _motif_hit(win_d, r"N[^P][ST]")
        chem += wc["removes_nglycan"] if has_nglycan else 0

        # deamidation: NG > NS > NH > ND (排除 N-glycan 重叠)
        # 如果是 N-glycan，跳过 NG/NS 避免双计
        if not has_nglycan:
            chem += wc["removes_deamidation_ng"] if _motif_hit(win_d, r"NG") else 0
            chem += wc["removes_deamidation_ns"] if _motif_hit(win_d, r"NS") else 0
        chem += wc["removes_deamidation_nh"] if _motif_hit(win_d, r"NH") else 0
        chem += wc["removes_deamidation_nd"] if _motif_hit(win_d, r"ND") else 0

        # isomerization: DG > DS = DT > DH
        chem += wc["removes_isomerization_dg"] if _motif_hit(win_d, r"DG") else 0
        chem += wc["removes_isomerization_ds"] if _motif_hit(win_d, r"DS") else 0
        chem += wc["removes_isomerization_dt"] if _motif_hit(win_d, r"DT") else 0
        chem += wc["removes_isomerization_dh"] if _motif_hit(win_d, r"DH") else 0

        # acid hydrolysis (D-X where X is small residue, 不含 G/S/T/H/D)
        # DD 单独计分，避免双计
        has_dd = _motif_hit(win_d, r"DD")
        if not has_dd:
            chem += wc["removes_acid_hydrolysis"] if _motif_hit(win_d, r"D[AVLIP]") else 0
        chem += wc["removes_acid_hydrolysis_dd"] if has_dd else 0

        # oxidation (M / W) - C 仅在非保守位点计分
        chem += wc["removes_oxidation"] if _motif_hit(win_d, r"[MW]") else 0
        if donor_aa == "C" and num not in {22, 92} if chain_type == "H" else {23, 88}:
            chem += wc["removes_oxidation"]  # 非保守 Cys

        # base hydrolysis (K-X where X is D/E)
        chem += wc["removes_base_hydrolysis"] if _motif_hit(win_d, r"K[DE]") else 0

        # metalloprotease cleavage (MK)
        chem += wc["removes_met_lyscleavage"] if _motif_hit(win_d, r"MK") else 0

        # ---- tier ----
        tier = _assign_tier(features, buried, is_vhh, donor_aa, pos)

        composite = round(100 * (
            WEIGHTS["blend"][0] * structural
            + WEIGHTS["blend"][1] * benefit
            + WEIGHTS["blend"][2] * max(0, min(1, chem))
        ), 1)

        # Gold-standard demotion: positions empirically shown to tolerate the
        # human residue (docs/backtest_report.md) are demoted to T3 UNLESS
        # structure data supports burial or CDR contact (AF3 mode overrides).
        demoted_note = ""
        if (tier in ("T1", "T2") and num in EMPIRICAL_NO_EFFECT[chain_type]
                and buried is not True and cdr_contact is not True):
            tier = "T3"
            composite = min(composite, 40.0)
            demoted_note = f"empirical: {EMPIRICAL_NO_EFFECT_NOTE.get(pos, '')}"

        # ---- empirical calibration (from experiment data) ----
        empirical_ddG = None
        empirical_n = 0
        empirical_note = ""
        if calibration:
            entry = calibration.get(pos)
            if entry:
                empirical_ddG = float(entry.get("ddG_kcal", 0.0))
                empirical_n = int(entry.get("n_variants", 0))
                adj = effect_thresholds(empirical_ddG, empirical_n)
                if adj == "keep_donor":
                    if tier in ("T1", "T2", "T3"):
                        tier = "T1"
                        empirical_note = (f"empirical: donor retains affinity "
                                          f"(ddG +{empirical_ddG:.2f}, n={empirical_n})")
                        composite = max(composite, 70.0)
                elif adj == "promote":
                    if tier == "T3":
                        tier = "T2"
                    empirical_note = (f"empirical: mild benefit of donor "
                                      f"(ddG +{empirical_ddG:.2f}, n={empirical_n})")
                elif adj == "demote":
                    if tier in ("T1", "T2"):
                        tier = "T3"
                    empirical_note = (f"empirical: human residue tolerated "
                                      f"(ddG {empirical_ddG:.2f}, n={empirical_n})")
                    composite = min(composite, 40.0)
                elif adj == "neutral" and empirical_n >= 2:
                    empirical_note = f"empirical: no effect (ddG {empirical_ddG:+.2f}, n={empirical_n})"

        rationale = _rationale(features, tier, buried, cdr_contact, ag_contact,
                               exposure, conservation_val, donor_aa, human_aa)
        if empirical_note:
            rationale.append(empirical_note)
        if demoted_note:
            rationale.append(demoted_note)

        candidates.append(BackMutationCandidate(
            position=pos,
            donor_aa=donor_aa,
            human_aa=human_aa,
            features=sorted(features),
            structural_score=round(structural, 2),
            benefit_score=round(benefit, 2),
            chemical_score=round(chem, 2),
            composite=composite,
            tier=tier,
            rationale=rationale,
            buried=buried,
            cdr_contact=cdr_contact,
            antigen_contact=ag_contact,
            empirical_ddG=empirical_ddG,
            empirical_n=empirical_n,
        ))

    return BackMutationResult(
        chain_type=chain_type,
        candidates=candidates,
        v_gene=v_gene,
        germline_conservation=conservation,
    )


def _assign_tier(features, buried, is_vhh, donor_aa, pos) -> str:
    if "vhh_hallmark" in features:
        return "KEEP_DONOR"
    if "disulfide_cys" in features:
        return "KEEP_DONOR"
    if "interface_core" in features:
        if buried is False:
            return "T2"
        return "T1"
    if "vernier" in features:
        if buried is True or "canonical" in features:
            return "T1"
        return "T2"
    if "canonical" in features:
        return "T2"
    if "cdr_contact" in features or "antigen_contact" in features:
        return "T2"
    if "buried" in features:
        return "T2"
    if "interface_extended" in features:
        return "T2"
    return "T3"


def _rationale(features, tier, buried, cdr_contact, ag_contact,
               exposure, conservation, donor_aa, human_aa) -> List[str]:
    out = []
    if tier == "KEEP_DONOR":
        if "vhh_hallmark" in features:
            out.append("VHH hallmark residue (FR2 hydrophilic patch); must stay donor for single-domain fold")
        if "disulfide_cys" in features:
            out.append("Framework Cys (potential CDR3 disulfide partner); must stay donor")
    if "interface_core" in features:
        out.append("Core VH/VL packing interface (Chothia)")
    if "vernier" in features:
        out.append("Vernier zone (Foote & Winter): supports CDR conformation")
    if "canonical" in features:
        out.append("Canonical-structure framework residue (Chothia & Lesk)")
    if "interface_extended" in features:
        out.append("Extended VH/VL interface position")
    if buried is True:
        out.append("Buried in domain core (structure)")
    if cdr_contact:
        out.append("Contacts CDR loop (<4.5 A, structure)")
    if ag_contact:
        out.append("Contacts antigen (<4.5 A, structure)")
    if exposure >= 0.7:
        out.append("Surface-exposed: direct immunogenicity risk")
    if conservation <= 0.2:
        out.append(f"Donor {donor_aa} rare in human germlines (conservation {conservation:.0%})")
    if not features and exposure < 0.5:
        out.append("No structural role identified; buried")
    return out
