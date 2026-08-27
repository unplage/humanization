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
      buried          {pos: Optional[bool]}  None = uncertain (relSASA 0.15-0.25)
      cdr_contact     {pos: bool}      framework residue contacts any CDR
      antigen_contact {pos: bool}      framework/CDR residue contacts antigen
      cdr_partners    {fr_pos: [cdr_pos, ...]}  which CDR residues a
                       framework residue contacts (heavy atom < 4.5 A)
      plddt           {pos: float}     AF3 confidence score (B-factor)
      rel_sasa        {pos: float}     relative SASA (0-1)
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

    def plddt(self, pos: str) -> Optional[float]:
        """AF3 confidence score for this position (0-100)."""
        d = self.data.get("plddt")
        return d.get(pos) if d else None

    def rel_sasa(self, pos: str) -> Optional[float]:
        """Relative SASA for this position (0-1)."""
        d = self.data.get("rel_sasa")
        return d.get(pos) if d else None

    def exposure(self, pos: str) -> float:
        b = self.buried("", pos)
        if b is None:
            return 0.5
        return 0.15 if b else 0.85


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

        # 方案1: pLDDT 加权 - 低置信度区域降低结构证据权重
        plddt_val = structure.plddt(pos)
        if plddt_val is not None and plddt_val < 50:
            # pLDDT < 50: 结构证据不可靠，大幅降权
            structural *= 0.3
        elif plddt_val is not None and plddt_val < 70:
            # pLDDT 50-70: 结构证据中等置信度，适度降权
            structural *= 0.7

        # ---- immunogenicity benefit ----
        exposure = structure.exposure(pos)
        conservation_val = conservation.get(pos, 0.5)
        # rare donor residue among germlines -> more human-like to revert
        rare = 1.0 - conservation_val
        benefit = 0.3 + 0.5 * exposure * rare
        # 无结构特征位点: 人源化收益低，但仍保留 conservation 相对梯度
        # （仅设上限，不抹平排序信息）
        if not features:
            benefit = min(benefit, 0.50)

        # ---- chemical score (developability) ----
        # WARNING: 仅基于序列模式检测，未考虑结构暴露状态；
        # 埋藏位点实际风险较低，表面暴露位点风险更高。
        # 已排除保守 Cys（VH 22/92, VL 23/88）和双计风险。
        # 每个 motif 均锚定当前回复位点（避免窗口内无关 motif 误报）。
        wc = WEIGHTS["chemical"]
        chem = 0.0
        rpos = donor.residue(pos)
        ridx = rpos.index if rpos is not None else max(0, donor.sequence.find(donor_aa))
        seq = donor.sequence

        def _aa(off: int) -> str:
            j = ridx + off
            return seq[j] if 0 <= j < len(seq) else ""

        a0, a1, a2 = _aa(0), _aa(1), _aa(2)

        # N-glycan (NxS/T): 当前位点是 N，N+1 非 P，N+2 是 S/T
        ngly = a0 == "N" and a1 not in ("P", "") and a2 in ("S", "T")
        if ngly:
            chem += wc["removes_nglycan"]

        # deamidation: 当前位点是 N（N-glycan 已计分时跳过 NS/NG 防双计）
        if a0 == "N" and not ngly:
            chem += wc["removes_deamidation_ng"] if a1 == "G" else 0
            chem += wc["removes_deamidation_ns"] if a1 == "S" else 0
            chem += wc["removes_deamidation_nh"] if a1 == "H" else 0
            chem += wc["removes_deamidation_nd"] if a1 == "D" else 0

        # isomerization / acid hydrolysis: 当前位点是 D
        if a0 == "D":
            chem += wc["removes_isomerization_dg"] if a1 == "G" else 0
            chem += wc["removes_isomerization_ds"] if a1 == "S" else 0
            chem += wc["removes_isomerization_dt"] if a1 == "T" else 0
            chem += wc["removes_isomerization_dh"] if a1 == "H" else 0
            # acid hydrolysis: D-X (X=A/V/L/I/P，不含 G/S/T/H/D 避免重叠)；DD 单独计分
            if a1 == "D":
                chem += wc["removes_acid_hydrolysis_dd"]
            elif a1 in ("A", "V", "L", "I", "P"):
                chem += wc["removes_acid_hydrolysis"]

        # oxidation: 当前位点是 M/W，或非保守 C
        if a0 in ("M", "W"):
            chem += wc["removes_oxidation"]
        # 保守二硫键 Cys：VH 22/92；kappa VL 23/88；lambda VL FR1 短一个
        # 残基，第一个保守 Cys 落在 22（与 developability.py 保持一致）。
        conserved_cys = ({22, 92} if chain_type == "H" else {22, 23, 88})
        if a0 == "C" and num not in conserved_cys:
            chem += wc["removes_oxidation"]  # 非保守 Cys（排除保守二硫键 Cys）

        # base hydrolysis: 当前位点是 K，K+1 是 D/E
        if a0 == "K" and a1 in ("D", "E"):
            chem += wc["removes_base_hydrolysis"]

        # metalloprotease cleavage: 当前位点是 M，M+1 是 K
        if a0 == "M" and a1 == "K":
            chem += wc["removes_met_lyscleavage"]

        # ---- introduced risk penalty (回复可能引入新风险) ----
        # 回复将当前位点改为 human_aa 后，检查是否新形成 N-glycan (N-X-S/T)。
        # 覆盖三种锚定方式（当前位点是 N / N+1 / N+2）。
        # 注：回复要求 donor != human，故 donor 的 N-glycan 已被移除（见上方
        # removes_nglycan），此处只惩罚"新引入"，不与移除奖励冲突。
        h0 = human_aa
        prev1 = seq[ridx - 1] if ridx - 1 >= 0 else ""
        prev2 = seq[ridx - 2] if ridx - 2 >= 0 else ""
        nxt1 = seq[ridx + 1] if ridx + 1 < len(seq) else ""
        # 情况1: 当前位点是 N，N+1 非 P，N+2 是 S/T
        case1 = h0 == "N" and a1 not in ("P", "") and a2 in ("S", "T")
        # 情况2: 当前位点是 N-glycan 的 X（N-X-S/T），X 由 human_aa 提供且非 P
        case2 = prev1 == "N" and h0 not in ("P", "") and nxt1 in ("S", "T")
        # 情况3: 当前位点是 N-glycan 的 S/T（N-X-S/T），S/T 由 human_aa 提供
        case3 = prev2 == "N" and prev1 not in ("P", "") and h0 in ("S", "T")
        if case1 or case2 or case3:
            chem += wc["introduces_nglycan"]  # -0.8 惩罚

        # ---- tier ----
        tier = _assign_tier(features, buried, is_vhh, donor_aa, pos)

        composite = round(100 * (
            WEIGHTS["blend"][0] * structural
            + WEIGHTS["blend"][1] * benefit
            # positive chemical rewards are capped at 1, but negative
            # penalties (introduces_nglycan) must pass through unclamped,
            # otherwise the penalty would silently vanish (regression-tested)
            + WEIGHTS["blend"][2] * min(1, chem)
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
    return out
