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
    FR4_STRUCTURAL_WEIGHTS,
    FR4_REVERSION_THRESHOLD,
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

# Backbone atom names (used to separate side-chain-mediated contacts from
# trivial main-chain contacts in the Ig beta-sandwich).
BACKBONE_ATOMS = {"N", "CA", "C", "O"}

# Residue side-chain volume (Zamyatnin 1972, A^3). Used to flag buried
# substitutions whose volume change is large enough to strain core packing.
RESIDUE_VOLUME = {
    "G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "T": 116.1, "P": 122.7,
    "V": 140.0, "D": 111.1, "E": 138.4, "N": 114.1, "Q": 143.8, "M": 162.9,
    "L": 166.7, "I": 166.7, "K": 168.6, "R": 173.4, "H": 153.2, "F": 189.9,
    "Y": 193.6, "W": 227.8,
}
# A buried donor<->human substitution changing side-chain volume by this much
# (A^3) is treated as a genuine packing perturbation and keeps its tier.
VOLUME_STRAIN_THRESHOLD = 40.0


def _volume_strain(donor_aa: str, human_aa: str) -> Optional[float]:
    """Absolute side-chain volume change for a donor->human substitution."""
    vd = RESIDUE_VOLUME.get(donor_aa)
    vh = RESIDUE_VOLUME.get(human_aa)
    if vd is None or vh is None:
        return None
    return abs(vh - vd)


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def _noisy_or(weights: List[float]) -> float:
    """Probabilistic OR of independent evidence weights (0-1)."""
    prod = 1.0
    for w in weights:
        prod *= (1.0 - max(0.0, min(1.0, w)))
    return 1.0 - prod


def _plddt_factor(plddt: Optional[float]) -> float:
    """Continuous confidence factor for structural evidence.

    None (no pLDDT available, e.g. experimental PDB) -> 1.0 (no change).
    Otherwise scales 0.2 (pLDDT <= 50) to 1.0 (pLDDT >= 90).
    """
    if plddt is None:
        return 1.0
    frac = (float(plddt) - 50.0) / 40.0
    frac = max(0.0, min(1.0, frac))
    return 0.2 + 0.8 * frac


_LIABILITY_PATTERNS: Optional[List[Tuple[re.Pattern, float]]] = None


def _liability_score(sequence: str) -> float:
    """Sum of developability liability weights present in ``sequence``."""
    global _LIABILITY_PATTERNS
    if _LIABILITY_PATTERNS is None:
        from .config import LIABILITY_MOTIFS
        _LIABILITY_PATTERNS = [
            (re.compile(pat), weight) for pat, weight in LIABILITY_MOTIFS.values()
        ]
    return sum(len(rx.findall(sequence)) * w for rx, w in _LIABILITY_PATTERNS)


def _calibration_matches(entry: dict, donor_aa: str, human_aa: str) -> bool:
    """True when a calibration entry applies to this donor/human pair.

    Older calibration files may omit the residue fields; treat those as
    position-only and accept them for backward compatibility.
    """
    ed = (entry.get("donor_aa") or "").upper()
    eh = (entry.get("human_aa") or "").upper()
    if not ed and not eh:
        return True
    return ed == donor_aa.upper() and eh == human_aa.upper()


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
    side_chain_contact: Optional[bool] = None  # P0: side-chain mediated CDR contact
    antigen_contact: Optional[bool] = None
    empirical_ddG: Optional[float] = None   # kcal/mol from experiments
    empirical_n: int = 0
    empirical_note: str = ""
    immunogenicity_score: Optional[float] = None  # optional MHC-II epitope proxy
    dev_risk_score: Optional[float] = None  # P0: developability risk score (0-1)


@dataclass
class BackMutationResult:
    chain_type: str
    candidates: List[BackMutationCandidate]
    v_gene: GermlineGene
    germline_conservation: Dict[str, float] = field(default_factory=dict)
    fr_indels: List = field(default_factory=list)  # List[FRIndel]

    def by_tier(self, tier: str) -> List[BackMutationCandidate]:
        return [c for c in self.candidates if c.tier == tier]

    def revert_positions(self, tiers=("T1", "T2")) -> List[str]:
        return [c.position for c in self.candidates if c.tier in tiers]

    @property
    def indel_insertion_positions(self) -> List[str]:
        """Positions of donor insertions (for V2 default inclusion)."""
        return [ind.position for ind in self.fr_indels
                if ind.indel_type == "insertion"]

    @property
    def indel_deletion_positions(self) -> List[str]:
        """Positions of donor deletions (germline residues to add)."""
        return [ind.position for ind in self.fr_indels
                if ind.indel_type == "deletion"]


class StructureHints:
    """Structural annotations from AF3 (optional). None = unknown.

    data keys:
      buried          {pos: Optional[bool]}  None = uncertain (relSASA 0.15-0.25)
      cdr_contact     {pos: bool}      framework residue contacts any CDR
      cdr_contact_sc  {pos: bool}      framework SIDE CHAIN contacts any CDR
                       (backbone-mediated contacts excluded: they are fixed
                       beta-sheet geometry, unchanged by side-chain swaps)
      antigen_contact {pos: bool}      framework/CDR residue contacts antigen
      cdr_partners    {fr_pos: [cdr_pos, ...]}  which CDR residues a
                       framework residue contacts (heavy atom < 4.5 A)
      cdr_partners_sc {fr_pos: [cdr_pos, ...]}  side-chain-mediated subset
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

    def cdr_contact_sc(self, chain: str, pos: str) -> Optional[bool]:
        """Side-chain-mediated CDR contact. None when the structure hints
        predate side-chain attribution (caller should fall back)."""
        d = self.data.get("cdr_contact_sc")
        return d.get(pos) if d else None

    def antigen_contact(self, chain: str, pos: str) -> Optional[bool]:
        d = self.data.get("antigen_contact")
        return d.get(pos) if d else None

    def cdr_partners(self, pos: str) -> set:
        d = self.data.get("cdr_partners") or {}
        v = d.get(pos) or []
        return set(v)

    def cdr_partners_sc(self, pos: str) -> set:
        d = self.data.get("cdr_partners_sc") or {}
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


def _evaluate_fr4_structural(
    pos: str,
    num: int,
    chain_type: str,
    dmap: dict,
    gmap: dict,
    structure: StructureHints,
) -> Optional[BackMutationCandidate]:
    """Evaluate if an FR4 position needs reversion to donor based on structure.

    FR4 comes from the human J gene by construction. However, if structure data
    reveals that an FR4 position contacts CDR3 or the antigen, reverting to the
    donor residue may be necessary to preserve binding.

    Returns BackMutationCandidate if reversion is recommended, None otherwise.
    """
    donor_aa = dmap.get(pos, "").upper()
    human_aa = gmap.get(pos, "").upper()

    # Skip if residues are identical
    if donor_aa == human_aa or donor_aa in ("", "X"):
        return None

    # Collect structural features for FR4
    features = []
    scores = []

    # 1. CDR3 contact (most important - FR4 follows CDR3)
    cdr3_contact = structure.cdr_contact(chain_type, pos)
    if cdr3_contact:
        features.append("cdr3_contact")
        scores.append(FR4_STRUCTURAL_WEIGHTS["cdr3_contact"])

    # 2. Antigen contact
    ag_contact = structure.antigen_contact(chain_type, pos)
    if ag_contact:
        features.append("antigen_contact")
        scores.append(FR4_STRUCTURAL_WEIGHTS["antigen_contact"])

    # 3. Buried status
    buried = structure.buried(chain_type, pos)
    if buried:
        features.append("buried")
        scores.append(FR4_STRUCTURAL_WEIGHTS["buried"])

    # 4. VH/VL interface
    if num in INTERFACE_CORE.get(chain_type, set()):
        features.append("interface_core")
        scores.append(FR4_STRUCTURAL_WEIGHTS["interface_core"])

    # No structural evidence - skip
    if not scores:
        return None

    # Calculate composite score
    composite = sum(scores)

    # Check threshold
    if composite < FR4_REVERSION_THRESHOLD:
        return None

    # Calculate individual scores for BackMutationCandidate
    structural = max(scores)
    benefit = 0.0  # FR4 reversion is structure-driven, not immunogenicity
    chem = 0.0

    # Build rationale
    rationale_parts = [f"FR4 structural reversion: {', '.join(features)}"]
    if cdr3_contact:
        rationale_parts.append("Contacts CDR3 loop (<4.5 Å)")
    if ag_contact:
        rationale_parts.append("Contacts antigen (<4.5 Å)")
    if buried:
        rationale_parts.append("Buried position")

    return BackMutationCandidate(
        position=pos,
        donor_aa=donor_aa,
        human_aa=human_aa,
        features=sorted(features),
        structural_score=round(structural, 2),
        benefit_score=round(benefit, 2),
        chemical_score=round(chem, 2),
        composite=round(100 * composite, 1),
        tier="T_FR4",
        rationale=rationale_parts,
        buried=buried,
        cdr_contact=cdr3_contact,
        antigen_contact=ag_contact,
    )


def analyze_backmutations(
    donor: NumberedChain,
    v_gene: GermlineGene,
    is_vhh: bool = False,
    structure: Optional[StructureHints] = None,
    top_germlines: Optional[List[Tuple[GermlineGene, dict]]] = None,
    calibration: Optional[Dict[str, dict]] = None,
    indel_overrides: Optional[Dict[str, str]] = None,
    j_gene: Optional[GermlineGene] = None,
    immunogenicity: Optional[Dict[str, float]] = None,
) -> BackMutationResult:
    """Score all framework positions where donor != chosen germline.

    Positions are walked in germline numbering and resolved to the aligned
    donor residue through the FR correspondence, so a donor framework insertion
    does not shift every downstream position.

    indel_overrides: dict mapping FR region to user-selected insertion position
                     (e.g. {"FR1": "H6A"})
    j_gene:          human J gene, enabling FR4 structural reversion analysis.
    """
    chain_type = donor.chain_type
    if v_gene.numbered is None:
        raise ValueError(f"[{chain_type}] germline gene without numbering")
    dmap = donor.posmap()
    gmap = v_gene.numbered.posmap()
    structure = structure or StructureHints()

    # FR indels (respecting user overrides) + donor<->germline correspondence
    from .fr_indel import (
        build_fr_correspondence,
        detect_fr_indels,
        update_indel_selection,
    )
    fr_indels = detect_fr_indels(donor, v_gene)
    if indel_overrides:
        for i, indel in enumerate(fr_indels):
            if indel.fr_region in indel_overrides:
                fr_indels[i] = update_indel_selection(
                    indel, indel_overrides[indel.fr_region])
    corr = build_fr_correspondence(donor, v_gene, fr_indels)
    insertion_set = set(corr.insertion_positions)

    conservation: Dict[str, float] = {}

    # Precompute the frequency-weighted residue distribution of the reference
    # panel per position (single pass over the repertoire), so each candidate
    # lookup is O(1). The panel is the full human repertoire (see
    # pipeline._top_homologous_germlines), so conservation is measured against
    # the whole repertoire rather than 20 close homologs.
    from .germline_frequency import get_frequency
    ref_freq: Dict[str, Dict[str, float]] = {}
    ref_total: Dict[str, float] = {}
    for g, _s in (top_germlines or []):
        gm = g.numbered.posmap() if g.numbered else {}
        if not gm:
            continue
        w = get_frequency(chain_type, g.gene_id) or 0.0
        for p, aa in gm.items():
            ref_freq.setdefault(p, {})
            ref_freq[p][aa] = ref_freq[p].get(aa, 0.0) + w
            ref_total[p] = ref_total.get(p, 0.0) + w

    def _conservation(dpos: str, gpos: str, donor_aa: str) -> float:
        total = ref_total.get(gpos, 0.0)
        if total <= 0:
            return 0.0
        return ref_freq.get(gpos, {}).get(donor_aa, 0.0) / total

    candidates: List[BackMutationCandidate] = []
    for gpos in sorted(gmap, key=lambda p: (_pos_num(p), p)):
        dpos = corr.germline_to_donor.get(gpos)
        if dpos is None or dpos in insertion_set:
            # donor deletion (no donor residue) or donor insertion
            # (scored as an indel candidate, not a substitution)
            continue
        # Structural feature membership and the FR/J boundaries are properties
        # of the DONOR Kabat numbering (the donor is the antibody being
        # humanized); the germline label can differ around insertions.
        num = _pos_num(dpos)
        if num >= J_ANCHOR[chain_type]:
            continue  # FR4 is handled separately (human J region)
        donor_res = donor.residue(dpos)
        if donor_res is None or donor_res.region not in FR_REGIONS:
            continue
        # H93/H94 carry the first two CDR3-loop residues (strict-Kabat FR3
        # labels, IMGT CDR3 105-106). They are grafted from the donor as part
        # of the loop and must never become back-mutation candidates.
        if chain_type == "H" and num in (93, 94):
            continue
        donor_aa = dmap.get(dpos, "").upper()
        human_aa = gmap.get(gpos, "").upper()
        if donor_aa == human_aa or donor_aa in ("", "X"):
            continue
        pos = dpos  # candidate label = donor position (graft resolves it)
        conservation_val = _conservation(dpos, gpos, donor_aa)
        conservation[pos] = conservation_val

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
        # P0: side-chain mediated CDR contact (higher priority than backbone)
        sc_contact = structure.cdr_contact_sc(chain_type, pos)
        if sc_contact is None:
            sc_contact = cdr_contact  # fallback if structure hints predate sc attribution
        if buried:
            features.append("buried")
        if cdr_contact:
            features.append("cdr_contact")
        if ag_contact:
            features.append("antigen_contact")

        # ---- structural score ----
        # Noisy-OR combination: independent pieces of structural evidence
        # accumulate (multiple weak features reinforce each other) instead of
        # being collapsed to the single strongest feature by max().
        w = WEIGHTS["structural"]
        structural = _noisy_or([w[f] for f in features if f in w])
        if features and buried is True:
            structural = max(structural, 0.7)
        if features and buried is False:
            structural = structural * 0.85

        # pLDDT: continuous confidence weighting (not a hard threshold). The
        # structure evidence is scaled smoothly from 0.2 at pLDDT <= 50 to 1.0
        # at pLDDT >= 90; unknown pLDDT leaves the evidence untouched.
        plddt_val = structure.plddt(pos)
        structural *= _plddt_factor(plddt_val)

        # ---- immunogenicity benefit ----
        # Default: surface-exposure x germline-rarity proxy. When a per-position
        # MHC-II epitope score is supplied (NetMHCIIpan or a precomputed map),
        # it dominates: a donor residue inside a strong predicted T-cell
        # epitope gives a high benefit for reverting to the human germline.
        exposure = structure.exposure(pos)
        rare = 1.0 - conservation_val
        epitope = immunogenicity.get(pos) if immunogenicity else None
        if epitope is not None:
            benefit = 0.3 + 0.5 * (0.4 * exposure * rare + 0.6 * float(epitope))
        else:
            benefit = 0.3 + 0.5 * exposure * rare
        benefit = max(0.3, min(0.8, benefit))
        # 无结构特征位点: 人源化收益低，但仍保留 conservation 相对梯度
        # （仅设上限，不抹平排序信息；表位证据存在时不压低）
        if not features and epitope is None:
            benefit = min(benefit, 0.50)

        # ---- chemical score (developability) ----
        # Symmetric liability delta at this position:
        #   chem > 0 -> reverting to the donor REMOVES a liability that the
        #               human (graft) state carries  -> reward reversion
        #   chem < 0 -> reverting INTRODUCES / retains a donor liability
        #               -> penalise reversion
        # This is a pure sequence-pattern scan (exposure-agnostic): buried
        # liabilities are less risky in practice, but the term stays
        # conservative. Conserved disulfide Cys are not in the motif table.
        d_res = donor.residue(pos)
        ridx = d_res.index if d_res is not None else max(0, donor.sequence.find(donor_aa))
        seq = donor.sequence
        human_state = seq[:ridx] + human_aa + seq[ridx + 1:]
        chem = _liability_score(human_state) - _liability_score(seq)

        # ---- tier ----
        # P0: Use structure-aware tier assignment
        strain = _volume_strain(donor_aa, human_aa)
        tier = _assign_tier_with_structure(
            features, buried, is_vhh, donor_aa, pos,
            side_chain_contact=sc_contact,
            cdr_contact=cdr_contact,
            antigen_contact=ag_contact,
            volume_strain=strain,
            dev_risk_score=None  # Will be calculated later
        )

        # ---- Step 3 structure-driven tier adjustment ----
        # Two distinct refinements, both gated on reliable structure:
        #
        # (A) "Must-revert" pillars (T1) require *functional* evidence: a
        #     SIDE-CHAIN contact to a CDR/antigen, or a substitution large
        #     enough to strain buried packing. A buried Vernier/canonical/
        #     interface-core residue whose only CDR proximity is main-chain
        #     mediated, and whose donor->human swap is near-isosteric, is a
        #     literature heuristic rather than a structural necessity: it is
        #     downgraded T1 -> T2 (still recommended, no longer mandatory).
        #
        # (B) An exposed, non-contacting framework residue has no structural
        #     reason to stay donor and carries real immunogenicity risk: it
        #     goes to T3.
        #
        # Backbone-only CDR "contacts" are pervasive in the Ig beta-sandwich
        # and are unchanged by a side-chain swap, so they never block (A).
        # Interface-core positions are still protected from (B) because the
        # contact detector cannot see the VH/VL partner chain.
        structural_demoted = False
        soft_demoted = False
        
        # P0: Calculate developability risk score
        # Based on sequence motifs that affect developability
        dev_risk_score = _calculate_dev_risk_score(donor_aa, pos, donor, chain_type)
        
        if structure.data and tier in ("T1", "T2"):
            sc_partners = structure.cdr_partners_sc(pos)
            functional_contact = (
                sc_contact is True or ag_contact is True
                or bool(sc_partners)
            )
            strained = strain is not None and strain >= VOLUME_STRAIN_THRESHOLD
            plddt_ok = (plddt_val is None or plddt_val >= 50)
            if plddt_ok and not functional_contact:
                if buried is False and "interface_core" not in features:
                    # exposed, no functional contact: no reason to stay donor
                    # (immunogenicity dominates; volume strain irrelevant)
                    tier = "T3"
                    structural_demoted = True
                elif tier == "T1" and not strained:
                    # buried core/interface pillar whose side chain does not
                    # touch the paratope and whose human swap is near-isosteric:
                    # recommended, but not structurally mandatory
                    tier = "T2"
                    soft_demoted = True
            
            # P0: Developability risk filtering for exposed positions
            # Exposed positions with high developability risk should be demoted
            if buried is False and dev_risk_score > 0.7 and tier in ("T1", "T2"):
                tier = "T3"
                structural_demoted = True
                rationale_item = f"dev_risk: high developability risk ({dev_risk_score:.2f})"
                # Will be added to rationale later

        composite = round(100 * (
            WEIGHTS["blend"][0] * structural
            + WEIGHTS["blend"][1] * benefit
            # positive chemical rewards are capped at 1, but negative
            # penalties (introduces_nglycan) must pass through unclamped,
            # otherwise the penalty would silently vanish (regression-tested)
            + WEIGHTS["blend"][2] * min(1, chem)
        ), 1)

        if structural_demoted:
            composite = min(composite, 40.0)
        elif soft_demoted:
            # keep it a strong T2 but below genuine functional T1s
            composite = min(composite, 65.0)

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
            # A calibration effect is only meaningful for the exact donor/human
            # substitution pair it was measured on; a different germline at the
            # same position is a different substitution.
            if entry and not _calibration_matches(entry, donor_aa, human_aa):
                entry = None
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

        rationale = _rationale_with_structure(features, tier, buried, cdr_contact, ag_contact,
                               exposure, conservation_val, donor_aa, human_aa,
                               side_chain_contact=sc_contact, dev_risk_score=dev_risk_score)
        if structural_demoted:
            rationale.append("structural: exposed and no CDR/antigen contact; "
                            "demoted from T1/T2 to T3")
        if soft_demoted:
            rationale.append(
                "structural: buried with no side-chain CDR/antigen contact "
                "and near-isosteric substitution; literature pillar downgraded "
                "from T1 (must revert) to T2 (recommended)")
        if empirical_note:
            rationale.append(empirical_note)
        if demoted_note:
            rationale.append(demoted_note)
        if epitope is not None and epitope >= 0.5:
            rationale.append(
                f"immunogenicity: donor residue in predicted MHC-II epitope "
                f"(score {epitope:.2f}); reverting lowers ADA risk")

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
            side_chain_contact=sc_contact,  # P0: side-chain mediated CDR contact
            antigen_contact=ag_contact,
            empirical_ddG=empirical_ddG,
            empirical_n=empirical_n,
            immunogenicity_score=epitope,
            dev_risk_score=dev_risk_score,  # P0: developability risk score
        ))

    # ---- FR4 structural reversion (only with J gene + structure data) ----
    # FR4 is human by construction, but a J-region residue that contacts CDR3
    # or the antigen may need to stay donor. This is the one sanctioned
    # exception to the "never back-mutate FR4" rule.
    if j_gene is not None and j_gene.numbered is not None and structure.data:
        jmap = j_gene.numbered.posmap()
        for jpos in sorted(jmap, key=lambda p: (_pos_num(p), p)):
            if _pos_num(jpos) < J_ANCHOR[chain_type]:
                continue
            donor_aa = dmap.get(jpos, "").upper()
            human_aa = jmap.get(jpos, "").upper()
            if donor_aa == human_aa or donor_aa in ("", "X"):
                continue
            fr4_candidate = _evaluate_fr4_structural(
                jpos, _pos_num(jpos), chain_type, dmap, jmap, structure)
            if fr4_candidate:
                candidates.append(fr4_candidate)

    # ---- FR indel candidates (donor insertions) ----
    for indel in fr_indels:
        if indel.indel_type == "insertion":
            # Donor insertion: score whether to keep donor (indel) or revert
            pos = indel.position
            donor_aa = indel.donor_aa
            # No germline counterpart — treat as "human_aa = missing"
            features = ["fr_indel", "donor_specific"]

            buried = structure.buried(chain_type, pos) if structure else None
            cdr_contact = structure.cdr_contact(chain_type, pos) if structure else None
            ag_contact = structure.antigen_contact(chain_type, pos) if structure else None

            if buried:
                features.append("buried")
            if cdr_contact:
                features.append("cdr_contact")

            # Structural score: noisy-OR of the indel evidence, pLDDT-weighted
            w = WEIGHTS["structural"]
            structural = _noisy_or([w[f] for f in features if f in w]) or 0.3
            if buried is True:
                structural = max(structural, 0.7)
            elif buried is False:
                structural = structural * 0.85
            structural *= _plddt_factor(structure.plddt(pos))

            # Benefit score: indel is donor-specific, high immunogenicity risk
            exposure = 0.0 if buried is True else (1.0 if buried is False else 0.5)
            benefit = 0.3 + 0.5 * exposure * 0.8  # moderate rarity

            # Chemical score: insertion liability delta (donor residue kept vs
            # absent); the insertion residue itself is scanned in context.
            d_res = donor.residue(pos)
            ridx = d_res.index if d_res is not None else max(0, donor.sequence.find(donor_aa))
            inserted_state = donor.sequence[:ridx] + donor_aa + donor.sequence[ridx + 1:]
            chem = _liability_score(inserted_state) - _liability_score(donor.sequence)

            composite = 100 * (WEIGHTS["blend"][0] * structural
                               + WEIGHTS["blend"][1] * benefit
                               + WEIGHTS["blend"][2] * min(1, chem))
            composite = round(composite, 1)

            # Tier: default T2 for insertions (need structure verification)
            tier = _assign_tier(features, buried, is_vhh, donor_aa, pos)

            rationale = ["FR insertion: donor-specific residue not in germline"]
            if buried is True:
                rationale.append("Buried: structurally tolerated")
            elif buried is False:
                rationale.append("Exposed: immunogenicity risk")

            candidates.append(BackMutationCandidate(
                position=pos,
                donor_aa=donor_aa,
                human_aa="-",
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
            ))

    return BackMutationResult(
        chain_type=chain_type,
        candidates=candidates,
        v_gene=v_gene,
        germline_conservation=conservation,
        fr_indels=fr_indels,
    )


def _annotate_with_functional_conservation(
    result: BackMutationResult,
    germline_db,
    donor,
    structure_data=None,
) -> BackMutationResult:
    """P1: Annotate candidates with functional conservation scores.
    
    Adds functional conservation, evolutionary rate, and structural constraint
    scores to each candidate's rationale.
    """
    from .conservation import calculate_functional_conservation
    
    for candidate in result.candidates:
        scores = calculate_functional_conservation(
            candidate.position,
            result.chain_type,
            germline_db,
            donor,
            structure_data,
        )
        
        # Store scores in candidate (add new fields if needed)
        # For now, add to rationale
        if scores["functional_conservation"] > 0.8:
            candidate.rationale.append(
                f"functional conservation: {scores['functional_conservation']:.2f} "
                f"(germline: {scores['germline_conservation']:.2f}, "
                f"evolution: {scores['evolutionary_rate']:.2f}, "
                f"struct: {scores['structural_constraint']:.2f})"
            )
        elif scores["functional_conservation"] < 0.3:
            candidate.rationale.append(
                f"low functional conservation ({scores['functional_conservation']:.2f})"
            )
    
    return result


def _assign_tier(features, buried, is_vhh, donor_aa, pos) -> str:
    """Assign tier based on structural features (pre-structure adjustment).
    
    P0 improvement: This is the initial tier assignment.
    Step 3 structure adjustment will further refine based on:
    - side_chain_contact (P0: higher priority than backbone contact)
    - dev_risk_score (P0: developability risk filtering)
    """
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


def _assign_tier_with_structure(features, buried, is_vhh, donor_aa, pos,
                                 side_chain_contact=None, cdr_contact=None,
                                 antigen_contact=None, volume_strain=None,
                                 dev_risk_score=None) -> str:
    """P0 improvement: Assign tier with structure and developability context.
    
    Priority order:
    1. KEEP_DONOR (VHH hallmark, disulfide)
    2. T1: side_chain_contact + buried (highest structural priority)
    3. T1: interface_core + buried OR vernier + buried
    4. T2: side_chain_contact + exposed OR backbone_contact + buried
    5. T2: other structural features
    6. T3: exposed + high dev_risk OR no structural support
    """
    if "vhh_hallmark" in features:
        return "KEEP_DONOR"
    if "disulfide_cys" in features:
        return "KEEP_DONOR"
    
    # P0: side-chain contact takes highest priority
    if side_chain_contact is True:
        if buried is True:
            return "T1"  # Side-chain contact + buried = highest priority
        elif buried is False:
            # Exposed but side-chain contacts CDR: still important
            return "T2"
    elif side_chain_contact is False and (cdr_contact is True or antigen_contact is True):
        # P0: backbone contact only (no side-chain contact)
        # Lower priority than side-chain contact
        if buried is True:
            return "T2"  # Backbone contact + buried = T2
        elif buried is False:
            # Exposed + backbone contact: check dev_risk first
            if dev_risk_score is not None and dev_risk_score > 0.7:
                return "T3"  # High dev_risk overrides backbone contact
            return "T2"
    
    # P0: Developability risk filtering (check early for exposed positions)
    # Exposed positions with high dev risk should be demoted to T3
    if buried is False and dev_risk_score is not None and dev_risk_score > 0.7:
        return "T3"
    
    # Interface core (VH/VL packing)
    if "interface_core" in features:
        if buried is False:
            return "T2"
        return "T1"
    
    # Vernier zone (supports CDR conformation)
    if "vernier" in features:
        if buried is True or "canonical" in features:
            return "T1"
        return "T2"
    
    # Canonical structure
    if "canonical" in features:
        return "T2"
    
    # CDR contact (backbone-mediated, lower priority)
    if cdr_contact is True or antigen_contact is True:
        return "T2"
    
    # Buried in core
    if buried is True:
        return "T2"
    
    # Extended interface
    if "interface_extended" in features:
        return "T2"
    
    return "T3"


def _calculate_dev_risk_score(donor_aa: str, pos: str, donor, chain_type: str) -> float:
    """P0: Calculate developability risk score for a position.
    
    Scores based on sequence motifs that affect developability:
    - N-glycosylation (N-X-S/T)
    - Deamidation (N-G, N-S, N-H)
    - Isomerization (D-G, D-S, D-T)
    - Oxidation (M, W)
    - Acid hydrolysis (D-X)
    
    Returns: 0.0 (low risk) to 1.0 (high risk)
    """
    from .config import RISK_MOTIFS
    import re
    
    risk_score = 0.0
    
    # Get the residue and surrounding context
    seq = donor.sequence
    idx = seq.find(donor_aa) if donor_aa else -1
    if idx < 0:
        return 0.0
    
    # Check each risk motif
    for motif_name, pattern in RISK_MOTIFS.items():
        # Check if this position is part of a risk motif
        for match in re.finditer(pattern, seq):
            if match.start() <= idx < match.end():
                # This position is in a risk motif
                if "N-glycan" in motif_name:
                    risk_score = max(risk_score, 0.9)
                elif "NG" in motif_name:
                    risk_score = max(risk_score, 0.8)
                elif "NS" in motif_name:
                    risk_score = max(risk_score, 0.7)
                elif "DG" in motif_name:
                    risk_score = max(risk_score, 0.8)
                elif "M" in motif_name and len(match.group()) == 1:
                    risk_score = max(risk_score, 0.5)
                else:
                    risk_score = max(risk_score, 0.6)
                break
    
    return risk_score


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


def _rationale_with_structure(features, tier, buried, cdr_contact, ag_contact,
                               exposure, conservation, donor_aa, human_aa,
                               side_chain_contact=None, dev_risk_score=None) -> List[str]:
    """P0: Enhanced rationale with structure and developability context."""
    out = _rationale(features, tier, buried, cdr_contact, ag_contact,
                     exposure, conservation, donor_aa, human_aa)
    
    # P0: Add side-chain contact rationale
    if side_chain_contact is True:
        out.append("Side-chain mediated CDR contact (<4.5 A, structure)")
    
    # P0: Add developability risk rationale
    if dev_risk_score is not None and dev_risk_score > 0.7:
        out.append(f"High developability risk (score {dev_risk_score:.2f})")
    
    return out
