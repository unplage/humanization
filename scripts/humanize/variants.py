"""Variant ladder assembly (V0 pure graft ... V3 + optional T3).

V0  pure graft (all-human framework, donor CDRs) — FR indel positions excluded
V1  V0 + Tier-1 back-mutations (structural pillars) — no indel positions
V2  V1 + Tier-2 back-mutations + donor FR insertions (recommended)
V3  V2 + selected Tier-3 (exposed, low-risk immunogenic positions)
SDR grafted variant is produced in structure mode (see sdr module).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .backmut import BackMutationResult
from .germline import GermlineGene
from .graft import GraftResult, graft_variant
from .numbering import NumberedChain


@dataclass
class Variant:
    name: str
    description: str
    graft: GraftResult
    backmutations: List[str] = field(default_factory=list)

    @property
    def sequence(self) -> str:
        return self.graft.sequence


def assemble_variants(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str,
    backmut: BackMutationResult,
    is_vhh: bool = False,
    extra_t3_max: int = 8,
    fr_indels: Optional[list] = None,
) -> List[Variant]:
    chain_type = donor.chain_type
    if fr_indels is None:
        fr_indels = backmut.fr_indels

    def _dedupe(positions):
        seen = set()
        out = []
        for p in positions:
            if p not in seen:
                seen.add(p)
                out.append(p)
        return out

    def build(name, desc, positions, exclude_indel=False):
        positions = _dedupe(positions)
        g = graft_variant(donor, v_gene, j_gene, scheme, positions,
                          is_vhh=is_vhh, exclude_indel=exclude_indel,
                          fr_indels=fr_indels)
        return Variant(name=name, description=desc, graft=g, backmutations=positions)

    # FR insertion positions (donor-only, not in germline)
    indel_ins = backmut.indel_insertion_positions
    # FR4 structural reversions (structure mode only; sanctioned exception)
    t_fr4 = backmut.revert_positions(("T_FR4",))

    # T2 sub-categorization by structural importance
    # T2a: buried + cdr_contact (highest priority - directly contacts CDR)
    # T2b: buried + (canonical/vernier/interface) OR exposed + cdr_contact
    # T2c: other T2 (buried only, low priority)
    t2a_positions = []
    t2b_positions = []
    t2c_positions = []
    
    for c in backmut.candidates:
        if c.tier != "T2":
            continue
        if c.buried is True and c.cdr_contact is True:
            # Highest priority: buried + CDR contact
            t2a_positions.append(c.position)
        elif c.buried is True and any(f in (c.features or "") for f in ["canonical", "vernier", "interface"]):
            # Medium priority: buried + structural feature
            t2b_positions.append(c.position)
        elif c.buried is False and c.cdr_contact is True:
            # Medium priority: exposed but contacts CDR
            t2b_positions.append(c.position)
        else:
            # Low priority: buried only or other
            t2c_positions.append(c.position)

    variants = []
    # V0: pure graft — no back-mutations, no indel positions
    # (germline lacks insertion positions, so they are absent from V0)
    variants.append(build(
        f"{chain_type}_V0",
        "pure graft: human framework + donor CDRs (FR indel excluded)",
        [], exclude_indel=True))
    # V1: T1 structural pillars only — no indel positions
    # (inherits V0's pure graft base, no donor insertions)
    t1 = backmut.revert_positions(("T1",))
    variants.append(build(
        f"{chain_type}_V1",
        "V0 + Tier-1 back-mutations (structural pillars)",
        t1, exclude_indel=True))
    
    # V2a: T1 + T2a (buried + CDR contact) - minimal functional set
    v2a_positions = _dedupe(t1 + t2a_positions)
    v2a_desc = f"V0 + T1 + T2a ({len(t2a_positions)} buried+CDR contact)"
    if indel_ins:
        v2a_desc += " + %d FR insertion(s)" % len(indel_ins)
    variants.append(build(
        f"{chain_type}_V2a", v2a_desc, v2a_positions))
    
    # V2b: T1 + T2a + T2b (buried + canonical/vernier/interface)
    v2b_positions = _dedupe(t1 + t2a_positions + t2b_positions)
    v2b_desc = f"V0 + T1 + T2a/T2b ({len(t2a_positions)}+{len(t2b_positions)} structural)"
    if indel_ins:
        v2b_desc += " + %d FR insertion(s)" % len(indel_ins)
    variants.append(build(
        f"{chain_type}_V2b", v2b_desc, v2b_positions))
    
    # V2: T1 + T2a + T2b + T2c + all donor FR insertions + FR4 structural (default include)
    # This is the current V2 logic (all T1+T2)
    t2 = backmut.revert_positions(("T1", "T2"))
    v2_positions = _dedupe(t2 + indel_ins + t_fr4)
    v2_desc = "V0 + Tier-1/2 back-mutations"
    if indel_ins:
        v2_desc += " + %d FR insertion(s): %s" % (len(indel_ins), ", ".join(indel_ins))
    if t_fr4:
        v2_desc += " + %d FR4 structural reversion(s)" % len(t_fr4)
    variants.append(build(
        f"{chain_type}_V2", v2_desc, v2_positions))
    
    # Tier 3: only exposed positions (immunogenicity drivers), capped
    t3_exposed = [
        c.position for c in sorted(
            (c for c in backmut.candidates
             if c.tier == "T3" and (c.buried is False or c.buried is None)),
            key=lambda c: (-c.composite, c.position),
        )
    ][:extra_t3_max]
    has_structure = any(c.buried is not None for c in backmut.candidates)
    v3_desc = ("V0 + Tier-1/2 + selected exposed Tier-3"
               if has_structure else
               "V0 + Tier-1/2 + top-composite Tier-3 (no structure data: "
               "exposure unknown, ranked by composite)")
    if indel_ins:
        v3_desc += " + %d FR insertion(s)" % len(indel_ins)
    variants.append(build(
        f"{chain_type}_V3", v3_desc, v2_positions + t3_exposed))
    return variants


def structure_guided_panel(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str,
    backmut: BackMutationResult,
    structure=None,
    is_vhh: bool = False,
    n_extra: int = 3,
) -> List[Variant]:
    """Graded multi-objective panel on top of V2.

    Starting from the V2 (T1+T2 + insertions + FR4) set, greedily add the
    candidate that recovers the most currently-uncovered framework->CDR
    contacts (from the donor structure); ties and the no-structure case fall
    back to the composite score. Each step yields one variant, giving a
    screening gradient of increasing back-mutation count / contact recovery.
    """
    chain_type = donor.chain_type

    def _dedupe(positions):
        seen = set()
        out = []
        for p in positions:
            if p not in seen:
                seen.add(p)
                out.append(p)
        return out

    base = _dedupe(
        list(backmut.revert_positions(("T1", "T2")))
        + backmut.indel_insertion_positions
        + backmut.revert_positions(("T_FR4",))
    )
    base_set = set(base)
    cand = {c.position: c for c in backmut.candidates}
    pool = [c for c in backmut.candidates if c.position not in base_set]
    if not pool:
        return []

    coverage = {}
    if structure is not None:
        for c in pool:
            partners = structure.cdr_partners(c.position)
            if partners:
                coverage[c.position] = set(partners)

    ordered: List[str] = []
    if coverage:
        uncovered = set().union(*coverage.values())
        for p in base:
            if structure is not None:
                uncovered -= structure.cdr_partners(p)
        remaining = {p: set(ps) for p, ps in coverage.items()}
        while len(ordered) < n_extra:
            best, best_gain = None, -1
            for p, pairs in remaining.items():
                if p in ordered:
                    continue
                gain = len(pairs & uncovered)
                if gain > best_gain or (
                    gain == best_gain and best is not None
                    and cand[p].composite > cand[best].composite
                ):
                    best, best_gain = p, gain
            if best is None or best_gain <= 0:
                break
            ordered.append(best)
            uncovered -= remaining[best]

    if len(ordered) < n_extra:
        rest = sorted((c for c in pool if c.position not in ordered),
                      key=lambda c: -c.composite)
        ordered += [c.position for c in rest[:n_extra - len(ordered)]]

    panel: List[Variant] = []
    for i in range(1, len(ordered) + 1):
        positions = _dedupe(base + ordered[:i])
        graft = graft_variant(
            donor, v_gene, j_gene, scheme, positions,
            is_vhh=is_vhh, fr_indels=backmut.fr_indels,
        )
        panel.append(Variant(
            name=f"{chain_type}_V_opt{i}",
            description=f"structure-guided panel step {i} (+{ordered[i - 1]})",
            graft=graft,
            backmutations=positions,
        ))
    return panel
