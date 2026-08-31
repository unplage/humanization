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
) -> List[Variant]:
    chain_type = donor.chain_type

    def build(name, desc, positions, exclude_indel=False):
        g = graft_variant(donor, v_gene, j_gene, scheme, positions,
                          is_vhh=is_vhh, exclude_indel=exclude_indel)
        return Variant(name=name, description=desc, graft=g, backmutations=positions)

    # FR insertion positions (donor-only, not in germline)
    indel_ins = backmut.indel_insertion_positions

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
    # V2: T1 + T2 + all donor FR insertions (user requirement: default include)
    t2 = backmut.revert_positions(("T1", "T2"))
    v2_positions = t2 + indel_ins
    v2_desc = "V0 + Tier-1/2 back-mutations"
    if indel_ins:
        v2_desc += " + %d FR insertion(s): %s" % (len(indel_ins), ", ".join(indel_ins))
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
