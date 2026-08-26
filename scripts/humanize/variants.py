"""Variant ladder assembly (V0 pure graft ... V3 + optional T3).

V0  pure graft (all-human framework, donor CDRs)
V1  V0 + Tier-1 back-mutations (structural pillars)
V2  V1 + Tier-2 back-mutations (recommended structural/CDR support)
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

    def build(name, desc, positions):
        g = graft_variant(donor, v_gene, j_gene, scheme, positions, is_vhh=is_vhh)
        return Variant(name=name, description=desc, graft=g, backmutations=positions)

    variants = []
    variants.append(build(
        f"{chain_type}_V0", "pure graft: human framework + donor CDRs", []))
    t1 = backmut.revert_positions(("T1",))
    variants.append(build(
        f"{chain_type}_V1", "V0 + Tier-1 back-mutations (structural pillars)", t1))
    t2 = backmut.revert_positions(("T1", "T2"))
    variants.append(build(
        f"{chain_type}_V2", "V0 + Tier-1/2 back-mutations", t2))
    # Tier 3: only exposed positions (immunogenicity drivers), capped
    # 按 composite 降序（收益优先），而非位置顺序。
    t3_exposed = [
        c.position for c in sorted(
            (c for c in backmut.candidates
             if c.tier == "T3" and (c.buried is False or c.buried is None)),
            key=lambda c: (-c.composite, c.position),
        )
    ][:extra_t3_max]
    # 无结构数据时 buried 恒为 None，"exposed" 筛选退化为按 composite 排序；
    # 描述中如实说明，避免误导。
    has_structure = any(c.buried is not None for c in backmut.candidates)
    v3_desc = ("V0 + Tier-1/2 + selected exposed Tier-3"
               if has_structure else
               "V0 + Tier-1/2 + top-composite Tier-3 (no structure data: "
               "exposure unknown, ranked by composite)")
    variants.append(build(
        f"{chain_type}_V3", v3_desc, t2 + t3_exposed))
    return variants
