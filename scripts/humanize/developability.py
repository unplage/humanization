"""Sequence-level developability checks (no structure required).

WARNING: 以下检测仅基于序列模式，未考虑三维结构暴露状态。
埋藏在蛋白核心的位点实际风险较低，表面暴露位点风险更高。
如需精确评估，请结合 AF3/结构数据判断 relSASA 或溶剂可及性。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .config import RISK_MOTIFS


@dataclass
class DevelopabilityIssue:
    position: str
    motif: str
    sequence_context: str
    risk_level: str = "medium"   # low / medium / high
    structure_note: str = ""     # 结构暴露提示（需 AF3 数据）


def _assess_risk_level(motif: str) -> str:
    """根据风险类型评估风险等级（序列水平）。

    等级: high > medium > low
    """
    # 高风险
    if "N-glycan" in motif:
        return "high"
    if "deamidation" in motif and "NG" in motif:
        return "high"
    if "isomerization" in motif and "DG" in motif:
        return "high"
    if "DD" in motif:  # DD 高风险：异构化 + 酸性水解
        return "high"
    # 中风险
    if "deamidation" in motif:  # NS/NH/ND
        return "medium"
    if "isomerization" in motif:  # DS/DT/DH
        return "medium"
    if "acid hydrolysis" in motif:  # D-X (非DD)
        return "medium"
    if "oxidation" in motif:
        return "medium"
    if "unpaired Cys" in motif:
        return "medium"
    # 低风险
    if "base hydrolysis" in motif:
        return "low"
    if "met-lyscleavage" in motif:
        return "low"
    return "medium"


def scan_sequence(chain) -> List[DevelopabilityIssue]:
    """Scan one numbered chain for sequence-level developability risks.

    WARNING: 仅检测序列模式，未考虑结构暴露；埋藏位点风险可能被高估。

    保守二硫键 Cys（VH 22/92、kappa VL 23/88、lambda VL 22/88）是结构必需的
    配对 Cys，不计入"unpaired Cys"/"oxidation (C)"风险。同一个 Cys 位置的
    oxidation (C) 和 unpaired Cys 只报告一次，避免双重计数。
    """
    issues: List[DevelopabilityIssue] = []
    seq = chain.sequence.upper()
    # Conserved intradomain disulfide Cys: VH 22/92; kappa VL 23/88;
    # lambda VL has a one-residue-shorter FR1, so its first Cys sits at 22.
    conserved_cys = {f"{chain.chain_type}{n}" for n in
                     (22, 92) if chain.chain_type == "H"} | \
                    {f"{chain.chain_type}{n}" for n in
                     (22, 23, 88) if chain.chain_type == "L"}
    cys_seen = set()   # positions already reported under a C-related motif
    for motif, pattern in RISK_MOTIFS.items():
        cys_related = "unpaired Cys" in motif or motif.startswith("oxidation (C)")
        for m in re.finditer(pattern, seq):
            context = seq[max(0, m.start() - 3): m.end() + 3]
            pos = None
            for r in chain.residues:
                if r.index == m.start():
                    pos = r.pos
                    break
            if cys_related and pos in conserved_cys:
                continue   # 结构必需配对 Cys，非风险
            if cys_related and pos in cys_seen:
                continue   # 同一 Cys 在另一个 C 类 motif 已报告过，去重
            if cys_related and pos is not None:
                cys_seen.add(pos)
            risk_level = _assess_risk_level(motif)
            issues.append(DevelopabilityIssue(
                position=pos or f"seq{m.start() + 1}",
                motif=motif,
                sequence_context=context,
                risk_level=risk_level,
                structure_note="（需 AF3 结构数据确认暴露状态）",
            ))
    return issues


def count_risks(chains: List) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for c in chains:
        for iss in scan_sequence(c):
            out[iss.motif] = out.get(iss.motif, 0) + 1
    return out
