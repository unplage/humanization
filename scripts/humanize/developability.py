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
    """根据风险类型评估风险等级（序列水平）。"""
    if "N-glycan" in motif:
        return "high"
    if "deamidation" in motif and "NG" in motif:
        return "high"
    if "isomerization" in motif and "DG" in motif:
        return "high"
    if "acid hydrolysis" in motif:
        return "medium"
    if "oxidation" in motif:
        return "medium"
    if "base hydrolysis" in motif:
        return "low"
    if "met-lyscleavage" in motif:
        return "low"
    if "unpaired Cys" in motif:
        return "medium"
    return "medium"


def scan_sequence(chain) -> List[DevelopabilityIssue]:
    """Scan one numbered chain for sequence-level developability risks.

    WARNING: 仅检测序列模式，未考虑结构暴露；埋藏位点风险可能被高估。

    保守二硫键 Cys（VH 22/92、VL 23/88）是结构必需的配对 Cys，不计入
    "unpaired Cys" / "oxidation (C)" 风险。
    """
    issues: List[DevelopabilityIssue] = []
    seq = chain.sequence.upper()
    conserved_cys = {f"{chain.chain_type}{n}" for n in
                     (22, 92) if chain.chain_type == "H"} | \
                    {f"{chain.chain_type}{n}" for n in
                     (23, 88) if chain.chain_type == "L"}
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
