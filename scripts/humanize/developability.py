"""Sequence-level developability checks (no structure required)."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List

from .config import RISK_MOTIFS


@dataclass
class DevelopabilityIssue:
    position: str
    motif: str
    sequence_context: str


def scan_sequence(chain) -> List[DevelopabilityIssue]:
    """Scan one numbered chain for sequence-level developability risks."""
    issues: List[DevelopabilityIssue] = []
    seq = chain.sequence.upper()
    for motif, pattern in RISK_MOTIFS.items():
        for m in re.finditer(pattern, seq):
            context = seq[max(0, m.start() - 3): m.end() + 3]
            pos = None
            for r in chain.residues:
                if r.index == m.start():
                    pos = r.pos
                    break
            issues.append(DevelopabilityIssue(
                position=pos or f"seq{m.start() + 1}",
                motif=motif,
                sequence_context=context,
            ))
    return issues


def count_risks(chains: List) -> Dict[str, int]:
    out: Dict[str, int] = {}
    for c in chains:
        for iss in scan_sequence(c):
            out[iss.motif] = out.get(iss.motif, 0) + 1
    return out
