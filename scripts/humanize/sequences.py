"""Input parsing and chain classification.

Handles FASTA with VH/VL (Fab) or single VHH chains, detects chain types via
the numbering engine, and validates that inputs are full V domains.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .numbering import NumberedChain, is_vhh, number_heavy, number_light


class InputError(ValueError):
    pass


@dataclass
class InputChain:
    name: str
    sequence: str
    chain_type: str            # "H", "L"
    numbered: Optional[NumberedChain] = None
    is_vhh: bool = False
    vhh_score: float = 0
    warnings: List[str] = field(default_factory=list)

    @property
    def clean(self) -> str:
        return re.sub(r"[^A-Za-z]", "", self.sequence).upper()


def parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Parse a FASTA file into (name, raw_sequence) pairs."""
    records = []
    name, lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(lines)))
                name = line[1:].split()[0]
                lines = []
            else:
                lines.append(line)
    if name is not None:
        records.append((name, "".join(lines)))
    return records


def classify_sequences(records: List[Tuple[str, str]]) -> List[InputChain]:
    """Number every record and classify as VH/VHH or VL."""
    chains: List[InputChain] = []
    seen_names: Dict[str, int] = {}
    for name, seq in records:
        # Auto-rename duplicates to prevent downstream collisions
        if name in seen_names:
            seen_names[name] += 1
            name = f"{name}_{seen_names[name]}"
        else:
            seen_names[name] = 0
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        if len(s) < 80:
            raise InputError(f"[{name}] sequence too short ({len(s)} aa) to be a V domain")
        # try heavy first (VHH has no light chain), then light
        chain, chain_type = None, None
        errors = []
        try:
            chain = number_heavy(s)
            chain_type = chain.chain_type
        except ValueError as e:
            errors.append(str(e))
            try:
                chain = number_light(s)
                chain_type = chain.chain_type
            except ValueError as e2:
                errors.append(str(e2))
                raise InputError(
                    f"[{name}] could not number as VH or VL:\n  VH: {errors[0]}\n  VL: {errors[1]}"
                ) from e2
        ic = InputChain(
            name=name,
            sequence=s,
            chain_type=chain_type,
            numbered=chain,
            warnings=list(chain.warnings),
        )
        if chain_type == "H":
            ic.is_vhh, _score, _ = is_vhh(chain)
            ic.vhh_score = _score
        chains.append(ic)
    return chains


def detect_format(chains: List[InputChain], explicit: Optional[str] = None) -> str:
    """Determine Fab (VH+VL) vs VHH (single VH with camelid hallmark)."""
    n_h = sum(1 for c in chains if c.chain_type == "H")
    n_l = sum(1 for c in chains if c.chain_type == "L")
    if explicit:
        if explicit.lower() not in ("fab", "vhh", "auto"):
            raise InputError(f"unknown format '{explicit}' (fab|vhh|auto)")
        if explicit.lower() != "auto":
            return explicit.lower()
    if n_h == 1 and n_l == 1:
        return "fab"
    if n_h == 1 and n_l == 0:
        if chains[0].is_vhh:
            return "vhh"
        return "vhh_suspect"   # single heavy without hallmark: warn downstream
    if n_h == 2 and n_l == 0:
        return "vhh_pair"
    raise InputError(
        f"could not infer format: found {n_h} heavy and {n_l} light chain(s); "
        "expected VH+VL (Fab) or a single VHH"
    )


def parse_input(path: str, explicit_format: Optional[str] = None) -> Tuple[str, List[InputChain]]:
    """Full input pipeline: parse FASTA, classify, detect format."""
    records = parse_fasta(path)
    if not records:
        raise InputError(f"no sequences found in {path}")
    chains = classify_sequences(records)
    fmt = detect_format(chains, explicit_format)
    return fmt, chains
