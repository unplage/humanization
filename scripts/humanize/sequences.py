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


def _anarci_chain_type(s: str) -> Optional[str]:
    """Return 'H'/'K'/'L' from ANARCI's top HMM hit, or None when ANARCI is
    unavailable, errors, or finds no confident hit (truncated/garbage seq).

    ANARCI's HMMs are the most reliable chain-type classifier available
    (verified: VH/VHH->H, kappa->K, lambda->L with e-values ~1e-50+; no
    false positives on fragments). Used to route numbering so that light
    chains always reach number_light instead of being hijacked by the
    try-heavy-first fallback chain.
    """
    try:
        from .anarci_adapter import is_anarci_available
        if not is_anarci_available():
            return None
        from anarci import anarci as _anarci_fn
        result = _anarci_fn([("query", s.upper().strip())], scheme="kabat",
                            output=False)
        hits = result[1][0] if result and result[1] else []
        if not hits:
            return None
        ct = hits[0].get("chain_type")
        return ct if ct in ("H", "K", "L") else None
    except Exception:
        return None


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
        chain, chain_type, errors = None, None, []
        # 1) ANARCI classification (when available): H -> number_heavy,
        #    K/L -> number_light (built-in engine; AbRSA VL numbering deviates
        #    from the Kabat L27A-F insertion convention and has no lambda
        #    support).
        ct = _anarci_chain_type(s)
        if ct in ("K", "L"):
            try:
                chain = number_light(s)
                chain_type = chain.chain_type
            except ValueError as e:
                errors.append(str(e))
        elif ct == "H":
            try:
                chain = number_heavy(s)
                chain_type = chain.chain_type
            except ValueError as e:
                errors.append(str(e))
        # 2) Fallback (no ANARCI or no hit): try heavy first (VHH has no
        #    light chain), then light.
        if chain is None:
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
            else:
                # The built-in heavy engine treats a failed W103 search as a
                # V-region-only germline (whole tail -> CDR3, no FR4), which
                # would swallow a lambda VL (its CDR3 'W' never matches
                # W[GST]). User FASTA always carries FR4: if the sequence
                # numbers cleanly as light with an F98-anchored FR4, it IS
                # a light chain.
                if chain is not None and not any(r.pos == "H103" for r in chain.residues):
                    try:
                        light = number_light(s)
                        if any(r.pos == "L98" for r in light.residues):
                            chain, chain_type = light, light.chain_type
                    except ValueError:
                        pass
        if chain is None or chain_type is None:
            raise InputError(
                f"[{name}] could not number as VH or VL:\n  " + "\n  ".join(errors)
            )
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
