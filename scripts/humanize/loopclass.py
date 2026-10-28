"""Per-antibody loop-class weighting and position-set correction (D3).

The literature position sets (Foote & Winter 1992 vernier; Chothia 1985
interface; Chothia & Lesk 1987 / Al-Lazikani 1997 canonical) are population
averages derived from a small number of structures. Once a structure of THIS
antibody exists, the measured loop geometry decides whether a literature
claim actually applies:

  * **Removal (false positive).** Canonical support residues are loop
    specific. A position the literature assigns to CDR2 but which is more than
    one shell (> 8 A) from every CDR2 atom in this antibody does not support
    that loop's actual class, so the canonical claim is corrected out. A
    literature "vernier" residue that is > 12 A from *every* CDR is not in any
    vernier zone and its claim is removed.

    Both cutoffs are calibrated on 86 real antibody chains (see
    ``tests/test_pipeline.py::test_loop_class_thresholds_on_real_structure``):
    canonical support residues sit 3.2 A (median) from the loop they support
    and only 3% exceed 8 A; vernier residues sit 3.1 A (median) from the
    nearest loop and only 0.4% exceed 12 A. The cutoffs therefore sit just
    outside the empirical support distribution and fire rarely on genuine
    claims.

  * **Addition (false negative).** A framework residue in the *second shell*
    of a loop (side chain 4.5-8 A away, no side-chain contact) that no
    literature set contains is added as a measured ``loop_shell`` residue.
    It is deliberately WEAK evidence (weight 0.35, no tier effect): a
    side-chain contact with a CDR is already reported as ``cdr_contact_sc``,
    and a mere 8 A proximity is a neighbourhood, not a packing contact. The
    shell evidence therefore refines the ranking only.

    To avoid double counting, ``loop_shell`` is added only to positions that
    carry no other structural evidence (including ``buried``): the shell rule
    is conditional on burial, so adding it on top of ``buried`` would count
    the same physical evidence twice in the noisy-OR.

  * **Class weighting.** Canonical classes are defined by main-chain
    conformation, so canonical support is only meaningful when the loop itself
    is reliably modelled. The weight is scaled by the loop's mean pLDDT:
    a poorly modelled loop has no dependable class to support.

The correction is structure-mode only and is a strict no-op when the per-loop
hints are absent, so the portable (sequence-only) path is unchanged. The
canonical/loop membership table is Kabat-space, so canonical set correction is
skipped under IMGT numbering; the measured shell evidence and the vernier
correction (which do not depend on the table) still apply.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set

from .config import CANONICAL_SUPPORT_BY_LOOP, LOOP_CLASS

# Weak, ranking-only evidence: the residue's side chain is in the second
# shell of a CDR loop (no side-chain contact). Never creates a tier.
SHELL_FEATURE = "loop_shell"

# Any structural evidence at all. A measured shell residue is only ADDED to
# positions that have none of these, so the noisy-OR never counts the same
# physical evidence (burial + loop proximity) twice.
ANY_STRUCTURAL = {
    "vernier", "canonical", "interface_core", "interface_extended",
    "buried", "cdr_contact", "antigen_contact", SHELL_FEATURE,
}


@dataclass
class LoopCorrection:
    """Result of the per-antibody loop-class correction for one position."""

    features: List[str]
    canonical_factor: float = 1.0
    notes: List[str] = field(default_factory=list)
    applied: bool = False


def canonical_loops(chain_type: str, num: int) -> Set[str]:
    """Loops this position canonically supports (Kabat-space table)."""
    out: Set[str] = set()
    for loop, positions in CANONICAL_SUPPORT_BY_LOOP.get(chain_type, {}).items():
        if num in positions:
            out.add(loop)
    return out


def _has_loop_geometry(hints) -> bool:
    if hints is None:
        return False
    data = getattr(hints, "data", None)
    if not data:
        return False
    return bool(data.get("loop_min_dist"))


def _floats(mapping) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for key, value in (mapping or {}).items():
        try:
            out[key] = float(value)
        except (TypeError, ValueError):
            continue
    return out


def shell_evidence(hints, chain_type: str, pos: str) -> Dict[str, str]:
    """Measured loop proximity of one residue: {loop: "contact"|"shell"}.

    ``contact`` = side chain within the contact cutoff (<= 4.5 A) of a loop
    atom; ``shell`` = side chain in the second shell (4.5-8 A). Only the
    shell is new information: a side-chain contact with a CDR loop is already
    reported by the generic ``cdr_contact_sc`` hint.
    """
    out: Dict[str, str] = {}
    if not _has_loop_geometry(hints):
        return out
    direct = LOOP_CLASS["direct"]
    shell = LOOP_CLASS["shell"]
    for loop, d in _floats(hints.loop_min_dist_sc(pos)).items():
        if d <= direct:
            out[loop] = "contact"
        elif d <= shell:
            out[loop] = "shell"
    return out


def is_buriedish(hints, chain_type: str, pos: str) -> bool:
    """True when the residue is buried enough to pack a loop."""
    buried = hints.buried(chain_type, pos)
    if buried is True:
        return True
    rel = hints.rel_sasa(pos)
    return rel is not None and rel <= LOOP_CLASS["shell_rel_sasa"]


def vote_threshold(n_models: int) -> int:
    """Minimum number of consensus models that must agree."""
    if n_models <= 1:
        return 1
    frac = LOOP_CLASS["vote_fraction"]
    return max(1, int(math.ceil(frac * n_models)))


def _canonical_factor(loop_plddt: Optional[float]) -> float:
    """Scale the canonical weight by the confidence of the loop it supports."""
    if loop_plddt is None:
        return 1.0
    full = LOOP_CLASS["canonical_plddt_full"]
    low = LOOP_CLASS["canonical_plddt_low"]
    floor = LOOP_CLASS["canonical_factor_floor"]
    if loop_plddt >= full:
        return 1.0
    if loop_plddt <= low:
        return floor
    span = max(full - low, 1e-6)
    return round(floor + (1.0 - floor) * (loop_plddt - low) / span, 3)


def apply(hints, chain_type: str, pos: str, num: int, features: List[str],
          numbering: str = "kabat") -> LoopCorrection:
    """Correct literature loop-role features with the measured loop geometry.

    Parameters
    ----------
    hints : StructureHints
        Structure annotations (needs the ``loop_*`` keys).
    chain_type : "H" | "L"
    pos, num : donor position label and its numeric part.
    features : list of literature-derived structural features for the position.
    numbering : coordinate space of ``pos``/``num`` ("kabat" | "imgt").

    Returns a :class:`LoopCorrection`; ``features`` is a new list (the input is
    never mutated).
    """
    res = LoopCorrection(features=list(features))
    if not _has_loop_geometry(hints):
        return res

    dists = _floats(hints.loop_min_dist(pos))
    if not dists:
        return res

    lengths = hints.loop_lengths() or {}
    if set(lengths) - set(dists):
        # A loop present in the numbering but absent from the model: distances
        # would be missing and every "far" call would be an artefact.
        return res

    near = {loop for loop, d in dists.items() if d <= LOOP_CLASS["shell"]}
    far_all = all(d > LOOP_CLASS["far"] for d in dists.values())
    classes = hints.loop_classes() or {}
    applied = False

    # --- canonical: loop-specific support, confidence-weighted --------------
    if numbering == "kabat" and "canonical" in res.features:
        predicted = canonical_loops(chain_type, num)
        if predicted:
            supported = predicted & near
            if not supported:
                nearest = min(dists.get(loop, 1e9) for loop in predicted)
                res.features.remove("canonical")
                res.notes.append(
                    f"loop-class: canonical support claimed for "
                    f"{'/'.join(sorted(predicted))} but nearest such loop is "
                    f"{nearest:.1f} Å away; canonical feature corrected out")
                applied = True
            else:
                # A canonical class only exists for a reliably modelled loop:
                # scale the claim by the loop's mean pLDDT.
                plddts = []
                for loop in supported:
                    v = (hints.loop_plddt() or {}).get(loop)
                    if isinstance(v, (int, float)):
                        plddts.append(float(v))
                if plddts:
                    reg = max(plddts)
                    factor = _canonical_factor(reg)
                    if factor < 1.0:
                        res.canonical_factor = factor
                        res.notes.append(
                            f"loop-class: "
                            f"{'/'.join(sorted(supported))} modelled at mean "
                            f"pLDDT {reg:.0f}; canonical support weight "
                            f"×{factor:.2f}")
                        applied = True

    # --- vernier: must be within the CDR neighbourhood at all --------------
    if "vernier" in res.features and far_all:
        nearest = min(dists.values())
        res.features.remove("vernier")
        res.notes.append(
            f"loop-class: literature vernier is {nearest:.1f} Å from every CDR "
            f"in this antibody; vernier feature corrected out")
        applied = True

    # --- addition: measured second-shell loop evidence --------------------
    # Require the residue's own SIDE CHAIN (not just its backbone) to be in
    # the loop's second shell: a backbone-only neighbour does not mean the
    # donor side chain is doing any work, and the substitution changes the
    # side chain. A side-chain CONTACT is already reported as cdr_contact_sc,
    # so only the shell band adds new information - and only for positions
    # without any other structural evidence, otherwise the burial and the
    # loop proximity would be counted twice in the noisy-OR.
    evidence = shell_evidence(hints, chain_type, pos)
    shell_loops = {loop for loop, kind in evidence.items() if kind == "shell"}
    if shell_loops and is_buriedish(hints, chain_type, pos) \
            and not (ANY_STRUCTURAL & set(res.features)):
        n_models = int((hints.data or {}).get("n_models") or 1)
        need = vote_threshold(n_models)
        votes = (hints.data or {}).get("loop_shell_votes")
        if votes is None:
            # single-model run: no votes to disagree with
            agreed = list(shell_loops)
        else:
            per_pos = votes.get(pos) or {}
            agreed = [loop for loop in shell_loops
                      if int(per_pos.get(loop, 0) or 0) >= need]
        if agreed:
            best = min(agreed, key=lambda loop: hints.loop_min_dist_sc(pos).get(loop, 1e9))
            res.features.append(SHELL_FEATURE)
            res.notes.append(
                f"loop-class: buried and side chain "
                f"{hints.loop_min_dist_sc(pos).get(best, 0.0):.1f} Å from "
                f"{best} (second shell, {need}/{n_models} models; weak "
                f"evidence, ranking only)")
            applied = True

    res.applied = applied
    return res
