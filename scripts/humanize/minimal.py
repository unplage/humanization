"""Affinity-preserving minimal-reversion design.

Implements the industry findings on "how to preserve affinity with the
fewest back-mutations" (see docs/improvement_roadmap.md):

  1. CVI homology metric  - Boehringer Ingelheim, JBC 300:105555 (2024):
     conservation of canonical + vernier + interface residues predicts both
     expression and affinity; reported per candidate germline.

  2. Minimal reversion set - structure mode: greedy set-cover over the
     framework<->CDR contacts lost when human residues replace donor
     residues. Only positions that actually support CDR loops are reverted;
     everything else stays human. Affinity-critical = maximal contact
     restoration with minimal mutations.

  3. Germline matrix variants - framework hopping: grafts on the runner-up
     germlines (top-N), each with its own back-mutation set. Empirically
     (BI 2024; CUMAb 2022) the best molecule is often NOT the most
     homologous graft.

  4. Paratope (SDR-style) grafting - structure mode: graft only the CDR
     residues that contact the antigen, plus structural pillars
     (XtalPi AI-guided paratope grafting, bioRxiv 2025).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .backmut import BackMutationResult, StructureHints
from .config import CANONICAL, INTERFACE_CORE, INTERFACE_EXTENDED, VERNIER_ZONE
from .germline import GermlineGene
from .graft import GraftResult, graft_chain, graft_variant
from .numbering import NumberedChain


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


CVI_SETS = {
    "H": (CANONICAL["H"] | VERNIER_ZONE["H"]
          | INTERFACE_CORE["H"] | INTERFACE_EXTENDED["H"]),
    "L": (CANONICAL["L"] | VERNIER_ZONE["L"]
          | INTERFACE_CORE["L"] | INTERFACE_EXTENDED["L"]),
}


def cvi_homology(donor: NumberedChain, germline: GermlineGene) -> float:
    """CVI (canonical/vernier/interface) residue identity donor vs germline.
    BI 2024: correlates with expression titer and retained affinity."""
    if germline.numbered is None:
        return 0.0
    ctype = donor.chain_type
    d = donor.posmap()
    g = germline.numbered.posmap()
    n = same = 0
    j_anchor = 103 if ctype == "H" else 98
    for pos, aa in d.items():
        num = _pos_num(pos)
        if num >= j_anchor:
            continue
        if num not in CVI_SETS[ctype]:
            continue
        if donor.region_of(pos) not in ("FR1", "FR2", "FR3"):
            continue
        if pos in g and g[pos]:
            n += 1
            same += 1 if g[pos] == aa else 0
    return round(same / n, 4) if n else 0.0


# ---------------------------------------------------------------------------
# Minimal reversion set (structure mode)
# ---------------------------------------------------------------------------

@dataclass
class MinimalReversion:
    positions: List[str] = field(default_factory=list)
    covered_contacts: int = 0
    total_contacts: int = 0
    method: str = "tier"          # "tier" | "set_cover"
    note: str = ""


def minimal_reversion_set(
    donor: NumberedChain,
    backmut: BackMutationResult,
    structure: Optional[StructureHints] = None,
    require_contact: bool = False,
) -> MinimalReversion:
    """Smallest back-mutation set that preserves CDR-supporting contacts.

    Set-cover formulation (structure mode):
      universe  = all (framework position, CDR residue) contact pairs of the
                  donor structure within 4.5 A, where the framework position
                  is a back-mutation candidate (donor != human germline);
      each candidate position covers its own contact pairs;
      greedy: repeatedly pick the candidate covering the most uncovered
      pairs; ties broken by higher composite score.
    Without structure data, falls back to Tier-1 positions (structural
    pillars by literature), which is the "minimal safe" set.
    """
    if structure is None or not (structure.data.get("cdr_partners") or {}):
        t1 = backmut.revert_positions(("T1",))
        return MinimalReversion(
            positions=t1,
            method="tier",
            note="no structure hints: Tier-1 (literature pillars) is the minimal safe set",
        )

    universe: set = set()
    coverage: Dict[str, set] = {}
    for c in backmut.candidates:
        partners = structure.cdr_partners(c.position)
        if partners:
            pairs = {(c.position, p) for p in partners}
            coverage[c.position] = pairs
            universe |= pairs
    if not universe:
        t1 = backmut.revert_positions(("T1",))
        return MinimalReversion(
            positions=t1,
            method="tier",
            note="no framework->CDR contacts found for candidates; fell back to Tier-1",
        )

    remaining = set(universe)
    chosen: List[str] = []
    candidates = {c.position: c for c in backmut.candidates}
    while remaining:
        best_pos, best_gain = None, -1
        for pos, pairs in coverage.items():
            if pos in chosen:
                continue
            gain = len(pairs & remaining)
            # tie-break: composite score, then structural score
            if gain > best_gain or (
                gain == best_gain and best_pos is not None
                and (candidates[pos].composite
                     > candidates[best_pos].composite)
            ):
                best_pos, best_gain = pos, gain
        if best_pos is None or best_gain <= 0:
            break
        chosen.append(best_pos)
        remaining -= coverage[best_pos]

    return MinimalReversion(
        positions=sorted(chosen, key=lambda p: (_pos_num(p), p)),
        covered_contacts=len(universe) - len(remaining),
        total_contacts=len(universe),
        method="set_cover",
        note="greedy set-cover over donor framework<->CDR contacts "
             "(4.5 A, AF3 structure)",
    )


# ---------------------------------------------------------------------------
# Paratope (SDR-style) grafting variant (structure mode)
# ---------------------------------------------------------------------------

def paratope_positions(donor: NumberedChain, structure: StructureHints) -> List[str]:
    """CDR residues that contact the antigen in the complex model."""
    out = []
    for r in donor.residues:
        if r.region.startswith("CDR") and structure.antigen_contact("", r.pos):
            out.append(r.pos)
    return out


def build_paratope_variant(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str,
    backmut: BackMutationResult,
    structure: StructureHints,
    is_vhh: bool = False,
) -> Optional[GraftResult]:
    """Precision graft: human framework + structural pillars (T1+T2) +
    only the antigen-contacting CDR residues (paratope). Requires the
    antigen complex structure; returns None otherwise."""
    if not (structure.data.get("antigen_contact") or {}):
        return None
    sdr = paratope_positions(donor, structure)
    if not sdr:
        return None

    base = graft_chain(donor, v_gene, j_gene, scheme, is_vhh=is_vhh)
    # revert T1+T2 pillars
    pillars = backmut.revert_positions(("T1", "T2"))
    # remove CDR residues that are NOT in the paratope (humanize them)
    force_human = []
    for pos in base.donor_positions:
        reg = donor.region_of(pos) or ""
        if reg.startswith("CDR") and pos not in sdr:
            force_human.append(pos)
    result = graft_variant(
        donor, v_gene, j_gene, scheme, pillars,
        is_vhh=is_vhh, force_human=force_human,
    )
    return result


# ---------------------------------------------------------------------------
# Germline matrix variants (framework hopping)
# ---------------------------------------------------------------------------

@dataclass
class MatrixEntry:
    germline: GermlineGene
    backmut: BackMutationResult
    cvi: float
    graft_v2: GraftResult


def matrix_alternatives(
    donor: NumberedChain,
    alternatives: List[Tuple[GermlineGene, dict]],
    j_gene: GermlineGene,
    scheme: str,
    is_vhh: bool = False,
    n: int = 3,
) -> List[MatrixEntry]:
    """Design V2-class variants on the top-N alternative germlines.
    Each alternative gets its own back-mutation set (analyzed on its own
    germline) — framework hopping, per BI/CUMAb evidence."""
    from .backmut import analyze_backmutations

    entries: List[MatrixEntry] = []
    for g, _s in alternatives[:n]:
        bm = analyze_backmutations(donor, g, is_vhh=is_vhh)
        positions = bm.revert_positions(("T1", "T2"))
        graft = graft_variant(donor, g, j_gene, scheme, positions, is_vhh=is_vhh)
        entries.append(MatrixEntry(
            germline=g,
            backmut=bm,
            cvi=cvi_homology(donor, g),
            graft_v2=graft,
        ))
    return entries
