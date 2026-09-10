"""CDR grafting onto human germline frameworks.

The grafted V domain is assembled position-wise:
    FR1-FR3  <- human germline V gene
    CDR1-3   <- donor (per chosen CDR definition)
    FR4      <- human J gene
    VHH      <- additionally keep donor residues at the FR2 hallmark
                positions (Kabat 37/44/45/47) and at structural Cys pairs
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .germline import GermlineDB, GermlineGene
from .germline import _region_of_pos
from .numbering import (
    NumberedChain,
    NumberedResidue,
    _cdr_segments,
)
from .config import VHH_HALLMARK
from .fr_indel import (
    FRCorrespondence,
    build_fr_correspondence,
    detect_fr_indels,
)

# CDR position sets per scheme (Kabat space). Framework-flanking positions
# included in Chothia/AbM/IMGT CDR1 (the 26-30 stem) are grafted as well.
#
# NOTE on CDR3-H: the strict-Kabat CDR3 is 95-102, but the first two residues
# of the CDR3 loop occupy Kabat positions 93/94 (they are labelled FR3 in the
# strict scheme; IMGT 105-106). The full loop must be grafted from the donor,
# so the GRAFT region for CDR3-H is 93-102 in every scheme. The Kabat table
# in the report still displays H93/H94 under FR3.
CDR_POS_SETS = {
    "kabat": {
        "H": {"CDR1": (31, 35), "CDR2": (50, 65), "CDR3": (93, 102)},
        "L": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
    },
    "chothia": {
        "H": {"CDR1": (26, 35), "CDR2": (52, 56), "CDR3": (93, 102)},
        "L": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
    },
    "abm": {
        "H": {"CDR1": (26, 35), "CDR2": (50, 58), "CDR3": (93, 102)},
        "L": {"CDR1": (24, 34), "CDR2": (50, 56), "CDR3": (89, 97)},
    },
    "imgt": {
        # IMGT CDR definitions (Lefranc 2003)
        # Note: These are Kabat position numbers that correspond to IMGT CDR boundaries
        # IMGT CDR1: IMGT 27-38 = Kabat 26-35
        # IMGT CDR2: IMGT 56-65 = Kabat 51-57
        # IMGT CDR3: IMGT 105-117 = Kabat 93-102
        "H": {"CDR1": (26, 35), "CDR2": (51, 57), "CDR3": (93, 102)},
        "L": {"CDR1": (24, 34), "CDR2": (51, 56), "CDR3": (89, 97)},
    },
}


def is_cdr_loop_position(chain_type: str, num: int, scheme: str = "kabat") -> bool:
    """True if the position is part of a grafted CDR loop (incl. CDR3-H 93/94).

    Used wherever loop membership matters for sequence-level logic (grafting,
    paratope/SDR, structure contact hints, back-mutation filtering)."""
    return num in cdr_positions_for(scheme, chain_type)


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def cdr_positions_for(scheme: str, chain_type: str) -> set:
    """Position numbers included in each CDR for the scheme (for grafting)."""
    out = set()
    for (lo, hi) in CDR_POS_SETS[scheme][chain_type].values():
        out.update(range(lo, hi + 1))
    return out


@dataclass
class GraftResult:
    scheme: str
    chain_type: str
    sequence: str
    numbered: NumberedChain
    origin: Dict[str, str]      # position -> "donor" | "germline" | "j"
    donor_positions: List[str]
    warnings: List[str] = field(default_factory=list)
    fr_indels: List = field(default_factory=list)  # List[FRIndel]
    fr_correspondence: Optional["FRCorrespondence"] = None


def chain_from_origin_map(
    seq: str,
    out: Dict[str, str],
    chain_type: str,
    donor: Optional[NumberedChain],
    gmap_src: Optional[NumberedChain],
    jmap_src: Optional[NumberedChain],
) -> NumberedChain:
    """Build a numbered chain directly from an assembled position map.

    Region labels come from the source chains (donor / germline / J) so the
    strict-Kabat annotations (e.g. H93/H94 = FR3) match the donor. Positions
    absent from all sources fall back to a numeric region rule.
    """
    sources = [s for s in (donor, gmap_src, jmap_src) if s is not None]
    res = []
    for i, (pos, aa) in enumerate(sorted(out.items(),
                                         key=lambda kv: (kv[0][0], _pos_num(kv[0]), kv[0]))):
        region = ""
        for s in sources:
            r = s.residue(pos)
            if r is not None:
                region = r.region
                break
        if not region:
            region = _region_of_pos(pos)
        res.append(NumberedResidue(pos, aa, region, i))
    chain = NumberedChain(chain_type, seq, res)
    chain.cdrs = _cdr_segments(chain)
    return chain


def _order_key(pos: str):
    return (pos[0], _pos_num(pos), pos)


def graft_chain(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str = "kabat",
    is_vhh: bool = False,
    exclude_indel: bool = False,
    fr_indels: Optional[list] = None,
) -> GraftResult:
    """Build one humanized V domain (CDR grafting).

    Framework positions are populated from the human germline using the
    donor<->germline correspondence (so a donor FR insertion never overwrites
    the germline residue). CDRs come from the donor.

    exclude_indel: if True, exclude donor-only FR positions (insertions) from
                   the graft. Used for V0 (pure graft without donor indels).
    fr_indels:     pre-computed FR indel list (keeps user-confirmed insertion
                   positions); detected internally when omitted.
    """
    chain_type = donor.chain_type
    if v_gene.numbered is None or j_gene.numbered is None:
        raise ValueError(f"[{chain_type}] germline gene without numbering: {v_gene.gene_id}")
    gmap_src: NumberedChain = v_gene.numbered
    jmap_src: NumberedChain = j_gene.numbered
    warnings: List[str] = []

    dmap = donor.posmap()                # position -> aa (donor)
    gmap = gmap_src.posmap()             # human germline V
    jmap = jmap_src.posmap()             # human J (FR4)

    if fr_indels is None:
        fr_indels = detect_fr_indels(donor, v_gene)
    corr = build_fr_correspondence(donor, v_gene, fr_indels)

    cdr_nums = cdr_positions_for(scheme, chain_type)
    j_anchor = 103 if chain_type == "H" else 98
    insertion_set = set(corr.insertion_positions)

    out: Dict[str, str] = {}
    origin: Dict[str, str] = {}

    # ---- 1. CDR loops from donor (donor defines the loop, incl. insertions)
    for pos, aa in dmap.items():
        num = _pos_num(pos)
        if num >= j_anchor:
            continue
        if num in cdr_nums and aa:
            out[pos] = aa
            origin[pos] = "donor"

    # ---- 2. Framework FR1-FR3: human germline template, aligned via corr
    for gpos, gaa in gmap.items():
        num = _pos_num(gpos)
        if num >= j_anchor or num in cdr_nums:
            continue
        dpos = corr.germline_to_donor.get(gpos)
        key = dpos if dpos is not None else gpos
        if (is_vhh and chain_type == "H" and num in VHH_HALLMARK
                and dpos is not None and dmap.get(dpos)):
            out[key] = dmap[dpos]
            origin[key] = "donor(vhh)"
        else:
            out[key] = gaa
            origin[key] = "germline"

    # ---- 3. Deletions (germline residues the donor lacks) stay human
    for gpos in corr.deletion_positions:
        if gpos not in out and gmap.get(gpos):
            out[gpos] = gmap[gpos]
            origin[gpos] = "germline"

    # ---- 4. Donor insertions (FR): kept unless pure graft (V0)
    for dpos in corr.insertion_positions:
        if exclude_indel:
            continue
        if dmap.get(dpos):
            out[dpos] = dmap[dpos]
            if is_vhh and chain_type == "H" and _pos_num(dpos) in VHH_HALLMARK:
                origin[dpos] = "donor(vhh)"
            else:
                origin[dpos] = "donor(indel)"

    # ---- 5. FR4 / J region (human J by construction; donor fallback)
    fr4_positions = sorted(
        (set(jmap) | {p for p in dmap if _pos_num(p) >= j_anchor}),
        key=_order_key,
    )
    for pos in fr4_positions:
        if _pos_num(pos) < j_anchor:
            continue
        if jmap.get(pos):
            out[pos] = jmap[pos]
            origin[pos] = "j"
        elif dmap.get(pos):
            # Position in donor but not covered by the J gene — keep the
            # donor residue so FR4 is not truncated (warn).
            out[pos] = dmap[pos]
            origin[pos] = "donor"
            warnings.append(
                f"[{chain_type}] FR4 position {pos} absent from J "
                f"gene; keeping donor residue to avoid truncation")

    seq = "".join(aa for pos, aa in sorted(out.items(), key=lambda kv: _order_key(kv[0])))

    # Build the numbered chain directly from the assembled position map.
    # Positions come from consistently-numbered donor/germline/J maps, so
    # framework lengths (e.g. VH3-family FR2 = 13 with gap at H49 vs a donor
    # that fills H49) are preserved by construction.
    numbered = chain_from_origin_map(seq, out, chain_type, donor, gmap_src, jmap_src)

    donor_positions = [pos for pos in out if origin[pos].startswith("donor")]

    return GraftResult(
        scheme=scheme,
        chain_type=chain_type,
        sequence=seq,
        numbered=numbered,
        origin=origin,
        donor_positions=donor_positions,
        warnings=warnings,
        fr_indels=fr_indels,
        fr_correspondence=corr,
    )


def _materialize(
    origin: Dict[str, str],
    corr: "FRCorrespondence",
    dmap: Dict[str, str],
    gmap: Dict[str, str],
    jmap: Dict[str, str],
) -> Tuple[str, Dict[str, str]]:
    """Turn an origin map into the amino-acid sequence.

    Donor labels are resolved directly from the donor map; germline labels are
    resolved through the correspondence (a donor label may correspond to a
    different germline position because of an upstream insertion).
    """
    aas: Dict[str, str] = {}
    for pos, src in origin.items():
        if src in ("donor", "donor(vhh)", "donor(indel)"):
            aas[pos] = dmap.get(pos, "")
        elif src == "germline":
            gpos = corr.donor_to_germline.get(pos, pos)
            aas[pos] = gmap.get(gpos) or dmap.get(pos, "")
        else:  # J
            aas[pos] = jmap.get(pos) or dmap.get(pos, "")
    seq = "".join(aa for pos, aa in sorted(aas.items(), key=lambda kv: _order_key(kv[0])))
    return seq, aas


def graft_variant(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str,
    backmutations: Optional[List[str]] = None,
    is_vhh: bool = False,
    force_human: Optional[List[str]] = None,
    exclude_indel: bool = False,
    fr_indels: Optional[list] = None,
) -> GraftResult:
    """Graft + apply a list of back-mutations (position labels like 'H67').

    backmutations: donor residues to restore at these positions (position ->
    donor aa). force_human: positions to keep human even if recommended.
    exclude_indel: if True, exclude donor FR insertions from base graft (V0).
    fr_indels:     pre-computed FR indels (keeps user-confirmed insertion pos).
    """
    base = graft_chain(donor, v_gene, j_gene, scheme, is_vhh,
                       exclude_indel=exclude_indel, fr_indels=fr_indels)
    if not backmutations and not force_human:
        return base
    if v_gene.numbered is None or j_gene.numbered is None:
        raise ValueError(f"[{base.chain_type}] germline gene without numbering: {v_gene.gene_id}")
    dmap = donor.posmap()
    gmap = v_gene.numbered.posmap()
    jmap = j_gene.numbered.posmap()
    corr = base.fr_correspondence or build_fr_correspondence(donor, v_gene, fr_indels)
    out = dict(base.origin)
    for pos in backmutations or []:
        if pos in dmap and dmap[pos]:
            out[pos] = "donor"
    for pos in force_human or []:
        if pos in out and out[pos] == "donor(vhh)":
            continue   # VHH hallmark must not be touched
        out[pos] = "germline"
    # rebuild sequence from origin map (germline labels resolve through corr)
    seq, seqs = _materialize(out, corr, dmap, gmap, jmap)
    # rebuild the numbered chain from the position map (see graft_chain)
    numbered = chain_from_origin_map(
        seq, seqs, base.chain_type, donor,
        v_gene.numbered, j_gene.numbered,
    )
    return GraftResult(
        scheme=scheme,
        chain_type=base.chain_type,
        sequence=seq,
        numbered=numbered,
        origin=out,
        donor_positions=[p for p in out if out[p].startswith("donor")],
        warnings=base.warnings,
        fr_indels=base.fr_indels,
        fr_correspondence=corr,
    )
