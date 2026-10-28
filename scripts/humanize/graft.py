"""CDR grafting onto human germline frameworks.

The grafted V domain is assembled position-wise:
    FR1-FR3  <- human germline V gene
    CDR1-3   <- donor (per chosen CDR definition)
    FR4      <- human J gene
    VHH      <- additionally keep donor residues at the FR2 hallmark
                positions (Kabat 37/44/45/47); conserved-Cys protection is
                handled by the back-mutation tier rules (backmut.py), not here
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .germline import GermlineDB, GermlineGene
from .germline import _region_of_pos, _imgt_region_of_pos
from .numbering import (
    NumberedChain,
    NumberedResidue,
    _cdr_segments,
)
from .schemespec import cdr_windows, get_spec
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
        # Chothia CDR-H1 is 26-32 (not 26-35, which is the AbM/Martin view).
        "H": {"CDR1": (26, 32), "CDR2": (52, 56), "CDR3": (93, 102)},
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


def is_cdr_loop_position(chain_type: str, num: int, scheme: str = "kabat",
                          numbering: str = "kabat") -> bool:
    """True if the position is part of a grafted CDR loop (incl. CDR3-H 93/94).

    Used wherever loop membership matters for sequence-level logic (grafting,
    paratope/SDR, structure contact hints, back-mutation filtering)."""
    return num in cdr_positions_for(scheme, chain_type, numbering)


def _pos_num(pos: str) -> int:
    return int("".join(c for c in pos if c.isdigit()))


def cdr_positions_for(scheme: str, chain_type: str, numbering: str = "kabat") -> set:
    """Position numbers included in each CDR for (numbering, scheme)."""
    out = set()
    for (lo, hi) in cdr_windows(numbering, scheme)[chain_type].values():
        out.update(range(lo, hi + 1))
    return out


def loop_of_position(chain_type: str, num: int, scheme: str = "kabat",
                     numbering: str = "kabat") -> Optional[str]:
    """CDR1/CDR2/CDR3 name for a CDR position number, else None.

    The per-loop name is what the loop-class correction needs: the literature
    canonical support residues are loop-specific, and only a structure can tell
    whether this antibody's residue actually supports the loop it is assigned
    to. ``cdr_windows`` returns disjoint windows ordered CDR1/2/3.
    """
    for name, (lo, hi) in cdr_windows(numbering, scheme)[chain_type].items():
        if lo <= num <= hi:
            return name
    return None


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


def _order_key(pos: str):
    return (pos[0], _pos_num(pos), pos)


def _sorted_positions(keys, donor: NumberedChain, numbering: str) -> List[str]:
    """Sequence-order the assembled position labels.

    Kabat labels sort correctly by (prefix, number, label). IMGT insertion
    letters at 111/112 are NOT emitted in sequence order by ANARCI
    (e.g. H112B, H112A, H112), so labels present in the donor are ordered by
    the donor's residue order; non-donor labels (FR deletions, J residues
    beyond the donor) are bare ints merged in numeric position."""
    keys = list(keys)
    if numbering != "imgt":
        return sorted(keys, key=_order_key)
    dindex = {r.pos: i for i, r in enumerate(donor.residues)}
    donor_keys = sorted((k for k in keys if k in dindex), key=lambda k: dindex[k])
    rest = sorted((k for k in keys if k not in dindex),
                  key=lambda k: (_pos_num(k), k))

    def _base(p: str) -> int:
        return _pos_num(p)

    out: List[str] = []
    di = ri = 0
    while di < len(donor_keys) or ri < len(rest):
        if di >= len(donor_keys):
            out.append(rest[ri]); ri += 1
        elif ri >= len(rest):
            out.append(donor_keys[di]); di += 1
        elif _base(rest[ri]) < _base(donor_keys[di]):
            out.append(rest[ri]); ri += 1
        else:
            out.append(donor_keys[di]); di += 1
    return out


def chain_from_origin_map(
    seq: str,
    out: Dict[str, str],
    chain_type: str,
    donor: Optional[NumberedChain],
    gmap_src: Optional[NumberedChain],
    jmap_src: Optional[NumberedChain],
    scheme: str = "kabat",
    numbering: str = "kabat",
) -> NumberedChain:
    """Build a numbered chain directly from an assembled position map.

    Region labels come from the source chains (donor / germline / J) so the
    strict-Kabat annotations (e.g. H93/H94 = FR3) match the donor. Positions
    absent from all sources fall back to a numeric region rule.
    """
    region_fn = _imgt_region_of_pos if numbering == "imgt" else _region_of_pos
    res = []
    if donor is not None:
        order = _sorted_positions(out.keys(), donor, numbering)
    else:
        order = sorted(out.keys(), key=_order_key)
    for i, pos in enumerate(order):
        aa = out[pos]
        region = ""
        for s in (x for x in (donor, gmap_src, jmap_src) if x is not None):
            r = s.residue(pos)
            if r is not None:
                region = r.region
                break
        if not region:
            region = region_fn(pos)
        res.append(NumberedResidue(pos, aa, region, i))
    chain = NumberedChain(chain_type, seq, res)
    chain.cdrs = _cdr_segments(chain, scheme, numbering)
    return chain


def graft_chain(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str = "kabat",
    is_vhh: bool = False,
    exclude_indel: bool = False,
    fr_indels: Optional[list] = None,
    numbering: str = "kabat",
) -> GraftResult:
    """Build one humanized V domain (CDR grafting).

    Framework positions are populated from the human germline using the
    donor<->germline correspondence (so a donor FR insertion never overwrites
    the germline residue). CDRs come from the donor.

    exclude_indel: if True, exclude donor-only FR positions (insertions) from
                   the graft. Used for V0 (pure graft without donor indels).
    fr_indels:     pre-computed FR indel list (keeps user-confirmed insertion
                   positions); detected internally when omitted.
    numbering:     coordinate space of donor/germline ("kabat" | "imgt").
    """
    chain_type = donor.chain_type
    if v_gene.numbered is None or j_gene.numbered is None:
        raise ValueError(f"[{chain_type}] germline gene without numbering: {v_gene.gene_id}")
    gmap_src: NumberedChain = v_gene.numbered
    jmap_src: NumberedChain = j_gene.numbered
    warnings: List[str] = []
    spec = get_spec(numbering)

    dmap = donor.posmap()                # position -> aa (donor)
    gmap = gmap_src.posmap()             # human germline V
    jmap = jmap_src.posmap()             # human J (FR4)

    if fr_indels is None:
        fr_indels = detect_fr_indels(donor, v_gene, numbering)
    corr = build_fr_correspondence(donor, v_gene, fr_indels, numbering)

    cdr_nums = cdr_positions_for(scheme, chain_type, numbering)
    j_anchor = spec.j_anchor_for(chain_type)
    hallmark = spec.vhh_hallmark
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
        if (is_vhh and chain_type == "H" and num in hallmark
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
            if is_vhh and chain_type == "H" and _pos_num(dpos) in hallmark:
                origin[dpos] = "donor(vhh)"
            else:
                origin[dpos] = "donor(indel)"

    # ---- 5. FR4 / J region (human J by construction; donor fallback)
    for pos in set(jmap) | {p for p in dmap if _pos_num(p) >= j_anchor}:
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

    seq = "".join(aa for pos in _sorted_positions(out.keys(), donor, numbering)
                  for aa in (out[pos],))

    # Build the numbered chain directly from the assembled position map.
    # Positions come from consistently-numbered donor/germline/J maps, so
    # framework lengths (e.g. VH3-family FR2 = 13 with gap at H49 vs a donor
    # that fills H49) are preserved by construction.
    numbered = chain_from_origin_map(seq, out, chain_type, donor, gmap_src,
                                     jmap_src, scheme, numbering)

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
    donor: Optional[NumberedChain] = None,
    numbering: str = "kabat",
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
    if donor is not None:
        order = _sorted_positions(aas.keys(), donor, numbering)
    else:
        order = sorted(aas.keys(), key=_order_key)
    seq = "".join(aas[pos] for pos in order)
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
    numbering: str = "kabat",
) -> GraftResult:
    """Graft + apply a list of back-mutations (position labels like 'H67').

    backmutations: donor residues to restore at these positions (position ->
    donor aa). force_human: positions to keep human even if recommended.
    exclude_indel: if True, exclude donor FR insertions from base graft (V0).
    fr_indels:     pre-computed FR indels (keeps user-confirmed insertion pos).
    """
    base = graft_chain(donor, v_gene, j_gene, scheme, is_vhh,
                       exclude_indel=exclude_indel, fr_indels=fr_indels,
                       numbering=numbering)
    if not backmutations and not force_human:
        return base
    if v_gene.numbered is None or j_gene.numbered is None:
        raise ValueError(f"[{base.chain_type}] germline gene without numbering: {v_gene.gene_id}")
    dmap = donor.posmap()
    gmap = v_gene.numbered.posmap()
    jmap = j_gene.numbered.posmap()
    corr = base.fr_correspondence or build_fr_correspondence(
        donor, v_gene, fr_indels, numbering)
    out = dict(base.origin)
    for pos in backmutations or []:
        if pos in dmap and dmap[pos]:
            out[pos] = "donor"
    for pos in force_human or []:
        if pos in out and out[pos] == "donor(vhh)":
            continue   # VHH hallmark must not be touched
        out[pos] = "germline"
    # rebuild sequence from origin map (germline labels resolve through corr)
    seq, seqs = _materialize(out, corr, dmap, gmap, jmap, donor, numbering)
    # rebuild the numbered chain from the position map (see graft_chain)
    numbered = chain_from_origin_map(
        seq, seqs, base.chain_type, donor,
        v_gene.numbered, j_gene.numbered, scheme, numbering,
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
