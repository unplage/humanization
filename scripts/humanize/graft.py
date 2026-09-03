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
from .fr_indel import detect_fr_indels

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


def _is_shifted_germline_residue(
    donor_pos: str,
    dmap: Dict[str, str],
    gmap: Dict[str, str],
    fr_indels: list,
    chain_type: str = "H",
) -> Optional[str]:
    """Check if a donor-only FR position is a shifted germline residue.
    
    When the donor has an upstream insertion (e.g. H6), downstream positions
    get shifted (donor H6A = germline H6). This function detects this by
    checking if the donor's downstream sequence matches the germline's
    sequence from the corresponding position.
    
    Returns the germline amino acid if this is a shifted residue, else None.
    """
    if not fr_indels:
        return None
    
    d_num = _pos_num(donor_pos)
    
    for indel in fr_indels:
        if indel.indel_type != "insertion":
            continue
        ins_num = _pos_num(indel.position)
        if d_num <= ins_num:
            continue  # This position is before or at the insertion
        
        # Check if downstream sequence matches germline
        # Donor H6A+ should match germline H6+ (shifted by 1)
        ins_base = indel.position.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        germ_pos = ins_base  # germline position corresponding to donor insertion + 1
        
        # Build sequences from donor position and germline position
        donor_seq = []
        germ_seq = []
        for offset in range(5):
            d_p = f"{chain_type}{d_num + offset}"
            g_p = f"{chain_type}{int(''.join(c for c in germ_pos if c.isdigit())) + offset}"
            if d_p in dmap:
                donor_seq.append(dmap[d_p])
            if g_p in gmap:
                germ_seq.append(gmap[g_p])
        
        if len(donor_seq) >= 3 and len(germ_seq) >= 3 and donor_seq[:3] == germ_seq[:3]:
            # Sequences match — this is a shifted germline residue
            # Donor position d_num corresponds to germline position d_num - 1
            # (shifted right by the insertion)
            germ_num = d_num - 1
            germ_key = f"{chain_type}{germ_num}"
            if germ_key in gmap:
                return gmap[germ_key]
    
    return None


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


def graft_chain(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str = "kabat",
    is_vhh: bool = False,
    exclude_indel: bool = False,
) -> GraftResult:
    """Build one humanized V domain (CDR grafting).
    
    exclude_indel: if True, exclude donor-only FR positions (insertions) from
                   the graft. Used for V0 (pure graft without donor indels).
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

    # Detect FR indels for shifted-residue detection
    fr_indels = detect_fr_indels(donor, v_gene)

    cdr_nums = cdr_positions_for(scheme, chain_type)
    j_anchor = 103 if chain_type == "H" else 98

    # build the grafted map
    out: Dict[str, str] = {}
    origin: Dict[str, str] = {}
    all_positions = sorted(
        set(dmap) | set(gmap) | set(jmap),
        key=lambda p: (p[0], _pos_num(p), p),
    )
    for pos in all_positions:
        num = _pos_num(pos)
        if num >= j_anchor:
            if pos in jmap and jmap[pos]:
                out[pos] = jmap[pos]
                origin[pos] = "j"
            elif pos in dmap and dmap[pos]:
                # Position in donor but not covered by the J gene — keep
                # the donor residue so the sequence stays complete (warn
                # but do NOT drop: a missing J residue truncates FR4).
                out[pos] = dmap[pos]
                origin[pos] = "donor"
                warnings.append(
                    f"[{chain_type}] FR4 position {pos} absent from J "
                    f"gene; keeping donor residue to avoid truncation")
            continue
        is_cdr = num in cdr_nums
        keep_donor = is_vhh and chain_type == "H" and num in VHH_HALLMARK
        if is_cdr:
            if pos in dmap and dmap[pos]:
                out[pos] = dmap[pos]
                origin[pos] = "donor"
        elif keep_donor and pos in dmap and dmap[pos]:
            out[pos] = dmap[pos]
            origin[pos] = "donor(vhh)"
        else:
            if pos in gmap and gmap[pos]:
                out[pos] = gmap[pos]
                origin[pos] = "germline"
            elif pos in dmap and dmap[pos]:
                # germline lacks this FR position. Two cases:
                # 1. True insertion (e.g. VH3-family H49 for VHH): keep donor
                # 2. Shifted germline residue due to upstream insertion:
                #    donor's H6A(Q) = germline's H6(Q) — use germline residue
                is_shifted = _is_shifted_germline_residue(
                    pos, dmap, gmap, fr_indels, chain_type)
                if is_shifted:
                    # This donor position corresponds to a germline position
                    # shifted by an upstream insertion — use germline residue
                    out[pos] = is_shifted
                    origin[pos] = "germline"
                elif exclude_indel:
                    # Pure graft mode (V0): skip true donor insertion positions
                    continue
                else:
                    out[pos] = dmap[pos]
                    origin[pos] = "donor(vhh)" if (is_vhh and chain_type == "H") else "donor(indel)"

    seq = "".join(aa for pos, aa in sorted(out.items(), key=lambda kv: (kv[0][0], _pos_num(kv[0]), kv[0])))
    # Build the numbered chain directly from the assembled position map.
    # Positions come from consistently-numbered donor/germline/J maps, so
    # framework lengths (e.g. VH3-family FR2 = 13 with gap at H49 vs a donor
    # that fills H49) are preserved by construction. Re-running the heuristic
    # anchor numbering on the assembled sequence would mis-derive FR2 length
    # from sequence content ("KG" -> 13) and shift CDR2 by one residue.
    numbered = chain_from_origin_map(seq, out, chain_type, donor, gmap_src, jmap_src)

    donor_positions = [pos for pos in out if origin[pos].startswith("donor")]

    # Detect FR indels between donor and germline
    fr_indels = detect_fr_indels(donor, v_gene)

    return GraftResult(
        scheme=scheme,
        chain_type=chain_type,
        sequence=seq,
        numbered=numbered,
        origin=origin,
        donor_positions=donor_positions,
        warnings=warnings,
        fr_indels=fr_indels,
    )


def graft_variant(
    donor: NumberedChain,
    v_gene: GermlineGene,
    j_gene: GermlineGene,
    scheme: str,
    backmutations: Optional[List[str]] = None,
    is_vhh: bool = False,
    force_human: Optional[List[str]] = None,
    exclude_indel: bool = False,
) -> GraftResult:
    """Graft + apply a list of back-mutations (position labels like 'H67').

    backmutations: donor residues to restore at these positions (position ->
    donor aa). force_human: positions to keep human even if recommended.
    exclude_indel: if True, exclude donor FR insertions from base graft (V0).
    """
    base = graft_chain(donor, v_gene, j_gene, scheme, is_vhh,
                       exclude_indel=exclude_indel)
    if not backmutations and not force_human:
        return base
    if v_gene.numbered is None or j_gene.numbered is None:
        raise ValueError(f"[{base.chain_type}] germline gene without numbering: {v_gene.gene_id}")
    dmap = donor.posmap()
    out = dict(base.origin)
    for pos in backmutations or []:
        if pos in dmap and dmap[pos]:
            out[pos] = "donor"
    for pos in force_human or []:
        if pos in out and out[pos] == "donor(vhh)":
            continue   # VHH hallmark must not be touched
        out[pos] = "germline"
    # rebuild sequence from origin map
    seqs: Dict[str, str] = {}
    for pos, src in out.items():
        if src in ("donor", "donor(vhh)", "donor(indel)"):
            seqs[pos] = dmap[pos]
        elif src == "germline":
            if v_gene.numbered is not None:
                seqs[pos] = v_gene.numbered.posmap()[pos]
        else:
            if j_gene.numbered is not None:
                seqs[pos] = j_gene.numbered.posmap()[pos]
    seq = "".join(aa for pos, aa in sorted(seqs.items(), key=lambda kv: (kv[0][0], _pos_num(kv[0]), kv[0])))
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
    )
