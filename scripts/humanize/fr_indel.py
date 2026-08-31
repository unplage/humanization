"""FR insertion/deletion detection between donor and germline.

Compares donor and germline FR positions (FR1/FR2/FR3) to find:
  - Donor insertions: positions in donor but not in germline
  - Donor deletions: positions in germline but not in donor

These indels represent structural differences between donor and germline
frameworks that must be handled explicitly during grafting and back-mutation
scoring.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Tuple

from .germline import GermlineGene
from .numbering import NumberedChain

# Kabat FR region boundaries (position number ranges, inclusive)
FR_REGIONS_H = {
    "FR1": (1, 30),
    "FR2": (36, 49),
    "FR3": (66, 92),
}
FR_REGIONS_L = {
    "FR1": (1, 23),
    "FR2": (35, 49),
    "FR3": (57, 88),
}

FR_REGIONS = {"H": FR_REGIONS_H, "L": FR_REGIONS_L}


def _pos_num(pos: str) -> int:
    """Extract integer from Kabat position (e.g. 'H6A' -> 6, 'H100B' -> 100)."""
    return int("".join(c for c in pos if c.isdigit()))


def _pos_ins(pos: str) -> str:
    """Extract insertion code from Kabat position (e.g. 'H6A' -> 'A')."""
    m = re.match(r'^[HL]\d+([A-Z]*)$', pos)
    return m.group(1) if m else ""


def _pos_base(pos: str) -> str:
    """Get base position without insertion (e.g. 'H6A' -> 'H6')."""
    m = re.match(r'^([HL]\d+)', pos)
    return m.group(1) if m else pos


def _region_of_num(chain_type: str, num: int) -> Optional[str]:
    """Get FR region name from position number, or None if in CDR/J."""
    for region, (lo, hi) in FR_REGIONS[chain_type].items():
        if lo <= num <= hi:
            return region
    return None


def _group_positions_by_region(
    posmap: Dict[str, str], chain_type: str
) -> Dict[str, Dict[str, str]]:
    """Group a posmap by FR region. Returns {region: {pos: aa}}."""
    regions: Dict[str, Dict[str, str]] = {}
    for pos, aa in posmap.items():
        num = _pos_num(pos)
        region = _region_of_num(chain_type, num)
        if region:
            regions.setdefault(region, {})[pos] = aa
    return regions


@dataclass
class FRIndel:
    """A detected insertion or deletion in a FR region."""
    chain_type: str          # 'H' or 'L'
    indel_type: str          # 'insertion' (donor has extra) or 'deletion' (donor missing)
    fr_region: str           # 'FR1', 'FR2', 'FR3'
    position: str            # Kabat position of the indel (e.g. 'H6A')
    donor_aa: str            # insertion: donor aa; deletion: '-'
    germline_aa: str         # insertion: '-'; deletion: germline aa
    donor_count: int         # number of residues in donor FR region
    germline_count: int      # number of residues in germline FR region
    nearby_positions: List[str] = field(default_factory=list)

    def __str__(self) -> str:
        if self.indel_type == "insertion":
            return (f"{self.chain_type} {self.fr_region}: donor insertion "
                    f"{self.position}({self.donor_aa}) — "
                    f"donor {self.donor_count} vs germline {self.germline_count} residues")
        else:
            return (f"{self.chain_type} {self.fr_region}: donor deletion "
                    f"{self.position}({self.germline_aa}) — "
                    f"donor {self.donor_count} vs germline {self.germline_count} residues")


def detect_fr_indels(
    donor: NumberedChain, v_gene: GermlineGene
) -> List[FRIndel]:
    """Detect FR insertions and deletions between donor and germline.

    Uses sequence alignment (not just position comparison) to find the
    exact insertion/deletion point. When donor FR has more residues than
    germline, we align the sequences to find where the extra residue is.

    Returns a list of FRIndel objects, one per detected indel.
    """
    if v_gene.numbered is None:
        return []

    chain_type = donor.chain_type
    dmap = donor.posmap()
    gmap = v_gene.numbered.posmap()

    d_regions = _group_positions_by_region(dmap, chain_type)
    g_regions = _group_positions_by_region(gmap, chain_type)

    indels: List[FRIndel] = []

    for region in ("FR1", "FR2", "FR3"):
        d_posset = d_regions.get(region, {})
        g_posset = g_regions.get(region, {})

        d_count = len(d_posset)
        g_count = len(g_posset)

        if d_count == g_count:
            continue

        # Get ordered position lists and AA sequences
        d_positions = sorted(d_posset.keys(), key=lambda p: (_pos_num(p), _pos_ins(p)))
        g_positions = sorted(g_posset.keys(), key=lambda p: (_pos_num(p), _pos_ins(p)))
        d_aas = [d_posset[p] for p in d_positions]
        g_aas = [g_posset[p] for p in g_positions]

        if d_count > g_count:
            # Donor has more residues → find insertion point by alignment
            ins_pos, ins_aa = _find_insertion_point(d_aas, g_aas, d_positions, g_positions)
            if ins_pos is not None:
                nearby = _find_nearby(ins_pos, set(d_positions), set(g_positions), chain_type)
                indels.append(FRIndel(
                    chain_type=chain_type,
                    indel_type="insertion",
                    fr_region=region,
                    position=ins_pos,
                    donor_aa=ins_aa,
                    germline_aa="-",
                    donor_count=d_count,
                    germline_count=g_count,
                    nearby_positions=nearby,
                ))
        else:
            # Germline has more residues → find deletion point
            del_pos, del_aa = _find_insertion_point(g_aas, d_aas, g_positions, d_positions)
            if del_pos is not None:
                nearby = _find_nearby(del_pos, set(d_positions), set(g_positions), chain_type)
                indels.append(FRIndel(
                    chain_type=chain_type,
                    indel_type="deletion",
                    fr_region=region,
                    position=del_pos,
                    donor_aa="-",
                    germline_aa=del_aa,
                    donor_count=d_count,
                    germline_count=g_count,
                    nearby_positions=nearby,
                ))

    return indels


def _find_insertion_point(
    long_aas: List[str], short_aas: List[str],
    long_positions: List[str], short_positions: List[str],
) -> Tuple[Optional[str], str]:
    """Find where the extra residue is in the longer sequence.

    Compares the two sequences and finds the position where the longer
    sequence has an extra residue (insertion).

    Returns (position_label, amino_acid) of the insertion, or (None, '') if
    not found.
    """
    if len(long_aas) <= len(short_aas):
        return None, ''

    # Try aligning by skipping one residue in the longer sequence
    # at each position, checking if the rest matches
    for skip_idx in range(len(long_aas)):
        # Build aligned sequences: long without skip_idx, vs short
        long_aligned = long_aas[:skip_idx] + long_aas[skip_idx + 1:]
        if long_aligned == short_aas:
            # Found the insertion: residue at skip_idx in long sequence
            return long_positions[skip_idx], long_aas[skip_idx]

    # Fallback: try simple left-to-right alignment
    di, gi = 0, 0
    while di < len(long_aas) and gi < len(short_aas):
        if long_aas[di] == short_aas[gi]:
            di += 1
            gi += 1
        else:
            # Mismatch: try skipping one in long
            if di + 1 < len(long_aas) and long_aas[di + 1] == short_aas[gi]:
                return long_positions[di], long_aas[di]
            elif gi + 1 < len(short_aas) and long_aas[di] == short_aas[gi + 1]:
                # Short has insertion? Shouldn't happen if long is longer
                gi += 1
            else:
                di += 1
                gi += 1

    return None, ''


def _find_nearby(
    target: str,
    source_set: set,
    other_set: set,
    chain_type: str,
    window: int = 3,
) -> List[str]:
    """Find positions near target that exist in both sets (alignment anchors)."""
    target_num = _pos_num(target)
    nearby = []
    for pos in sorted(source_set | other_set, key=lambda p: (_pos_num(p), _pos_ins(p))):
        num = _pos_num(pos)
        if abs(num - target_num) <= window and pos != target:
            nearby.append(pos)
    return nearby[:5]
