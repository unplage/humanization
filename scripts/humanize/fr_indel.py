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
from typing import Dict, List, Optional, Tuple

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
class FRIndelCandidate:
    """A single candidate insertion/deletion position."""
    position: str          # Kabat position (e.g. 'H4', 'H6A')
    donor_aa: str          # amino acid at this position in donor
    confidence: float      # confidence score 0-1
    alignment_score: int   # alignment score
    context_match: float   # context match percentage (0-1)
    is_recommended: bool = False  # recommended by algorithm

    def __str__(self) -> str:
        rec = " ← 推荐" if self.is_recommended else ""
        return f"{self.position}({self.donor_aa}) 置信度={self.confidence:.2f} 上下文={self.context_match:.0%}{rec}"


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
    candidates: List[FRIndelCandidate] = field(default_factory=list)  # all candidates
    selected_candidate: Optional[str] = None  # user-selected position

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

    Returns a list of FRIndel objects, one per detected indel. Each FRIndel
    contains a list of candidates with confidence scores for interactive selection.
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
            # Donor has more residues → find all insertion candidates
            candidates = _find_all_insertion_candidates(d_aas, g_aas, d_positions, g_positions)
            if candidates:
                # Use the best candidate as default
                best = candidates[0]
                nearby = _find_nearby(best.position, set(d_positions), set(g_positions), chain_type)
                indels.append(FRIndel(
                    chain_type=chain_type,
                    indel_type="insertion",
                    fr_region=region,
                    position=best.position,
                    donor_aa=best.donor_aa,
                    germline_aa="-",
                    donor_count=d_count,
                    germline_count=g_count,
                    nearby_positions=nearby,
                    candidates=candidates,
                ))
        else:
            # Germline has more residues → find deletion point
            del_candidates = _find_all_insertion_candidates(g_aas, d_aas, g_positions, d_positions)
            if del_candidates:
                best = del_candidates[0]
                nearby = _find_nearby(best.position, set(d_positions), set(g_positions), chain_type)
                indels.append(FRIndel(
                    chain_type=chain_type,
                    indel_type="deletion",
                    fr_region=region,
                    position=best.position,
                    donor_aa="-",
                    germline_aa=best.donor_aa,
                    donor_count=d_count,
                    germline_count=g_count,
                    nearby_positions=nearby,
                    candidates=del_candidates,
                ))

    return indels


def _needleman_wunsch(seq1: List[str], seq2: List[str]) -> List[Tuple[str, str]]:
    """Needleman-Wunsch global alignment. Returns list of (char_from_seq1, char_from_seq2)
    where '-' indicates a gap. Gap penalty = -2, match = +1, mismatch = -1."""
    n, m = len(seq1), len(seq2)
    # Score matrix
    score = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = -2 * i
    for j in range(1, m + 1):
        score[0][j] = -2 * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = 1 if seq1[i - 1] == seq2[j - 1] else -1
            score[i][j] = max(score[i - 1][j] - 2,      # gap in seq2
                              score[i][j - 1] - 2,      # gap in seq1
                              score[i - 1][j - 1] + s)  # match/mismatch
    # Traceback
    alignment = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = 1 if seq1[i - 1] == seq2[j - 1] else -1
            if score[i][j] == score[i - 1][j - 1] + s:
                alignment.append((seq1[i - 1], seq2[j - 1]))
                i -= 1
                j -= 1
                continue
        if i > 0 and score[i][j] == score[i - 1][j] - 2:
            alignment.append((seq1[i - 1], '-'))
            i -= 1
        else:
            alignment.append(('-', seq2[j - 1]))
            j -= 1
    alignment.reverse()
    return alignment


def _find_insertion_point(
    long_aas: List[str], short_aas: List[str],
    long_positions: List[str], short_positions: List[str],
) -> Tuple[Optional[str], str]:
    """Find where the extra residue is in the longer sequence.

    Uses Needleman-Wunsch global alignment to find the insertion point.
    Returns (position_label, amino_acid) of the insertion, or (None, '') if
    not found.
    """
    if len(long_aas) <= len(short_aas):
        return None, ''

    alignment = _needleman_wunsch(long_aas, short_aas)

    # Find gaps in seq2 (short) = insertions in seq1 (long)
    # Collect insertion positions and their context scores
    best_pos = None
    best_aa = ''
    best_score = -1

    li = 0  # index into long_positions
    for al, as_ in alignment:
        if al != '-' and as_ != '-':
            # Both match/mismatch - count as alignment evidence
            pass
        if al != '-' and as_ == '-':
            # Gap in short = insertion in long
            pos_label = long_positions[li] if li < len(long_positions) else None
            # Score by counting matches in the alignment
            match_count = sum(1 for a, b in alignment if a == b and a != '-')
            if match_count > best_score:
                best_score = match_count
                best_pos = pos_label
                best_aa = al
        if al != '-':
            li += 1

    # Also check: gap in long = deletion from short (not needed here but
    # useful for completeness in detect_fr_indels)

    if best_pos is not None:
        # Verify alignment quality: at least 40% identity
        total = min(len(long_aas), len(short_aas))
        if best_score >= total * 0.3:
            return best_pos, best_aa

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


def select_insertion_interactive(indel: FRIndel) -> Optional[str]:
    """Interactive mode: let user select the correct insertion position.
    
    Args:
        indel: FRIndel object with candidates list
    
    Returns:
        Selected position string, or None if user skips
    """
    if not indel.candidates or len(indel.candidates) <= 1:
        return indel.position if indel.candidates else None
    
    print(f"\n{'='*60}")
    print(f"检测到 {indel.fr_region} {indel.indel_type}")
    print(f"供体 {indel.donor_count} 残基 vs 种系 {indel.germline_count} 残基")
    print(f"{'='*60}")
    print(f"\n{'#':<4} {'位置':<8} {'残基':<6} {'置信度':<10} {'上下文匹配':<12} {'备注':<8}")
    print("-" * 60)
    
    for i, c in enumerate(indel.candidates):
        rec = "← 推荐" if c.is_recommended else ""
        print(f"{i+1:<4} {c.position:<8} {c.donor_aa:<6} {c.confidence:<10.2f} {c.context_match:<12.0%} {rec}")
    
    # Get user input
    default_idx = 0  # recommended candidate
    for i, c in enumerate(indel.candidates):
        if c.is_recommended:
            default_idx = i
            break
    
    while True:
        try:
            user_input = input(f"\n请选择插入位置 (1-{len(indel.candidates)}) [默认: {default_idx + 1}]: ").strip()
            if not user_input:
                selected_idx = default_idx
                break
            selected_idx = int(user_input) - 1
            if 0 <= selected_idx < len(indel.candidates):
                break
            print(f"请输入 1-{len(indel.candidates)} 之间的数字")
        except ValueError:
            print("请输入有效的数字")
        except EOFError:
            # Non-interactive mode, use default
            selected_idx = default_idx
            break
    
    selected = indel.candidates[selected_idx]
    print(f"已选择: {selected.position} ({selected.donor_aa})")
    
    return selected.position


def update_indel_selection(indel: FRIndel, selected_position: str) -> FRIndel:
    """Update FRIndel with user-selected position.
    
    Args:
        indel: Original FRIndel object
        selected_position: User-selected position string
    
    Returns:
        Updated FRIndel with new position and selected_candidate set
    """
    # Find the candidate with matching position
    selected_candidate = None
    for c in indel.candidates:
        if c.position == selected_position:
            selected_candidate = c
            break
    
    if selected_candidate is None:
        # Position not in candidates, return original
        return indel
    
    # Create new FRIndel with updated position
    return FRIndel(
        chain_type=indel.chain_type,
        indel_type=indel.indel_type,
        fr_region=indel.fr_region,
        position=selected_candidate.position,
        donor_aa=selected_candidate.donor_aa if indel.indel_type == "insertion" else indel.donor_aa,
        germline_aa=indel.germline_aa if indel.indel_type == "insertion" else selected_candidate.donor_aa,
        donor_count=indel.donor_count,
        germline_count=indel.germline_count,
        nearby_positions=indel.nearby_positions,
        candidates=indel.candidates,
        selected_candidate=selected_candidate.position,
    )


def _score_insertion_context(
    pos_index: int,
    long_aas: List[str],
    short_aas: List[str],
    window: int = 3,
) -> float:
    """Score how well the context around an insertion point matches.
    
    This function evaluates how well the surrounding residues match when we
    assume the insertion is at pos_index. The key insight is that after the
    insertion, the downstream residues in the long sequence should match the
    downstream residues in the short sequence (shifted by 1).
    
    We give higher weight to downstream context (0.6) than upstream (0.4)
    because downstream matches are more informative about the insertion point.
    
    Args:
        pos_index: Index in long_aas where the insertion is located
        long_aas: Full amino acid list of longer sequence
        short_aas: Full amino acid list of shorter sequence
        window: Number of residues to check on each side
    
    Returns:
        Context match score (0-1)
    """
    # Check upstream context (before the insertion point)
    # Upstream residues should match directly (no shift)
    up_matches = 0
    up_total = 0
    for i in range(1, window + 1):
        long_idx = pos_index - i
        short_idx = pos_index - i
        if long_idx >= 0 and short_idx >= 0 and short_idx < len(short_aas):
            up_total += 1
            if long_aas[long_idx] == short_aas[short_idx]:
                up_matches += 1
    
    # Check downstream context (after the insertion point)
    # Downstream residues should match with a shift of 1
    down_matches = 0
    down_total = 0
    for i in range(1, window + 1):
        long_idx = pos_index + i
        short_idx = pos_index + i - 1  # -1 because insertion shifts downstream
        if long_idx < len(long_aas) and short_idx < len(short_aas):
            down_total += 1
            if long_aas[long_idx] == short_aas[short_idx]:
                down_matches += 1
    
    up_rate = up_matches / up_total if up_total > 0 else 0.0
    down_rate = down_matches / down_total if down_total > 0 else 0.0
    
    # Weighted combination: downstream is more informative
    return 0.4 * up_rate + 0.6 * down_rate


def _find_all_insertion_candidates(
    long_aas: List[str],
    short_aas: List[str],
    long_positions: List[str],
    short_positions: List[str],
) -> List[FRIndelCandidate]:
    """Find all candidate insertion positions with confidence scores.
    
    This function evaluates all possible insertion positions in the longer
    sequence, not just the gaps found by Needleman-Wunsch. For each position,
    it calculates a confidence score based on how well the context matches
    when we assume the insertion is at that position.
    
    Positions with Kabat insertion codes (e.g. H6A, H100B) receive a bonus
    because the numbering algorithm already identified them as insertions.
    
    Returns a list of FRIndelCandidate objects sorted by confidence (highest first).
    """
    if len(long_aas) <= len(short_aas):
        return []
    
    candidates = []
    
    # Evaluate each position in the longer sequence as a potential insertion point
    for i in range(len(long_aas)):
        pos_label = long_positions[i] if i < len(long_positions) else None
        if not pos_label:
            continue
        
        # Calculate context match for this position
        context_match = _score_insertion_context(i, long_aas, short_aas)
        
        # Calculate how many residues match when we assume insertion at position i
        matches = 0
        total = 0
        
        # Upstream: long[0:i] should match short[0:i]
        for j in range(i):
            if j < len(short_aas):
                total += 1
                if long_aas[j] == short_aas[j]:
                    matches += 1
        
        # Downstream: long[i+1:] should match short[i:] (shifted by 1)
        for j in range(i + 1, len(long_aas)):
            short_idx = j - 1
            if short_idx < len(short_aas):
                total += 1
                if long_aas[j] == short_aas[short_idx]:
                    matches += 1
        
        # Calculate confidence
        alignment_ratio = matches / total if total > 0 else 0
        confidence = 0.4 * alignment_ratio + 0.6 * context_match
        
        candidates.append(FRIndelCandidate(
            position=pos_label,
            donor_aa=long_aas[i],
            confidence=confidence,
            alignment_score=matches,
            context_match=context_match,
        ))
    
    # Sort by confidence (highest first)
    candidates.sort(key=lambda c: (-c.confidence, -c.context_match))
    
    # Mark the best candidate as recommended
    if candidates:
        candidates[0].is_recommended = True
    
    return candidates
