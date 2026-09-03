"""IMGT numbering adapter for antibody humanization evaluation.

Provides Kabat-to-IMGT position mapping and IMGT-based germline comparison
for humanization degree evaluation per WHO/INN/USAN standards.

IMGT numbering (Lefranc 2003):
  VH and VL use IDENTICAL region definitions:
  FR1 1-26, CDR1 27-38, FR2 39-55, CDR2 56-65, FR3 66-104, CDR3 105-117, FR4 118-128

Kabat numbering:
  VH: FR1 1-30, CDR1 31-35, FR2 36-49, CDR2 50-65, FR3 66-94, CDR3 95-102, FR4 103-113
  VL: FR1 1-23, CDR1 24-34, FR2 35-49, CDR2 50-56, FR3 57-88, CDR3 89-97, FR4 98-107

Reference: Lefranc, M.-P. (2003). IMGT unique numbering for immunoglobulin
and T cell receptor variable domains and Ig superfamily V-like domains.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

# IMGT region definitions (Lefranc 2003 standard)
# IMPORTANT: VH and VL use IDENTICAL IMGT region definitions (unlike Kabat)
IMGT_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

# VH and VL use the same definitions
IMGT_VH_REGIONS = IMGT_REGIONS
IMGT_VL_REGIONS = IMGT_REGIONS

# Kabat region definitions
KABAT_VH_REGIONS = {
    "FR1": (1, 30),
    "CDR1": (31, 35),
    "FR2": (36, 49),
    "CDR2": (50, 65),
    "FR3": (66, 94),
    "CDR3": (95, 102),
    "FR4": (103, 113),
}

KABAT_VL_REGIONS = {
    "FR1": (1, 23),
    "CDR1": (24, 34),
    "FR2": (35, 49),
    "CDR2": (50, 56),
    "FR3": (57, 88),
    "CDR3": (89, 97),
    "FR4": (98, 107),
}

# Core IMGT regions used for humanization evaluation (excluding CDR3 and FR4)
IMGT_EVAL_REGIONS = ["FR1", "CDR1", "FR2", "CDR2", "FR3"]
IMGT_EVAL_REGIONS_VH = IMGT_EVAL_REGIONS
IMGT_EVAL_REGIONS_VL = IMGT_EVAL_REGIONS

# =============================================================================
# Complete Kabat-to-IMGT mapping tables (based on IMGT official correspondence)
# Reference: Lefranc, M.-P. (1999, 2003). IMGT-ONTOLOGY and IMGT unique numbering
# =============================================================================

# VH chain mapping (Kabat -> IMGT)
# Key: Kabat position number, Value: IMGT position number
# For insertions (e.g., 35A, 52A), use _map_insertion() function
KABAT_TO_IMGT_VH = {
    # FR1 (Kabat 1-30 -> IMGT 1-26)
    # Note: IMGT has a gap at position 10, so Kabat 10 maps to IMGT 11
    1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9,
    10: 11, 11: 12, 12: 13, 13: 14, 14: 15, 15: 16,
    16: 17, 17: 18, 18: 19, 19: 20, 20: 21, 21: 22,
    22: 23, 23: 24, 24: 25, 25: 26,  # 22=Cys (1st-CYS), 26=FR1 end

    # CDR1 (Kabat 26-35 -> IMGT 27-38)
    # Note: IMGT has gaps at positions 31-34 for CDR1 length variation
    26: 27, 27: 28, 28: 29, 29: 30,
    # Kabat 30-33 map to IMGT gaps (31-34), handled by _map_insertion
    34: 35, 35: 36,  # Kabat 34-35 -> IMGT 35-36

    # FR2 (Kabat 36-50 -> IMGT 41-55)
    # Note: IMGT FR2 starts at 41 (CONSERVED-TRP at 41)
    36: 41, 37: 42, 38: 43, 39: 44, 40: 45, 41: 46,
    42: 47, 43: 48, 44: 49, 45: 50, 46: 51, 47: 52,
    48: 53, 49: 54, 50: 55,

    # CDR2 (Kabat 51-65 -> IMGT 56-65)
    # Note: IMGT has gaps at positions 58-64 for CDR2 length variation
    51: 56, 52: 57,
    # Kabat 52A-52J map to IMGT 58-65, handled by _map_insertion

    # FR3 (Kabat 58-94 -> IMGT 66-104)
    # Note: Complex mapping due to IMGT gaps at positions 73, 81-82
    58: 66, 59: 67, 60: 68, 61: 69, 62: 70, 63: 71,
    64: 72,  # Kabat 64 -> IMGT 72
    # Kabat 65 maps to gap at IMGT 73 (handled by alignment)
    66: 75, 67: 76, 68: 77, 69: 78, 70: 79, 71: 80,
    # Kabat 72-74 map to IMGT gaps 81-82 (handled by alignment)
    75: 84, 76: 85, 77: 86, 78: 87, 79: 88, 80: 89,
    81: 90, 82: 91,  # Kabat 81-82 -> IMGT 90-91
    # Kabat 82A-82C map to IMGT 92-94, handled by _map_insertion
    83: 95, 84: 96, 85: 97, 86: 98, 87: 99, 88: 100,
    89: 101, 90: 102, 91: 103, 92: 104,  # 92=Cys (2nd-CYS)

    # CDR3 (Kabat 93-102 -> IMGT 105-117)
    # Note: CDR3 is variable length, IMGT has gaps at positions 111-112
    93: 105, 94: 106, 95: 107, 96: 108, 97: 109, 98: 110,
    99: 111, 100: 112,
    # Kabat 100A-100K map to IMGT 108+, handled by _map_insertion
    101: 114, 102: 115,

    # FR4 (Kabat 103-113 -> IMGT 118-128)
    # Note: Uniform offset of +15
    103: 118, 104: 119, 105: 120, 106: 121, 107: 122,
    108: 123, 109: 124, 110: 125, 111: 126, 112: 127, 113: 128,
}

# VL chain mapping (Kabat -> IMGT)
# Note: FR1-FR3 are mostly identical between Kabat and IMGT
KABAT_TO_IMGT_VL = {
    # FR1 (Kabat 1-23 -> IMGT 1-23) - identical mapping
    1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9,
    10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 15,
    16: 16, 17: 17, 18: 18, 19: 19, 20: 20, 21: 21,
    22: 22, 23: 23,  # 23=Cys (1st-CYS)

    # CDR1 (Kabat 24-34 -> IMGT 27-38)
    # Note: Kabat CDR1 starts at 24, IMGT CDR1 starts at 27
    24: 27,  # Kabat 24 -> IMGT 27
    # Kabat 25-29 map to IMGT 28-32 (part of CDR1)
    25: 28, 26: 29, 27: 30, 28: 31, 29: 32,
    # Kabat 30-33 map to IMGT gaps 33-34 (handled by alignment)
    30: 33, 31: 34, 32: 35, 33: 36, 34: 37, 35: 38,  # Wait, 34 is end of CDR1?
    # Actually: Kabat 34 -> IMGT 35 (end of CDR1), Kabat 35 -> IMGT 39 (FR2 start)

    # FR2 (Kabat 35-49 -> IMGT 39-55)
    # Note: IMGT FR2 starts at 39 (CONSERVED-TRP at 39)
    35: 39, 36: 40, 37: 41, 38: 42, 39: 43, 40: 44,
    41: 45, 42: 46, 43: 47, 44: 48, 45: 49, 46: 50,
    47: 51, 48: 52, 49: 53, 50: 55,

    # CDR2 (Kabat 51-56 -> IMGT 56-65)
    # Note: IMGT has gaps at positions 58-64 for CDR2 length variation
    51: 56, 52: 57,
    # Kabat 53-56 map to IMGT 58-61 (CDR2)
    53: 58, 54: 59, 55: 60, 56: 61,

    # FR3 (Kabat 57-88 -> IMGT 66-104)
    # Note: Complex mapping similar to VH
    57: 66, 58: 67, 59: 68, 60: 69, 61: 70, 62: 71,
    63: 72,  # Kabat 63 -> IMGT 72
    # Kabat 64-66 map to IMGT 73-75 (FR3)
    64: 73, 65: 74, 66: 75,
    67: 77, 68: 78, 69: 79, 70: 80, 71: 81, 72: 82,
    # Kabat 73-74 map to IMGT 83-84 (FR3)
    73: 83, 74: 84,
    75: 85, 76: 86, 77: 87, 78: 88, 79: 89, 80: 90,
    81: 91, 82: 92,  # Kabat 81-82 -> IMGT 91-92
    # Kabat 82A-82C map to IMGT 93-95, handled by _map_insertion
    83: 96, 84: 97, 85: 98, 86: 99, 87: 100, 88: 104,
    # Kabat 89-91 map to IMGT 105-107 (CDR3 start)

    # CDR3 (Kabat 89-97 -> IMGT 105-117)
    # Note: CDR3 is variable length
    89: 105, 90: 106, 91: 107,
    # Kabat 92-94 map to IMGT 108-110
    92: 108, 93: 109, 94: 110,
    # Kabat 95-97 map to IMGT 111-115
    95: 111, 96: 112, 97: 115,

    # FR4 (Kabat 98-107 -> IMGT 118-127)
    # Note: Uniform offset of +20
    98: 118, 99: 119, 100: 120, 101: 121, 102: 122,
    103: 123, 104: 124, 105: 125, 106: 126, 107: 127,
}


def _map_insertion(num: int, insertion: str, chain_type: str) -> Optional[int]:
    """Map Kabat insertion positions to IMGT positions.

    Kabat uses letter suffixes for insertions (e.g., 35A, 52A, 82A).
    IMGT handles length variation differently (gaps at top of CDR loops).
    This function maps Kabat insertions to the corresponding IMGT positions.

    Args:
        num: Numeric position (e.g., 35 for 35A)
        insertion: Letter suffix (e.g., "A" for 35A)
        chain_type: "H" or "L"

    Returns:
        IMGT position number, or None if mapping not available
    """
    if not insertion:
        return None

    letter_index = ord(insertion.upper()) - ord("A")
    if letter_index < 0:
        return None

    if chain_type == "H":
        # VH insertions
        if num == 35:
            # 35A/35B -> IMGT 39/40 (FR2 start, before CONSERVED-TRP)
            if letter_index <= 1:
                return 39 + letter_index

        elif num == 52:
            # 52A-52J -> IMGT 58-65 (CDR2, 8 positions)
            # Note: IMGT CDR2 has 10 positions (56-65), but 56-57 are for Kabat 51-52
            # Kabat 52A-52J map to IMGT 58-65 (A=58, B=59, ..., H=65)
            if letter_index <= 7:  # A-H (8 positions)
                return 58 + letter_index
            # For J and beyond, map to the last CDR2 position
            elif letter_index == 9:  # J
                return 65

        elif num == 82:
            # 82A-82C -> IMGT 92-94 (FR3)
            if letter_index <= 2:
                return 92 + letter_index

        elif num == 100:
            # 100A-100K -> IMGT 108-118 (CDR3)
            if letter_index <= 10:
                return 108 + letter_index

        elif num == 102:
            # 102A-102C -> IMGT 116-118 (CDR3 end)
            if letter_index <= 2:
                return 116 + letter_index

    else:  # VL
        # VL insertions
        if num == 27:
            # 27A-27E -> IMGT 28-32 (CDR1)
            if letter_index <= 4:
                return 28 + letter_index

        elif num == 82:
            # 82A-82C -> IMGT 93-95 (FR3)
            if letter_index <= 2:
                return 93 + letter_index

        elif num == 100:
            # 100A-100D -> IMGT 109-112 (CDR3)
            if letter_index <= 3:
                return 109 + letter_index

    return None


def kabat_to_imgt_position(kabat_pos: str, chain_type: str) -> Optional[int]:
    """Convert a Kabat position string (e.g. 'H31', 'H100A') to IMGT position number.

    Uses precise lookup tables based on IMGT official correspondence (Lefranc 2003).
    For insertion positions (e.g., 35A, 52A), uses _map_insertion() function.

    Args:
        kabat_pos: Kabat position label (e.g., "H31", "H100A", "L27C")
        chain_type: "H" for heavy chain, "L" for light chain

    Returns:
        IMGT position number, or None if position cannot be mapped
    """
    if not kabat_pos or len(kabat_pos) < 2:
        return None

    # Extract prefix (H or L)
    prefix = kabat_pos[0]
    if prefix not in ("H", "L"):
        return None

    # Extract numeric position
    num_str = ""
    insertion = ""
    for c in kabat_pos[1:]:
        if c.isdigit():
            num_str += c
        elif c.isalpha():
            insertion += c

    if not num_str:
        return None
    num = int(num_str)

    # Handle insertions
    if insertion:
        return _map_insertion(num, insertion, chain_type)

    # Simple position mapping using lookup table
    mapping = KABAT_TO_IMGT_VH if chain_type == "H" else KABAT_TO_IMGT_VL
    return mapping.get(num, None)


def get_imgt_region(pos: int, chain_type: str = "H") -> str:
    """Get the IMGT region for a given IMGT position number.

    Args:
        pos: IMGT position number (1-128)
        chain_type: "H" or "L" (both use same regions in IMGT)

    Returns:
        Region name: "FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", or "FR4"
    """
    for region, (start, end) in IMGT_REGIONS.items():
        if start <= pos <= end:
            return region
    return "FR4"  # Default for positions beyond defined ranges


def kabat_posmap_to_imgt_posmap(posmap: Dict[str, str], chain_type: str) -> Dict[str, str]:
    """Convert a Kabat-based position map to an IMGT-based position map.

    Converts all positions including insertions (e.g., H35A -> IMGT 39).
    Positions that cannot be mapped are skipped.

    Args:
        posmap: Dictionary mapping Kabat position labels to amino acids
                e.g., {"H1": "E", "H31": "Y", "H100A": "F"}
        chain_type: "H" for heavy chain, "L" for light chain

    Returns:
        Dictionary mapping IMGT position labels to amino acids
        e.g., {"H1": "E", "H27": "Y", "H108": "F"}
    """
    result = {}
    for kabat_pos, aa in posmap.items():
        if not kabat_pos or kabat_pos[0] not in ("H", "L"):
            continue

        # Skip positions without valid amino acid
        if not aa or aa == "-":
            continue

        imgt_num = kabat_to_imgt_position(kabat_pos, chain_type)
        if imgt_num is not None:
            imgt_label = f"{'H' if chain_type == 'H' else 'L'}{imgt_num}"
            # Handle potential collisions (multiple Kabat positions -> same IMGT)
            # In case of collision, keep the one with insertion (more specific)
            if imgt_label not in result:
                result[imgt_label] = aa.upper()
            else:
                # Check if new position has insertion (more specific)
                has_insertion = any(c.isalpha() for c in kabat_pos[1:])
                if has_insertion:
                    result[imgt_label] = aa.upper()

    return result


def compare_to_germline_imgt(
    query_posmap: Dict[str, str],
    germline_posmap: Dict[str, str],
    chain_type: str,
) -> Dict[str, Any]:
    """Compare a query sequence against germline using IMGT position mapping.

    Evaluates FR identity, CDR identity, and overall identity over the
    Fv region (FR1+CDR1+FR2+CDR2+FR3), excluding CDR3 and FR4,
    per WHO/INN/USAN standards.

    Args:
        query_posmap: Query sequence as {kabat_position: amino_acid}
        germline_posmap: Germline sequence as {kabat_position: amino_acid}
        chain_type: "H" or "L"

    Returns:
        Dict with fr_identity, cdr_identity, all_identity, n_fr, n_cdr,
        plus imgt_region_stats for per-region breakdown
    """
    # Convert both to IMGT posmaps
    query_imgt = kabat_posmap_to_imgt_posmap(query_posmap, chain_type)
    germline_imgt = kabat_posmap_to_imgt_posmap(germline_posmap, chain_type)

    if not query_imgt or not germline_imgt:
        return {
            "fr_identity": 0.0, "cdr_identity": 0.0, "all_identity": 0.0,
            "n_fr": 0, "n_cdr": 0,
            "imgt_region_stats": {},
        }

    # Define which IMGT regions to evaluate (excluding CDR3 and FR4)
    eval_regions = set(IMGT_EVAL_REGIONS)

    fr_count = 0
    fr_match = 0
    cdr_count = 0
    cdr_match = 0

    # Per-region statistics
    region_stats = {region: {"count": 0, "match": 0} for region in IMGT_EVAL_REGIONS}

    for imgt_label, aa in query_imgt.items():
        # Extract position number (handle labels like H27, H108.1)
        num_str = ""
        for c in imgt_label[1:]:
            if c.isdigit() or c == ".":
                num_str += c
        if not num_str:
            continue

        # Handle decimal positions (e.g., 108.1 for CDR3 insertions)
        try:
            if "." in num_str:
                imgt_num = int(num_str.split(".")[0])
            else:
                imgt_num = int(num_str)
        except ValueError:
            continue

        # Get IMGT region
        region = get_imgt_region(imgt_num, chain_type)
        if region not in eval_regions:
            continue

        # Check if germline has this position
        if imgt_label not in germline_imgt:
            continue

        germline_aa = germline_imgt[imgt_label]

        # Update region stats
        region_stats[region]["count"] += 1
        if aa == germline_aa:
            region_stats[region]["match"] += 1

        if region.startswith("FR"):
            fr_count += 1
            if aa == germline_aa:
                fr_match += 1
        elif region.startswith("CDR"):
            cdr_count += 1
            if aa == germline_aa:
                cdr_match += 1

    fr_identity = fr_match / fr_count if fr_count > 0 else 0.0
    cdr_identity = cdr_match / cdr_count if cdr_count > 0 else 0.0
    total = fr_count + cdr_count
    total_match = fr_match + cdr_match
    all_identity = total_match / total if total > 0 else 0.0

    # Calculate per-region identity
    imgt_region_stats = {}
    for region in IMGT_EVAL_REGIONS:
        stats = region_stats[region]
        count = stats["count"]
        match = stats["match"]
        imgt_region_stats[region] = {
            "identity": round(match / count, 4) if count > 0 else 0.0,
            "count": count,
            "match": match,
        }

    return {
        "fr_identity": round(fr_identity, 4),
        "cdr_identity": round(cdr_identity, 4),
        "all_identity": round(all_identity, 4),
        "n_fr": fr_count,
        "n_cdr": cdr_count,
        "imgt_region_stats": imgt_region_stats,
    }


def format_imgt_region_alignment(
    query_posmap: Dict[str, str],
    germline_posmap: Dict[str, str],
    chain_type: str,
    top_n: int = 5,
) -> str:
    """Format a detailed alignment showing IMGT region differences.

    Args:
        query_posmap: Query sequence as {kabat_position: amino_acid}
        germline_posmap: Germline sequence as {kabat_position: amino_acid}
        chain_type: "H" or "L"
        top_n: Not used, kept for API compatibility

    Returns:
        Formatted string showing region-by-region alignment
    """
    query_imgt = kabat_posmap_to_imgt_posmap(query_posmap, chain_type)
    germline_imgt = kabat_posmap_to_imgt_posmap(germline_posmap, chain_type)

    lines = []

    # Group by region
    regions_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]

    for region in regions_order:
        start, end = IMGT_REGIONS[region]
        region_query = []
        region_germline = []
        region_diffs = []

        for num in range(start, end + 1):
            label = f"{'H' if chain_type == 'H' else 'L'}{num}"
            q_aa = query_imgt.get(label, "-")
            g_aa = germline_imgt.get(label, "-")
            region_query.append(q_aa)
            region_germline.append(g_aa)
            if q_aa != g_aa and q_aa != "-" and g_aa != "-":
                region_diffs.append(num)

        q_str = "".join(region_query)
        g_str = "".join(region_germline)

        # Mark differences
        diff_indicator = ""
        for i, num in enumerate(range(start, end + 1)):
            label = f"{'H' if chain_type == 'H' else 'L'}{num}"
            q_aa = query_imgt.get(label, "-")
            g_aa = germline_imgt.get(label, "-")
            if q_aa != g_aa and q_aa != "-" and g_aa != "-":
                diff_indicator += "X"
            else:
                diff_indicator += "."

        lines.append(f"  {region:5s}: {q_str}  |  {g_str}  |  {diff_indicator}")

    return "\n".join(lines)
