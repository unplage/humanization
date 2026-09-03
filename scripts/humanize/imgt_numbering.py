"""IMGT numbering adapter for antibody humanization evaluation.

Provides Kabat-to-IMGT position mapping and IMGT-based germline comparison
for humanization degree evaluation per WHO/INN/USAN standards.

IMGT numbering (Lefranc 2003):
  VH: FR1 1-26, CDR1 27-38, FR2 39-55, CDR2 56-65, FR3 66-104, CDR3 105-117, FR4 118-128
  VL: FR1 1-23, CDR1 24-34, FR2 35-55, CDR2 56-56, FR3 57-88, CDR3 89-104, FR4 105-118

Kabat numbering:
  VH: FR1 1-30, CDR1 31-35, FR2 36-49, CDR2 50-65, FR3 66-94, CDR3 95-102, FR4 103-113
  VL: FR1 1-23, CDR1 24-34, FR2 35-49, CDR2 50-56, FR3 57-88, CDR3 89-97, FR4 98-107

Reference: Lefranc, M.-P. (2003). IMGT unique numbering for immunoglobulin
and T cell receptor variable domains and Ig superfamily V-like domains.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

# IMGT region definitions (numbering ranges)
IMGT_VH_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}

IMGT_VL_REGIONS = {
    "FR1": (1, 23),
    "CDR1": (24, 34),
    "FR2": (35, 55),
    "CDR2": (56, 65),
    "FR3": (66, 88),
    "CDR3": (89, 104),
    "FR4": (105, 118),
}

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
IMGT_EVAL_REGIONS_VH = ["FR1", "CDR1", "FR2", "CDR2", "FR3"]
IMGT_EVAL_REGIONS_VL = ["FR1", "CDR1", "FR2", "CDR2", "FR3"]


def kabat_to_imgt_position(kabat_pos: str, chain_type: str) -> Optional[int]:
    """Convert a Kabat position string (e.g. 'H31', 'H100A') to IMGT position number.

    For positions that don't map cleanly (e.g. insertions), returns None.
    The caller should handle these cases separately.
    """
    # Extract numeric position from Kabat label
    prefix = kabat_pos[0]  # H or L
    num_str = "".join(c for c in kabat_pos[1:] if c.isdigit())
    if not num_str:
        return None
    num = int(num_str)

    # Check for insertions (letter suffix like H100A, H52A)
    has_insertion = any(c.isalpha() for c in kabat_pos[1:])

    if chain_type == "H":
        if has_insertion:
            # IMGT doesn't use insertions the same way; skip these positions
            return None

        # Kabat VH -> IMGT VH mapping (simplified, covers standard positions)
        if num <= 30:
            return num  # FR1: Kabat 1-30 -> IMGT 1-26 (roughly aligned)
        elif num <= 35:
            return num + (-4)  # CDR1: Kabat 31-35 -> IMGT 27-31
        elif num <= 49:
            return num + (-3)  # FR2: Kabat 36-49 -> IMGT 33-46
        elif num <= 65:
            return num + (-6)  # CDR2: Kabat 50-65 -> IMGT 44-59
        elif num <= 94:
            return num + (-1)  # FR3: Kabat 66-94 -> IMGT 65-93
        elif num <= 102:
            return num + (10)  # CDR3: Kabat 95-102 -> IMGT 105-112
        elif num <= 113:
            return num + (13)  # FR4: Kabat 103-113 -> IMGT 116-126
        else:
            return None
    else:  # VL
        if has_insertion:
            return None

        if num <= 23:
            return num  # FR1: Kabat 1-23 -> IMGT 1-23
        elif num <= 34:
            return num  # CDR1: Kabat 24-34 -> IMGT 24-34
        elif num <= 49:
            return num  # FR2: Kabat 35-49 -> IMGT 35-49
        elif num <= 56:
            return num  # CDR2: Kabat 50-56 -> IMGT 50-56
        elif num <= 88:
            return num  # FR3: Kabat 57-88 -> IMGT 57-88
        elif num <= 97:
            return num + (7)  # CDR3: Kabat 89-97 -> IMGT 96-104
        elif num <= 107:
            return num + (8)  # FR4: Kabat 98-107 -> IMGT 106-115
        else:
            return None


def get_imgt_region(pos: int, chain_type: str) -> str:
    """Get the IMGT region for a given IMGT position number."""
    regions = IMGT_VH_REGIONS if chain_type == "H" else IMGT_VL_REGIONS
    for region, (start, end) in regions.items():
        if start <= pos <= end:
            return region
    return "FR4"  # Default for positions beyond defined ranges


def kabat_posmap_to_imgt_posmap(posmap: Dict[str, str], chain_type: str) -> Dict[str, str]:
    """Convert a Kabat-based position map to an IMGT-based position map.

    Returns {imgt_position_label: amino_acid} dict.
    Positions that cannot be mapped are skipped.
    """
    result = {}
    for kabat_pos, aa in posmap.items():
        if kabat_pos[0] not in ("H", "L"):
            continue
        imgt_num = kabat_to_imgt_position(kabat_pos, chain_type)
        if imgt_num is not None:
            imgt_label = f"{'H' if chain_type == 'H' else 'L'}{imgt_num}"
            result[imgt_label] = aa
    return result


def compare_to_germline_imgt(
    query_posmap: Dict[str, str],
    germline_posmap: Dict[str, str],
    chain_type: str,
) -> Dict[str, float]:
    """Compare a query sequence against germline using IMGT position mapping.

    Evaluates FR identity, CDR identity, and overall identity over the
    Fv region (FR1+CDR1+FR2+CDR2+FR3), excluding CDR3 and FR4,
    per WHO/INN/USAN standards.

    Args:
        query_posmap: Query sequence as {kabat_position: amino_acid}
        germline_posmap: Germline sequence as {kabat_position: amino_acid}
        chain_type: "H" or "L"

    Returns:
        Dict with fr_identity, cdr_identity, all_identity, n_fr, n_cdr
    """
    # Convert both to IMGT posmaps
    query_imgt = kabat_posmap_to_imgt_posmap(query_posmap, chain_type)
    germline_imgt = kabat_posmap_to_imgt_posmap(germline_posmap, chain_type)

    if not query_imgt or not germline_imgt:
        return {"fr_identity": 0.0, "cdr_identity": 0.0, "all_identity": 0.0, "n_fr": 0, "n_cdr": 0}

    # Define which IMGT regions to evaluate (excluding CDR3 and FR4)
    eval_regions = set(IMGT_EVAL_REGIONS_VH if chain_type == "H" else IMGT_EVAL_REGIONS_VL)

    fr_count = 0
    fr_match = 0
    cdr_count = 0
    cdr_match = 0

    for imgt_label, aa in query_imgt.items():
        # Extract position number
        num_str = "".join(c for c in imgt_label[1:] if c.isdigit())
        if not num_str:
            continue
        imgt_num = int(num_str)

        # Get IMGT region
        region = get_imgt_region(imgt_num, chain_type)
        if region not in eval_regions:
            continue

        # Check if germline has this position
        if imgt_label not in germline_imgt:
            continue

        germline_aa = germline_imgt[imgt_label]

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

    return {
        "fr_identity": round(fr_identity, 4),
        "cdr_identity": round(cdr_identity, 4),
        "all_identity": round(all_identity, 4),
        "n_fr": fr_count,
        "n_cdr": cdr_count,
    }


def format_imgt_region_alignment(
    query_posmap: Dict[str, str],
    germline_posmap: Dict[str, str],
    chain_type: str,
    top_n: int = 5,
) -> str:
    """Format a detailed alignment showing IMGT region differences."""
    query_imgt = kabat_posmap_to_imgt_posmap(query_posmap, chain_type)
    germline_imgt = kabat_posmap_to_imgt_posmap(germline_posmap, chain_type)

    lines = []

    # Group by region
    regions_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
    regions = IMGT_VH_REGIONS if chain_type == "H" else IMGT_VL_REGIONS

    for region in regions_order:
        start, end = regions[region]
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
