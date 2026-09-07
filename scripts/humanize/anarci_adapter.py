"""ANARCI adapter for IMGT numbering.

Wraps the ANARCI library to provide IMGT-numbered antibody sequences
in the same NumberedChain format used by the pipeline.

ANARCI (Antibody Numbering and Antigen Receptor ClassIfication) uses
HMMER profile searches for robust IMGT numbering. It is more reliable
than AbRSA for IMGT numbering, especially for VL chains.

Reference:
    Dunbar J, Deane CM. ANARCI: antigen receptor numbering and annotation
    tool. Nucleic Acids Res. 2016;44(W1):W27-W32.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from .numbering import NumberedChain

# IMGT region definitions for assignment
_IMGT_REGIONS = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}


def is_anarci_available() -> bool:
    """Check if ANARCI is importable."""
    try:
        from anarci import anarci  # noqa: F401
        return True
    except ImportError:
        return False


def _get_imgt_region(pos_num: int) -> str:
    """Get IMGT region from position number."""
    for region, (start, end) in _IMGT_REGIONS.items():
        if start <= pos_num <= end:
            return region
    return "FR4"


def _detect_lambda_chain(sequence: str) -> bool:
    """Detect if a VL sequence is a lambda chain.

    Note: ANARCI correctly handles lambda chains natively. This function
    is kept for compatibility but returns False since ANARCI doesn't need
    special lambda handling (unlike AbRSA).
    """
    return False


def _anarci_to_numbered_chain(result_anarci, chain_type: str, sequence: str):
    """Convert ANARCI output to NumberedChain for any scheme (Kabat/IMGT).

    Args:
        result_anarci: Raw ANARCI output tuple
        chain_type: "H" or "L"
        sequence: Original amino acid sequence

    Returns:
        NumberedChain or None if conversion fails
    """
    from .numbering import NumberedChain, NumberedResidue  # noqa: F811

    try:
        numbering_list = result_anarci[0][0][0][0]
    except (IndexError, TypeError):
        return None

    if not numbering_list:
        return None

    seq = sequence.upper().strip()
    prefix = chain_type

    residues = []
    seq_idx = 0
    for (pos_num, ins_flag), aa in numbering_list:
        if aa == '-':
            continue

        # Build position label
        if ins_flag and ins_flag != ' ':
            pos_label = f"{prefix}{pos_num}{ins_flag}"
        else:
            pos_label = f"{prefix}{pos_num}"

        # Determine region from position number (works for both Kabat and IMGT)
        region = _get_region_from_position(pos_num, chain_type)

        residues.append(NumberedResidue(
            pos=pos_label,
            aa=aa.upper(),
            region=region,
            index=seq_idx,
        ))
        seq_idx += 1

    if not residues:
        return None

    return NumberedChain(
        chain_type=chain_type,
        sequence=seq,
        residues=residues,
        species_hint="unknown",
        warnings=[],
        cdrs={},
    )


def _get_region_from_position(pos_num: int, chain_type: str) -> str:
    """Get region from position number (Kabat convention)."""
    if chain_type == 'H':
        if pos_num <= 30:
            return "FR1"
        elif pos_num <= 35:
            return "CDR1"
        elif pos_num <= 49:
            return "FR2"
        elif pos_num <= 65:
            return "CDR2"
        elif pos_num <= 94:
            return "FR3"
        elif pos_num <= 102:
            return "CDR3"
        else:
            return "FR4"
    else:  # VL
        if pos_num <= 23:
            return "FR1"
        elif pos_num <= 34:
            return "CDR1"
        elif pos_num <= 49:
            return "FR2"
        elif pos_num <= 56:
            return "CDR2"
        elif pos_num <= 88:
            return "FR3"
        elif pos_num <= 97:
            return "CDR3"
        else:
            return "FR4"


def number_with_anarci_imgt(sequence: str, chain_type: str) -> Optional[NumberedChain]:
    """Number a sequence using ANARCI with IMGT scheme.

    Args:
        sequence: Amino acid sequence
        chain_type: "H" for heavy/VHH, "L" for light

    Returns:
        NumberedChain with IMGT numbering, or None if ANARCI fails.
    """
    from .numbering import NumberedChain, NumberedResidue  # noqa: F811

    if not is_anarci_available():
        return None

    try:
        from anarci import anarci
    except ImportError:
        return None

    seq = sequence.upper().strip()

    # Lambda chain detection: ANARCI misclassifies lambda as kappa
    # We detect lambda and adjust chain_type accordingly
    is_lambda = False
    if chain_type == 'L' and _detect_lambda_chain(seq):
        is_lambda = True

    # Run ANARCI
    try:
        result = anarci(
            [('query', seq)],
            scheme='imgt',
            output=False,
        )
    except Exception:
        return None

    # Parse result structure: (numbering, alignments, germline)
    # numbering[0] = per-chain result
    # numbering[0][0] = single chain result
    # numbering[0][0][0] = (numbering_list, confidence, chain_length)
    try:
        chain_result = result[0][0][0]
        numbering_list = chain_result[0]
    except (IndexError, TypeError):
        return None

    if not numbering_list:
        return None

    # Determine actual chain type from ANARCI alignment
    actual_chain_type = chain_type
    if result[1] and result[1][0]:
        for align in result[1][0]:
            if isinstance(align, dict):
                detected = align.get('chain_type', '')
                if detected in ('H', 'K', 'L'):
                    # K = kappa, L = lambda, but ANARCI may misclassify
                    if is_lambda:
                        actual_chain_type = 'L'
                    elif detected == 'H':
                        actual_chain_type = 'H'
                    else:
                        actual_chain_type = 'L'
                    break

    # Convert ANARCI tuples to NumberedResidue objects
    residues = []
    seq_idx = 0
    for (pos_num, ins_flag), aa in numbering_list:
        if aa == '-':
            continue  # Skip gaps

        # Build position label: "L56" or "H100A"
        prefix = actual_chain_type
        if ins_flag and ins_flag != ' ':
            pos_label = f"{prefix}{pos_num}{ins_flag}"
        else:
            pos_label = f"{prefix}{pos_num}"

        # Get IMGT region
        region = _get_imgt_region(pos_num)

        residues.append(NumberedResidue(
            pos=pos_label,
            aa=aa.upper(),
            region=region,
            index=seq_idx,
        ))
        seq_idx += 1

    if not residues:
        return None

    # Build CDR segments from ANARCI numbering
    cdrs = _build_cdrs_from_numbering(numbering_list, actual_chain_type)

    return NumberedChain(
        chain_type=actual_chain_type,
        sequence=seq,
        residues=residues,
        species_hint="unknown",
        warnings=[],
        cdrs=cdrs,
    )


def _build_cdrs_from_numbering(numbering_list, chain_type: str):
    """Derive CDR boundaries from ANARCI IMGT numbering.

    Finds first and last non-gap positions in each CDR region.
    """
    cdrs = {}
    prefix = chain_type

    for cdr_name, (start, end) in [("CDR1", (27, 38)), ("CDR2", (56, 65)), ("CDR3", (105, 117))]:
        first_pos = None
        last_pos = None
        for (pos_num, ins_flag), aa in numbering_list:
            if start <= pos_num <= end and aa != '-':
                if first_pos is None:
                    first_pos = pos_num
                last_pos = pos_num

        if first_pos is not None and last_pos is not None:
            # Use last insertion letter for end position if present
            end_label = f"{prefix}{last_pos}"
            for (pos_num, ins_flag), aa in numbering_list:
                if pos_num == last_pos and ins_flag and ins_flag != ' ':
                    end_label = f"{prefix}{last_pos}{ins_flag}"
                    break

            start_label = f"{prefix}{first_pos}"
            cdrs[cdr_name] = (start_label, end_label)

    return cdrs
