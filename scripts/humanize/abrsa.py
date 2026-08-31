"""AbRSA wrapper for antibody numbering.

Calls the AbRSA binary and parses output into NumberedChain objects.
Falls back to built-in numbering engine if AbRSA is unavailable.

Reference:
    L Li, S Chen, Z Miao, Y Liu, X Liu, ZX Xiao and Y Cao.
    AbRSA: a Robust Tool for Antibody Numbering.
    Protein Science. 28, 1524-1531, 2019
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from .numbering import NumberedChain

_ABRSA_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'AbRSA', 'ABRSA')
_ABRSA_BIN = os.path.join(_ABRSA_DIR, 'AbRSA')


def is_abrsa_available() -> bool:
    """Check if AbRSA binary is available and executable."""
    return os.path.isfile(_ABRSA_BIN) and os.access(_ABRSA_BIN, os.X_OK)


def _parse_numbering_line(line: str) -> Optional[Tuple[str, str, str]]:
    """Parse a single numbering line like 'H   1  E' or 'H  6A Q' or 'H 100A P'.
    
    Returns (chain_type, position, amino_acid) or None if not parseable.
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None
    
    parts = line.split()
    if len(parts) < 3:
        return None
    
    chain_type = parts[0]
    pos = parts[1]
    aa = parts[2]
    
    # Validate chain type
    if chain_type not in ('H', 'L'):
        return None
    
    # Validate position format (e.g., "1", "6A", "100A", "100B")
    if not re.match(r'^[HL]\d+[A-Z]*$', f"{chain_type}{pos}"):
        return None
    
    # Validate amino acid (allow any single letter A-Z, not just standard amino acids)
    # AbRSA may output non-standard characters in test sequences
    if not re.match(r'^[A-Z]$', aa.upper()):
        return None
    
    return (chain_type, pos, aa.upper())


def _pos_to_region(pos: str, chain_type: str) -> str:
    """Convert a Kabat position to a region (FR1, CDR1, etc.)."""
    # Extract numeric position and any insertion code
    m = re.match(r'^[HL](\d+)([A-Z]*)$', f"{chain_type}{pos}")
    if not m:
        return "FR1"  # fallback
    
    num = int(m.group(1))
    ins = m.group(2)
    
    if chain_type == 'H':
        # VH regions (Kabat scheme)
        if num <= 30:
            return "FR1"
        elif num <= 35 or (num == 35 and ins):
            return "CDR1"
        elif num <= 49:
            return "FR2"
        elif num <= 65 or (num == 65 and ins):
            return "CDR2"
        elif num <= 92:
            return "FR3"
        elif num <= 102 or (num == 102 and ins):
            return "CDR3"
        else:
            return "FR4"
    else:
        # VL regions (Kabat scheme)
        if num <= 23:
            return "FR1"
        elif num <= 34 or (num == 34 and ins):
            return "CDR1"
        elif num <= 49:
            return "FR2"
        elif num <= 56:
            return "CDR2"
        elif num <= 88:
            return "FR3"
        elif num <= 97:
            return "CDR3"
        else:
            return "FR4"


def _build_region_map_from_summary(cdr_regions: dict, seq: str, chain_type: str) -> dict:
    """Build seq-index → region map from AbRSA stdout CDR summary dict.
    
    The stdout summary gives correct CDR boundaries as contiguous sequence
    segments. We parse these to build an authoritative region map that
    overrides the hardcoded Kabat boundaries in _pos_to_region.
    
    cdr_regions: dict like {'H_FR1': 'EVQL...', 'H_CDR1': 'NYWLG', ...}
    Returns dict: {seq_index: 'FR1'|'CDR1'|'FR2'|'CDR2'|'FR3'|'CDR3'|'FR4'}
    """
    region_map = {}
    if not cdr_regions:
        return region_map
    
    # Ordered region keys
    prefix = chain_type
    region_order = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    
    pos = 0
    for region in region_order:
        key = f"{prefix}_{region}"
        seg = cdr_regions.get(key, '')
        if seg:
            for j in range(len(seg)):
                if pos + j < len(seq):
                    region_map[pos + j] = region
            pos += len(seg)
    
    return region_map


def _fix_vl_cdr3_numbering(residues: List[Tuple[str, str, str]]) -> List[Tuple[str, str, str]]:
    """Fix AbRSA VL CDR3 numbering bug where positions are skipped.
    
    AbRSA sometimes skips positions in VL CDR3 (e.g., L95 is skipped and L96
    is assigned to what should be L95). This function detects and corrects
    such cases by renumbering consecutive positions.
    """
    if not residues:
        return residues
    
    # Group residues by chain type
    chains = {}
    for ct, pos, aa in residues:
        if ct not in chains:
            chains[ct] = []
        chains[ct].append((ct, pos, aa))
    
    fixed = []
    for ct, chain_residues in chains.items():
        if ct != 'L':  # Only fix VL
            fixed.extend(chain_residues)
            continue
        
        # Parse positions to detect gaps
        parsed = []
        for ct2, pos, aa in chain_residues:
            m = re.match(r'^L(\d+)([A-Z]*)$', f"L{pos}")
            if m:
                num = int(m.group(1))
                ins = m.group(2)
                parsed.append((num, ins, aa, pos))
            else:
                parsed.append((0, '', aa, pos))
        
        # Check for gaps in CDR3 region (positions 89-97)
        # Find CDR3 residues
        cdr3_start = None
        cdr3_end = None
        for i, (num, ins, aa, orig_pos) in enumerate(parsed):
            if 89 <= num <= 97 and not ins:
                if cdr3_start is None:
                    cdr3_start = i
                cdr3_end = i
        
        if cdr3_start is None or cdr3_end is None:
            fixed.extend(chain_residues)
            continue
        
        # Check for gaps in CDR3
        has_gap = False
        for i in range(cdr3_start, cdr3_end):
            num1 = parsed[i][0]
            num2 = parsed[i + 1][0]
            if num2 - num1 > 1 and not parsed[i + 1][1]:  # Gap and not an insertion
                has_gap = True
                break
        
        if not has_gap:
            fixed.extend(chain_residues)
            continue
        
        # Fix the gap by renumbering
        # Find the gap position and renumber
        for i in range(cdr3_start, cdr3_end):
            num1 = parsed[i][0]
            num2 = parsed[i + 1][0]
            if num2 - num1 > 1 and not parsed[i + 1][1]:
                # Found a gap, renumber from here
                # The residue at i+1 should be num1+1
                for j in range(i + 1, cdr3_end + 1):
                    old_num, old_ins, old_aa, old_pos = parsed[j]
                    new_num = num1 + (j - i)
                    parsed[j] = (new_num, old_ins, old_aa, f"L{new_num}{old_ins}")
                break
        
        # Convert back to tuple format
        for num, ins, aa, orig_pos in parsed:
            fixed.append((ct, f"{num}{ins}" if ins else str(num), aa))
    
    return fixed


def _parse_abrsa_stdout(stdout: str) -> List[Tuple[str, str, str]]:
    """Parse AbRSA stdout (detailed numbering)."""
    residues = []
    for line in stdout.split('\n'):
        parsed = _parse_numbering_line(line)
        if parsed:
            residues.append(parsed)
    return residues


def _parse_abrsa_stderr(stderr: str) -> Dict[str, str]:
    """Parse AbRSA stderr (CDR region summary)."""
    regions = {}
    for line in stderr.split('\n'):
        line = line.strip()
        if ':' not in line:
            continue
        
        # Match patterns like "H_FR1 : EVQLLEQ..." or "H_CDR1: NYWLG"
        m = re.match(r'^([HL]_(?:FR|CDR)\d+)\s*:\s*(.+)$', line)
        if m:
            key = m.group(1)
            seq = m.group(2).strip()
            regions[key] = seq
    
    return regions


def _build_cdr_segments(residues: List[Tuple[str, str, str]], 
                        chain_type: str,
                        cdr_regions: Dict[str, str]) -> Dict[str, Tuple[str, str]]:
    """Build CDR segment dictionary from AbRSA output."""
    cdrs = {}
    
    for key, seq in cdr_regions.items():
        # key format: "H_CDR1", "L_CDR2", etc.
        if f'{chain_type}_CDR' not in key:
            continue
        
        cdr_name = key.split('_')[1]  # e.g., "CDR1"
        
        # Find start and end positions from residues
        start_pos = None
        end_pos = None
        
        for ct, pos, aa in residues:
            if ct != chain_type:
                continue
            
            region = _pos_to_region(pos, chain_type)
            if region == cdr_name:
                if start_pos is None:
                    start_pos = pos
                end_pos = pos
        
        if start_pos and end_pos:
            cdrs[cdr_name] = (f"{chain_type}{start_pos}", f"{chain_type}{end_pos}")
    
    return cdrs


def number_with_abrsa(seq: str, chain_type: str = "H",
                      scheme: str = "kabat") -> Optional[NumberedChain]:
    """Number a sequence using AbRSA.
    
    Args:
        seq: Amino acid sequence
        chain_type: "H" for heavy/VHH, "L" for light
        scheme: "kabat", "chothia", or "imgt"
    
    Returns:
        NumberedChain if successful, None if AbRSA fails or unavailable.
    """
    # Lazy import to avoid circular dependency
    from .numbering import NumberedChain, NumberedResidue  # noqa: F811
    
    if not is_abrsa_available():
        return None
    
    # AbRSA has known bugs with lambda chains (positions shifted by 1)
    # Fall back to built-in engine for lambda chains
    # Lambda chains are detected by having Cys at position 22 (not 23)
    s = seq.upper().strip()
    if chain_type == 'L':
        # Check if this is a lambda chain by looking for Cys at index 21 (L22)
        # Lambda chains have Cys at L22, kappa chains have Cys at L23
        cys_at_22 = False
        for i in range(20, min(27, len(s))):
            if s[i] == 'C':
                if i == 21:  # Lambda chain (Cys at index 21 = L22)
                    cys_at_22 = True
                break
        
        if cys_at_22:
            # This is a lambda chain, fall back to built-in engine
            return None
    
    # Create temporary files
    fasta_path = None
    output_path = None
    
    try:
        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">query\n{seq}\n")
            fasta_path = f.name
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            output_path = f.name
        
        # Map scheme to flag
        scheme_flags = {'kabat': '-k', 'chothia': '-c', 'imgt': '-g'}
        flag = scheme_flags.get(scheme, '-k')
        
        # Run AbRSA with -o flag to get detailed numbering
        result = subprocess.run(
            [_ABRSA_BIN, "-i", fasta_path, flag, "-o", output_path],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode != 0:
            return None
        
        # Read output file
        with open(output_path, 'r') as f:
            output = f.read()
        
        # Parse output file for detailed numbering
        residues_data = _parse_abrsa_stdout(output)
        
        # Parse stdout for CDR region summary
        cdr_regions = _parse_abrsa_stderr(result.stdout)
        
        if not residues_data:
            return None
        
        # Detect actual chain type from AbRSA output
        actual_chain_type = residues_data[0][0] if residues_data else chain_type
        
        # Fix VL CDR3 numbering bug if present
        if actual_chain_type == 'L':
            residues_data = _fix_vl_cdr3_numbering(residues_data)
        
        # Build region map from stdout CDR summary (authoritative source)
        region_map = _build_region_map_from_summary(cdr_regions, seq.upper().strip(), actual_chain_type)
        
        # Build NumberedResidue list using region_map
        residues = []
        for i, (ct, pos, aa) in enumerate(residues_data):
            region = region_map.get(i, _pos_to_region(pos, actual_chain_type))
            residues.append(NumberedResidue(
                pos=f"{actual_chain_type}{pos}",
                aa=aa,
                region=region,
                index=i
            ))
        
        # Build CDR segments
        cdrs = _build_cdr_segments(residues_data, actual_chain_type, cdr_regions)
        
        # Create NumberedChain
        chain = NumberedChain(
            chain_type=actual_chain_type,
            sequence=seq.upper().strip(),
            residues=residues,
            species_hint="unknown",
            warnings=[],
            cdrs=cdrs
        )
        
        return chain
        
    except (subprocess.TimeoutExpired, Exception) as e:
        return None
    finally:
        for path in (fasta_path, output_path):
            if path:
                try:
                    os.unlink(path)
                except OSError:
                    pass
