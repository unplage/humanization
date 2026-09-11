#!/usr/bin/env python3
"""TAP (Therapeutic Antibody Profiler) Analysis Tool

Calculates 5 developability metrics based on TAP guidelines:
1. Total CDR Length (Ltot) - Sum of all 6 CDR lengths
2. Patches of Surface Hydrophobicity (PSH) - Hydrophobic patches on CDR vicinity
3. Patches of Positive Charge (PPC) - Positive charge patches  
4. Patches of Negative Charge (PNC) - Negative charge patches
5. Structural Fv Charge Symmetry Parameter (SFvCSP) - Charge balance VH vs VL

Usage:
    python3 tools/tap_profiler/tap_analyzer.py --pdb <structure.pdb>
    python3 tools/tap_profiler/tap_analyzer.py --vh <seq> --vl <seq>
    python3 tools/tap_profiler/tap_analyzer.py --pdb-dir <af3_dir>

Requirements:
    - freesasa (pip install freesasa)
    - python-docx (optional, for Word reports)
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import freesasa
except ImportError:
    print("ERROR: freesasa not installed. Run: pip install freesasa", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# TAP Thresholds (from TAP2 2025 update)
# ---------------------------------------------------------------------------

TAP_THRESHOLDS = {
    "Ltot": {
        "green": (43, 54),  # Acceptable range
        "amber": [(37, 42), (55, 65)],  # Caution zones
        "red": [(0, 36), (66, 200)],  # Danger zones
    },
    "PSH": {
        "green": (111.41, 167.63),  # Acceptable range
        "amber": [(95.77, 111.40), (167.64, 211.65)],  # Caution zones
        "red": [(0, 95.76), (211.66, 1000)],  # Danger zones
    },
    "PPC": {
        "green": (0, 1.33),  # Acceptable range (low is good)
        "amber": [(1.34, 4.20)],  # Caution zone
        "red": [(4.21, 100)],  # Danger zone
    },
    "PNC": {
        "green": (0, 1.98),  # Acceptable range (low is good)
        "amber": [(1.99, 4.43)],  # Caution zone
        "red": [(4.44, 100)],  # Danger zone
    },
    "SFvCSP": {
        "green": (-5.99, 100),  # Acceptable range (less negative is good)
        "amber": [(-30.60, -6.00)],  # Caution zone
        "red": [(-100, -30.61)],  # Danger zone (very negative)
    },
}

# Amino acid classifications for TAP metrics
# Charged residues (positive: K, R, H; negative: D, E)
# Hydrophobic residues: A, V, L, I, M, F, W, Y
CHARGED_POS = {'K': 1, 'R': 1, 'H': 0.5}  # His is partially charged at physiological pH
CHARGED_NEG = {'D': -1, 'E': -1}
HYDROPHOBIC = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'Y'}

# Kabat CDR definitions (H and L chains)
KABAT_CDRS = {
    "H": {
        "CDR1": (31, 35),   # H31-H35
        "CDR2": (50, 65),   # H50-H65
        "CDR3": (95, 102),  # H95-H102
    },
    "L": {
        "CDR1": (24, 34),   # L24-L34
        "CDR2": (50, 56),   # L50-L56
        "CDR3": (89, 97),   # L89-L97
    },
}


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class TAPMetrics:
    """Container for TAP metric values."""
    Ltot: int  # Total CDR length
    PSH: float  # Patches of Surface Hydrophobicity
    PPC: float  # Patches of Positive Charge
    PNC: float  # Patches of Negative Charge
    SFvCSP: float  # Structural Fv Charge Symmetry Parameter
    
    # Individual CDR lengths
    cdr1_h_len: int = 0
    cdr2_h_len: int = 0
    cdr3_h_len: int = 0
    cdr1_l_len: int = 0
    cdr2_l_len: int = 0
    cdr3_l_len: int = 0
    
    # Flags
    Ltot_flag: str = "green"
    PSH_flag: str = "green"
    PPC_flag: str = "green"
    PNC_flag: str = "green"
    SFvCSP_flag: str = "green"


@dataclass
class TAPResult:
    """Complete TAP analysis result."""
    name: str
    vh_sequence: str
    vl_sequence: str
    metrics: TAPMetrics
    flags: Dict[str, str]
    risk_summary: str  # Overall risk assessment
    details: Dict


# ---------------------------------------------------------------------------
# PDB Parsing and SASA Calculation
# ---------------------------------------------------------------------------

def parse_pdb(pdb_path: str) -> Dict[str, Dict[int, List[dict]]]:
    """Parse PDB file and return atoms grouped by chain and residue number."""
    atoms = defaultdict(lambda: defaultdict(list))
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    chain = line[21].strip() or ' '
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms[chain][res_num].append({
                        'name': atom_name,
                        'resname': res_name,
                        'resnum': res_num,
                        'x': x, 'y': y, 'z': z
                    })
                except (ValueError, IndexError):
                    continue
    return dict(atoms)


def get_chain_sequence(atoms: Dict[str, Dict[int, List[dict]]]) -> Dict[str, str]:
    """Extract amino acid sequence from parsed atoms."""
    aa_3to1 = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    }
    
    sequences = {}
    for chain, residues in atoms.items():
        seq = []
        for res_num in sorted(residues.keys()):
            res_name = residues[res_num][0]['resname']
            if res_name in aa_3to1:
                seq.append(aa_3to1[res_name])
            else:
                seq.append('X')  # Unknown residue
        sequences[chain] = ''.join(seq)
    return sequences


def calculate_sasa(pdb_path: str) -> Dict[str, Dict[int, float]]:
    """Calculate per-residue SASA using FreeSASA."""
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    
    sasa_by_residue = defaultdict(lambda: defaultdict(float))
    n = structure.nAtoms()
    
    for i in range(n):
        chain = structure.chainLabel(i)
        res_num = int(str(structure.residueNumber(i)).strip())
        area = result.atomArea(i)
        sasa_by_residue[chain][res_num] += area
    
    return dict(sasa_by_residue)


# ---------------------------------------------------------------------------
# Sequence Analysis (No Structure Required)
# ---------------------------------------------------------------------------

def calculate_cdr_length_from_sequence(sequence: str, chain_type: str) -> Dict[str, int]:
    """Calculate CDR lengths from sequence using Kabat definition.
    
    Note: This is an approximation. For accurate CDR lengths,
    structural numbering is required.
    """
    # For sequence-only analysis, we use typical CDR length ranges
    # This is less accurate than structure-based analysis
    cdr_lengths = {}
    
    if chain_type == "H":
        # VH is approximately 120 residues
        # CDR lengths vary, but typical ranges:
        # CDR1: 5-8, CDR2: 16-19, CDR3: 6-20+
        seq_len = len(sequence)
        if seq_len >= 110:
            # Estimate CDR lengths based on position
            cdr_lengths["CDR1"] = 6  # Typical
            cdr_lengths["CDR2"] = 17  # Typical
            cdr_lengths["CDR3"] = min(max(seq_len - 95, 6), 20)  # Variable
        else:
            cdr_lengths["CDR1"] = 5
            cdr_lengths["CDR2"] = 16
            cdr_lengths["CDR3"] = 8
    else:  # VL
        # VL is approximately 107 residues
        seq_len = len(sequence)
        if seq_len >= 100:
            cdr_lengths["CDR1"] = 11  # Typical kappa
            cdr_lengths["CDR2"] = 7   # Typical
            cdr_lengths["CDR3"] = 9   # Typical
        else:
            cdr_lengths["CDR1"] = 10
            cdr_lengths["CDR2"] = 6
            cdr_lengths["CDR3"] = 8
    
    return cdr_lengths


def calculate_Ltot_sequence(vh_sequence: str, vl_sequence: str) -> Tuple[int, Dict[str, int]]:
    """Calculate total CDR length from sequences only (approximation)."""
    vh_cdrs = calculate_cdr_length_from_sequence(vh_sequence, "H")
    vl_cdrs = calculate_cdr_length_from_sequence(vl_sequence, "L")
    
    all_cdrs = {**vh_cdrs, **vl_cdrs}
    Ltot = sum(all_cdrs.values())
    
    return Ltot, all_cdrs


# ---------------------------------------------------------------------------
# Structure-Based Analysis
# ---------------------------------------------------------------------------

def calculate_cdr_lengths_from_structure(atoms: Dict[str, Dict[int, List[dict]]]) -> Dict[str, Dict[str, int]]:
    """Calculate CDR lengths from structure using Kabat numbering.
    
    Determines chain type based on:
    1. Chain name (H = heavy, L = light, A/B = first/second chain)
    2. Residue count (VH ~120, VL ~107)
    """
    cdr_lengths = {}
    
    # Sort chains by residue count to identify VH vs VL
    chain_info = []
    for chain, residues in atoms.items():
        n_residues = len(residues)
        chain_info.append((chain, n_residues, residues))
    
    # Sort by residue count (VH is typically longer)
    chain_info.sort(key=lambda x: x[1], reverse=True)
    
    # Assign chain types
    for i, (chain, n_residues, residues) in enumerate(chain_info):
        # Direct chain name mapping
        if chain in ("H", "H1", "H2"):
            chain_type = "H"
        elif chain in ("L", "L1", "L2"):
            chain_type = "L"
        # For A/B chains: A=VH, B=VL (common in AF3 output)
        elif chain == "A" and i == 0:
            chain_type = "H"  # First chain (longer) = VH
        elif chain == "B" and i == 1:
            chain_type = "L"  # Second chain (shorter) = VL
        # Fallback: use residue count
        elif n_residues > 115:
            chain_type = "H"
        else:
            chain_type = "L"
        
        cdr_lengths[chain] = {}
        for cdr_name, (start, end) in KABAT_CDRS[chain_type].items():
            # Count residues in CDR region
            count = sum(1 for res_num in residues.keys() if start <= res_num <= end)
            cdr_lengths[chain][cdr_name] = max(count, 1)  # Minimum 1
    
    return cdr_lengths


def calculate_SASA_metrics(sasa_data: Dict[str, Dict[int, float]],
                           sequences: Dict[str, str]) -> Tuple[float, float, float]:
    """Calculate PSH, PPC, PNC from SASA data.
    
    PSH: Number of hydrophobic patches on CDR vicinity surface
    PPC: Number of positive charge patches
    PNC: Number of negative charge patches
    
    TAP uses patch counting, not raw SASA summation.
    A patch is defined as a contiguous cluster of residues with specific properties.
    """
    # Maximum SASA values for normalization (Tien et al. 2013)
    MAX_SASA = {
        'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0,
        'C': 167.0, 'E': 223.0, 'Q': 225.0, 'G': 104.0,
        'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0,
        'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0,
        'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0,
    }
    
    # Collect all surface-exposed residues with their properties
    exposed_residues = []
    
    for chain, res_sasa in sasa_data.items():
        seq = sequences.get(chain, "")
        
        for res_num, sasa in res_sasa.items():
            if res_num > len(seq):
                continue
            
            aa = seq[res_num - 1]  # Convert to 0-based
            
            # Calculate relative SASA
            max_sasa = MAX_SASA.get(aa, 150.0)
            rel_sasa = sasa / max_sasa if max_sasa > 0 else 0
            
            # Only count surface-exposed residues (relSASA > 0.2)
            if rel_sasa > 0.2:
                exposed_residues.append({
                    'chain': chain,
                    'res_num': res_num,
                    'aa': aa,
                    'sasa': sasa,
                    'rel_sasa': rel_sasa,
                    'is_hydrophobic': aa in HYDROPHOBIC,
                    'charge': CHARGED_POS.get(aa, 0) + CHARGED_NEG.get(aa, 0),
                })
    
    # Count patches using simple clustering
    # A patch is a group of 3+ consecutive residues with same property
    
    def count_patches(residues, property_check):
        """Count patches of residues satisfying property_check."""
        patches = 0
        current_patch = 0
        
        for i, res in enumerate(residues):
            if property_check(res):
                current_patch += 1
            else:
                if current_patch >= 3:  # Minimum patch size
                    patches += 1
                current_patch = 0
        
        # Check last patch
        if current_patch >= 3:
            patches += 1
        
        return patches
    
    # Sort by chain and position for patch detection
    exposed_residues.sort(key=lambda x: (x['chain'], x['res_num']))
    
    # Count hydrophobic patches (PSH)
    psh = count_patches(exposed_residues, lambda r: r['is_hydrophobic'])
    
    # Count positive charge patches (PPC)
    ppc = count_patches(exposed_residues, lambda r: r['charge'] > 0)
    
    # Count negative charge patches (PNC)
    pnc = count_patches(exposed_residues, lambda r: r['charge'] < 0)
    
    return psh, ppc, pnc


def calculate_SFvCSP(sasa_data: Dict[str, Dict[int, float]],
                     sequences: Dict[str, str]) -> float:
    """Calculate Structural Fv Charge Symmetry Parameter.
    
    SFvCSP = net_charge_VH - net_charge_VL
    
    More negative values indicate greater charge asymmetry.
    TAP uses this to assess whether VH and VL have balanced charges.
    """
    # Maximum SASA values for normalization
    MAX_SASA = {
        'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0,
        'C': 167.0, 'E': 223.0, 'Q': 225.0, 'G': 104.0,
        'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0,
        'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0,
        'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0,
    }
    
    vh_charge = 0.0
    vl_charge = 0.0
    
    for chain, res_sasa in sasa_data.items():
        seq = sequences.get(chain, "")
        
        # Determine if this is VH or VL based on chain name or length
        # In PDB, VH is typically chain H or A, VL is chain L or B
        # Or we can check residue count
        is_vh = chain in ("H", "H1", "H2", "A") or len(res_sasa) > 115
        
        for res_num, sasa in res_sasa.items():
            if res_num > len(seq):
                continue
            
            aa = seq[res_num - 1]  # Convert to 0-based
            
            # Calculate relative SASA
            max_sasa = MAX_SASA.get(aa, 150.0)
            rel_sasa = sasa / max_sasa if max_sasa > 0 else 0
            
            # Only count surface-exposed residues (relSASA > 0.2)
            if rel_sasa < 0.2:
                continue
            
            # Calculate charge contribution
            if aa in CHARGED_POS:
                charge = CHARGED_POS[aa]
            elif aa in CHARGED_NEG:
                charge = CHARGED_NEG[aa]
            else:
                charge = 0
            
            if is_vh:
                vh_charge += charge
            else:
                vl_charge += charge
    
    # SFvCSP = VH charge - VL charge
    # More negative = more asymmetric = higher risk
    sfdcsp = vh_charge - vl_charge
    
    return sfdcsp


# ---------------------------------------------------------------------------
# Flag Calculation
# ---------------------------------------------------------------------------

def calculate_flag(value: float, thresholds: Dict) -> str:
    """Calculate flag (green/amber/red) for a metric value."""
    # Check green (acceptable)
    green_low, green_high = thresholds["green"]
    if green_low <= value <= green_high:
        return "green"
    
    # Check amber (caution)
    for amber_low, amber_high in thresholds["amber"]:
        if amber_low <= value <= amber_high:
            return "amber"
    
    # Check red (danger)
    for red_low, red_high in thresholds["red"]:
        if red_low <= value <= red_high:
            return "red"
    
    # Default to red if outside all ranges
    return "red"


def calculate_risk_summary(flags: Dict[str, str]) -> str:
    """Calculate overall risk summary from individual flags."""
    red_count = sum(1 for f in flags.values() if f == "red")
    amber_count = sum(1 for f in flags.values() if f == "amber")
    
    if red_count >= 2:
        return "HIGH RISK"
    elif red_count == 1:
        return "MEDIUM-HIGH RISK"
    elif amber_count >= 2:
        return "MEDIUM RISK"
    elif amber_count == 1:
        return "LOW-MEDIUM RISK"
    else:
        return "LOW RISK"


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------

def analyze_tap(pdb_path: str, name: str = "") -> TAPResult:
    """Perform complete TAP analysis on a structure."""
    if not name:
        name = Path(pdb_path).stem
    
    # Parse structure
    atoms = parse_pdb(pdb_path)
    sequences = get_chain_sequence(atoms)
    
    # Calculate SASA
    sasa_data = calculate_sasa(pdb_path)
    
    # Calculate CDR lengths
    cdr_lengths = calculate_cdr_lengths_from_structure(atoms)
    
    # Calculate Ltot
    Ltot = 0
    cdr_details = {}
    for chain, cdrs in cdr_lengths.items():
        for cdr_name, length in cdrs.items():
            key = f"{cdr_name}_{chain}"
            cdr_details[key] = length
            Ltot += length
    
    # Calculate SASA-based metrics
    psh, ppc, pnc = calculate_SASA_metrics(sasa_data, sequences)
    
    # Calculate SFvCSP
    sfdcsp = calculate_SFvCSP(sasa_data, sequences)
    
    # Determine which chain is VH and which is VL
    # Sort by residue count (VH is typically longer)
    chain_sorted = sorted(cdr_lengths.keys(), 
                          key=lambda c: len(atoms.get(c, {})), 
                          reverse=True)
    
    vh_chain = chain_sorted[0] if len(chain_sorted) > 0 else "H"
    vl_chain = chain_sorted[1] if len(chain_sorted) > 1 else "L"
    
    # Create metrics object
    metrics = TAPMetrics(
        Ltot=Ltot,
        PSH=psh,
        PPC=ppc,
        PNC=pnc,
        SFvCSP=sfdcsp,
        cdr1_h_len=cdr_lengths.get(vh_chain, {}).get("CDR1", 0),
        cdr2_h_len=cdr_lengths.get(vh_chain, {}).get("CDR2", 0),
        cdr3_h_len=cdr_lengths.get(vh_chain, {}).get("CDR3", 0),
        cdr1_l_len=cdr_lengths.get(vl_chain, {}).get("CDR1", 0),
        cdr2_l_len=cdr_lengths.get(vl_chain, {}).get("CDR2", 0),
        cdr3_l_len=cdr_lengths.get(vl_chain, {}).get("CDR3", 0),
    )
    
    # Calculate flags
    flags = {
        "Ltot": calculate_flag(Ltot, TAP_THRESHOLDS["Ltot"]),
        "PSH": calculate_flag(psh, TAP_THRESHOLDS["PSH"]),
        "PPC": calculate_flag(ppc, TAP_THRESHOLDS["PPC"]),
        "PNC": calculate_flag(pnc, TAP_THRESHOLDS["PNC"]),
        "SFvCSP": calculate_flag(sfdcsp, TAP_THRESHOLDS["SFvCSP"]),
    }
    
    # Update metrics with flags
    metrics.Ltot_flag = flags["Ltot"]
    metrics.PSH_flag = flags["PSH"]
    metrics.PPC_flag = flags["PPC"]
    metrics.PNC_flag = flags["PNC"]
    metrics.SFvCSP_flag = flags["SFvCSP"]
    
    # Risk summary
    risk_summary = calculate_risk_summary(flags)
    
    # Detailed information
    details = {
        "pdb_file": pdb_path,
        "vh_sequence": sequences.get("H", ""),
        "vl_sequence": sequences.get("L", ""),
        "vh_length": len(sequences.get("H", "")),
        "vl_length": len(sequences.get("L", "")),
        "cdr_lengths": cdr_details,
        "sasa_per_chain": {
            chain: sum(res_sasa.values()) 
            for chain, res_sasa in sasa_data.items()
        },
    }
    
    return TAPResult(
        name=name,
        vh_sequence=sequences.get("H", ""),
        vl_sequence=sequences.get("L", ""),
        metrics=metrics,
        flags=flags,
        risk_summary=risk_summary,
        details=details,
    )


def analyze_tap_from_sequences(vh_sequence: str, vl_sequence: str, 
                                name: str = "sequence") -> TAPResult:
    """Perform TAP analysis from sequences only (approximation)."""
    # Calculate Ltot from sequences
    Ltot, cdr_details = calculate_Ltot_sequence(vh_sequence, vl_sequence)
    
    # For sequence-only analysis, we cannot calculate PSH/PPC/PNC/SFvCSP
    # These require structural information
    psh = 0.0
    ppc = 0.0
    pnc = 0.0
    sfdcsp = 0.0
    
    # Create metrics object
    metrics = TAPMetrics(
        Ltot=Ltot,
        PSH=psh,
        PPC=ppc,
        PNC=pnc,
        SFvCSP=sfdcsp,
        cdr1_h_len=cdr_details.get("CDR1", 0),
        cdr2_h_len=cdr_details.get("CDR2", 0),
        cdr3_h_len=cdr_details.get("CDR3", 0),
        cdr1_l_len=cdr_details.get("CDR1", 0),
        cdr2_l_len=cdr_details.get("CDR2", 0),
        cdr3_l_len=cdr_details.get("CDR3", 0),
    )
    
    # Calculate flags (only Ltot is available for sequence-only)
    flags = {
        "Ltot": calculate_flag(Ltot, TAP_THRESHOLDS["Ltot"]),
        "PSH": "unknown",  # Requires structure
        "PPC": "unknown",  # Requires structure
        "PNC": "unknown",  # Requires structure
        "SFvCSP": "unknown",  # Requires structure
    }
    
    # Risk summary (based only on available metrics)
    risk_summary = "UNKNOWN (sequence-only analysis)"
    
    # Detailed information
    details = {
        "analysis_type": "sequence_only",
        "vh_sequence": vh_sequence,
        "vl_sequence": vl_sequence,
        "vh_length": len(vh_sequence),
        "vl_length": len(vl_sequence),
        "cdr_lengths": cdr_details,
        "note": "PSH/PPC/PNC/SFvCSP require structural analysis",
    }
    
    return TAPResult(
        name=name,
        vh_sequence=vh_sequence,
        vl_sequence=vl_sequence,
        metrics=metrics,
        flags=flags,
        risk_summary=risk_summary,
        details=details,
    )


# ---------------------------------------------------------------------------
# Output Formatting
# ---------------------------------------------------------------------------

def format_flag(flag: str) -> str:
    """Format flag with color indicator."""
    if flag == "green":
        return "✓ GREEN"
    elif flag == "amber":
        return "⚠ AMBER"
    elif flag == "red":
        return "✗ RED"
    else:
        return "? UNKNOWN"


def format_text_report(result: TAPResult, title: str = "TAP Analysis") -> str:
    """Format results as text report."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"  {title}")
    lines.append("=" * 80)
    lines.append(f"\n  Name: {result.name}")
    lines.append(f"  VH length: {len(result.vh_sequence)} residues")
    lines.append(f"  VL length: {len(result.vl_sequence)} residues")
    
    lines.append(f"\n{'='*80}")
    lines.append("  TAP METRICS")
    lines.append(f"{'='*80}")
    
    m = result.metrics
    lines.append(f"\n  {'Metric':<15} {'Value':<12} {'Flag':<15} {'Threshold (Green)'}")
    lines.append(f"  {'-'*15} {'-'*12} {'-'*15} {'-'*30}")
    
    lines.append(f"  {'Ltot':<15} {m.Ltot:<12} {format_flag(m.Ltot_flag):<15} {'43-54'}")
    lines.append(f"  {'PSH':<15} {m.PSH:<12.2f} {format_flag(m.PSH_flag):<15} {'111-168'}")
    lines.append(f"  {'PPC':<15} {m.PPC:<12.2f} {format_flag(m.PPC_flag):<15} {'0-1.33'}")
    lines.append(f"  {'PNC':<15} {m.PNC:<12.2f} {format_flag(m.PNC_flag):<15} {'0-1.98'}")
    lines.append(f"  {'SFvCSP':<15} {m.SFvCSP:<12.2f} {format_flag(m.SFvCSP_flag):<15} {'>-6.00'}")
    
    lines.append(f"\n{'='*80}")
    lines.append("  CDR LENGTHS")
    lines.append(f"{'='*80}")
    lines.append(f"\n  {'CDR':<10} {'VH':<8} {'VL':<8}")
    lines.append(f"  {'-'*10} {'-'*8} {'-'*8}")
    lines.append(f"  {'CDR1':<10} {m.cdr1_h_len:<8} {m.cdr1_l_len:<8}")
    lines.append(f"  {'CDR2':<10} {m.cdr2_h_len:<8} {m.cdr2_l_len:<8}")
    lines.append(f"  {'CDR3':<10} {m.cdr3_h_len:<8} {m.cdr3_l_len:<8}")
    lines.append(f"  {'Total':<10} {m.cdr1_h_len+m.cdr2_h_len+m.cdr3_h_len:<8} "
                f"{m.cdr1_l_len+m.cdr2_l_len+m.cdr3_l_len:<8}")
    
    lines.append(f"\n{'='*80}")
    lines.append("  RISK ASSESSMENT")
    lines.append(f"{'='*80}")
    lines.append(f"\n  Overall: {result.risk_summary}")
    
    # Flag summary
    red_flags = [k for k, v in result.flags.items() if v == "red"]
    amber_flags = [k for k, v in result.flags.items() if v == "amber"]
    
    if red_flags:
        lines.append(f"  Red flags: {', '.join(red_flags)}")
    if amber_flags:
        lines.append(f"  Amber flags: {', '.join(amber_flags)}")
    if not red_flags and not amber_flags:
        lines.append(f"  No concerning flags detected")
    
    lines.append("")
    
    return "\n".join(lines)


def save_json_report(result: TAPResult, output_path: str):
    """Save results as JSON file."""
    data = {
        "name": result.name,
        "vh_sequence": result.vh_sequence,
        "vl_sequence": result.vl_sequence,
        "metrics": asdict(result.metrics),
        "flags": result.flags,
        "risk_summary": result.risk_summary,
        "details": result.details,
    }
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def save_markdown_report(result: TAPResult, output_path: str, title: str = "TAP Analysis"):
    """Save results as Markdown file."""
    lines = []
    lines.append(f"# {title}\n")
    lines.append(f"**Name:** {result.name}\n")
    lines.append(f"**VH length:** {len(result.vh_sequence)} residues")
    lines.append(f"**VL length:** {len(result.vl_sequence)} residues\n")
    
    lines.append(f"\n## TAP Metrics\n")
    lines.append("| Metric | Value | Flag | Threshold (Green) |")
    lines.append("|--------|-------|------|-------------------|")
    
    m = result.metrics
    lines.append(f"| Ltot | {m.Ltot} | {format_flag(m.Ltot_flag)} | 43-54 |")
    lines.append(f"| PSH | {m.PSH:.2f} | {format_flag(m.PSH_flag)} | 111-168 |")
    lines.append(f"| PPC | {m.PPC:.2f} | {format_flag(m.PPC_flag)} | 0-1.33 |")
    lines.append(f"| PNC | {m.PNC:.2f} | {format_flag(m.PNC_flag)} | 0-1.98 |")
    lines.append(f"| SFvCSP | {m.SFvCSP:.2f} | {format_flag(m.SFvCSP_flag)} | >-6.00 |")
    
    lines.append(f"\n## CDR Lengths\n")
    lines.append("| CDR | VH | VL |")
    lines.append("|-----|-----|-----|")
    lines.append(f"| CDR1 | {m.cdr1_h_len} | {m.cdr1_l_len} |")
    lines.append(f"| CDR2 | {m.cdr2_h_len} | {m.cdr2_l_len} |")
    lines.append(f"| CDR3 | {m.cdr3_h_len} | {m.cdr3_l_len} |")
    lines.append(f"| **Total** | **{m.cdr1_h_len+m.cdr2_h_len+m.cdr3_h_len}** | "
                f"**{m.cdr1_l_len+m.cdr2_l_len+m.cdr3_l_len}** |")
    
    lines.append(f"\n## Risk Assessment\n")
    lines.append(f"**Overall:** {result.risk_summary}\n")
    
    red_flags = [k for k, v in result.flags.items() if v == "red"]
    amber_flags = [k for k, v in result.flags.items() if v == "amber"]
    
    if red_flags:
        lines.append(f"- **Red flags:** {', '.join(red_flags)}")
    if amber_flags:
        lines.append(f"- **Amber flags:** {', '.join(amber_flags)}")
    if not red_flags and not amber_flags:
        lines.append(f"- No concerning flags detected")
    
    with open(output_path, 'w') as f:
        f.write("\n".join(lines))


def save_word_report(result: TAPResult, output_path: str, title: str = "TAP Analysis"):
    """Save results as Word document."""
    try:
        from docx import Document
        from docx.shared import Inches, Pt, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.table import WD_TABLE_ALIGNMENT
    except ImportError:
        print("  WARNING: python-docx not installed. Skipping Word report.", file=sys.stderr)
        print("  Install with: pip install python-docx", file=sys.stderr)
        return
    
    doc = Document()
    
    # Title
    heading = doc.add_heading(title, 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER
    
    # Basic info
    doc.add_paragraph(f"Name: {result.name}", style='Intense Quote')
    doc.add_paragraph(f"VH length: {len(result.vh_sequence)} residues")
    doc.add_paragraph(f"VL length: {len(result.vl_sequence)} residues")
    
    # TAP Metrics table
    doc.add_heading('TAP Metrics', level=1)
    
    table = doc.add_table(rows=6, cols=4)
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    
    # Header
    headers = ['Metric', 'Value', 'Flag', 'Threshold (Green)']
    for i, header in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = header
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.bold = True
                run.font.size = Pt(9)
    
    # Data rows
    m = result.metrics
    rows_data = [
        ('Ltot', str(m.Ltot), format_flag(m.Ltot_flag), '43-54'),
        ('PSH', f'{m.PSH:.2f}', format_flag(m.PSH_flag), '111-168'),
        ('PPC', f'{m.PPC:.2f}', format_flag(m.PPC_flag), '0-1.33'),
        ('PNC', f'{m.PNC:.2f}', format_flag(m.PNC_flag), '0-1.98'),
        ('SFvCSP', f'{m.SFvCSP:.2f}', format_flag(m.SFvCSP_flag), '>-6.00'),
    ]
    
    for i, (metric, value, flag, threshold) in enumerate(rows_data):
        row = table.rows[i + 1]
        row.cells[0].text = metric
        row.cells[1].text = value
        row.cells[2].text = flag
        row.cells[3].text = threshold
        
        # Color code flag
        flag_cell = row.cells[2]
        for paragraph in flag_cell.paragraphs:
            for run in paragraph.runs:
                if "RED" in flag:
                    run.font.color.rgb = RGBColor(255, 0, 0)
                    run.font.bold = True
                elif "AMBER" in flag:
                    run.font.color.rgb = RGBColor(255, 165, 0)
                    run.font.bold = True
                elif "GREEN" in flag:
                    run.font.color.rgb = RGBColor(0, 128, 0)
        
        # Set font size
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(9)
    
    # CDR Lengths table
    doc.add_heading('CDR Lengths', level=1)
    
    cdr_table = doc.add_table(rows=5, cols=3)
    cdr_table.style = 'Light Grid Accent 1'
    cdr_table.alignment = WD_TABLE_ALIGNMENT.CENTER
    
    cdr_headers = ['CDR', 'VH', 'VL']
    for i, header in enumerate(cdr_headers):
        cell = cdr_table.rows[0].cells[i]
        cell.text = header
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.bold = True
                run.font.size = Pt(9)
    
    cdr_rows = [
        ('CDR1', str(m.cdr1_h_len), str(m.cdr1_l_len)),
        ('CDR2', str(m.cdr2_h_len), str(m.cdr2_l_len)),
        ('CDR3', str(m.cdr3_h_len), str(m.cdr3_l_len)),
        ('Total', str(m.cdr1_h_len+m.cdr2_h_len+m.cdr3_h_len),
         str(m.cdr1_l_len+m.cdr2_l_len+m.cdr3_l_len)),
    ]
    
    for i, (cdr, vh_len, vl_len) in enumerate(cdr_rows):
        row = cdr_table.rows[i + 1]
        row.cells[0].text = cdr
        row.cells[1].text = vh_len
        row.cells[2].text = vl_len
        
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(9)
    
    # Risk Assessment
    doc.add_heading('Risk Assessment', level=1)
    
    risk_para = doc.add_paragraph()
    risk_para.add_run('Overall: ').bold = True
    risk_para.add_run(result.risk_summary)
    
    red_flags = [k for k, v in result.flags.items() if v == "red"]
    amber_flags = [k for k, v in result.flags.items() if v == "amber"]
    
    if red_flags:
        p = doc.add_paragraph()
        p.add_run('Red flags: ').bold = True
        p.add_run(', '.join(red_flags))
        for run in p.runs:
            if ', '.join(red_flags) in run.text:
                run.font.color.rgb = RGBColor(255, 0, 0)
    
    if amber_flags:
        p = doc.add_paragraph()
        p.add_run('Amber flags: ').bold = True
        p.add_run(', '.join(amber_flags))
        for run in p.runs:
            if ', '.join(amber_flags) in run.text:
                run.font.color.rgb = RGBColor(255, 165, 0)
    
    if not red_flags and not amber_flags:
        doc.add_paragraph('No concerning flags detected')
    
    doc.save(output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="TAP (Therapeutic Antibody Profiler) Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze structure
  python3 tools/tap_profiler/tap_analyzer.py --pdb af3/C45H4_coVL/Structure\\ Prediction\\ \\(Boltz-2\\)/rank_1.pdb

  # Analyze sequences only (approximation)
  python3 tools/tap_profiler/tap_analyzer.py --vh <VH_SEQUENCE> --vl <VL_SEQUENCE>

  # Batch analysis
  python3 tools/tap_profiler/tap_analyzer.py --pdb-dir af3/ --output outputs/tap/
        """
    )
    parser.add_argument("--pdb", help="Path to PDB/CIF structure file")
    parser.add_argument("--pdb-dir", help="Directory containing PDB files (recursive)")
    parser.add_argument("--vh", help="VH amino acid sequence")
    parser.add_argument("--vl", help="VL amino acid sequence")
    parser.add_argument("--name", default="", help="Name for the analysis")
    parser.add_argument("--output", "-o", default="outputs/tap_profiler",
                        help="Output directory (default: outputs/tap_profiler)")
    parser.add_argument("--json", action="store_true", help="Save JSON report")
    parser.add_argument("--markdown", action="store_true", help="Save Markdown report")
    parser.add_argument("--word", action="store_true", help="Save Word report")
    parser.add_argument("--all-formats", action="store_true",
                        help="Save all formats (text + JSON + Markdown + Word)")
    
    args = parser.parse_args()
    
    if not args.pdb and not args.pdb_dir and not args.vh:
        parser.error("Either --pdb, --pdb-dir, or --vh/--vl must be specified")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Collect results
    results = []
    
    if args.pdb:
        # Analyze single structure
        print(f"\nAnalyzing structure: {args.pdb}")
        result = analyze_tap(args.pdb, args.name or Path(args.pdb).stem)
        results.append(result)
        
        # Print text report
        text_report = format_text_report(result, f"TAP Analysis: {result.name}")
        print(text_report)
        
    elif args.pdb_dir:
        # Batch analysis
        pdb_dir = Path(args.pdb_dir)
        pdb_files = list(pdb_dir.rglob("*.pdb"))
        
        if not pdb_files:
            print("No PDB files found.", file=sys.stderr)
            sys.exit(1)
        
        print(f"\nFound {len(pdb_files)} PDB files")
        
        for pdb_path in pdb_files:
            print(f"\nAnalyzing: {pdb_path}")
            try:
                result = analyze_tap(str(pdb_path), pdb_path.stem)
                results.append(result)
                
                # Print text report
                text_report = format_text_report(result, f"TAP Analysis: {result.name}")
                print(text_report)
                
            except Exception as e:
                print(f"  ERROR: {e}", file=sys.stderr)
                continue
    
    elif args.vh:
        # Sequence-only analysis
        if not args.vl:
            parser.error("--vl is required when using --vh")
        
        print(f"\nAnalyzing sequences (approximation)")
        result = analyze_tap_from_sequences(args.vh, args.vl, args.name or "sequence")
        results.append(result)
        
        # Print text report
        text_report = format_text_report(result, f"TAP Analysis: {result.name}")
        print(text_report)
    
    # Save reports
    for result in results:
        # Use parent directory name for unique filenames
        if hasattr(result, 'details') and 'pdb_file' in result.details:
            pdb_path = Path(result.details['pdb_file'])
            parent_name = pdb_path.parent.parent.name.replace("Structure Prediction (Boltz-2)_", "")
            stem = f"{parent_name}_{pdb_path.stem}" if parent_name else pdb_path.stem
        else:
            stem = result.name
        
        if args.json or args.all_formats:
            json_path = os.path.join(args.output, f"{stem}_tap.json")
            save_json_report(result, json_path)
            print(f"  JSON report saved: {json_path}")
        
        if args.markdown or args.all_formats:
            md_path = os.path.join(args.output, f"{stem}_tap.md")
            save_markdown_report(result, md_path, f"TAP Analysis: {result.name}")
            print(f"  Markdown report saved: {md_path}")
        
        if args.word or args.all_formats:
            docx_path = os.path.join(args.output, f"{stem}_tap.docx")
            save_word_report(result, docx_path, f"TAP Analysis: {result.name}")
            print(f"  Word report saved: {docx_path}")
    
    # Save combined JSON if multiple files
    if len(results) > 1 and (args.json or args.all_formats):
        combined_path = os.path.join(args.output, "combined_tap_results.json")
        combined_data = [asdict(r) for r in results]
        with open(combined_path, 'w') as f:
            json.dump(combined_data, f, indent=2, ensure_ascii=False)
        print(f"\nCombined JSON report saved: {combined_path}")


if __name__ == "__main__":
    main()
