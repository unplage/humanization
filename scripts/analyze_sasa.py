#!/usr/bin/env python3
"""Analyze FreeSASA for antibody structure to determine buried/exposed positions."""

import sys
import freesasa
from collections import defaultdict

def parse_pdb_atoms(pdb_path):
    """Parse PDB file and return atoms grouped by chain and residue."""
    atoms = defaultdict(lambda: defaultdict(list))
    
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
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
    return atoms

def calculate_sasa_per_residue(pdb_path):
    """Calculate per-residue SASA using FreeSASA."""
    # Parse structure
    structure = freesasa.Structure(pdb_path)
    
    # Calculate SASA
    result = freesasa.calc(structure)
    
    # Group by chain and residue
    sasa_by_residue = defaultdict(lambda: defaultdict(float))
    n = structure.nAtoms()
    
    for i in range(n):
        chain = structure.chainLabel(i)
        res_num = structure.residueNumber(i)
        area = result.atomArea(i)
        sasa_by_residue[chain][res_num] += area
    
    return sasa_by_residue, structure

def get_reference_sasa():
    """Get maximum SASA values for amino acids (Tien et al. 2013)."""
    return {
        'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
        'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0,
        'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
        'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
        'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
    }

def main():
    pdb_path = sys.argv[1] if len(sys.argv) > 1 else "af3/c81_af3/Structure Prediction (Boltz-2)/rank_1.pdb"
    
    print("=" * 70)
    print("FreeSASA Analysis for C81 Antibody")
    print("=" * 70)
    
    # Calculate SASA
    sasa_by_residue, structure = calculate_sasa_per_residue(pdb_path)
    ref_sasa = get_reference_sasa()
    
    # Analyze each chain
    for chain in ['A', 'B']:
        chain_type = 'VH' if chain == 'A' else 'VL'
        print(f"\n{'='*70}")
        print(f"  Chain {chain} ({chain_type})")
        print(f"{'='*70}")
        
        residues = sasa_by_residue[chain]
        print(f"{'Kabat':<10} {'Residue':<8} {'AbsSASA':<12} {'RelSASA':<10} {'Classification':<15}")
        print("-" * 60)
        
        for res_num in sorted(residues.keys()):
            abs_sasa = residues[res_num]
            # Get residue name from structure
            resname = None
            for i in range(structure.nAtoms()):
                if structure.chainLabel(i) == chain and structure.residueNumber(i) == res_num:
                    resname = structure.residueName(i)
                    break
            
            if resname and resname in ref_sasa:
                rel_sasa = abs_sasa / ref_sasa[resname]
            else:
                rel_sasa = 0.0
            
            # Classification
            if rel_sasa < 0.20:
                classification = "BURIED"
            elif rel_sasa < 0.50:
                classification = "INTERMEDIATE"
            else:
                classification = "EXPOSED"
            
            print(f"{res_num:<10} {resname:<8} {abs_sasa:<12.1f} {rel_sasa:<10.3f} {classification:<15}")

if __name__ == "__main__":
    main()
