#!/usr/bin/env python3
"""Set up directory structure for structural risk assessment."""

import os
import shutil

# Directories
IGFOLD_DIR = "/home/jiemiaoxing/work/opencode_task/humanization/af3/amg110_rerun_combinatorial"
STRUCTURE_DIR = "/home/jiemiaoxing/work/opencode_task/humanization/outputs/immunogenicity_49/structures"
IMMUNOGENICITY_JSON = "/home/jiemiaoxing/work/opencode_task/humanization/outputs/immunogenicity_49/all_49_variants_immunogenicity.json"

# H and L variants
H_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]
L_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]

def setup_structures():
    """Create symlinks for PDB files with expected naming patterns."""
    os.makedirs(STRUCTURE_DIR, exist_ok=True)
    
    for h_var in H_VARIANTS:
        for l_var in L_VARIANTS:
            # Source PDB file
            src_pdb = os.path.join(IGFOLD_DIR, f"H_{h_var}_L_{l_var}.pdb")
            if not os.path.exists(src_pdb):
                print(f"WARNING: {src_pdb} not found")
                continue
            
            # Create symlinks with various naming patterns
            # Pattern 1: HV0_LV0_VH.pdb and HV0_LV0_VL.pdb (matching immunogenicity variant names)
            for chain_suffix in ["_VH", "_VL"]:
                link_name = f"H{h_var}_L{l_var}{chain_suffix}.pdb"
                link_path = os.path.join(STRUCTURE_DIR, link_name)
                if not os.path.exists(link_path):
                    os.symlink(src_pdb, link_path)
                    print(f"Created symlink: {link_name} -> H_{h_var}_L_{l_var}.pdb")
    
    print(f"\nStructure directory: {STRUCTURE_DIR}")
    print(f"Files created: {len(os.listdir(STRUCTURE_DIR))}")

if __name__ == "__main__":
    setup_structures()
