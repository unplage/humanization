#!/usr/bin/env python3
"""Create combined FASTA for immunogenicity analysis of 49 IgFold variants."""

import os

COMBINATORIAL_DIR = "/home/jiemiaoxing/work/opencode_task/humanization/af3/amg110_rerun_combinatorial"
OUTPUT_DIR = "/home/jiemiaoxing/work/opencode_task/humanization/outputs/immunogenicity_49"

# H variants (for naming)
H_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]
L_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]

def create_combined_fasta():
    """Create combined FASTA with all 49 variants."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    combined_fasta = os.path.join(OUTPUT_DIR, "all_49_variants.fasta")
    
    with open(combined_fasta, 'w') as out:
        for h_var in H_VARIANTS:
            for l_var in L_VARIANTS:
                fasta_file = os.path.join(COMBINATORIAL_DIR, f"H_{h_var}_L_{l_var}.fasta")
                if os.path.exists(fasta_file):
                    with open(fasta_file) as f:
                        lines = f.readlines()
                        # Read VH and VL sequences
                        vh_seq = ""
                        vl_seq = ""
                        current = None
                        for line in lines:
                            if line.startswith(">H"):
                                current = "H"
                            elif line.startswith(">L"):
                                current = "L"
                            elif current == "H":
                                vh_seq = line.strip()
                            elif current == "L":
                                vl_seq = line.strip()
                        
                        # Write combined entry
                        pair_name = f"H{h_var}_L{l_var}"
                        out.write(f">{pair_name}_VH\n{vh_seq}\n")
                        out.write(f">{pair_name}_VL\n{vl_seq}\n")
    
    print(f"Created combined FASTA: {combined_fasta}")
    return combined_fasta

if __name__ == "__main__":
    create_combined_fasta()
