#!/usr/bin/env python3
"""Generate combinatorial VH+VL pairs for IgFold structure prediction.

Creates all permutations of H and L chain variants for comprehensive
structural validation of humanized antibody candidates.

Usage:
    python3 scripts/humanize/combinatorial_igfold.py \
        --input outputs/amg110_step3/step3/variants.fasta \
        --output af3/amg110_combinatorial/ \
        --device cpu
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def parse_variants_fasta(fasta_path: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Parse variants FASTA and extract H and L chain sequences.
    
    Returns:
        (h_chains, l_chains): Dicts mapping variant name to sequence
    """
    h_chains = {}
    l_chains = {}
    
    current_name = None
    current_seq = []
    current_chain = None
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_name and current_seq:
                    seq = ''.join(current_seq)
                    if current_chain == 'H':
                        h_chains[current_name] = seq
                    else:
                        l_chains[current_name] = seq
                
                # Parse header: >AMG110.H|H_V0|description
                parts = line[1:].split('|')
                if len(parts) >= 2:
                    chain_type = parts[0].split('.')[-1]  # H or L
                    variant_name = parts[1]  # H_V0, L_V0, etc.
                    current_chain = chain_type
                    current_name = variant_name
                    current_seq = []
            elif line and not line.startswith('#'):
                current_seq.append(line)
        
        # Save last sequence
        if current_name and current_seq:
            seq = ''.join(current_seq)
            if current_chain == 'H':
                h_chains[current_name] = seq
            else:
                l_chains[current_name] = seq
    
    return h_chains, l_chains


def generate_combinations(
    h_chains: Dict[str, str],
    l_chains: Dict[str, str],
    include_original: bool = True
) -> List[Tuple[str, str, str, str]]:
    """Generate all H x L combinations.
    
    Returns:
        List of (combo_name, h_name, l_name, description) tuples
    """
    combinations = []
    
    # Define priority order for H chains
    h_priority = ['H_V0', 'H_V1', 'H_V2a', 'H_V2b', 'H_V2', 'H_V3', 'H_Vmin']
    l_priority = ['L_V0', 'L_V1', 'L_V2a', 'L_V2b', 'L_V2', 'L_V3', 'L_Vmin']
    
    # Filter to available variants
    h_available = [v for v in h_priority if v in h_chains]
    l_available = [v for v in l_priority if v in l_chains]
    
    # Generate all combinations
    for h_name in h_available:
        for l_name in l_available:
            # Extract tier from names
            h_tier = h_name.replace('H_', '')
            l_tier = l_name.replace('L_', '')
            
            # Skip original pairings if not requested
            if not include_original and h_tier == l_tier:
                continue
            
            # Create combination name
            if h_tier == l_tier:
                combo_name = f"{h_tier}"
                desc = f"Original pairing"
            else:
                combo_name = f"{h_tier}_{l_tier}"
                desc = f"H_{h_tier} + L_{l_tier}"
            
            combinations.append((combo_name, h_name, l_name, desc))
    
    return combinations


def create_combo_fasta(
    combo_name: str,
    h_name: str,
    h_seq: str,
    l_name: str,
    l_seq: str,
    output_dir: str
) -> str:
    """Create FASTA file for a combination.
    
    Returns:
        Path to created FASTA file
    """
    fasta_path = os.path.join(output_dir, f"combo_{combo_name}.fasta")
    
    with open(fasta_path, 'w') as f:
        f.write(f">H|{h_name}|{combo_name} heavy chain\n")
        f.write(f"{h_seq}\n")
        f.write(f">L|{l_name}|{combo_name} light chain\n")
        f.write(f"{l_seq}\n")
    
    return fasta_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate combinatorial VH+VL pairs for IgFold"
    )
    parser.add_argument("--input", required=True,
                        help="Path to variants FASTA file")
    parser.add_argument("--output", required=True,
                        help="Output directory for FASTA files")
    parser.add_argument("--device", default="cpu",
                        help="Device for IgFold (cpu/cuda:0)")
    parser.add_argument("--skip-original", action="store_true",
                        help="Skip original pairings (V0+V0, V1+V1, etc.)")
    parser.add_argument("--list-only", action="store_true",
                        help="Only list combinations, don't create files")
    
    args = parser.parse_args()
    
    # Parse input
    print(f"Parsing variants from: {args.input}")
    h_chains, l_chains = parse_variants_fasta(args.input)
    
    print(f"Found {len(h_chains)} H chain variants: {list(h_chains.keys())}")
    print(f"Found {len(l_chains)} L chain variants: {list(l_chains.keys())}")
    
    # Generate combinations
    combinations = generate_combinations(
        h_chains, l_chains,
        include_original=not args.skip_original
    )
    
    print(f"\nGenerated {len(combinations)} combinations:")
    print(f"{'#':<4} {'Combo':<12} {'H Chain':<10} {'L Chain':<10} {'Description'}")
    print("-" * 60)
    
    for i, (combo_name, h_name, l_name, desc) in enumerate(combinations, 1):
        print(f"{i:<4} {combo_name:<12} {h_name:<10} {l_name:<10} {desc}")
    
    if args.list_only:
        return
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Create FASTA files
    print(f"\nCreating FASTA files in: {args.output}")
    fasta_files = []
    
    for combo_name, h_name, l_name, desc in combinations:
        h_seq = h_chains[h_name]
        l_seq = l_chains[l_name]
        
        fasta_path = create_combo_fasta(
            combo_name, h_name, h_seq, l_name, l_seq, args.output
        )
        fasta_files.append(fasta_path)
        print(f"  Created: {os.path.basename(fasta_path)}")
    
    # Create batch script
    batch_script = os.path.join(args.output, "run_igfold_batch.sh")
    with open(batch_script, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# Batch IgFold prediction for combinatorial variants\n\n")
        f.write(f"OUTPUT_DIR={args.output}\n")
        f.write(f"DEVICE={args.device}\n\n")
        
        for fasta_path in fasta_files:
            fasta_name = Path(fasta_path).stem
            f.write(f"echo 'Predicting {fasta_name}...'\n")
            f.write(f"conda run -n igfold python tools/igfold/igfold_predict.py \\\n")
            f.write(f"  --input {fasta_path} \\\n")
            f.write(f"  --output $OUTPUT_DIR \\\n")
            f.write(f"  --device $DEVICE\n")
            f.write(f"\n")
    
    os.chmod(batch_script, 0o755)
    print(f"\nBatch script created: {batch_script}")
    
    # Create README
    readme_path = os.path.join(args.output, "README.md")
    with open(readme_path, 'w') as f:
        f.write("# Combinatorial IgFold Structure Prediction\n\n")
        f.write("## Overview\n\n")
        f.write(f"This directory contains {len(combinations)} combinatorial VH+VL pairs\n")
        f.write("for comprehensive structural validation of humanized AMG110 variants.\n\n")
        f.write("## Combinations\n\n")
        f.write("| # | Combo | H Chain | L Chain | Description |\n")
        f.write("|---|-------|---------|---------|-------------|\n")
        
        for i, (combo_name, h_name, l_name, desc) in enumerate(combinations, 1):
            f.write(f"| {i} | {combo_name} | {h_name} | {l_name} | {desc} |\n")
        
        f.write("\n## Usage\n\n")
        f.write("```bash\n")
        f.write("# Run all predictions\n")
        f.write("bash run_igfold_batch.sh\n\n")
        f.write("# Or run individual predictions\n")
        f.write("conda run -n igfold python tools/igfold/igfold_predict.py \\\n")
        f.write("  --input combo_V2a_V2b.fasta \\\n")
        f.write("  --output . \\\n")
        f.write("  --device cpu\n")
        f.write("```\n\n")
        f.write("## Expected Runtime\n\n")
        f.write("- CPU: ~6-10 seconds per structure\n")
        f.write(f"- Total: ~{len(combinations) * 8} seconds (~{len(combinations) * 8 // 60} minutes)\n")
    
    print(f"README created: {readme_path}")
    print(f"\nTotal combinations: {len(combinations)}")
    print(f"Estimated runtime: ~{len(combinations) * 8} seconds")


if __name__ == "__main__":
    main()
