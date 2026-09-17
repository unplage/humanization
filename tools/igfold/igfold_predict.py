#!/usr/bin/env python3
"""IgFold Structure Prediction Script.

Predicts antibody structures using IgFold from FASTA sequences.
Outputs PDB files for RMSD evaluation with the humanization pipeline.

Requirements:
    - conda environment: igfold (Python >= 3.11, PyTorch >= 2.0)
    - Install: conda create -n igfold python=3.11 && conda run -n igfold pip install igfold
    - AntiBERTy weights: First run requires internet to download ~100MB weights from Hugging Face
      If offline, set ANTIbertY_WEIGHTS_DIR to a local copy

Usage:
    # Predict a single antibody (VH + VL)
    conda run -n igfold python tools/igfold_predict.py --input af3/amg110_variants/amg110_V1.fasta --output af3/amg110_variants/

    # Predict all variants in a directory
    conda run -n igfold python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/

    # Predict with refinement and renumbering
    conda run -n igfold python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/ --refine openmm --renumber

    # Use GPU
    conda run -n igfold python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/ --device cuda:0

    # Offline mode (if weights already downloaded)
    conda run -n igfold python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/

Author: Humanization Pipeline
Date: 2026
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse FASTA file and return {chain_id: sequence} dict.
    
    Supports two formats:
    1. Standard FASTA with chain IDs in headers:
       >H
       EVQLVQSG...
       >L
       DVVMTQ...
    
    2. FASTA with descriptive headers (extracts first word after >):
       >V1_H Heavy chain
       EVQLVQSG...
       >V1_L Light chain
       DVVMTQ...
    
    Returns:
        Dict mapping chain_id (H/L) to sequence
    """
    chains = {}
    current_chain = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous chain if exists
                if current_chain and current_seq:
                    chains[current_chain] = ''.join(current_seq)
                
                # Parse header
                header = line[1:].strip()
                # Try to extract chain ID (H or L)
                parts = header.split()
                if parts:
                    first_part = parts[0]
                    # Check if it's a single letter H or L
                    if first_part in ('H', 'L'):
                        current_chain = first_part
                    # Check if it contains _H or _L
                    elif '_H' in first_part or first_part.endswith('H'):
                        current_chain = 'H'
                    elif '_L' in first_part or first_part.endswith('L'):
                        current_chain = 'L'
                    else:
                        # Default: first chain is H, second is L
                        current_chain = 'H' if not chains else 'L'
                else:
                    current_chain = 'H' if not chains else 'L'
                
                current_seq = []
            elif line and not line.startswith('#'):
                # Remove any whitespace and convert to uppercase
                current_seq.append(line.upper().replace(' ', ''))
        
        # Save last chain
        if current_chain and current_seq:
            chains[current_chain] = ''.join(current_seq)
    
    return chains


def check_network_connectivity() -> bool:
    """Check if network is available for downloading weights."""
    try:
        import urllib.request
        import urllib.error
        urllib.request.urlopen("https://huggingface.co", timeout=5)
        return True
    except (urllib.error.URLError, OSError, TimeoutError):
        return False


def check_weights_available() -> bool:
    """Check if AntiBERTy weights are already available locally.

    Returns True if weights exist in:
    1. ANTIbertY_WEIGHTS_DIR environment variable
    2. Packaged trained_models/AntiBERTy_md_smooth directory
    3. Hugging Face cache
    """
    import os
    # Check environment variable
    env_dir = os.environ.get("ANTIBERTY_WEIGHTS_DIR")
    if env_dir and os.path.isdir(env_dir):
        if os.path.isfile(os.path.join(env_dir, "config.json")) and \
           (os.path.isfile(os.path.join(env_dir, "model.safetensors")) or
            os.path.isfile(os.path.join(env_dir, "pytorch_model.bin"))):
            return True

    # Check packaged weights
    try:
        import antiberty
        pkg_dir = os.path.dirname(os.path.realpath(antiberty.__file__))
        packaged_dir = os.path.join(pkg_dir, "trained_models", "AntiBERTy_md_smooth")
        if os.path.isdir(packaged_dir):
            if os.path.isfile(os.path.join(packaged_dir, "config.json")) and \
               (os.path.isfile(os.path.join(packaged_dir, "model.safetensors")) or
                os.path.isfile(os.path.join(packaged_dir, "pytorch_model.bin"))):
                return True
    except ImportError:
        pass

    # Check Hugging Face cache
    hf_cache = os.path.expanduser("~/.cache/huggingface/hub")
    if os.path.isdir(hf_cache):
        for entry in os.listdir(hf_cache):
            if "AntiBERTy" in entry or "antiberty" in entry.lower():
                return True

    return False


def run_igfold_cli(
    sequences: Dict[str, str],
    output_path: str,
    refine: str = "none",
    renumber: bool = False,
    num_models: int = 4,
    device: str = "cpu",
    quiet: bool = True,
) -> bool:
    """Run IgFold using CLI interface.
    
    Args:
        sequences: Dict mapping chain_id to sequence (e.g., {"H": "...", "L": "..."})
        output_path: Path to output PDB file
        refine: Refinement method ("none", "openmm", "pyrosetta")
        renumber: Whether to renumber with Chothia scheme
        num_models: Number of ensemble models (1-4)
        device: Device to use ("cpu", "cuda:0", etc.)
        quiet: Suppress progress messages
    
    Returns:
        True if successful, False otherwise
    """
    # Build command
    cmd = [
        "igfold", "fold",
        "-o", output_path,
        "--num-models", str(num_models),
        "--device", device,
        "--refine", refine,
    ]
    
    if renumber:
        cmd.append("--renumber")
    
    if quiet:
        cmd.append("-q")
    
    # Add sequences
    if "H" in sequences:
        cmd.extend(["-H", sequences["H"]])
    if "L" in sequences:
        cmd.extend(["-L", sequences["L"]])
    
    # For single chain (nanobody)
    if len(sequences) == 1:
        chain_id = list(sequences.keys())[0]
        seq = sequences[chain_id]
        cmd = [c for c in cmd if c not in ("-H", "-L", sequences.get("H", ""), sequences.get("L", ""))]
        cmd.extend(["--chain", chain_id, seq])
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minutes timeout
        )
        
        if result.returncode != 0:
            print(f"  [ERROR] IgFold failed: {result.stderr[:500]}", file=sys.stderr)
            return False
        
        return True
    
    except subprocess.TimeoutExpired:
        print(f"  [ERROR] IgFold timed out after 5 minutes", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  [ERROR] IgFold failed: {e}", file=sys.stderr)
        return False


def run_igfold_python(
    sequences: Dict[str, str],
    output_path: str,
    do_refine: bool = False,
    do_renum: bool = False,
    device: str = "cpu",
) -> bool:
    """Run IgFold using Python API.
    
    Args:
        sequences: Dict mapping chain_id to sequence
        output_path: Path to output PDB file
        do_refine: Whether to refine with OpenMM
        do_renum: Whether to renumber with Chothia scheme
        device: Device to use
    
    Returns:
        True if successful, False otherwise
    """
    try:
        import torch
        from igfold import IgFoldRunner
        
        # Resolve device
        if device.startswith("cuda") and torch.cuda.is_available():
            device_id = int(device.split(":")[-1]) if ":" in device else 0
            resolved_device = torch.device(f"cuda:{device_id}")
        else:
            resolved_device = torch.device("cpu")
        
        igfold = IgFoldRunner(device=resolved_device)
        
        out = igfold.fold(
            output_path,
            sequences=sequences,
            do_refine=do_refine,
            do_renum=do_renum,
        )
        
        return True
    
    except ImportError as e:
        print(f"  [ERROR] Python API not available: {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  [ERROR] Python API failed: {e}", file=sys.stderr)
        return False


def validate_pdb(pdb_path: str) -> bool:
    """Validate that a PDB file is not empty and contains ATOM records.
    
    Args:
        pdb_path: Path to PDB file
    
    Returns:
        True if valid, False otherwise
    """
    if not os.path.isfile(pdb_path):
        return False
    
    try:
        with open(pdb_path, 'r') as f:
            content = f.read()
        
        # Check file is not empty
        if not content.strip():
            return False
        
        # Check for ATOM records
        if 'ATOM' not in content:
            return False
        
        return True
    except (IOError, OSError):
        return False


def predict_single_fasta(
    fasta_path: str,
    output_dir: str,
    refine: str = "none",
    renumber: bool = False,
    num_models: int = 4,
    device: str = "cpu",
    use_python_api: bool = False,
) -> Optional[str]:
    """Predict structure from a single FASTA file.
    
    Args:
        fasta_path: Path to input FASTA file
        output_dir: Directory to output PDB file
        refine: Refinement method
        renumber: Whether to renumber
        num_models: Number of ensemble models
        device: Device to use
        use_python_api: Use Python API instead of CLI
    
    Returns:
        Path to output PDB file if successful, None otherwise
    """
    # Parse FASTA
    sequences = parse_fasta(fasta_path)
    
    if not sequences:
        print(f"  [ERROR] No sequences found in {fasta_path}", file=sys.stderr)
        return None
    
    # Determine output filename
    fasta_name = Path(fasta_path).stem
    output_path = os.path.join(output_dir, f"{fasta_name}.pdb")
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Run prediction
    if use_python_api:
        success = run_igfold_python(
            sequences,
            output_path,
            do_refine=(refine != "none"),
            do_renum=renumber,
            device=device,
        )
    else:
        success = run_igfold_cli(
            sequences,
            output_path,
            refine=refine,
            renumber=renumber,
            num_models=num_models,
            device=device,
        )
    
    if success and os.path.exists(output_path):
        if validate_pdb(output_path):
            return output_path
        else:
            print(f"  [WARNING] Output PDB file is invalid or empty: {output_path}", file=sys.stderr)
            return None
    
    return None


def predict_batch(
    input_dir: str,
    output_dir: str,
    refine: str = "none",
    renumber: bool = False,
    num_models: int = 4,
    device: str = "cpu",
    use_python_api: bool = False,
) -> List[str]:
    """Predict structures for all FASTA files in a directory.
    
    Args:
        input_dir: Directory containing FASTA files
        output_dir: Directory to output PDB files
        refine: Refinement method
        renumber: Whether to renumber
        num_models: Number of ensemble models
        device: Device to use
        use_python_api: Use Python API instead of CLI
    
    Returns:
        List of output PDB file paths
    """
    input_path = Path(input_dir)
    fasta_files = sorted(input_path.glob("*.fasta")) + sorted(input_path.glob("*.fa"))
    
    if not fasta_files:
        print(f"[WARNING] No FASTA files found in {input_dir}", file=sys.stderr)
        return []
    
    print(f"[IgFold] Found {len(fasta_files)} FASTA files to predict")
    print(f"[IgFold] Output directory: {output_dir}")
    print(f"[IgFold] Refinement: {refine}")
    print(f"[IgFold] Renumber: {renumber}")
    print(f"[IgFold] Device: {device}")
    print()
    
    output_files = []
    total_time = 0
    
    for i, fasta_file in enumerate(fasta_files, 1):
        print(f"[{i}/{len(fasta_files)}] Predicting {fasta_file.name}...")
        
        start_time = time.time()
        
        output_path = predict_single_fasta(
            str(fasta_file),
            output_dir,
            refine=refine,
            renumber=renumber,
            num_models=num_models,
            device=device,
            use_python_api=use_python_api,
        )
        
        elapsed = time.time() - start_time
        total_time += elapsed
        
        if output_path:
            output_files.append(output_path)
            print(f"  ✓ Output: {os.path.basename(output_path)} ({elapsed:.1f}s)")
        else:
            print(f"  ✗ Failed to predict {fasta_file.name}")
        
        print()
    
    print(f"[IgFold] Completed: {len(output_files)}/{len(fasta_files)} structures predicted")
    print(f"[IgFold] Total time: {total_time:.1f}s ({total_time/len(fasta_files):.1f}s per structure)")
    
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description="IgFold Structure Prediction for Humanization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Predict a single antibody
    python tools/igfold_predict.py --input af3/amg110_variants/amg110_V1.fasta --output af3/amg110_variants/

    # Predict all variants in a directory
    python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/

    # Predict with refinement and renumbering
    python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/ --refine openmm --renumber

    # Use GPU
    python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/ --device cuda:0

    # Use Python API
    python tools/igfold_predict.py --input-dir af3/amg110_variants/ --output af3/amg110_variants/ --python-api

Requirements:
    - conda environment: igfold (Python >= 3.11, PyTorch >= 2.0)
    - Install: conda create -n igfold python=3.11 && conda run -n igfold pip install igfold
        """,
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        type=str,
        help="Input FASTA file (single antibody prediction)",
    )
    input_group.add_argument(
        "--input-dir",
        type=str,
        help="Input directory containing FASTA files (batch prediction)",
    )
    
    # Output
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Output directory for PDB files",
    )
    
    # Prediction options
    parser.add_argument(
        "--refine",
        choices=["none", "openmm", "pyrosetta"],
        default="none",
        help="Refinement method (default: none)",
    )
    parser.add_argument(
        "--renumber",
        action="store_true",
        help="Renumber residues with Chothia scheme",
    )
    parser.add_argument(
        "--num-models",
        type=int,
        default=4,
        choices=[1, 2, 3, 4],
        help="Number of ensemble models (1-4, default: 4)",
    )
    parser.add_argument(
        "--device",
        type=str,
        default="cpu",
        help="Device to use: cpu, cuda:0, cuda:1, mps (default: cpu)",
    )
    parser.add_argument(
        "--python-api",
        action="store_true",
        help="Use Python API instead of CLI (experimental)",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Skip network connectivity check (use when weights are already downloaded)",
    )
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)
    
    # Auto-detect offline mode if weights are already available
    if not args.offline and check_weights_available():
        print("[IgFold] AntiBERTy weights found locally - using offline mode")
        args.offline = True
    
    # Check network connectivity for first-time setup (unless --offline)
    if not args.offline:
        print("[IgFold] Checking network connectivity...")
        if not check_network_connectivity():
            print("[WARNING] Network is unreachable. IgFold requires internet for first-time setup", file=sys.stderr)
            print("          to download AntiBERTy weights (~100MB) from Hugging Face.", file=sys.stderr)
            print()
            print("Options:", file=sys.stderr)
            print("  1. Run on a machine with internet access", file=sys.stderr)
            print("  2. Set ANTIbertY_WEIGHTS_DIR to a local copy of the weights", file=sys.stderr)
            print("  3. Manually download weights from: https://huggingface.co/jeffruffolo/AntiBERTy", file=sys.stderr)
            print("  4. Use --offline flag if weights are already downloaded", file=sys.stderr)
            print()
            sys.exit(1)
    else:
        print("[IgFold] Offline mode - skipping network check")
    
    # Run prediction
    if args.input:
        # Single file prediction
        output_path = predict_single_fasta(
            args.input,
            args.output,
            refine=args.refine,
            renumber=args.renumber,
            num_models=args.num_models,
            device=args.device,
            use_python_api=args.python_api,
        )
        
        if output_path:
            print(f"\n✓ Structure predicted: {output_path}")
        else:
            print(f"\n✗ Failed to predict structure", file=sys.stderr)
            sys.exit(1)
    
    elif args.input_dir:
        # Batch prediction
        output_files = predict_batch(
            args.input_dir,
            args.output,
            refine=args.refine,
            renumber=args.renumber,
            num_models=args.num_models,
            device=args.device,
            use_python_api=args.python_api,
        )
        
        if not output_files:
            sys.exit(1)


if __name__ == "__main__":
    main()
