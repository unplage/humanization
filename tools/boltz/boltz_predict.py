#!/usr/bin/env python3
"""Boltz-2 Structure Prediction Script.

Predicts antibody/protein structures using Boltz-2 from YAML or FASTA sequences.
Outputs PDB files with pLDDT in B-factor column for the humanization pipeline.

Requirements:
    - conda environment: boltz (Python >= 3.11, PyTorch >= 2.2)
    - Install: conda create -n boltz python=3.11 && conda run -n boltz pip install "torch==2.5.1+cu121" -f https://download.pytorch.org/whl/cu121 && conda run -n boltz pip install boltz==2.2.1
    - Weights: Place boltz2_conf.ckpt, boltz2_aff.ckpt, mols.tar in ~/.boltz/ (or $BOLTZ_CACHE)
    - MSA: Uses ColabFold server by default (--use-msa-server); for offline use add 'msa: empty' to YAML
    - Kernels: Disabled by default because cuequivariance_torch is optional. Pass --kernels only when it is installed.

Usage:
    # Predict a single antibody (H + L chains) from YAML
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/

    # Predict from FASTA (auto-generates YAML with chains H/L)
    conda run -n boltz python tools/boltz/boltz_predict.py --input antibody.fasta --output predictions/

    # Predict all FASTA files in a directory
    conda run -n boltz python tools/boltz/boltz_predict.py --input-dir variants/ --output predictions/

    # Single chain (VHH/nanobody)
    conda run -n boltz python tools/boltz/boltz_predict.py --input vhh.fasta --output predictions/ --single-chain

    # Offline mode (no MSA server, requires pre-computed MSA or single-sequence)
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --offline

    # Use multiple devices
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --devices 2

    # Multiple diffusion samples
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --diffusion-samples 3

    # Enable cuequivariance kernels (requires cuequivariance_torch)
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --kernels

    # Change the per-structure timeout (default 600s)
    conda run -n boltz python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --timeout 1800

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
from typing import Dict, List, Optional, Tuple, Any
import json


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """Parse FASTA file and return {chain_id: sequence} dict.
    
    Supports standard FASTA with chain IDs in headers:
        >H
        EVQLVQSG...
        >L
        DVVMTQ...
    
    Or descriptive headers:
        >V1_H Heavy chain
        EVQLVQSG...
        >V1_L Light chain
        DVVMTQ...
    
    Returns:
        Dict mapping chain_id (H/L/A) to sequence
    """
    chains = {}
    current_chain = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_chain and current_seq:
                    chains[current_chain] = ''.join(current_seq)
                
                header = line[1:].strip()
                parts = header.split()
                if parts:
                    first_part = parts[0]
                    if first_part in ('H', 'L', 'A'):
                        current_chain = first_part
                    elif '_H' in first_part or first_part.endswith('H'):
                        current_chain = 'H'
                    elif '_L' in first_part or first_part.endswith('L'):
                        current_chain = 'L'
                    elif '_A' in first_part or first_part.endswith('A'):
                        current_chain = 'A'
                    else:
                        current_chain = 'H' if not chains else ('L' if 'H' in chains else 'A')
                else:
                    current_chain = 'H' if not chains else ('L' if 'H' in chains else 'A')
                current_seq = []
            elif line and not line.startswith('#'):
                current_seq.append(line.upper().replace(' ', ''))
    
    if current_chain and current_seq:
        chains[current_chain] = ''.join(current_seq)
    
    return chains


def parse_yaml_input(yaml_path: str) -> Dict[str, Any]:
    """Parse Boltz YAML input file."""
    import yaml
    with open(yaml_path, 'r') as f:
        return yaml.safe_load(f)


def build_yaml_from_sequences(
    sequences: Dict[str, str],
    output_path: str,
    use_msa_server: bool = True,
    msa_dir: Optional[str] = None,
) -> str:
    """Build Boltz YAML input from sequences dict."""
    import yaml
    
    seq_list = []
    for chain_id, sequence in sequences.items():
        entry = {
            "protein": {
                "id": chain_id,
                "sequence": sequence,
            }
        }
        if use_msa_server:
            entry["protein"]["msa"] = "empty"  # Will use --use-msa-server flag
        elif msa_dir:
            msa_file = os.path.join(msa_dir, f"{chain_id}.a3m")
            if os.path.exists(msa_file):
                entry["protein"]["msa"] = msa_file
            else:
                entry["protein"]["msa"] = "empty"
        else:
            entry["protein"]["msa"] = "empty"
        seq_list.append(entry)
    
    yaml_content = {
        "version": 1,
        "sequences": seq_list,
    }
    
    with open(output_path, 'w') as f:
        yaml.dump(yaml_content, f, default_flow_style=False)
    
    return output_path


def check_boltz_available() -> bool:
    """Check if Boltz is importable."""
    try:
        import boltz
        return True
    except ImportError:
        return False


def check_weights_available() -> bool:
    """Check if Boltz weights are available locally."""
    boltz_cache = os.environ.get("BOLTZ_CACHE", os.path.expanduser("~/.boltz"))
    required = ["boltz2_conf.ckpt", "boltz2_aff.ckpt", "mols"]
    for req in required:
        if not os.path.exists(os.path.join(boltz_cache, req)):
            return False
    return True


def validate_pdb(pdb_path: str) -> bool:
    """Validate that a PDB file is not empty and contains ATOM records."""
    if not os.path.isfile(pdb_path):
        return False
    try:
        with open(pdb_path, 'r') as f:
            content = f.read()
        if not content.strip():
            return False
        if 'ATOM' not in content:
            return False
        return True
    except (IOError, OSError):
        return False


def run_boltz_cli(
    yaml_path: str,
    output_dir: str,
    use_msa_server: bool = True,
    devices: int = 1,
    diffusion_samples: int = 1,
    no_kernels: bool = True,
    timeout: int = 600,
    extra_args: Optional[List[str]] = None,
) -> bool:
    """Run Boltz using CLI interface.
    
    Args:
        yaml_path: Path to input YAML file
        output_dir: Output directory for predictions
        use_msa_server: Whether to use ColabFold MSA server
        devices: Number of GPU devices
        diffusion_samples: Number of diffusion samples
        no_kernels: Disable cuequivariance kernels
        timeout: Per-structure timeout in seconds
        extra_args: Additional arguments to pass to boltz predict
    
    Returns:
        True if successful, False otherwise
    """
    cmd = [
        "boltz", "predict",
        yaml_path,
        "--out_dir", output_dir,
        "--output_format", "pdb",
        "--devices", str(devices),
        "--diffusion_samples", str(diffusion_samples),
    ]
    
    if use_msa_server:
        cmd.append("--use_msa_server")
    
    if no_kernels:
        cmd.append("--no_kernels")
    
    if extra_args:
        cmd.extend(extra_args)
    
    try:
        print(f"  [Boltz] Running: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        
        if result.returncode != 0:
            detail = (result.stderr or "").strip() or (result.stdout or "").strip()
            print(f"  [ERROR] Boltz failed (exit {result.returncode})", file=sys.stderr)
            if detail:
                print(detail, file=sys.stderr)
            return False
        
        return True
    
    except subprocess.TimeoutExpired:
        print(f"  [ERROR] Boltz timed out after {timeout} seconds", file=sys.stderr)
        return False
    except Exception as e:
        print(f"  [ERROR] Boltz failed: {e}", file=sys.stderr)
        return False


def find_output_pdb(output_dir: str, yaml_stem: str) -> Optional[str]:
    """Find the output PDB file from Boltz prediction."""
    # Boltz output structure: output_dir/boltz_results_<stem>/predictions/<stem>/<stem>_model_0.pdb
    boltz_results = Path(output_dir) / f"boltz_results_{yaml_stem}"
    if not boltz_results.exists():
        # Try to find any boltz_results dir
        for d in Path(output_dir).glob("boltz_results_*"):
            boltz_results = d
            break
    
    if not boltz_results or not boltz_results.exists():
        return None
    
    # Look for model_0.pdb (best confidence)
    for pdb in boltz_results.glob(f"predictions/*/*_model_0.pdb"):
        return str(pdb)
    
    # Fallback: any model
    for pdb in boltz_results.glob(f"predictions/*/*_model_*.pdb"):
        return str(pdb)
    
    return None


def predict_single(
    input_path: str,
    output_dir: str,
    use_msa_server: bool = True,
    devices: int = 1,
    diffusion_samples: int = 1,
    no_kernels: bool = True,
    single_chain: bool = False,
    timeout: int = 600,
    extra_args: Optional[List[str]] = None,
) -> Optional[str]:
    """Predict structure from a single input file (FASTA or YAML).
    
    Args:
        input_path: Path to input FASTA or YAML file
        output_dir: Directory to output PDB file
        use_msa_server: Whether to use ColabFold MSA server
        devices: Number of GPU devices
        diffusion_samples: Number of diffusion samples
        no_kernels: Disable cuequivariance kernels
        single_chain: Treat input as single chain (VHH)
        timeout: Per-structure timeout in seconds
        extra_args: Additional arguments to pass to boltz predict
    
    Returns:
        Path to output PDB file if successful, None otherwise
    """
    input_path_obj = Path(input_path)
    input_suffix = input_path_obj.suffix.lower()
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine sequences
    if input_suffix in ('.yaml', '.yml'):
        yaml_input = parse_yaml_input(str(input_path_obj))
        # Extract sequences from YAML
        sequences = {}
        for entry in yaml_input.get("sequences", []):
            if "protein" in entry:
                pid = entry["protein"]["id"]
                seq = entry["protein"]["sequence"]
                sequences[pid] = seq
    elif input_suffix in ('.fasta', '.fa'):
        sequences = parse_fasta(str(input_path_obj))
        if single_chain and len(sequences) > 1:
            # Keep only first chain
            first_key = list(sequences.keys())[0]
            sequences = {first_key: sequences[first_key]}
    else:
        print(f"  [ERROR] Unsupported input format: {input_suffix}", file=sys.stderr)
        return None
    
    if not sequences:
        print(f"  [ERROR] No sequences found in {input_path}", file=sys.stderr)
        return None
    
    # Build YAML input
    yaml_stem = input_path_obj.stem
    yaml_path = os.path.join(output_dir, f"{yaml_stem}.yaml")
    build_yaml_from_sequences(sequences, yaml_path, use_msa_server)
    
    print(f"  [Boltz] Input: {input_path_obj.name} -> {len(sequences)} chain(s): {list(sequences.keys())}")
    print(f"  [Boltz] YAML: {yaml_path}")
    
    # Run Boltz
    success = run_boltz_cli(
        yaml_path,
        output_dir,
        use_msa_server=use_msa_server,
        devices=devices,
        diffusion_samples=diffusion_samples,
        no_kernels=no_kernels,
        timeout=timeout,
        extra_args=extra_args,
    )
    
    if not success:
        return None
    
    # Find output PDB
    pdb_path = find_output_pdb(output_dir, yaml_stem)
    
    if pdb_path and validate_pdb(pdb_path):
        return pdb_path
    
    return None


def predict_batch(
    input_dir: str,
    output_dir: str,
    use_msa_server: bool = True,
    devices: int = 1,
    diffusion_samples: int = 1,
    no_kernels: bool = True,
    single_chain: bool = False,
    timeout: int = 600,
    extra_args: Optional[List[str]] = None,
) -> List[str]:
    """Predict structures for all FASTA/YAML files in a directory."""
    input_path = Path(input_dir)
    input_files = sorted(input_path.glob("*.fasta")) + \
                  sorted(input_path.glob("*.fa")) + \
                  sorted(input_path.glob("*.yaml")) + \
                  sorted(input_path.glob("*.yml"))
    
    if not input_files:
        print(f"[WARNING] No input files found in {input_dir}", file=sys.stderr)
        return []
    
    print(f"[Boltz] Found {len(input_files)} input files to predict")
    print(f"[Boltz] Output directory: {output_dir}")
    print(f"[Boltz] MSA server: {use_msa_server}")
    print(f"[Boltz] Devices: {devices}")
    print(f"[Boltz] Diffusion samples: {diffusion_samples}")
    print()
    
    output_files = []
    total_time = 0
    
    for i, input_file in enumerate(input_files, 1):
        print(f"[{i}/{len(input_files)}] Predicting {input_file.name}...")
        
        start_time = time.time()
        
        output_path = predict_single(
            str(input_file),
            output_dir,
            use_msa_server=use_msa_server,
            devices=devices,
            diffusion_samples=diffusion_samples,
            no_kernels=no_kernels,
            single_chain=single_chain,
            timeout=timeout,
            extra_args=extra_args,
        )
        
        elapsed = time.time() - start_time
        total_time += elapsed
        
        if output_path:
            output_files.append(output_path)
            print(f"  ✓ Output: {os.path.basename(output_path)} ({elapsed:.1f}s)")
        else:
            print(f"  ✗ Failed to predict {input_file.name}")
        
        print()
    
    print(f"[Boltz] Completed: {len(output_files)}/{len(input_files)} structures predicted")
    print(f"[Boltz] Total time: {total_time:.1f}s ({total_time/len(input_files):.1f}s per structure)")
    
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description="Boltz-2 Structure Prediction for Humanization Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Predict a single antibody from YAML
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/

    # Predict from FASTA (auto-generates YAML with H/L chains)
    python tools/boltz/boltz_predict.py --input antibody.fasta --output predictions/

    # Predict all variants in a directory
    python tools/boltz/boltz_predict.py --input-dir variants/ --output predictions/

    # Single chain (VHH/nanobody)
    python tools/boltz/boltz_predict.py --input vhh.fasta --output predictions/ --single-chain

    # Offline mode (no MSA server, single-sequence)
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --offline

    # Use multiple devices
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --devices 2

    # Multiple diffusion samples
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --diffusion-samples 3

    # Enable cuequivariance kernels (requires cuequivariance_torch)
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --kernels

    # Change the per-structure timeout (default 600s)
    python tools/boltz/boltz_predict.py --input fv.yaml --output predictions/ --timeout 1800

Requirements:
    - conda environment: boltz (Python >= 3.11, PyTorch >= 2.2)
    - Install: conda create -n boltz python=3.11 && conda run -n boltz pip install "torch==2.5.1+cu121" -f https://download.pytorch.org/whl/cu121 && conda run -n boltz pip install boltz==2.2.1
    - Weights: boltz2_conf.ckpt, boltz2_aff.ckpt, mols.tar in ~/.boltz/ (or $BOLTZ_CACHE)
    - Kernels are disabled by default; install cuequivariance_torch and pass --kernels to enable them.
        """,
    )
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--input",
        type=str,
        help="Input FASTA or YAML file (single prediction)",
    )
    input_group.add_argument(
        "--input-dir",
        type=str,
        help="Input directory containing FASTA/YAML files (batch prediction)",
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
        "--use-msa-server",
        action="store_true",
        default=True,
        help="Use ColabFold MSA server (default: True)",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Disable MSA server (requires pre-computed MSA or single-sequence mode)",
    )
    parser.add_argument(
        "--devices",
        type=int,
        default=1,
        help="Number of GPU devices (default: 1)",
    )
    parser.add_argument(
        "--diffusion-samples",
        type=int,
        default=1,
        help="Number of diffusion samples (default: 1)",
    )
    parser.add_argument(
        "--kernels",
        dest="no_kernels",
        action="store_false",
        help="Enable cuequivariance kernels (requires cuequivariance_torch)",
    )
    parser.add_argument(
        "--no-kernels",
        dest="no_kernels",
        action="store_true",
        help="Disable cuequivariance kernels (default)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="Per-structure timeout in seconds (default: 600)",
    )
    parser.add_argument(
        "--single-chain",
        action="store_true",
        help="Treat input as single chain (VHH/nanobody)",
    )
    parser.add_argument(
        "--extra-args",
        type=str,
        default="",
        help="Extra arguments to pass to 'boltz predict' (quoted string)",
    )
    
    args = parser.parse_args()
    
    # Handle offline flag
    use_msa_server = not args.offline and args.use_msa_server
    if args.offline:
        print("[Boltz] Offline mode - MSA server disabled")
    
    extra_args = args.extra_args.split() if args.extra_args else []
    
    # Ensure output directory exists
    os.makedirs(args.output, exist_ok=True)
    
    # Check Boltz availability
    if not check_boltz_available():
        print("[ERROR] Boltz not installed in current environment", file=sys.stderr)
        print("        Activate boltz conda env: conda activate boltz", file=sys.stderr)
        sys.exit(1)
    
    # Check weights
    if not check_weights_available():
        print("[WARNING] Boltz weights not found in ~/.boltz/ or $BOLTZ_CACHE", file=sys.stderr)
        print("          Required: boltz2_conf.ckpt, boltz2_aff.ckpt, mols/ directory", file=sys.stderr)
    
    # Run prediction
    if args.input:
        output_path = predict_single(
            args.input,
            args.output,
            use_msa_server=use_msa_server,
            devices=args.devices,
            diffusion_samples=args.diffusion_samples,
            no_kernels=args.no_kernels,
            single_chain=args.single_chain,
            timeout=args.timeout,
            extra_args=extra_args,
        )
        
        if output_path:
            print(f"\n✓ Structure predicted: {output_path}")
        else:
            print(f"\n✗ Failed to predict structure", file=sys.stderr)
            sys.exit(1)
    
    elif args.input_dir:
        output_files = predict_batch(
            args.input_dir,
            args.output,
            use_msa_server=use_msa_server,
            devices=args.devices,
            diffusion_samples=args.diffusion_samples,
            no_kernels=args.no_kernels,
            single_chain=args.single_chain,
            timeout=args.timeout,
            extra_args=extra_args,
        )
        
        if not output_files:
            sys.exit(1)


if __name__ == "__main__":
    main()