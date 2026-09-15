#!/usr/bin/env python3
"""Protein Stability Predictor - ΔG and ΔΔG analysis for antibody humanization.

Predicts absolute folding stability (ΔG) and mutational effects (ΔΔG) on
antibody variants. Uses ESM3ΔG and SaProtΔG models (MIT license) for
accurate, free stability prediction.

Features:
  - Absolute stability ΔG prediction for antibody structures
  - Mutational effect ΔΔG prediction for back-mutations
  - Batch analysis for V0-V3 variant panels
  - Integration with humanization pipeline

Usage:
    # Single structure analysis
    python3 tools/stability_predictor/stability_analyzer.py --pdb structure.pdb

    # Compare wildtype vs mutant
    python3 tools/stability_predictor/stability_analyzer.py --wt wildtype.pdb --mutant mutant.pdb

    # Batch analysis for variant panel
    python3 tools/stability_predictor/stability_analyzer.py --pdb-dir variants/ --output results/

    # From humanization pipeline output
    python3 tools/stability_predictor/stability_analyzer.py --pipeline-dir outputs/step3/
"""

import argparse
import json
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


@dataclass
class StabilityResult:
    """Result of stability prediction."""
    pdb_file: str
    chain_id: str = "A"
    dg: Optional[float] = None  # Absolute stability (kcal/mol)
    ddg: Optional[float] = None  # Mutational effect (kcal/mol)
    model: str = "ESM3dG"
    error: Optional[str] = None
    
    @property
    def is_stable(self) -> bool:
        """Check if structure is predicted to be stable."""
        if self.dg is not None:
            return self.dg < 0  # Negative ΔG = stable
        return True
    
    @property
    def is_stabilizing(self) -> bool:
        """Check if mutation is stabilizing."""
        if self.ddg is not None:
            return self.ddg < 0  # Negative ΔΔG = stabilizing
        return True


@dataclass
class BatchStabilityResult:
    """Batch results for multiple structures."""
    results: List[StabilityResult] = field(default_factory=list)
    summary: Dict = field(default_factory=dict)
    
    def add_result(self, result: StabilityResult):
        self.results.append(result)
        self._update_summary()
    
    def _update_summary(self):
        dg_values = [r.dg for r in self.results if r.dg is not None]
        ddg_values = [r.ddg for r in self.results if r.ddg is not None]
        
        self.summary = {
            "n_structures": len(self.results),
            "n_successful": sum(1 for r in self.results if r.error is None),
            "n_failed": sum(1 for r in self.results if r.error is not None),
            "dg_mean": sum(dg_values) / len(dg_values) if dg_values else None,
            "dg_min": min(dg_values) if dg_values else None,
            "dg_max": max(dg_values) if dg_values else None,
            "ddg_mean": sum(ddg_values) / len(ddg_values) if ddg_values else None,
            "ddg_min": min(ddg_values) if ddg_values else None,
            "ddg_max": max(ddg_values) if ddg_values else None,
            "n_stable": sum(1 for r in self.results if r.is_stable and r.error is None),
            "n_stabilizing": sum(1 for r in self.results if r.is_stabilizing and r.error is None),
        }


def check_esm3dg_available() -> bool:
    """Check if ESM3dG is available."""
    try:
        import ESM3dG
        return True
    except ImportError:
        return False


def check_saprodg_available() -> bool:
    """Check if SaProtΔG is available."""
    try:
        import SaProtdG
        return True
    except ImportError:
        return False


def predict_dg_esm3dg(pdb_file: str, chain_id: str = "A") -> StabilityResult:
    """Predict absolute stability ΔG using ESM3dG.
    
    Args:
        pdb_file: Path to PDB file
        chain_id: Chain to analyze (default: A)
    
    Returns:
        StabilityResult with ΔG value
    """
    result = StabilityResult(pdb_file=pdb_file, chain_id=chain_id, model="ESM3dG")
    
    try:
        import ESM3dG
        dg = ESM3dG.ESM3dG_predict(pdb_file)
        result.dg = dg
    except ImportError:
        result.error = "ESM3dG not installed. Run: pip install git+https://github.com/yehlincho/absolute-stability-predictor.git"
    except Exception as e:
        result.error = f"ESM3dG prediction failed: {str(e)}"
    
    return result


def predict_dg_saprodg(pdb_file: str, chain_id: str = "A") -> StabilityResult:
    """Predict absolute stability ΔG using SaProtΔG.
    
    Args:
        pdb_file: Path to PDB file
        chain_id: Chain to analyze (default: A)
    
    Returns:
        StabilityResult with ΔG value
    """
    result = StabilityResult(pdb_file=pdb_file, chain_id=chain_id, model="SaProtΔG")
    
    try:
        import SaProtdG
        dg = SaProtdG.SaProtdG_predict(pdb_file)
        result.dg = dg
    except ImportError:
        result.error = "SaProtΔG not installed. Run: pip install git+https://github.com/yehlincho/absolute-stability-predictor.git"
    except Exception as e:
        result.error = f"SaProtΔG prediction failed: {str(e)}"
    
    return result


def predict_ddg(wt_pdb: str, mutant_pdb: str, model: str = "ESM3dG") -> StabilityResult:
    """Predict mutational effect ΔΔG.
    
    Args:
        wt_pdb: Path to wildtype PDB file
        mutant_pdb: Path to mutant PDB file
        model: Model to use ("ESM3dG" or "SaProtΔG")
    
    Returns:
        StabilityResult with ΔΔG value
    """
    result = StabilityResult(
        pdb_file=mutant_pdb,
        model=model,
        ddg=None
    )
    
    try:
        if model == "ESM3dG":
            import ESM3dG
            ddg = ESM3dG.ESM3dG_predict(wt_pdb, mutant_pdb)
            result.ddg = ddg
        elif model == "SaProtΔG":
            import SaProtdG
            ddg = SaProtdG.SaProtdG_predict(wt_pdb, mutant_pdb)
            result.ddg = ddg
        else:
            result.error = f"Unknown model: {model}"
    except ImportError:
        result.error = f"{model} not installed. Run: pip install git+https://github.com/yehlincho/absolute-stability-predictor.git"
    except Exception as e:
        result.error = f"{model} prediction failed: {str(e)}"
    
    return result


def analyze_batch(pdb_dir: str, model: str = "ESM3dG") -> BatchStabilityResult:
    """Analyze all PDB files in a directory.
    
    Args:
        pdb_dir: Directory containing PDB files
        model: Model to use
    
    Returns:
        BatchStabilityResult with all predictions
    """
    batch_result = BatchStabilityResult()
    
    if not os.path.exists(pdb_dir):
        batch_result.results = [StabilityResult(
            pdb_file=pdb_dir,
            error=f"Directory not found: {pdb_dir}"
        )]
        return batch_result
    
    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    
    if not pdb_files:
        batch_result.results = [StabilityResult(
            pdb_file=pdb_dir,
            error="No PDB files found in directory"
        )]
        return batch_result
    
    for pdb_file in sorted(pdb_files):
        pdb_path = os.path.join(pdb_dir, pdb_file)
        
        if model == "ESM3dG":
            result = predict_dg_esm3dg(pdb_path)
        else:
            result = predict_dg_saprodg(pdb_path)
        
        batch_result.add_result(result)
    
    return batch_result


def format_result_text(result: StabilityResult) -> str:
    """Format single result as text."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  Stability Analysis: {os.path.basename(result.pdb_file)}")
    lines.append(f"{'='*70}")
    lines.append(f"\n  Model: {result.model}")
    lines.append(f"  Chain: {result.chain_id}")
    
    if result.error:
        lines.append(f"\n  ERROR: {result.error}")
    else:
        if result.dg is not None:
            lines.append(f"\n  Absolute Stability (ΔG): {result.dg:.2f} kcal/mol")
            lines.append(f"  Status: {'STABLE' if result.is_stable else 'UNSTABLE'}")
        
        if result.ddg is not None:
            lines.append(f"\n  Mutational Effect (ΔΔG): {result.ddg:.2f} kcal/mol")
            lines.append(f"  Effect: {'STABILIZING' if result.is_stabilizing else 'DESTABILIZING'}")
    
    return "\n".join(lines)


def format_batch_text(batch_result: BatchStabilityResult) -> str:
    """Format batch results as text."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  Batch Stability Analysis Summary")
    lines.append(f"{'='*70}")
    
    s = batch_result.summary
    lines.append(f"\n  Total structures: {s.get('n_structures', 0)}")
    lines.append(f"  Successful: {s.get('n_successful', 0)}")
    lines.append(f"  Failed: {s.get('n_failed', 0)}")
    
    if s.get('dg_mean') is not None:
        lines.append(f"\n  ΔG Statistics:")
        lines.append(f"    Mean: {s['dg_mean']:.2f} kcal/mol")
        lines.append(f"    Min:  {s['dg_min']:.2f} kcal/mol")
        lines.append(f"    Max:  {s['dg_max']:.2f} kcal/mol")
        lines.append(f"    Stable: {s.get('n_stable', 0)}/{s.get('n_structures', 0)}")
    
    if s.get('ddg_mean') is not None:
        lines.append(f"\n  ΔΔG Statistics:")
        lines.append(f"    Mean: {s['ddg_mean']:.2f} kcal/mol")
        lines.append(f"    Min:  {s['ddg_min']:.2f} kcal/mol")
        lines.append(f"    Max:  {s['ddg_max']:.2f} kcal/mol")
        lines.append(f"    Stabilizing: {s.get('n_stabilizing', 0)}/{s.get('n_structures', 0)}")
    
    lines.append(f"\n{'='*70}")
    lines.append(f"  Individual Results")
    lines.append(f"{'='*70}")
    
    for result in batch_result.results:
        lines.append(format_result_text(result))
    
    return "\n".join(lines)


def format_result_json(result: StabilityResult) -> Dict:
    """Format result as JSON-compatible dict."""
    return {
        "pdb_file": result.pdb_file,
        "chain_id": result.chain_id,
        "dg": result.dg,
        "ddg": result.ddg,
        "model": result.model,
        "error": result.error,
        "is_stable": result.is_stable if result.error is None else None,
        "is_stabilizing": result.is_stabilizing if result.error is None else None,
    }


def format_batch_json(batch_result: BatchStabilityResult) -> Dict:
    """Format batch results as JSON-compatible dict."""
    return {
        "summary": batch_result.summary,
        "results": [format_result_json(r) for r in batch_result.results],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Protein Stability Predictor - ΔG and ΔΔG analysis"
    )
    parser.add_argument("--pdb", help="Single PDB file to analyze")
    parser.add_argument("--pdb-dir", help="Directory of PDB files for batch analysis")
    parser.add_argument("--wt", help="Wildtype PDB file (for ΔΔG prediction)")
    parser.add_argument("--mutant", help="Mutant PDB file (for ΔΔG prediction)")
    parser.add_argument("--model", default="ESM3dG", choices=["ESM3dG", "SaProtΔG"],
                        help="Prediction model (default: ESM3dG)")
    parser.add_argument("--chain", default="A", help="Chain ID to analyze (default: A)")
    parser.add_argument("--output", help="Output directory for results")
    parser.add_argument("--json", action="store_true", help="Output JSON format")
    parser.add_argument("--markdown", action="store_true", help="Output Markdown format")
    parser.add_argument("--word", action="store_true", help="Output Word format")
    parser.add_argument("--all-formats", action="store_true", help="Output all formats")
    
    args = parser.parse_args()
    
    if not args.pdb and not args.pdb_dir and not (args.wt and args.mutant):
        parser.error("Provide --pdb, --pdb-dir, or --wt/--mutant")
    
    # Check if models are available
    print("[stability] Checking model availability...")
    esm3d_available = check_esm3dg_available()
    saprodg_available = check_saprodg_available()
    
    print(f"  ESM3dG: {'available' if esm3d_available else 'NOT INSTALLED'}")
    print(f"  SaProtΔG: {'available' if saprodg_available else 'NOT INSTALLED'}")
    
    if not esm3d_available and not saprodg_available:
        print("\n[stability] ERROR: No prediction models available!")
        print("  Install ESM3dG (recommended):")
        print("    pip install git+https://github.com/yehlincho/absolute-stability-predictor.git")
        print("\n  Or install SaProtΔG:")
        print("    pip install git+https://github.com/yehlincho/absolute-stability-predictor.git")
        return 1
    
    # Determine which model to use
    model = args.model
    if model == "ESM3dG" and not esm3d_available:
        if saprodg_available:
            print("[stability] ESM3dG not available, using SaProtΔG")
            model = "SaProtΔG"
        else:
            print("[stability] ERROR: No models available")
            return 1
    elif model == "SaProtΔG" and not saprodg_available:
        if esm3d_available:
            print("[stability] SaProtΔG not available, using ESM3dG")
            model = "ESM3dG"
        else:
            print("[stability] ERROR: No models available")
            return 1
    
    # Run analysis
    if args.wt and args.mutant:
        # ΔΔG prediction
        print(f"\n[stability] Predicting ΔΔG: {os.path.basename(args.wt)} -> {os.path.basename(args.mutant)}")
        result = predict_ddg(args.wt, args.mutant, model)
        
        if args.json:
            print(json.dumps(format_result_json(result), indent=2))
        else:
            print(format_result_text(result))
    
    elif args.pdb:
        # Single structure ΔG prediction
        print(f"\n[stability] Analyzing: {os.path.basename(args.pdb)}")
        
        if model == "ESM3dG":
            result = predict_dg_esm3dg(args.pdb, args.chain)
        else:
            result = predict_dg_saprodg(args.pdb, args.chain)
        
        if args.json:
            print(json.dumps(format_result_json(result), indent=2))
        else:
            print(format_result_text(result))
    
    elif args.pdb_dir:
        # Batch analysis
        print(f"\n[stability] Batch analyzing: {args.pdb_dir}")
        batch_result = analyze_batch(args.pdb_dir, model)
        
        if args.json:
            print(json.dumps(format_batch_json(batch_result), indent=2))
        else:
            print(format_batch_text(batch_result))
    
    # Save results if output directory specified
    if args.output:
        os.makedirs(args.output, exist_ok=True)
        
        if args.pdb:
            result = predict_dg_esm3dg(args.pdb, args.chain) if model == "ESM3dG" else predict_dg_saprodg(args.pdb, args.chain)
            data = format_result_json(result)
            output_file = os.path.join(args.output, "stability_result.json")
        elif args.wt and args.mutant:
            result = predict_ddg(args.wt, args.mutant, model)
            data = format_result_json(result)
            output_file = os.path.join(args.output, "stability_result.json")
        else:
            batch_result = analyze_batch(args.pdb_dir, model)
            data = format_batch_json(batch_result)
            output_file = os.path.join(args.output, "stability_batch_results.json")
        
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)
        print(f"\n[stability] Results saved to: {output_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
