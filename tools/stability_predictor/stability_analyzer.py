#!/usr/bin/env python3
"""Protein Stability Predictor - ΔΔG analysis for antibody humanization using ThermoMPNN.

Predicts mutational effects (ΔΔG) on antibody variants using ThermoMPNN
(Kuhlman Lab, MIT license). Focuses on back-mutation stability assessment
for antibody humanization pipelines.

Features:
  - Mutational effect (ΔΔG) prediction for back-mutations
  - Batch analysis for V0-V3 variant panels
  - Integration with humanization pipeline
  - Kabat/IMGT region annotation for each mutation
  - Multiple output formats (text, JSON, Markdown, Word)

Usage:
    # Single structure analysis with Kabat positions
    python3 tools/stability_predictor/stability_analyzer.py \
        --pdb structure.pdb --chain H \
        --positions 37,44,45,47,48,49,67,69,71,73,76,78,80,82,87,94 \
        --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python

    # Batch analysis for variant panel
    python3 tools/stability_predictor/stability_analyzer.py \
        --pdb-dir variants/ --chain H \
        --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
        --output results/

    # From humanization pipeline output
    python3 tools/stability_predictor/stability_analyzer.py \
        --pdb-dir outputs/step3/af3/ --chain H \
        --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
        --output outputs/stability/ --json

Requirements:
    - ThermoMPNN repo cloned locally
    - Conda environment with ThermoMPNN dependencies (Python 3.10+)
    - Model checkpoint (.ckpt) in ThermoMPNN/models/
"""

import argparse
import json
import os
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


@dataclass
class StabilityResult:
    """Result of stability prediction for a single mutation."""
    pdb_file: str
    chain_id: str = "A"
    position_pdb: int = 0          # Sequential position in PDB chain
    position_kabat: str = ""       # Kabat position (e.g., "H37")
    wildtype: str = ""             # WT amino acid
    mutation: str = ""             # Mutant amino acid
    ddg: Optional[float] = None    # Predicted ΔΔG (kcal/mol)
    region_kabat: str = ""         # Kabat region
    region_imgt: str = ""          # IMGT region
    model: str = "ThermoMPNN"
    error: Optional[str] = None
    
    @property
    def is_stabilizing(self) -> bool:
        """Check if mutation is stabilizing (ΔΔG < 0)."""
        if self.ddg is not None:
            return self.ddg < 0
        return False
    
    @property
    def is_significant(self) -> bool:
        """Check if |ΔΔG| >= 1.0 kcal/mol (significant effect)."""
        if self.ddg is not None:
            return abs(self.ddg) >= 1.0
        return False


@dataclass
class BatchStabilityResult:
    """Batch results for multiple structures."""
    results: List[StabilityResult] = field(default_factory=list)
    summary: Dict = field(default_factory=dict)
    
    def add_result(self, result: StabilityResult):
        self.results.append(result)
        self._update_summary()
    
    def _update_summary(self):
        ddg_values = [r.ddg for r in self.results if r.ddg is not None]
        
        self.summary = {
            "n_mutations": len(self.results),
            "n_successful": sum(1 for r in self.results if r.error is None),
            "n_failed": sum(1 for r in self.results if r.error is not None),
            "ddg_mean": sum(ddg_values) / len(ddg_values) if ddg_values else None,
            "ddg_min": min(ddg_values) if ddg_values else None,
            "ddg_max": max(ddg_values) if ddg_values else None,
            "n_stabilizing": sum(1 for r in self.results if r.is_stabilizing and r.error is None),
            "n_destabilizing": sum(1 for r in self.results if not r.is_stabilizing and r.error is None and r.ddg is not None),
            "n_significant": sum(1 for r in self.results if r.is_significant and r.error is None),
        }


# Import ThermoMPNN analyzer functions
try:
    from .thermompnn_analyzer import (
        analyze_structure,
        analyze_batch,
        StabilityReport,
    )
    THERMOMPNN_AVAILABLE = True
except ImportError:
    analyze_structure = None
    analyze_batch = None
    StabilityReport = None
    THERMOMPNN_AVAILABLE = False


def format_result_text(result: StabilityResult) -> str:
    """Format single result as text."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  Stability Analysis: {os.path.basename(result.pdb_file)} (chain {result.chain_id})")
    lines.append(f"{'='*70}")
    lines.append(f"\n  Model: {result.model}")
    lines.append(f"  Position: {result.position_kabat} (PDB pos {result.position_pdb})")
    lines.append(f"  Mutation: {result.wildtype} -> {result.mutation}")
    lines.append(f"  Kabat Region: {result.region_kabat or '-'}")
    lines.append(f"  IMGT Region: {result.region_imgt or '-'}")
    
    if result.error:
        lines.append(f"\n  ERROR: {result.error}")
    else:
        if result.ddg is not None:
            lines.append(f"\n  Mutational Effect (ΔΔG): {result.ddg:.2f} kcal/mol")
            effect = "STABILIZING" if result.is_stabilizing else "DESTABILIZING"
            signif = " (significant, |ΔΔG| ≥ 1.0)" if result.is_significant else ""
            lines.append(f"  Effect: {effect}{signif}")
    
    return "\n".join(lines)


def format_batch_text(batch_result: BatchStabilityResult) -> str:
    """Format batch results as text."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  Batch Stability Analysis Summary (ThermoMPNN)")
    lines.append(f"{'='*70}")
    
    s = batch_result.summary
    lines.append(f"\n  Total mutations: {s.get('n_mutations', 0)}")
    lines.append(f"  Successful: {s.get('n_successful', 0)}")
    lines.append(f"  Failed: {s.get('n_failed', 0)}")
    
    if s.get('ddg_mean') is not None:
        lines.append(f"\n  ΔΔG Statistics:")
        lines.append(f"    Mean: {s['ddg_mean']:.2f} kcal/mol")
        lines.append(f"    Min:  {s['ddg_min']:.2f} kcal/mol")
        lines.append(f"    Max:  {s['ddg_max']:.2f} kcal/mol")
        lines.append(f"    Stabilizing (ΔΔG < 0):   {s.get('n_stabilizing', 0)}")
        lines.append(f"    Destabilizing (ΔΔG > 0): {s.get('n_destabilizing', 0)}")
        lines.append(f"    Significant (|ΔΔG| ≥ 1.0): {s.get('n_significant', 0)}")
    
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
        "position_pdb": result.position_pdb,
        "position_kabat": result.position_kabat,
        "wildtype": result.wildtype,
        "mutation": result.mutation,
        "ddg": result.ddg,
        "region_kabat": result.region_kabat,
        "region_imgt": result.region_imgt,
        "model": result.model,
        "error": result.error,
        "is_stabilizing": result.is_stabilizing if result.error is None else None,
        "is_significant": result.is_significant if result.error is None else None,
    }


def format_batch_json(batch_result: BatchStabilityResult) -> Dict:
    """Format batch results as JSON-compatible dict."""
    return {
        "summary": batch_result.summary,
        "results": [format_result_json(r) for r in batch_result.results],
    }


def main():
    parser = argparse.ArgumentParser(
        description="Protein Stability Predictor - ΔΔG analysis with ThermoMPNN"
    )
    parser.add_argument("--pdb", help="Single PDB/CIF file to analyze")
    parser.add_argument("--pdb-dir", help="Directory of PDB/CIF files for batch analysis")
    parser.add_argument("--chain", default="H", help="Chain ID to analyze (default: H)")
    parser.add_argument("--positions", help="Comma-separated Kabat positions (e.g., '37,44,45' or 'H37,H44')")
    parser.add_argument("--thermompnn-dir", required=True, help="Path to ThermoMPNN repo (required)")
    parser.add_argument("--python", required=True, 
                        help="Python binary for ThermoMPNN env (e.g., '~/miniconda3/envs/thermoMPNN/bin/python' or 'conda run -n thermoMPNN python')")
    parser.add_argument("--output", help="Output directory for results")
    parser.add_argument("--json", action="store_true", help="Output JSON format")
    parser.add_argument("--markdown", action="store_true", help="Output Markdown format")
    parser.add_argument("--word", action="store_true", help="Output Word format")
    parser.add_argument("--all-formats", action="store_true", help="Output all formats")
    parser.add_argument("--verbose", action="store_true", help="Verbose ThermoMPNN output")
    
    args = parser.parse_args()
    
    if not args.pdb and not args.pdb_dir:
        parser.error("Provide --pdb or --pdb-dir")
    
    if not THERMOMPNN_AVAILABLE or analyze_structure is None or analyze_batch is None:
        print("[stability] ERROR: ThermoMPNN analyzer not available")
        return 1
    
    # Parse positions
    positions = None
    if args.positions:
        for p in args.positions.split(','):
            p = p.strip().upper()
            if not p:
                continue
            if not p.startswith(('H', 'L')):
                p = 'H' + p
            if positions is None:
                positions = []
            positions.append(p)
    
    # Handle conda run syntax
    python_bin = args.python
    
    # Output directory
    if args.output:
        os.makedirs(args.output, exist_ok=True)
    
    # Run analysis
    if args.pdb:
        struct_out_dir = args.output or "outputs/stability"
        report = analyze_structure(
            args.pdb, args.chain, args.thermompnn_dir, python_bin,
            struct_out_dir, positions, args.verbose
        )
        
        # Convert to StabilityResult list for batch formatting
        batch_result = BatchStabilityResult()
        for mut in report.mutations:
            batch_result.add_result(StabilityResult(
                pdb_file=report.pdb_file,
                chain_id=args.chain,
                position_pdb=mut.position_pdb,
                position_kabat=mut.position_kabat,
                wildtype=mut.wildtype,
                mutation=mut.mutation,
                ddg=mut.ddG,
                region_kabat=mut.region_kabat,
                region_imgt=mut.region_imgt,
                model="ThermoMPNN",
                error=report.error,
            ))
        
        if args.json:
            print(json.dumps(format_batch_json(batch_result), indent=2))
        else:
            print(format_batch_text(batch_result))
        
        # Save results
        if args.output:
            data = format_batch_json(batch_result)
            output_file = os.path.join(args.output, "stability_batch_results.json")
            with open(output_file, "w") as f:
                json.dump(data, f, indent=2)
            print(f"\n[stability] Results saved to: {output_file}")
            
            # Also save markdown/word if requested
            if args.markdown or args.all_formats:
                from .thermompnn_analyzer import save_markdown_report
                save_markdown_report(report, os.path.join(args.output, f"{Path(args.pdb).stem}_stability.md"))
            if args.word or args.all_formats:
                from .thermompnn_analyzer import save_word_report
                save_word_report(report, os.path.join(args.output, f"{Path(args.pdb).stem}_stability.docx"))
    
    elif args.pdb_dir:
        struct_out_dir = args.output or "outputs/stability"
        reports = analyze_batch(
            args.pdb_dir, args.chain, args.thermompnn_dir, python_bin,
            struct_out_dir, positions, args.verbose
        )
        
        # Convert to batch result
        batch_result = BatchStabilityResult()
        for report in reports:
            for mut in report.mutations:
                batch_result.add_result(StabilityResult(
                    pdb_file=report.pdb_file,
                    chain_id=args.chain,
                    position_pdb=mut.position_pdb,
                    position_kabat=mut.position_kabat,
                    wildtype=mut.wildtype,
                    mutation=mut.mutation,
                    ddg=mut.ddG,
                    region_kabat=mut.region_kabat,
                    region_imgt=mut.region_imgt,
                    model="ThermoMPNN",
                    error=report.error,
                ))
        
        if args.json:
            print(json.dumps(format_batch_json(batch_result), indent=2))
        else:
            print(format_batch_text(batch_result))
        
        # Save results
        if args.output:
            data = format_batch_json(batch_result)
            output_file = os.path.join(args.output, "stability_batch_results.json")
            with open(output_file, "w") as f:
                json.dump(data, f, indent=2)
            print(f"\n[stability] Results saved to: {output_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())