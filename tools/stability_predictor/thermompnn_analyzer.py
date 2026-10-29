#!/usr/bin/env python3
"""ThermoMPNN Stability Analyzer - ΔΔG prediction for antibody back-mutations.

Uses ThermoMPNN (Kuhlman Lab, MIT license) to predict mutational effects
on protein stability. Designed for evaluating back-mutation candidates
in antibody humanization pipelines.

Requirements:
    - ThermoMPNN repo cloned locally
    - Conda environment with ThermoMPNN dependencies (Python 3.10+)
    - Model checkpoint (.ckpt) available

Usage:
    # Single structure with specific Kabat positions
    python3 tools/stability_predictor/thermompnn_analyzer.py \
        --pdb structure.pdb --chain H \
        --positions 37,44,45,47,48,49,67,69,71,73,76,78,80,82,87,94 \
        --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
        --output outputs/stability/

    # Batch analysis of variant directory
    python3 tools/stability_predictor/thermompnn_analyzer.py \
        --pdb-dir variants/ --chain H \
        --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
        --output outputs/stability/ --json

Output:
    - JSON: detailed per-mutation ΔΔG with Kabat position, region, effect
    - Text: formatted report with stabilizing/destabilizing classification
    - Markdown: table format for documentation
    - Word: formatted document with color-coded risk levels
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Repo root for optional Kabat region annotation
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class MutationResult:
    """Result for a single mutation."""
    pdb_file: str
    chain: str
    position_pdb: int          # 1-based sequential position in PDB chain
    position_kabat: str        # Kabat position (e.g., "H37", "L50")
    wildtype: str              # WT amino acid (1-letter)
    mutation: str              # Mutant amino acid (1-letter)
    ddG: float                 # Predicted ΔΔG (kcal/mol)
    region_kabat: str          # Kabat region (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4)
    region_imgt: str           # IMGT region
    is_stabilizing: bool       # True if ddG < 0
    is_significant: bool       # True if |ddG| >= 1.0 kcal/mol
    error: Optional[str] = None


@dataclass
class StabilityReport:
    """Aggregated stability report for a structure."""
    pdb_file: str
    chain: str
    mutations: List[MutationResult] = field(default_factory=list)
    summary: Dict = field(default_factory=dict)
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# ThermoMPNN Integration
# ---------------------------------------------------------------------------

def check_thermompnn_installation(thermompnn_dir: str, python_bin: str) -> Tuple[bool, Optional[str]]:
    """Verify ThermoMPNN installation has required files."""
    tm_dir = Path(thermompnn_dir).expanduser().resolve()
    
    # Check repo structure
    required_files = [
        "analysis/custom_inference.py",
        "local.yaml",
    ]
    
    for f in required_files:
        if not (tm_dir / f).exists():
            return False, f"Missing {f} in {tm_dir}"
    
    # Check for checkpoint
    ckpt_candidates = list(tm_dir.glob("models/*.ckpt")) + list(tm_dir.glob("models/*.pt")) + \
                       list(tm_dir.glob("vanilla_model_weights/*.pt"))
    if not ckpt_candidates:
        return False, f"No model checkpoint found in {tm_dir}/models/ or {tm_dir}/vanilla_model_weights/"
    
    # Check Python environment
    try:
        result = subprocess.run(
            [python_bin, "-c", "import torch; import pytorch_lightning; import omegaconf; print('OK')"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return False, f"Python env check failed: {result.stderr}"
    except subprocess.TimeoutExpired:
        return False, "Python env check timed out"
    except FileNotFoundError:
        return False, f"Python binary not found: {python_bin}"
    
    return True, None


def run_thermompnn_inference(
    pdb_path: str,
    chain: str,
    thermompnn_dir: str,
    python_bin: str,
    output_dir: str,
    positions_kabat: Optional[List[str]] = None,
    verbose: bool = False,
) -> Tuple[List[MutationResult], Optional[str]]:
    """Run ThermoMPNN custom_inference.py and parse results.
    
    Args:
        pdb_path: Path to input PDB file
        chain: Chain ID to analyze (e.g., 'H', 'L', 'A', 'B')
        thermompnn_dir: Path to ThermoMPNN repo
        python_bin: Python binary for ThermoMPNN conda env
        output_dir: Output directory for ThermoMPNN CSV
        positions_kabat: Optional list of Kabat positions to filter (e.g., ['H37', 'H44'])
        verbose: Print subprocess output
    
    Returns:
        List of MutationResult objects, error message if any
    """
    tm_dir = Path(thermompnn_dir).expanduser().resolve()
    pdb_path_obj = Path(pdb_path).resolve()
    output_dir_obj = Path(output_dir).resolve()
    output_dir_obj.mkdir(parents=True, exist_ok=True)
    
    # Find checkpoint
    ckpt_candidates = list(tm_dir.glob("models/*.ckpt")) + list(tm_dir.glob("models/*.pt")) + \
                       list(tm_dir.glob("vanilla_model_weights/*.pt"))
    model_path = str(ckpt_candidates[0])
    
    # Run ThermoMPNN custom_inference.py
    # It generates a CSV with all single-point mutations
    script_path = tm_dir / "analysis" / "custom_inference.py"
    
    # Split python_bin if it's a command like "conda run -n env python"
    python_parts = python_bin.split() if " " in python_bin else [python_bin]
    cmd = python_parts + [
        str(script_path),
        "--pdb", str(pdb_path),
        "--chain", chain,
        "--model_path", model_path,
        "--out_dir", str(output_dir),
    ]
    
    if verbose:
        print(f"[ThermoMPNN] Running: {' '.join(cmd)}", file=sys.stderr)
    
    try:
        result = subprocess.run(cmd, capture_output=not verbose, text=True, timeout=300)
    except subprocess.TimeoutExpired:
        return [], "ThermoMPNN inference timed out (>5 min)"
    except Exception as e:
        return [], f"Failed to run ThermoMPNN: {e}"
    
    if result.returncode != 0:
        return [], f"ThermoMPNN failed (exit {result.returncode}): {result.stderr}"
    
    # Find output CSV
    pdb_stem = pdb_path_obj.stem
    csv_pattern = f"ThermoMPNN_inference_{pdb_stem}.csv"
    csv_files = list(output_dir_obj.glob(csv_pattern))
    if not csv_files:
        # Try glob
        csv_files = list(output_dir_obj.glob("*.csv"))
    if not csv_files:
        return [], f"No output CSV found in {output_dir_obj}"
    
    csv_path = csv_files[0]
    
    # Parse CSV
    mutations = parse_thermompnn_csv(csv_path, pdb_path_obj.name, chain, positions_kabat)
    
    return mutations, None


def parse_thermompnn_csv(
    csv_path: Path,
    pdb_name: str,
    chain: str,
    positions_kabat: Optional[List[str]] = None,
) -> List[MutationResult]:
    """Parse ThermoMPNN output CSV into MutationResult objects.
    
    CSV columns: Model, Dataset, ddG_pred, position, wildtype, mutation, pdb, chain
    """
    mutations = []
    position_filter = set(positions_kabat) if positions_kabat else None
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                pos_pdb = int(row['position'])
                wt = row['wildtype']
                mut = row['mutation']
                ddg = float(row['ddG_pred'])
                
                # Build Kabat position string (e.g., "H37")
                pos_kabat = f"{chain}{pos_pdb}"
                
                # Filter by requested Kabat positions
                if position_filter and pos_kabat not in position_filter:
                    continue
                
                mutations.append(MutationResult(
                    pdb_file=pdb_name,
                    chain=chain,
                    position_pdb=pos_pdb,
                    position_kabat=pos_kabat,
                    wildtype=wt,
                    mutation=mut,
                    ddG=ddg,
                    region_kabat="",  # Will be filled by annotate_regions
                    region_imgt="",
                    is_stabilizing=ddg < 0,
                    is_significant=abs(ddg) >= 1.0,
                ))
            except (KeyError, ValueError) as e:
                continue  # Skip malformed rows
    
    return mutations


def annotate_regions(mutations: List[MutationResult], sequence: str, chain_type: str) -> List[MutationResult]:
    """Add Kabat and IMGT region annotations to mutations using pipeline numbering."""
    if not mutations:
        return mutations
    
    try:
        from scripts.humanize.numbering import number_heavy, number_light
        from scripts.humanize.anarci_adapter import number_with_anarci_imgt
    except ImportError:
        return mutations
    
    # Kabat numbering
    kabat_chain = None
    try:
        if chain_type == "H":
            kabat_chain = number_heavy(sequence)
        else:
            kabat_chain = number_light(sequence)
    except Exception:
        pass
    
    # IMGT numbering
    imgt_chain = None
    try:
        imgt_chain = number_with_anarci_imgt(sequence, chain_type)
    except Exception:
        pass
    
    # Build position -> region maps
    kabat_regions = {}
    if kabat_chain:
        for r in kabat_chain.residues:
            if 0 <= r.index < len(sequence):
                kabat_regions[r.index + 1] = r.region
    
    imgt_regions = {}
    if imgt_chain:
        for r in imgt_chain.residues:
            if 0 <= r.index < len(sequence):
                imgt_regions[r.index + 1] = r.region
    
    # Annotate mutations (ThermoMPNN uses 1-based sequential PDB position)
    for mut in mutations:
        mut.region_kabat = kabat_regions.get(mut.position_pdb, "-")
        mut.region_imgt = imgt_regions.get(mut.position_pdb, "-")
    
    return mutations


def extract_chain_sequence(pdb_path: str, chain: str) -> str:
    """Extract amino acid sequence for a chain from PDB."""
    from Bio.PDB.PDBParser import PDBParser
    
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V',
        'TRP': 'W', 'TYR': 'Y', 'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
        'XLE': 'J', 'XAA': 'X'
    }
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('input', pdb_path)
    
    for model in structure:
        for ch in model:
            if ch.id == chain:
                seq = []
                for res in ch:
                    if res.id[0] == ' ' and res.resname in three_to_one:  # Standard residue
                        seq.append(three_to_one[res.resname])
                return ''.join(seq)
    return ""


# ---------------------------------------------------------------------------
# Output Formatting
# ---------------------------------------------------------------------------

def format_text_report(report: StabilityReport) -> str:
    """Format stability report as text."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"  ThermoMPNN Stability Analysis: {report.pdb_file} (chain {report.chain})")
    lines.append("=" * 80)
    
    if report.error:
        lines.append(f"\n  ERROR: {report.error}")
        return "\n".join(lines)
    
    if not report.mutations:
        lines.append("\n  No mutations analyzed (check --positions filter or structure).")
        return "\n".join(lines)
    
    # Table header
    lines.append(f"\n  {'Kabat':<8} {'WT':<4} {'Mut':<4} {'ΔΔG':>8} {'Effect':<14} {'Signif.':<10} {'Kabat Reg':<10} {'IMGT Reg':<10}")
    lines.append("  " + "-" * 76)
    
    for mut in sorted(report.mutations, key=lambda m: m.position_pdb):
        effect = "STABILIZING" if mut.is_stabilizing else "DESTABILIZING"
        signif = "YES" if mut.is_significant else "no"
        marker = " *" if mut.is_significant else ""
        lines.append(
            f"  {mut.position_kabat:<8} {mut.wildtype:<4} {mut.mutation:<4} "
            f"{mut.ddG:>8.2f} {effect:<14} {signif:<10} {mut.region_kabat:<10} {mut.region_imgt:<10}{marker}"
        )
    
    # Summary
    s = report.summary
    lines.append(f"\n{'='*80}")
    lines.append("  SUMMARY")
    lines.append(f"{'='*80}")
    lines.append(f"  Total mutations analyzed:  {s.get('total', 0)}")
    lines.append(f"  Stabilizing (ΔΔG < 0):     {s.get('stabilizing', 0)}")
    lines.append(f"  Destabilizing (ΔΔG > 0):   {s.get('destabilizing', 0)}")
    lines.append(f"  Significant (|ΔΔG| ≥ 1.0): {s.get('significant', 0)}")
    if s.get('mean_ddg') is not None:
        lines.append(f"  Mean ΔΔG:                  {s['mean_ddg']:.2f} kcal/mol")
        lines.append(f"  Min ΔΔG:                   {s['min_ddg']:.2f} kcal/mol")
        lines.append(f"  Max ΔΔG:                   {s['max_ddg']:.2f} kcal/mol")
    
    lines.append(f"\n  Thresholds: |ΔΔG| ≥ 1.0 kcal/mol = significant; ΔΔG < 0 = stabilizing")
    lines.append(f"  Model: ThermoMPNN (Kuhlman Lab, MIT license)")
    
    return "\n".join(lines)


def save_json_report(report: StabilityReport, output_path: str):
    """Save report as JSON."""
    data = {
        "pdb_file": report.pdb_file,
        "chain": report.chain,
        "mutations": [asdict(m) for m in report.mutations],
        "summary": report.summary,
        "error": report.error,
    }
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def save_markdown_report(report: StabilityReport, output_path: str):
    """Save report as Markdown."""
    lines = []
    lines.append(f"# ThermoMPNN Stability Analysis: {report.pdb_file}\n")
    lines.append(f"**Chain:** {report.chain}  \n")
    
    if report.error:
        lines.append(f"\n**ERROR:** {report.error}\n")
    elif not report.mutations:
        lines.append("\nNo mutations analyzed.\n")
    else:
        lines.append("| Kabat | WT | Mut | ΔΔG (kcal/mol) | Effect | Significant | Kabat Region | IMGT Region |")
        lines.append("|-------|-----|-----|----------------|--------|-------------|--------------|-------------|")
        
        for mut in sorted(report.mutations, key=lambda m: m.position_pdb):
            effect = "**STABILIZING**" if mut.is_stabilizing else "DESTABILIZING"
            signif = "✓" if mut.is_significant else ""
            lines.append(
                f"| {mut.position_kabat} | {mut.wildtype} | {mut.mutation} | "
                f"{mut.ddG:.2f} | {effect} | {signif} | {mut.region_kabat} | {mut.region_imgt} |"
            )
    
    # Summary
    s = report.summary
    lines.append(f"\n## Summary\n")
    lines.append(f"- **Total mutations:** {s.get('total', 0)}")
    lines.append(f"- **Stabilizing (ΔΔG < 0):** {s.get('stabilizing', 0)}")
    lines.append(f"- **Destabilizing (ΔΔG > 0):** {s.get('destabilizing', 0)}")
    lines.append(f"- **Significant (|ΔΔG| ≥ 1.0):** {s.get('significant', 0)}")
    if s.get('mean_ddg') is not None:
        lines.append(f"- **Mean ΔΔG:** {s['mean_ddg']:.2f} kcal/mol")
        lines.append(f"- **Min ΔΔG:** {s['min_ddg']:.2f} kcal/mol")
        lines.append(f"- **Max ΔΔG:** {s['max_ddg']:.2f} kcal/mol")
    
    with open(output_path, 'w') as f:
        f.write("\n".join(lines))


def save_word_report(report: StabilityReport, output_path: str):
    """Save report as Word document."""
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
    heading = doc.add_heading(f'ThermoMPNN Stability Analysis: {report.pdb_file}', 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER
    
    doc.add_paragraph(f"Chain: {report.chain}", style='Intense Quote')
    
    if report.error:
        doc.add_paragraph(f"ERROR: {report.error}", style='Intense Quote')
        doc.save(output_path)
        return
    
    if not report.mutations:
        doc.add_paragraph("No mutations analyzed.")
        doc.save(output_path)
        return
    
    # Mutations table
    table = doc.add_table(rows=1, cols=8)
    table.style = 'Light Grid Accent 1'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    
    headers = ['Kabat', 'WT', 'Mut', 'ΔΔG (kcal/mol)', 'Effect', 'Significant', 'Kabat Region', 'IMGT Region']
    for i, header in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = header
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.bold = True
                run.font.size = Pt(9)
    
    for mut in sorted(report.mutations, key=lambda m: m.position_pdb):
        row = table.add_row()
        row.cells[0].text = mut.position_kabat
        row.cells[1].text = mut.wildtype
        row.cells[2].text = mut.mutation
        row.cells[3].text = f"{mut.ddG:.2f}"
        row.cells[4].text = "STABILIZING" if mut.is_stabilizing else "DESTABILIZING"
        row.cells[5].text = "✓" if mut.is_significant else ""
        row.cells[6].text = mut.region_kabat
        row.cells[7].text = mut.region_imgt
        
        # Color code
        risk_cell = row.cells[4]
        for paragraph in risk_cell.paragraphs:
            for run in paragraph.runs:
                if mut.is_significant and mut.is_stabilizing:
                    run.font.color.rgb = RGBColor(0, 128, 0)  # Green
                    run.font.bold = True
                elif mut.is_significant and not mut.is_stabilizing:
                    run.font.color.rgb = RGBColor(255, 0, 0)  # Red
                    run.font.bold = True
                run.font.size = Pt(9)
        
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.size = Pt(9)
    
    # Summary
    doc.add_heading('Summary', level=1)
    s = report.summary
    summary_table = doc.add_table(rows=5, cols=2)
    summary_table.style = 'Light Grid Accent 1'
    
    summary_data = [
        ("Total mutations", s.get('total', 0)),
        ("Stabilizing (ΔΔG < 0)", s.get('stabilizing', 0)),
        ("Destabilizing (ΔΔG > 0)", s.get('destabilizing', 0)),
        ("Significant (|ΔΔG| ≥ 1.0)", s.get('significant', 0)),
    ]
    
    for i, (label, value) in enumerate(summary_data):
        summary_table.rows[i].cells[0].text = label
        summary_table.rows[i].cells[1].text = str(value)
    
    if s.get('mean_ddg') is not None:
        row = summary_table.add_row()
        row.cells[0].text = "Mean ΔΔG"
        row.cells[1].text = f"{s['mean_ddg']:.2f} kcal/mol"
        row = summary_table.add_row()
        row.cells[0].text = "Min ΔΔG"
        row.cells[1].text = f"{s['min_ddg']:.2f} kcal/mol"
        row = summary_table.add_row()
        row.cells[0].text = "Max ΔΔG"
        row.cells[1].text = f"{s['max_ddg']:.2f} kcal/mol"
    
    doc.add_paragraph()
    doc.add_paragraph("Thresholds: |ΔΔG| ≥ 1.0 kcal/mol = significant; ΔΔG < 0 = stabilizing", style='Intense Quote')
    doc.add_paragraph("Model: ThermoMPNN (Kuhlman Lab, MIT license)", style='Intense Quote')
    
    doc.save(output_path)


# ---------------------------------------------------------------------------
# Main Analysis Function
# ---------------------------------------------------------------------------

def analyze_structure(
    pdb_path: str,
    chain: str,
    thermompnn_dir: str,
    python_bin: str,
    output_dir: str,
    positions: Optional[List[str]] = None,
    verbose: bool = False,
) -> StabilityReport:
    """Analyze a single structure for back-mutation stability effects."""
    report = StabilityReport(pdb_file=os.path.basename(pdb_path), chain=chain)
    
    # Check installation
    ok, err = check_thermompnn_installation(thermompnn_dir, python_bin)
    if not ok:
        report.error = f"ThermoMPNN not properly installed: {err}"
        report.summary = {"total": 0, "stabilizing": 0, "destabilizing": 0, "significant": 0}
        return report
    
    # Extract sequence for region annotation
    sequence = extract_chain_sequence(pdb_path, chain)
    chain_type = chain.upper() if chain.upper() in ("H", "L") else "H"
    
    # Run ThermoMPNN
    mutations, err = run_thermompnn_inference(
        pdb_path, chain, thermompnn_dir, python_bin, output_dir,
        positions_kabat=positions, verbose=verbose
    )
    
    if err:
        report.error = err
        report.summary = {"total": 0, "stabilizing": 0, "destabilizing": 0, "significant": 0}
        return report
    
    # Annotate regions
    if sequence:
        mutations = annotate_regions(mutations, sequence, chain_type)
    
    report.mutations = mutations
    
    # Build summary
    total = len(mutations)
    stabilizing = sum(1 for m in mutations if m.is_stabilizing)
    destabilizing = total - stabilizing
    significant = sum(1 for m in mutations if m.is_significant)
    ddg_values = [m.ddG for m in mutations]
    
    report.summary = {
        "total": total,
        "stabilizing": stabilizing,
        "destabilizing": destabilizing,
        "significant": significant,
        "mean_ddg": sum(ddg_values) / len(ddg_values) if ddg_values else None,
        "min_ddg": min(ddg_values) if ddg_values else None,
        "max_ddg": max(ddg_values) if ddg_values else None,
    }
    
    return report


def analyze_batch(
    pdb_dir: str,
    chain: str,
    thermompnn_dir: str,
    python_bin: str,
    output_dir: str,
    positions: Optional[List[str]] = None,
    verbose: bool = False,
) -> List[StabilityReport]:
    """Analyze all PDB files in a directory."""
    reports = []
    pdb_dir_path = Path(pdb_dir)
    
    if not pdb_dir_path.exists():
        return [StabilityReport(pdb_file=pdb_dir, chain=chain, error=f"Directory not found: {pdb_dir}")]
    
    pdb_files = sorted(list(pdb_dir_path.rglob("*.pdb")) + list(pdb_dir_path.rglob("*.cif")))
    
    if not pdb_files:
        return [StabilityReport(pdb_file=str(pdb_dir), chain=chain, error="No PDB/CIF files found")]
    
    for pdb_file in pdb_files:
        print(f"\n[ThermoMPNN] Analyzing: {pdb_file}")
        # Create subdirectory per structure for ThermoMPNN output
        struct_out_dir = Path(output_dir) / "thermompnn_raw" / pdb_file.stem
        report = analyze_structure(
            str(pdb_file), chain, thermompnn_dir, python_bin,
            str(struct_out_dir), positions, verbose
        )
        reports.append(report)
    
    return reports


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_positions(pos_str: str) -> List[str]:
    """Parse comma-separated position string (e.g., 'H37,H44,H45' or '37,44,45')."""
    positions = []
    for p in pos_str.split(','):
        p = p.strip().upper()
        if not p:
            continue
        # If no chain prefix, assume H (heavy chain)
        if not re.match(r'^[HL]', p):
            p = 'H' + p
        positions.append(p)
    return positions


def main():
    parser = argparse.ArgumentParser(
        description="ThermoMPNN Stability Analyzer - ΔΔG prediction for antibody back-mutations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single structure with Kabat positions
  python3 tools/stability_predictor/thermompnn_analyzer.py \
      --pdb structure.pdb --chain H \
      --positions 37,44,45,47,48,49,67,69,71,73,76,78,80,82,87,94 \
      --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
      --output outputs/stability/

  # Batch analysis of AF3 variant directory
  python3 tools/stability_predictor/thermompnn_analyzer.py \
      --pdb-dir outputs/step3/af3/ --chain H \
      --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
      --output outputs/stability/ --json

  # Using conda run directly
  python3 tools/stability_predictor/thermompnn_analyzer.py \
      --pdb structure.pdb --chain H \
      --thermompnn-dir ~/ThermoMPNN --python "conda run -n thermoMPNN python" \
      --output outputs/stability/
        """
    )
    parser.add_argument("--pdb", help="Single PDB/CIF file to analyze")
    parser.add_argument("--pdb-dir", help="Directory of PDB/CIF files for batch analysis")
    parser.add_argument("--chain", default="H", help="Chain ID to analyze (default: H)")
    parser.add_argument("--positions", help="Comma-separated Kabat positions (e.g., '37,44,45' or 'H37,H44')")
    parser.add_argument("--thermompnn-dir", required=True, help="Path to ThermoMPNN repo (required)")
    parser.add_argument("--python", required=True, help="Python binary for ThermoMPNN env (e.g., '~/miniconda3/envs/thermoMPNN/bin/python' or 'conda run -n thermoMPNN python')")
    parser.add_argument("--output", help="Output directory for reports")
    parser.add_argument("--json", action="store_true", help="Save JSON report")
    parser.add_argument("--markdown", action="store_true", help="Save Markdown report")
    parser.add_argument("--word", action="store_true", help="Save Word report")
    parser.add_argument("--all-formats", action="store_true", help="Save all formats")
    parser.add_argument("--verbose", action="store_true", help="Verbose ThermoMPNN output")
    
    args = parser.parse_args()
    
    if not args.pdb and not args.pdb_dir:
        parser.error("Either --pdb or --pdb-dir must be specified")
    
    # Parse positions
    positions = parse_positions(args.positions) if args.positions else None
    
    # Handle conda run syntax
    python_bin = args.python
    if python_bin.startswith("conda run"):
        # Use as-is, will be split in subprocess
        pass
    
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
        
        # Print text report
        print(format_text_report(report))
        
        # Save reports
        if args.output and (args.json or args.markdown or args.word or args.all_formats):
            stem = Path(args.pdb).stem
            if args.json or args.all_formats:
                save_json_report(report, os.path.join(args.output, f"{stem}_stability.json"))
            if args.markdown or args.all_formats:
                save_markdown_report(report, os.path.join(args.output, f"{stem}_stability.md"))
            if args.word or args.all_formats:
                save_word_report(report, os.path.join(args.output, f"{stem}_stability.docx"))
    
    elif args.pdb_dir:
        reports = analyze_batch(
            args.pdb_dir, args.chain, args.thermompnn_dir, python_bin,
            args.output or "outputs/stability", positions, args.verbose
        )
        
        # Print and save each
        for report in reports:
            print(format_text_report(report))
            
            if args.output and (args.json or args.markdown or args.word or args.all_formats):
                stem = Path(report.pdb_file).stem
                if args.json or args.all_formats:
                    save_json_report(report, os.path.join(args.output, f"{stem}_stability.json"))
                if args.markdown or args.all_formats:
                    save_markdown_report(report, os.path.join(args.output, f"{stem}_stability.md"))
                if args.word or args.all_formats:
                    save_word_report(report, os.path.join(args.output, f"{stem}_stability.docx"))
        
        # Combined JSON
        if len(reports) > 1 and (args.json or args.all_formats):
            combined = {
                "reports": [
                    {
                        "pdb_file": r.pdb_file,
                        "chain": r.chain,
                        "mutations": [asdict(m) for m in r.mutations],
                        "summary": r.summary,
                        "error": r.error,
                    }
                    for r in reports
                ]
            }
            combined_path = os.path.join(args.output, "combined_stability_results.json")
            with open(combined_path, 'w') as f:
                json.dump(combined, f, indent=2, ensure_ascii=False)
            print(f"\nCombined JSON saved: {combined_path}")


if __name__ == "__main__":
    main()