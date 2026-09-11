#!/usr/bin/env python3
"""PTM Exposure Analysis Tool

Analyze high-risk PTM (Post-Translational Modification) sites in antibody structures
based on AF3/PDB files. Calculates actual solvent exposure (SASA) to assess true
risk levels, going beyond sequence-only predictions.

Usage:
    python3 tools/ptm_exposure/analyze_ptm.py --pdb <structure.pdb> [--output <dir>]
    python3 tools/ptm_exposure/analyze_ptm.py --pdb-dir <af3_dir> [--output <dir>]

Requirements:
    - freesasa (pip install freesasa)
    - python-docx (optional, for Word reports)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import freesasa
except ImportError:
    print("ERROR: freesasa not installed. Run: pip install freesasa", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# PTM Risk Patterns (from scripts/humanize/config.py)
# ---------------------------------------------------------------------------

PTM_PATTERNS = {
    # Deamidation (Asn -> Asp/isoAsp)
    "deamidation (NG)": {"pattern": r"NG", "risk": "high", "exposure_threshold": 0.30},
    "deamidation (NS)": {"pattern": r"NS", "risk": "medium", "exposure_threshold": 0.40},
    "deamidation (NH)": {"pattern": r"NH", "risk": "medium", "exposure_threshold": 0.40},
    "deamidation (ND)": {"pattern": r"ND", "risk": "low", "exposure_threshold": 0.50},

    # Isomerization (Asp -> isoAsp)
    "isomerization (DG)": {"pattern": r"DG", "risk": "high", "exposure_threshold": 0.30},
    "isomerization (DS)": {"pattern": r"DS", "risk": "medium", "exposure_threshold": 0.40},
    "isomerization (DT)": {"pattern": r"DT", "risk": "medium", "exposure_threshold": 0.40},
    "isomerization (DH)": {"pattern": r"DH", "risk": "low", "exposure_threshold": 0.50},

    # Acid hydrolysis
    "acid hydrolysis (DD)": {"pattern": r"DD", "risk": "high", "exposure_threshold": 0.30},
    "acid hydrolysis (D-X)": {"pattern": r"D[AVLIP]", "risk": "medium", "exposure_threshold": 0.40},

    # Oxidation
    "oxidation (M)": {"pattern": r"M", "risk": "medium", "exposure_threshold": 0.20},
    "oxidation (W)": {"pattern": r"W", "risk": "medium", "exposure_threshold": 0.20},
    "oxidation (C)": {"pattern": r"C", "risk": "medium", "exposure_threshold": 0.20},

    # Unpaired Cys
    "unpaired Cys": {"pattern": r"C", "risk": "medium", "exposure_threshold": 0.20},

    # N-glycosylation
    "N-glycan (NxS/T)": {"pattern": r"N[^P][ST]", "risk": "high", "exposure_threshold": 0.30},

    # Met-Lys cleavage
    "met-lyscleavage (MK)": {"pattern": r"MK", "risk": "low", "exposure_threshold": 0.50},
}

# Conserved disulfide Cys positions (should not be flagged as unpaired)
CONSERVED_CYS = {
    "H": {22, 92},      # VH
    "L": {22, 23, 88},  # VL (kappa and lambda)
}

# Maximum SASA values for amino acids (Tien et al. 2013)
MAX_SASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0,
    'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
    'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
    'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
}

# 3-letter to 1-letter amino acid mapping
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


# ---------------------------------------------------------------------------
# Data Classes
# ---------------------------------------------------------------------------

@dataclass
class ResidueInfo:
    chain: str
    res_num: int
    res_name: str
    res_name_1l: str
    abs_sasa: float
    rel_sasa: float
    classification: str  # BURIED, INTERMEDIATE, EXPOSED


@dataclass
class PTMSite:
    chain: str
    res_num: int
    res_name: str
    motif: str
    motif_risk: str  # high, medium, low
    sequence_context: str
    abs_sasa: float
    rel_sasa: float
    sasa_classification: str
    exposure_threshold: float
    is_exposed: bool
    adjusted_risk: str  # Risk after exposure adjustment


# ---------------------------------------------------------------------------
# PDB Parsing
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
    sequences = {}
    for chain, residues in atoms.items():
        seq = []
        for res_num in sorted(residues.keys()):
            res_name = residues[res_num][0]['resname']
            if res_name in AA_3TO1:
                seq.append(AA_3TO1[res_name])
            else:
                seq.append('X')  # Unknown residue
        sequences[chain] = ''.join(seq)
    return sequences


# ---------------------------------------------------------------------------
# SASA Calculation
# ---------------------------------------------------------------------------

def calculate_sasa(pdb_path: str) -> Dict[str, Dict[int, float]]:
    """Calculate per-residue SASA using FreeSASA."""
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)

    sasa_by_residue = defaultdict(lambda: defaultdict(float))
    n = structure.nAtoms()

    for i in range(n):
        chain = structure.chainLabel(i)
        # Convert residue number to int (FreeSASA may return string with spaces)
        res_num = int(str(structure.residueNumber(i)).strip())
        area = result.atomArea(i)
        sasa_by_residue[chain][res_num] += area

    return dict(sasa_by_residue)


def get_residue_info(chain: str, res_num: int, res_name: str,
                     abs_sasa: float) -> ResidueInfo:
    """Get detailed residue information including exposure classification."""
    if res_name in MAX_SASA:
        max_sasa = MAX_SASA[res_name]
        rel_sasa = abs_sasa / max_sasa if max_sasa > 0 else 0.0
    else:
        rel_sasa = 0.0

    if rel_sasa < 0.20:
        classification = "BURIED"
    elif rel_sasa < 0.50:
        classification = "INTERMEDIATE"
    else:
        classification = "EXPOSED"

    res_name_1l = AA_3TO1.get(res_name, 'X')

    return ResidueInfo(
        chain=chain,
        res_num=res_num,
        res_name=res_name,
        res_name_1l=res_name_1l,
        abs_sasa=abs_sasa,
        rel_sasa=rel_sasa,
        classification=classification
    )


# ---------------------------------------------------------------------------
# PTM Detection
# ---------------------------------------------------------------------------

def detect_ptm_sites(sequence: str, chain: str,
                     sasa_data: Dict[int, float],
                     residue_names: Dict[int, str]) -> List[PTMSite]:
    """Detect PTM sites in sequence and assess exposure."""
    sites = []
    conserved_cys = CONSERVED_CYS.get(chain, set())

    for motif_name, motif_info in PTM_PATTERNS.items():
        pattern = motif_info["pattern"]
        risk = motif_info["risk"]
        threshold = motif_info["exposure_threshold"]

        for m in re.finditer(pattern, sequence):
            start = m.start()
            # Get the position of the first residue in the match
            res_num = start + 1  # 1-based indexing

            # Skip if this is a conserved Cys and motif involves C
            if "Cys" in motif_name or "oxidation (C)" in motif_name:
                if res_num in conserved_cys:
                    continue

            # Get sequence context (3 residues before and after)
            ctx_start = max(0, start - 3)
            ctx_end = min(len(sequence), m.end() + 3)
            context = sequence[ctx_start:ctx_end]

            # Get SASA data
            abs_sasa = sasa_data.get(res_num, 0.0)
            res_name = residue_names.get(res_num, 'UNK')

            if res_name in MAX_SASA:
                max_sasa = MAX_SASA[res_name]
                rel_sasa = abs_sasa / max_sasa if max_sasa > 0 else 0.0
            else:
                rel_sasa = 0.0

            # Classification
            if rel_sasa < 0.20:
                sasa_class = "BURIED"
            elif rel_sasa < 0.50:
                sasa_class = "INTERMEDIATE"
            else:
                sasa_class = "EXPOSED"

            # Exposure assessment
            is_exposed = rel_sasa >= threshold

            # Adjusted risk based on exposure
            if is_exposed:
                adjusted_risk = risk  # Keep original risk
            else:
                # Downgrade risk if buried
                risk_levels = {"high": "medium", "medium": "low", "low": "low"}
                adjusted_risk = risk_levels.get(risk, "low")

            # Create 1-letter amino acid name
            res_name_1l = AA_3TO1.get(res_name, 'X')

            sites.append(PTMSite(
                chain=chain,
                res_num=res_num,
                res_name=res_name_1l,
                motif=motif_name,
                motif_risk=risk,
                sequence_context=context,
                abs_sasa=abs_sasa,
                rel_sasa=rel_sasa,
                sasa_classification=sasa_class,
                exposure_threshold=threshold,
                is_exposed=is_exposed,
                adjusted_risk=adjusted_risk
            ))

    return sites


# ---------------------------------------------------------------------------
# Main Analysis
# ---------------------------------------------------------------------------

def analyze_structure(pdb_path: str) -> Dict:
    """Analyze a single structure for PTM exposure."""
    # Parse PDB
    atoms = parse_pdb(pdb_path)
    sequences = get_chain_sequence(atoms)

    # Calculate SASA
    sasa_data = calculate_sasa(pdb_path)

    # Analyze each chain
    results = {
        "pdb_file": pdb_path,
        "chains": {},
        "summary": {
            "total_ptm_sites": 0,
            "high_risk": 0,
            "medium_risk": 0,
            "low_risk": 0,
            "exposed_high_risk": 0,
            "buried_high_risk": 0,
        }
    }

    for chain, sequence in sequences.items():
        if chain not in sasa_data:
            continue

        # Get residue names
        residue_names = {}
        for res_num in sasa_data[chain]:
            for atom in atoms[chain].get(res_num, []):
                residue_names[res_num] = atom['resname']
                break

        # Detect PTM sites
        ptm_sites = detect_ptm_sites(sequence, chain, sasa_data[chain], residue_names)

        # Get all residue info
        residue_info = []
        for res_num in sorted(sasa_data[chain].keys()):
            res_name = residue_names.get(res_num, 'UNK')
            abs_sasa = sasa_data[chain][res_num]
            info = get_residue_info(chain, res_num, res_name, abs_sasa)
            residue_info.append(info)

        # Update summary
        for site in ptm_sites:
            results["summary"]["total_ptm_sites"] += 1
            if site.adjusted_risk == "high":
                results["summary"]["high_risk"] += 1
                if site.is_exposed:
                    results["summary"]["exposed_high_risk"] += 1
                else:
                    results["summary"]["buried_high_risk"] += 1
            elif site.adjusted_risk == "medium":
                results["summary"]["medium_risk"] += 1
            else:
                results["summary"]["low_risk"] += 1

        results["chains"][chain] = {
            "sequence": sequence,
            "ptm_sites": [asdict(site) for site in ptm_sites],
            "residue_count": len(residue_info),
        }

    return results


# ---------------------------------------------------------------------------
# Output Formatting
# ---------------------------------------------------------------------------

def format_text_report(results: Dict, title: str = "PTM Exposure Analysis") -> str:
    """Format results as text report."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"  {title}")
    lines.append("=" * 80)
    lines.append(f"\n  File: {results['pdb_file']}")
    lines.append("")

    for chain, data in results["chains"].items():
        lines.append(f"\n  Chain {chain}:")
        lines.append(f"  {'='*70}")

        if not data["ptm_sites"]:
            lines.append("    No PTM sites detected.")
            continue

        lines.append(f"  {'Pos':<8} {'AA':<5} {'Motif':<25} {'AbsSASA':<10} {'RelSASA':<10} {'Class':<15} {'Risk':<10}")
        lines.append("  " + "-" * 85)

        for site in sorted(data["ptm_sites"], key=lambda x: x["res_num"]):
            risk_marker = "!" if site["adjusted_risk"] == "high" else " " if site["adjusted_risk"] == "medium" else " "
            lines.append(
                f"  {chain}{site['res_num']:<5} {site['res_name']:<5} {site['motif']:<25} "
                f"{site['abs_sasa']:<10.1f} {site['rel_sasa']:<10.3f} {site['sasa_classification']:<15} "
                f"{site['adjusted_risk'].upper():<10}{risk_marker}"
            )
            lines.append(f"         Context: ...{site['sequence_context']}...")

    # Summary
    lines.append(f"\n{'='*80}")
    lines.append("  SUMMARY")
    lines.append(f"{'='*80}")
    summary = results["summary"]
    lines.append(f"  Total PTM sites:      {summary['total_ptm_sites']}")
    lines.append(f"  High risk (exposed):  {summary['exposed_high_risk']}")
    lines.append(f"  High risk (buried):   {summary['buried_high_risk']}")
    lines.append(f"  Medium risk:          {summary['medium_risk']}")
    lines.append(f"  Low risk:             {summary['low_risk']}")
    lines.append("")

    return "\n".join(lines)


def save_json_report(results: Dict, output_path: str):
    """Save results as JSON file."""
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)


def save_markdown_report(results: Dict, output_path: str, title: str = "PTM Exposure Analysis"):
    """Save results as Markdown file."""
    lines = []
    lines.append(f"# {title}\n")
    lines.append(f"**File:** `{results['pdb_file']}`\n")

    for chain, data in results["chains"].items():
        lines.append(f"\n## Chain {chain}\n")

        if not data["ptm_sites"]:
            lines.append("No PTM sites detected.\n")
            continue

        lines.append("| Position | AA | Motif | AbsSASA | RelSASA | Classification | Risk |")
        lines.append("|----------|-----|-------|---------|---------|----------------|------|")

        for site in sorted(data["ptm_sites"], key=lambda x: x["res_num"]):
            risk_marker = "**" if site["adjusted_risk"] == "high" else ""
            lines.append(
                f"| {chain}{site['res_num']} | {site['res_name']} | {site['motif']} | "
                f"{site['abs_sasa']:.1f} | {site['rel_sasa']:.3f} | "
                f"{site['sasa_classification']} | {risk_marker}{site['adjusted_risk'].upper()}{risk_marker} |"
            )

    # Summary
    lines.append(f"\n## Summary\n")
    summary = results["summary"]
    lines.append(f"- **Total PTM sites:** {summary['total_ptm_sites']}")
    lines.append(f"- **High risk (exposed):** {summary['exposed_high_risk']}")
    lines.append(f"- **High risk (buried):** {summary['buried_high_risk']}")
    lines.append(f"- **Medium risk:** {summary['medium_risk']}")
    lines.append(f"- **Low risk:** {summary['low_risk']}")
    lines.append("")

    with open(output_path, 'w') as f:
        f.write("\n".join(lines))


def save_word_report(results: Dict, output_path: str, title: str = "PTM Exposure Analysis"):
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

    # File info
    doc.add_paragraph(f"File: {results['pdb_file']}", style='Intense Quote')

    # Process each chain
    for chain, data in results["chains"].items():
        doc.add_heading(f'Chain {chain}', level=1)

        if not data["ptm_sites"]:
            doc.add_paragraph("No PTM sites detected.")
            continue

        # Add sequence info
        seq = data["sequence"]
        doc.add_paragraph(f"Sequence length: {len(seq)} residues")
        
        # Create table for PTM sites
        table = doc.add_table(rows=1, cols=7)
        table.style = 'Light Grid Accent 1'
        table.alignment = WD_TABLE_ALIGNMENT.CENTER

        # Header row
        headers = ['Position', 'AA', 'Motif', 'AbsSASA', 'RelSASA', 'Classification', 'Risk']
        for i, header in enumerate(headers):
            cell = table.rows[0].cells[i]
            cell.text = header
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.bold = True
                    run.font.size = Pt(9)

        # Data rows
        for site in sorted(data["ptm_sites"], key=lambda x: x["res_num"]):
            row = table.add_row()
            row.cells[0].text = f"{chain}{site['res_num']}"
            row.cells[1].text = site["res_name"]
            row.cells[2].text = site["motif"]
            row.cells[3].text = f"{site['abs_sasa']:.1f}"
            row.cells[4].text = f"{site['rel_sasa']:.3f}"
            row.cells[5].text = site["sasa_classification"]
            row.cells[6].text = site["adjusted_risk"].upper()

            # Color code risk level
            risk_cell = row.cells[6]
            for paragraph in risk_cell.paragraphs:
                for run in paragraph.runs:
                    if site["adjusted_risk"] == "high":
                        run.font.color.rgb = RGBColor(255, 0, 0)  # Red
                        run.font.bold = True
                    elif site["adjusted_risk"] == "medium":
                        run.font.color.rgb = RGBColor(255, 165, 0)  # Orange

            # Set font size for all cells
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    for run in paragraph.runs:
                        run.font.size = Pt(9)

        # Add context for high/medium risk sites
        high_med_sites = [s for s in data["ptm_sites"] 
                         if s["adjusted_risk"] in ("high", "medium")]
        if high_med_sites:
            doc.add_heading('High/Medium Risk Details', level=2)
            for site in high_med_sites:
                p = doc.add_paragraph()
                p.add_run(f"{chain}{site['res_num']} ({site['res_name']})").bold = True
                p.add_run(f" - {site['motif']}")
                p.add_run(f"\nContext: ...{site['sequence_context']}...")
                p.add_run(f"\nRelSASA: {site['rel_sasa']:.3f} ({site['sasa_classification']})")

    # Summary section
    doc.add_heading('Summary', level=1)
    summary = results["summary"]
    
    summary_table = doc.add_table(rows=6, cols=2)
    summary_table.style = 'Light Grid Accent 1'
    
    summary_data = [
        ("Total PTM sites", summary['total_ptm_sites']),
        ("High risk (exposed)", summary['exposed_high_risk']),
        ("High risk (buried)", summary['buried_high_risk']),
        ("Medium risk", summary['medium_risk']),
        ("Low risk", summary['low_risk']),
    ]
    
    for i, (label, value) in enumerate(summary_data):
        summary_table.rows[i].cells[0].text = label
        summary_table.rows[i].cells[1].text = str(value)

    doc.save(output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyze PTM exposure in antibody structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze single structure
  python3 tools/ptm_exposure/analyze_ptm.py --pdb af3/E51H2_coVL/Structure\\ Prediction\\ \\(Boltz-2\\)/rank_1.pdb

  # Analyze with output directory
  python3 tools/ptm_exposure/analyze_ptm.py --pdb structure.pdb --output outputs/ptm/

  # Analyze all structures in a directory
  python3 tools/ptm_exposure/analyze_ptm.py --pdb-dir af3/ --output outputs/ptm_batch/
        """
    )
    parser.add_argument("--pdb", help="Path to PDB/CIF structure file")
    parser.add_argument("--pdb-dir", help="Directory containing PDB files (recursive)")
    parser.add_argument("--output", "-o", default="outputs/ptm_exposure",
                        help="Output directory (default: outputs/ptm_exposure)")
    parser.add_argument("--json", action="store_true", help="Save JSON report")
    parser.add_argument("--markdown", action="store_true", help="Save Markdown report")
    parser.add_argument("--word", action="store_true", help="Save Word report")
    parser.add_argument("--all-formats", action="store_true",
                        help="Save all formats (text + JSON + Markdown + Word)")

    args = parser.parse_args()

    if not args.pdb and not args.pdb_dir:
        parser.error("Either --pdb or --pdb-dir must be specified")

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Collect PDB files
    pdb_files = []
    if args.pdb:
        pdb_files.append(Path(args.pdb))
    elif args.pdb_dir:
        pdb_dir = Path(args.pdb_dir)
        for pdb_file in pdb_dir.rglob("*.pdb"):
            pdb_files.append(pdb_file)

    if not pdb_files:
        print("No PDB files found.", file=sys.stderr)
        sys.exit(1)

    # Process each file
    all_results = []
    for pdb_path in pdb_files:
        print(f"\nAnalyzing: {pdb_path}")

        try:
            results = analyze_structure(str(pdb_path))
            all_results.append(results)

            # Print text report
            title = f"PTM Exposure Analysis: {pdb_path.stem}"
            text_report = format_text_report(results, title)
            print(text_report)

            # Save reports - use parent directory name for unique filenames
            # e.g., af3/C45H4_coVL/Structure Prediction (Boltz-2)/rank_1.pdb -> C45H4_coVL_rank_1
            parent_name = pdb_path.parent.parent.name.replace("Structure Prediction (Boltz-2)_", "")
            stem = f"{parent_name}_{pdb_path.stem}" if parent_name else pdb_path.stem
            
            if args.json or args.all_formats:
                json_path = os.path.join(args.output, f"{stem}_ptm.json")
                save_json_report(results, json_path)
                print(f"  JSON report saved: {json_path}")

            if args.markdown or args.all_formats:
                md_path = os.path.join(args.output, f"{stem}_ptm.md")
                save_markdown_report(results, md_path, title)
                print(f"  Markdown report saved: {md_path}")

            if args.word or args.all_formats:
                docx_path = os.path.join(args.output, f"{stem}_ptm.docx")
                save_word_report(results, docx_path, title)
                print(f"  Word report saved: {docx_path}")

        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            continue

    # Save combined JSON if multiple files
    if len(all_results) > 1 and (args.json or args.all_formats):
        combined_path = os.path.join(args.output, "combined_ptm_results.json")
        with open(combined_path, 'w') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False)
        print(f"\nCombined JSON report saved: {combined_path}")


if __name__ == "__main__":
    main()
