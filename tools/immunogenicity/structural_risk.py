#!/usr/bin/env python3
"""Structural risk adjustment for immunogenicity prediction.

Two-step analysis workflow:
1. Sequence-based analysis (HLAIIPred/heuristic) → initial risk scores
2. Structure-based confirmation (AF3/IgFold) → adjust risk based on buriedness

Key insight: Buried residues are less accessible to MHC-II antigen processing,
so sequence-based predictions may overestimate risk for buried positions.

Reference:
- Tien et al. (2013) consensus 3D protein structure
- FreeSASA: Beer et al. (2016)
"""
from __future__ import annotations

import csv
import json
import os
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Try to import freesasa
try:
    import freesasa
    HAS_FREESASA = True
except ImportError:
    HAS_FREESASA = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reference SASA values from Tien et al. 2013 (A^2)
REFERENCE_SASA = {
    'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0, 'CYS': 140.0,
    'GLU': 183.0, 'GLN': 179.0, 'GLY': 85.0, 'HIS': 194.0, 'ILE': 182.0,
    'LEU': 180.0, 'LYS': 211.0, 'MET': 203.0, 'PHE': 218.0, 'PRO': 143.0,
    'SER': 122.0, 'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0,
}

# Buriedness thresholds (Tien et al. 2013)
BURIED_THRESHOLD = 0.20    # relSASA < 0.20 = buried
EXPOSED_THRESHOLD = 0.25   # relSASA > 0.25 = exposed
# 0.15 <= relSASA <= 0.25 = uncertain zone (linear interpolation)

# Risk adjustment factors
RISK_ADJUSTMENT = {
    'buried': 0.1,         # 90% risk reduction
    'exposed': 1.0,        # no reduction
    'uncertain': 0.5,      # 50% reduction
    'no_structure': 0.5,   # conservative: assume 50% risk
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PositionBuriedness:
    """Buriedness data for a single residue position."""
    position: str           # e.g. "H5"
    chain: str              # "H" or "L"
    seq_index: int          # 0-based index in sequence
    resname: str            # 3-letter residue name
    abs_sasa: float         # absolute SASA (A^2)
    rel_sasa: float         # relative SASA (0-1)
    buried: Optional[bool]  # True=buried, False=exposed, None=uncertain
    is_backmutation: bool   # True if this is a backmutation position
    donor_aa: Optional[str] = None
    human_aa: Optional[str] = None
    tier: Optional[str] = None  # T1/T2/T3
    cdr_contact: Optional[bool] = None
    antigen_contact: Optional[bool] = None

    @property
    def risk_class(self) -> str:
        """Risk classification based on buriedness."""
        if self.buried is True:
            return 'buried'
        elif self.buried is False:
            return 'exposed'
        else:
            return 'uncertain'

    @property
    def risk_factor(self) -> float:
        """Risk adjustment factor (0-1)."""
        return RISK_ADJUSTMENT.get(self.risk_class, 0.5)


@dataclass
class PeptideRiskAssessment:
    """Structural risk assessment for a single peptide."""
    peptide: str
    start_index: int        # 0-based start position in sequence
    end_index: int          # 0-based end position (exclusive)
    raw_score: float        # HLAIIPred/heuristic score
    positions: List[PositionBuriedness] = field(default_factory=list)
    adjustment_factor: float = 1.0
    adjusted_score: float = 0.0
    n_buried: int = 0
    n_exposed: int = 0
    n_uncertain: int = 0
    n_backmutation: int = 0
    rationale: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.adjusted_score = self.raw_score * self.adjustment_factor


@dataclass
class VariantStructuralRisk:
    """Structural risk assessment for an entire variant."""
    variant_name: str
    chain_type: str
    sequence: str
    peptides: List[PeptideRiskAssessment] = field(default_factory=list)
    overall_adjustment: float = 1.0
    structure_source: str = ""  # e.g. "V2.pdb", "donor.pdb"


# ---------------------------------------------------------------------------
# SASA computation
# ---------------------------------------------------------------------------

def compute_rel_sasa(pdb_path: str, chain_id: str = 'H') -> Dict[str, PositionBuriedness]:
    """Compute relative SASA for all residues in a PDB file.

    Args:
        pdb_path: Path to PDB structure file
        chain_id: Chain identifier to analyze

    Returns:
        Dict mapping position string (e.g. "H5") to PositionBuriedness
    """
    if not HAS_FREESASA:
        raise ImportError("FreeSASA is required. Install with: pip install freesasa")

    import freesasa as _freesasa
    structure = _freesasa.Structure(pdb_path)
    result = _freesasa.calc(structure)

    residue_areas = result.residueAreas()
    chain_areas = residue_areas.get(chain_id, {})

    rel_sasa = {}
    for resseq_str, area_obj in chain_areas.items():
        resseq = int(resseq_str)
        total_area = area_obj.total
        rel_total = area_obj.relativeTotal
        resname = area_obj.residueType

        # Determine buriedness
        if rel_total < BURIED_THRESHOLD:
            buried = True
        elif rel_total > EXPOSED_THRESHOLD:
            buried = False
        else:
            buried = None  # uncertain zone

        position = f'{chain_id}{resseq}'
        rel_sasa[position] = PositionBuriedness(
            position=position,
            chain=chain_id,
            seq_index=resseq - 1,  # convert to 0-based
            resname=resname,
            abs_sasa=total_area,
            rel_sasa=rel_total,
            buried=buried,
            is_backmutation=False,  # will be updated later
        )

    return rel_sasa


def compute_all_chains_sasa(pdb_path: str) -> Dict[str, Dict[str, PositionBuriedness]]:
    """Compute relSASA for all chains in a PDB file.

    Returns:
        Dict mapping chain_id -> {position -> PositionBuriedness}
    """
    if not HAS_FREESASA:
        raise ImportError("FreeSASA is required")

    import freesasa as _freesasa
    structure = _freesasa.Structure(pdb_path)
    result = _freesasa.calc(structure)

    residue_areas = result.residueAreas()
    all_chains = {}

    for chain_id, chain_areas in residue_areas.items():
        chain_data = {}
        for resseq_str, area_obj in chain_areas.items():
            resseq = int(resseq_str)
            total_area = area_obj.total
            rel_total = area_obj.relativeTotal
            resname = area_obj.residueType

            if rel_total < BURIED_THRESHOLD:
                buried = True
            elif rel_total > EXPOSED_THRESHOLD:
                buried = False
            else:
                buried = None

            position = f'{chain_id}{resseq}'
            chain_data[position] = PositionBuriedness(
                position=position,
                chain=chain_id,
                seq_index=resseq - 1,
                resname=resname,
                abs_sasa=total_area,
                rel_sasa=rel_total,
                buried=buried,
                is_backmutation=False,
            )
        all_chains[chain_id] = chain_data

    return all_chains


# ---------------------------------------------------------------------------
# Position mapping
# ---------------------------------------------------------------------------

def map_peptide_positions(
    peptide_start: int,
    peptide_end: int,
    sequence: str,
    chain_id: str,
    structure_data: Dict[str, PositionBuriedness],
    backmutations: Optional[List[Dict]] = None,
) -> List[PositionBuriedness]:
    """Map a peptide range to its constituent positions with buriedness data.

    Args:
        peptide_start: 0-based start index of peptide in sequence
        peptide_end: 0-based end index (exclusive)
        sequence: full chain sequence
        chain_id: chain identifier ("H" or "L")
        structure_data: dict from compute_rel_sasa
        backmutations: list of backmutation dicts from pipeline output

    Returns:
        List of PositionBuriedness for positions within the peptide
    """
    # Build backmutation lookup
    backmut_set = {}
    if backmutations:
        for bm in backmutations:
            pos_str = bm.get('position', '')
            if pos_str.startswith(chain_id):
                pos_num = int(pos_str[len(chain_id):])
                backmut_set[pos_num] = bm

    positions = []
    for seq_idx in range(peptide_start, min(peptide_end, len(sequence))):
        pos_num = seq_idx + 1  # 1-based
        position = f'{chain_id}{pos_num}'

        # Get structure data
        struct_pos = structure_data.get(position)

        # Check if this is a backmutation
        bm_info = backmut_set.get(pos_num)

        if struct_pos:
            # Update with backmutation info
            struct_pos.is_backmutation = bm_info is not None
            if bm_info:
                struct_pos.donor_aa = bm_info.get('donor_aa')
                struct_pos.human_aa = bm_info.get('human_aa')
                struct_pos.tier = bm_info.get('tier')
                struct_pos.cdr_contact = bm_info.get('cdr_contact')
                struct_pos.antigen_contact = bm_info.get('antigen_contact')
            positions.append(struct_pos)
        else:
            # No structure data for this position
            positions.append(PositionBuriedness(
                position=position,
                chain=chain_id,
                seq_index=seq_idx,
                resname='',
                abs_sasa=0.0,
                rel_sasa=0.0,
                buried=None,
                is_backmutation=bm_info is not None,
                donor_aa=bm_info.get('donor_aa') if bm_info else None,
                human_aa=bm_info.get('human_aa') if bm_info else None,
                tier=bm_info.get('tier') if bm_info else None,
                cdr_contact=bm_info.get('cdr_contact') if bm_info else None,
                antigen_contact=bm_info.get('antigen_contact') if bm_info else None,
            ))

    return positions


# ---------------------------------------------------------------------------
# Risk adjustment calculation
# ---------------------------------------------------------------------------

def calculate_adjustment_factor(
    positions: List[PositionBuriedness],
    method: str = "relSASA",
) -> Tuple[float, List[str]]:
    """Calculate structural risk adjustment factor for a peptide.

    Args:
        positions: list of PositionBuriedness for positions in peptide
        method: adjustment method ("relSASA", "binary", "tier")

    Returns:
        (adjustment_factor, rationale_list)
    """
    backmut_positions = [p for p in positions if p.is_backmutation]
    rationale = []

    if not backmut_positions:
        rationale.append("No backmutation positions in peptide: false positive")
        return 0.0, rationale

    n_total = len(backmut_positions)
    n_buried = sum(1 for p in backmut_positions if p.buried is True)
    n_exposed = sum(1 for p in backmut_positions if p.buried is False)
    n_uncertain = sum(1 for p in backmut_positions if p.buried is None)

    if method == "binary":
        if n_buried == n_total:
            rationale.append(f"All {n_total} backmutation positions buried (relSASA < {BURIED_THRESHOLD})")
            return RISK_ADJUSTMENT['buried'], rationale
        elif n_exposed > 0:
            factor = n_exposed / n_total
            rationale.append(f"{n_exposed}/{n_total} backmutation positions exposed (relSASA > {EXPOSED_THRESHOLD})")
            return factor, rationale
        else:
            rationale.append(f"All {n_total} backmutation positions uncertain")
            return RISK_ADJUSTMENT['uncertain'], rationale

    elif method == "relSASA":
        # Continuous weighting based on relSASA values
        weighted_sum = 0.0
        weight_sum = 0.0

        for pos in backmut_positions:
            # Base weight: importance of this backmutation position
            tier_weights = {
                'T1': 1.0,
                'T2': 0.8,
                'T3': 0.5,
                'T_FR4': 0.9,
                'KEEP_DONOR': 1.0,
            }
            base_weight = tier_weights.get(pos.tier or '', 0.6)

            # Structural weight based on relSASA
            if pos.rel_sasa > 0:
                # Use continuous relSASA weighting
                if pos.rel_sasa < 0.10:
                    struct_weight = 0.05
                elif pos.rel_sasa < 0.15:
                    struct_weight = 0.05 + 0.25 * (pos.rel_sasa - 0.10) / 0.05
                elif pos.rel_sasa < 0.25:
                    struct_weight = 0.30 + 0.40 * (pos.rel_sasa - 0.15) / 0.10
                elif pos.rel_sasa < 0.40:
                    struct_weight = 0.70 + 0.20 * (pos.rel_sasa - 0.25) / 0.15
                else:
                    struct_weight = 1.0
            elif pos.buried is not None:
                struct_weight = RISK_ADJUSTMENT['buried'] if pos.buried else RISK_ADJUSTMENT['exposed']
            else:
                struct_weight = RISK_ADJUSTMENT['no_structure']

            # CDR contact adjustment
            if pos.cdr_contact:
                struct_weight = max(struct_weight, 0.6)
                rationale.append(f"{pos.position}: CDR contact (structurally relevant)")

            weighted_sum += base_weight * struct_weight
            weight_sum += base_weight

            rationale.append(
                f"{pos.position}: relSASA={pos.rel_sasa:.3f}, "
                f"weight={struct_weight:.3f}"
            )

        factor = weighted_sum / weight_sum if weight_sum > 0 else 0.5
        return factor, rationale

    else:
        raise ValueError(f"Unknown method: {method}")


# ---------------------------------------------------------------------------
# Main structural risk assessment
# ---------------------------------------------------------------------------

def assess_variant_structural_risk(
    variant_name: str,
    chain_type: str,
    sequence: str,
    per_peptide_hits: List[Dict],
    structure_path: str,
    backmutations: Optional[List[Dict]] = None,
    method: str = "relSASA",
) -> VariantStructuralRisk:
    """Assess structural risk for a variant using its actual structure.

    This is the two-step approach:
    1. Use sequence-based scores (from HLAIIPred/heuristic)
    2. Adjust based on structural buriedness from actual variant structure

    Args:
        variant_name: variant identifier (e.g. "H_V2")
        chain_type: "H" or "L"
        sequence: full chain sequence
        per_peptide_hits: list of peptide dicts from HLAIIPred
        structure_path: path to variant PDB file
        backmutations: list of backmutation dicts from pipeline
        method: adjustment method

    Returns:
        VariantStructuralRisk with adjusted scores
    """
    # Compute relSASA from variant structure
    structure_data = compute_rel_sasa(structure_path, chain_type)

    # Assess each peptide
    peptides = []
    for hit in per_peptide_hits:
        peptide = hit.get('peptide', '')
        score = hit.get('score', 0.0)
        start_idx = hit.get('core_pos', 0) if 'core_pos' in hit else 0

        # Find peptide position in sequence
        # Try to find exact match
        seq_idx = sequence.find(peptide)
        if seq_idx >= 0:
            start_idx = seq_idx
        else:
            # Use provided position or approximate
            start_idx = hit.get('start_index', start_idx)

        end_idx = start_idx + len(peptide)

        # Map peptide to positions
        positions = map_peptide_positions(
            start_idx, end_idx, sequence, chain_type,
            structure_data, backmutations
        )

        # Calculate adjustment factor
        factor, rationale = calculate_adjustment_factor(positions, method)

        # Count position categories
        backmut_positions = [p for p in positions if p.is_backmutation]
        n_buried = sum(1 for p in backmut_positions if p.buried is True)
        n_exposed = sum(1 for p in backmut_positions if p.buried is False)
        n_uncertain = sum(1 for p in backmut_positions if p.buried is None)

        adjusted_score = score * factor
        print(f"DEBUG肽: {peptide}: score={score:.3f}, factor={factor:.3f}, adjusted={adjusted_score:.3f}")

        peptides.append(PeptideRiskAssessment(
            peptide=peptide,
            start_index=start_idx,
            end_index=end_idx,
            raw_score=score,
            positions=positions,
            adjustment_factor=factor,
            adjusted_score=adjusted_score,
            n_buried=n_buried,
            n_exposed=n_exposed,
            n_uncertain=n_uncertain,
            n_backmutation=len(backmut_positions),
            rationale=rationale,
        ))

    # Calculate overall adjustment (weighted by raw score)
    total_raw = sum(p.raw_score for p in peptides)
    total_adjusted = sum(p.adjusted_score for p in peptides)
    overall_adjustment = total_adjusted / total_raw if total_raw > 0 else 1.0

    print(f"DEBUG: total_raw={total_raw:.3f}, total_adjusted={total_adjusted:.3f}, overall={overall_adjustment:.3f}")

    return VariantStructuralRisk(
        variant_name=variant_name,
        chain_type=chain_type,
        sequence=sequence,
        peptides=peptides,
        overall_adjustment=overall_adjustment,
        structure_source=os.path.basename(structure_path),
    )


def load_backmutations_from_csv(csv_path: str) -> List[Dict]:
    """Load backmutation data from pipeline CSV output."""
    backmutations = []
    if not os.path.isfile(csv_path):
        return backmutations

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            backmutations.append({
                'position': row.get('position', ''),
                'donor_aa': row.get('donor_aa', ''),
                'human_aa': row.get('human_aa', ''),
                'tier': row.get('tier', ''),
                'buried': row.get('buried', '') == 'True',
                'cdr_contact': row.get('cdr_contact', '') == 'True',
                'antigen_contact': row.get('antigen_contact', '') == 'True',
                'composite': float(row.get('composite', 0)),
            })
    return backmutations


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_structural_risk_report(
    assessments: List[VariantStructuralRisk],
    output_path: str,
    title: str = "Structural Risk Assessment Report",
) -> None:
    """Generate detailed Markdown report of structural risk assessment."""
    lines = [
        f"# {title}",
        "",
        "## Overview",
        "",
        "This report applies structural risk adjustment to immunogenicity predictions.",
        "",
        "**Two-step workflow:**",
        "1. **Sequence-based analysis**: HLAIIPred/heuristic predicts MHC-II binding from sequence",
        "2. **Structure-based confirmation**: FreeSASA computes relSASA from actual variant structure",
        "",
        "**Key principle**: Buried residues (relSASA < 0.20) are less accessible to MHC-II",
        "antigen processing, so sequence-based predictions may overestimate risk.",
        "",
        "---",
        "",
    ]

    for assessment in assessments:
        lines.append(f"## {assessment.variant_name} ({assessment.chain_type} chain)")
        lines.append("")
        lines.append(f"- **Structure**: {assessment.structure_source}")
        lines.append(f"- **Overall adjustment**: {assessment.overall_adjustment:.3f}")
        lines.append("")

        # Peptide summary table
        lines.append("### Peptide Risk Summary")
        lines.append("")
        lines.append("| Peptide | Raw Score | Adj Factor | Adj Score | Buried | Exposed | Risk Level |")
        lines.append("|---------|-----------|------------|-----------|--------|---------|------------|")

        for p in sorted(assessment.peptides, key=lambda x: x.adjusted_score, reverse=True):
            if p.raw_score < 0.3:  # skip low-confidence peptides
                continue
            risk_level = "HIGH" if p.adjusted_score > 0.5 else ("MEDIUM" if p.adjusted_score > 0.2 else "LOW")
            lines.append(
                f"| `{p.peptide[:20]}...` | {p.raw_score:.3f} | "
                f"{p.adjustment_factor:.3f} | {p.adjusted_score:.3f} | "
                f"{p.n_buried} | {p.n_exposed} | {risk_level} |"
            )
        lines.append("")

        # Detailed position analysis for high-risk peptides
        high_risk = [p for p in assessment.peptides if p.adjusted_score > 0.3]
        if high_risk:
            lines.append("### Detailed Position Analysis")
            lines.append("")
            for p in high_risk[:5]:  # top 5
                lines.append(f"#### Peptide: `{p.peptide}`")
                lines.append(f"- Raw score: {p.raw_score:.3f}")
                lines.append(f"- Adjusted score: {p.adjusted_score:.3f}")
                lines.append("")
                lines.append("| Position | relSASA | Buried | CDR Contact | Tier | Risk Factor |")
                lines.append("|----------|---------|--------|-------------|------|-------------|")
                for pos in p.positions:
                    if pos.is_backmutation:
                        buried_str = "Yes" if pos.buried else ("No" if pos.buried is False else "Uncertain")
                        cdr_str = "Yes" if pos.cdr_contact else "No"
                        lines.append(
                            f"| {pos.position} | {pos.rel_sasa:.3f} | {buried_str} | "
                            f"{cdr_str} | {pos.tier} | {pos.risk_factor:.3f} |"
                        )
                lines.append("")

        # Rationale
        all_rationale = []
        for p in assessment.peptides:
            all_rationale.extend(p.rationale)
        if all_rationale:
            lines.append("### Rationale")
            lines.append("")
            for r in set(all_rationale):  # unique rationale
                lines.append(f"- {r}")
            lines.append("")

        lines.append("---")
        lines.append("")

    # Summary comparison
    if len(assessments) > 1:
        lines.append("## Cross-Variant Comparison")
        lines.append("")
        lines.append("| Variant | Structure | Overall Adj | High-Risk Peptides |")
        lines.append("|---------|-----------|-------------|-------------------|")
        for a in assessments:
            n_high = sum(1 for p in a.peptides if p.adjusted_score > 0.5)
            lines.append(
                f"| {a.variant_name} | {a.structure_source} | "
                f"{a.overall_adjustment:.3f} | {n_high} |"
            )
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def generate_structural_risk_word(
    assessments: List[VariantStructuralRisk],
    output_path: str,
    title: str = "Structural Risk Assessment Report",
) -> None:
    """Generate Word document report of structural risk assessment."""
    try:
        from docx import Document
        from docx.shared import Inches, Pt, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        print("  WARNING: python-docx not installed. Skipping Word report.",
              file=sys.stderr)
        return

    doc = Document()
    heading = doc.add_heading(title, 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # Overview
    doc.add_heading("Overview", level=1)
    doc.add_paragraph(
        "This report applies structural risk adjustment to immunogenicity predictions."
    )
    doc.add_paragraph(
        "Two-step workflow: (1) Sequence-based analysis (HLAIIPred/heuristic) "
        "predicts MHC-II binding; (2) Structure-based confirmation (FreeSASA) "
        "adjusts risk based on actual buriedness from variant structures."
    )

    for assessment in assessments:
        doc.add_heading(f"{assessment.variant_name} ({assessment.chain_type} chain)", level=1)
        doc.add_paragraph(f"Structure: {assessment.structure_source}")
        doc.add_paragraph(f"Overall adjustment: {assessment.overall_adjustment:.3f}")

        # Peptide summary table
        doc.add_heading("Peptide Risk Summary", level=2)
        high_risk_peptides = [p for p in assessment.peptides if p.raw_score >= 0.3]
        if high_risk_peptides:
            table = doc.add_table(rows=1 + len(high_risk_peptides), cols=7)
            table.style = "Table Grid"
            headers = ["Peptide", "Raw Score", "Adj Factor", "Adj Score",
                       "Buried", "Exposed", "Risk Level"]
            for i, h in enumerate(headers):
                cell = table.rows[0].cells[i]
                cell.text = h
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.bold = True

            for ri, pep in enumerate(sorted(high_risk_peptides,
                                            key=lambda x: x.adjusted_score,
                                            reverse=True), 1):
                row = table.rows[ri]
                risk_level = "HIGH" if pep.adjusted_score > 0.5 else (
                    "MEDIUM" if pep.adjusted_score > 0.2 else "LOW")
                vals = [
                    pep.peptide[:25] + ("..." if len(pep.peptide) > 25 else ""),
                    f"{pep.raw_score:.3f}",
                    f"{pep.adjustment_factor:.3f}",
                    f"{pep.adjusted_score:.3f}",
                    str(pep.n_buried),
                    str(pep.n_exposed),
                    risk_level,
                ]
                for ci, val in enumerate(vals):
                    row.cells[ci].text = val

        # Detailed position analysis
        detailed_peptides = [p for p in assessment.peptides if p.adjusted_score > 0.3]
        if detailed_peptides:
            doc.add_heading("Detailed Position Analysis", level=2)
            for pep in detailed_peptides[:3]:
                doc.add_heading(f"Peptide: {pep.peptide}", level=3)
                doc.add_paragraph(f"Raw score: {pep.raw_score:.3f}")
                doc.add_paragraph(f"Adjusted score: {pep.adjusted_score:.3f}")

                backmut_positions = [p for p in pep.positions if p.is_backmutation]
                if backmut_positions:
                    table2 = doc.add_table(rows=1 + len(backmut_positions), cols=5)
                    table2.style = "Table Grid"
                    for i, h in enumerate(["Position", "relSASA", "Buried",
                                           "CDR Contact", "Tier"]):
                        cell = table2.rows[0].cells[i]
                        cell.text = h
                        for p in cell.paragraphs:
                            for run in p.runs:
                                run.bold = True

                    for ri, pos in enumerate(backmut_positions, 1):
                        row = table2.rows[ri]
                        buried_str = "Yes" if pos.buried else (
                            "No" if pos.buried is False else "Uncertain")
                        vals = [
                            pos.position,
                            f"{pos.rel_sasa:.3f}",
                            buried_str,
                            "Yes" if pos.cdr_contact else "No",
                            pos.tier or "",
                        ]
                        for ci, val in enumerate(vals):
                            row.cells[ci].text = val

                        # Color-code based on buriedness
                        if pos.buried is True:
                            row.cells[1].paragraphs[0].runs[0].font.color.rgb = (
                                RGBColor(0, 128, 0))  # green
                        elif pos.buried is False:
                            row.cells[1].paragraphs[0].runs[0].font.color.rgb = (
                                RGBColor(255, 0, 0))  # red

    # Cross-variant comparison
    if len(assessments) > 1:
        doc.add_heading("Cross-Variant Comparison", level=1)
        table3 = doc.add_table(rows=1 + len(assessments), cols=4)
        table3.style = "Table Grid"
        for i, h in enumerate(["Variant", "Structure", "Overall Adj",
                               "High-Risk Peptides"]):
            cell = table3.rows[0].cells[i]
            cell.text = h
            for p in cell.paragraphs:
                for run in p.runs:
                    run.bold = True

        for ri, a in enumerate(assessments, 1):
            row = table3.rows[ri]
            n_high = sum(1 for p in a.peptides if p.adjusted_score > 0.5)
            vals = [
                a.variant_name,
                a.structure_source,
                f"{a.overall_adjustment:.3f}",
                str(n_high),
            ]
            for ci, val in enumerate(vals):
                row.cells[ci].text = val

    doc.save(output_path)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def run_structural_risk_assessment(
    immunogenicity_results: Dict,
    structure_dir: str,
    backmutation_csv_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    method: str = "relSASA",
    verbose: bool = False,
) -> Dict:
    """Run structural risk assessment on immunogenicity results.

    This is the main entry point for the two-step workflow.

    Args:
        immunogenicity_results: output from immunogenicity_analyzer
        structure_dir: directory containing variant PDB files
        backmutation_csv_dir: directory containing backmutation CSVs
        output_dir: directory for output files
        method: adjustment method
        verbose: print progress

    Returns:
        Updated results with structural risk adjustments
    """
    if output_dir is None:
        output_dir = os.path.dirname(immunogenicity_results.get('input_fasta', '.'))

    os.makedirs(output_dir, exist_ok=True)

    assessments = []

    for chain_type, cdata in immunogenicity_results.get('chains', {}).items():
        variants = cdata.get('variants', {})

        # Load backmutations if available
        backmutations = []
        if backmutation_csv_dir:
            # Try various naming patterns for backmutation CSV files
            possible_csv_names = [
                f'backmutations_{chain_type}.csv',
                f'backmutations_*.{chain_type}.csv',
            ]
            for pattern in possible_csv_names:
                if '*' in pattern:
                    # Use glob to find matching files
                    import glob
                    matches = glob.glob(os.path.join(backmutation_csv_dir, pattern))
                    if matches:
                        csv_path = matches[0]  # Use the first match
                        if os.path.isfile(csv_path):
                            backmutations = load_backmutations_from_csv(csv_path)
                            if verbose:
                                print(f"  Loaded {len(backmutations)} backmutations from {csv_path}")
                            break
                else:
                    csv_path = os.path.join(backmutation_csv_dir, pattern)
                    if os.path.isfile(csv_path):
                        backmutations = load_backmutations_from_csv(csv_path)
                        if verbose:
                            print(f"  Loaded {len(backmutations)} backmutations from {csv_path}")
                        break

        print(f"DEBUG: chain={chain_type}, backmutations={len(backmutations)}")

        for vname, vdata in variants.items():
            # Find structure file
            # Try various naming patterns
            structure_path = None

            # Extract variant name without chain prefix (e.g., "H_V0" -> "V0")
            variant_part = vname.split('_', 1)[1] if '_' in vname else vname

            possible_names = [
                f"igfold_amg110_{variant_part}.pdb",
                f"igfold_{variant_part}.pdb",
                f"igfold_amg110_{variant_part.lower()}.pdb",
                f"igfold_{variant_part.lower()}.pdb",
                f"amg110_{variant_part}.pdb",
                f"amg110_{variant_part.lower()}.pdb",
                f"{variant_part}.pdb",
                f"{variant_part.lower()}.pdb",
                f"igfold_{vname.replace('_', '')}.pdb",
                f"igfold_{vname.lower().replace('_', '')}.pdb",
                f"amg110_{vname.replace('_', '')}.pdb",
                f"amg110_{vname.lower().replace('_', '')}.pdb",
                f"{vname}.pdb",
                f"{vname.lower()}.pdb",
            ]
            for pname in possible_names:
                candidate = os.path.join(structure_dir, pname)
                if os.path.isfile(candidate):
                    structure_path = candidate
                    break

            # Also check subdirectories
            if structure_path is None:
                for root, dirs, files in os.walk(structure_dir):
                    for pname in possible_names:
                        candidate = os.path.join(root, pname)
                        if os.path.isfile(candidate):
                            structure_path = candidate
                            break
                    if structure_path:
                        break

            if structure_path is None:
                if verbose:
                    print(f"  WARNING: No structure found for {vname}")
                continue

            if verbose:
                print(f"  Assessing {vname} with {os.path.basename(structure_path)}")

            # Run assessment
            assessment = assess_variant_structural_risk(
                variant_name=vname,
                chain_type=chain_type,
                sequence=vdata.get('sequence', ''),
                per_peptide_hits=vdata.get('per_peptide_hits', []),
                structure_path=structure_path,
                backmutations=backmutations,
                method=method,
            )
            assessments.append(assessment)

            if verbose:
                print(f"    Assessment result: adjustment={assessment.overall_adjustment:.3f}, "
                      f"peptides={len(assessment.peptides)}")

            # Update results with adjusted scores
            vdata['structural_adjustment'] = assessment.overall_adjustment
            vdata['adjusted_peptides'] = [
                {
                    'peptide': p.peptide,
                    'raw_score': p.raw_score,
                    'adjusted_score': p.adjusted_score,
                    'adjustment_factor': p.adjustment_factor,
                    'n_buried': p.n_buried,
                    'n_exposed': p.n_exposed,
                }
                for p in assessment.peptides
                if p.raw_score >= 0.3
            ]
            if verbose:
                print(f"    Updated {vname}: adjustment={assessment.overall_adjustment:.3f}, "
                      f"peptides={len(vdata['adjusted_peptides'])}")

    # Generate reports
    if assessments:
        # Markdown report
        md_path = os.path.join(output_dir, 'structural_risk_assessment.md')
        generate_structural_risk_report(assessments, md_path)
        if verbose:
            print(f"  Generated Markdown report: {md_path}")

        # Word report
        docx_path = os.path.join(output_dir, 'structural_risk_assessment.docx')
        generate_structural_risk_word(assessments, docx_path)
        if verbose:
            print(f"  Generated Word report: {docx_path}")

    return immunogenicity_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Structural risk adjustment for immunogenicity prediction"
    )
    parser.add_argument("--immunogenicity-json", required=True,
                        help="Path to immunogenicity results JSON")
    parser.add_argument("--structure-dir", required=True,
                        help="Directory containing variant PDB files")
    parser.add_argument("--backmutation-dir",
                        help="Directory containing backmutation CSVs")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory for reports")
    parser.add_argument("--method", default="relSASA",
                        choices=["relSASA", "binary"],
                        help="Adjustment method")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    # Load immunogenicity results
    with open(args.immunogenicity_json) as f:
        results = json.load(f)

    # Run assessment
    updated = run_structural_risk_assessment(
        results,
        args.structure_dir,
        args.backmutation_dir,
        args.output_dir,
        method=args.method,
        verbose=args.verbose,
    )

    # Save updated results
    output_json = os.path.join(args.output_dir, 'immunogenicity_with_structural_risk.json')
    with open(output_json, 'w') as f:
        json.dump(updated, f, indent=2)
    print(f"Saved updated results to {output_json}")
