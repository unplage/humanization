# PTM Exposure Analysis Tool

Analyze high-risk PTM (Post-Translational Modification) sites in antibody structures based on AF3/PDB files. Calculates actual solvent exposure (SASA) to assess true risk levels, going beyond sequence-only predictions.

## Features

- **PTM Detection**: Identifies deamidation, isomerization, oxidation, glycosylation, and other PTM risks
- **Structural Exposure**: Uses FreeSASA to calculate actual solvent accessibility
- **Risk Adjustment**: Downgrades risk for buried (non-exposed) PTM sites
- **Multiple Output Formats**: Text, JSON, Markdown, and Word reports
- **Batch Processing**: Analyze single files or entire directories

## Requirements

```bash
pip install freesasa
pip install python-docx  # Optional, for Word reports
```

## Usage

### Single Structure Analysis

```bash
# Basic analysis (text output)
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb af3/E51H2_coVL/Structure\ Prediction\ \(Boltz-2\)/rank_1.pdb

# With JSON output
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb structure.pdb \
    --json

# With Markdown output
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb structure.pdb \
    --markdown

# With Word output
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb structure.pdb \
    --word

# All formats (text + JSON + Markdown + Word)
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb structure.pdb \
    --all-formats \
    --output outputs/ptm/
```

### Batch Processing

```bash
# Analyze all PDB files in a directory
python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb-dir af3/ \
    --output outputs/ptm_batch/ \
    --all-formats
```

## PTM Patterns Detected

| Risk Type | Pattern | Risk Level | Exposure Threshold |
|-----------|---------|------------|-------------------|
| Deamidation | NG | High | relSASA ≥ 0.30 |
| Deamidation | NS | Medium | relSASA ≥ 0.40 |
| Deamidation | NH | Medium | relSASA ≥ 0.40 |
| Deamidation | ND | Low | relSASA ≥ 0.50 |
| Isomerization | DG | High | relSASA ≥ 0.30 |
| Isomerization | DS | Medium | relSASA ≥ 0.40 |
| Isomerization | DT | Medium | relSASA ≥ 0.40 |
| Isomerization | DH | Low | relSASA ≥ 0.50 |
| Acid Hydrolysis | DD | High | relSASA ≥ 0.30 |
| Acid Hydrolysis | D-X | Medium | relSASA ≥ 0.40 |
| Oxidation | M | Medium | relSASA ≥ 0.20 |
| Oxidation | W | Medium | relSASA ≥ 0.20 |
| Oxidation | C | Medium | relSASA ≥ 0.20 |
| Unpaired Cys | C | Medium | relSASA ≥ 0.20 |
| N-glycosylation | N[^P][ST] | High | relSASA ≥ 0.30 |
| Met-Lys Cleavage | MK | Low | relSASA ≥ 0.50 |

## Exposure Classification

- **BURIED**: relSASA < 0.20 (low risk)
- **INTERMEDIATE**: 0.20 ≤ relSASA < 0.50 (medium risk)
- **EXPOSED**: relSASA ≥ 0.50 (high risk)

## Risk Adjustment

PTM sites are classified as "adjusted risk" based on both sequence pattern and structural exposure:

- **Exposed** (relSASA ≥ threshold): Keep original risk level
- **Buried** (relSASA < threshold): Downgrade risk level (high→medium, medium→low)

## Output Example

```
================================================================================
  PTM Exposure Analysis: E51H2-coVL
================================================================================

  Chain A (VH):
  ========================================================================

  Pos     AA    Motif                     AbsSASA    RelSASA    Class           Risk
  ---------------------------------------------------------------------------------------
  H26     D     isomerization (DG)        45.2       0.234      INTERMEDIATE    HIGH!
         Context: ...SDGTVF...
  H27     G     isomerization (DG)        12.1       0.116      BURIED          MEDIUM
         Context: ...DGTVFP...
  H33     W     oxidation (W)             89.3       0.313      INTERMEDIATE    MEDIUM
         Context: ...NYWLGW...
  H50     M     oxidation (M)             23.4       0.104      BURIED          LOW

================================================================================
  SUMMARY
================================================================================
  Total PTM sites:      4
  High risk (exposed):  1
  High risk (buried):   0
  Medium risk:          2
  Low risk:             1
```

## Integration with Pipeline

This tool is independent of the main humanization pipeline. It can be used:

1. **Before humanization**: Assess PTM risks in donor sequences
2. **After humanization**: Validate that optimized variants have acceptable PTM profiles
3. **For external structures**: Analyze any antibody PDB/CIF file

## Notes

- Conserved disulfide Cys (VH 22/92, VL 22/23/88) are excluded from "unpaired Cys" and "oxidation (C)" risks
- SASA is calculated using the FreeSASA algorithm (Lee & Richards)
- Maximum SASA values are from Tien et al. 2013
