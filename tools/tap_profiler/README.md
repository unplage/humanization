# TAP (Therapeutic Antibody Profiler) Analysis Tool

Calculate 5 developability metrics based on TAP guidelines to assess antibody developability risk.

## Overview

TAP (Therapeutic Antibody Profiler) compares antibody properties against clinical-stage therapeutics (CSTs) to identify potential developability issues. This tool implements the 5 key TAP metrics:

1. **Total CDR Length (Ltot)** - Sum of all 6 CDR lengths (IMGT numbering via ANARCI)
2. **Patches of Surface Hydrophobicity (PSH)** - Hydrophobic patches on CDR vicinity
3. **Patches of Positive Charge (PPC)** - Positive charge patches
4. **Patches of Negative Charge (PNC)** - Negative charge patches
5. **Structural Fv Charge Symmetry Parameter (SFvCSP)** - Charge balance VH vs VL

## Features

- **ANARCI IMGT numbering**: Uses HMM-based alignment for accurate CDR boundary detection
- **MSMS support**: Uses MSMS for molecular surface calculation (preferred)
- **FreeSASA fallback**: Falls back to FreeSASA if MSMS is not available
- **Sequence-only mode**: Analysis when structure is unavailable (Ltot only)
- **Multiple output formats**: Text, JSON, Markdown, and Word reports
- **Batch processing**: Analyze single files or entire directories
- **TAP2 thresholds**: Updated thresholds from 2025 TAP2 guidelines

## Requirements

```bash
# Required for IMGT numbering
pip install ANARCI

# For molecular surface calculation (choose one):
conda install bioconda::msms  # MSMS (preferred, Linux/macOS only)
pip install freesasa  # FreeSASA (fallback, cross-platform)

# Optional
pip install python-docx  # For Word reports
```

## Usage

### Check Available Tools

```bash
# Check which tools are available
python3 tools/tap_profiler/tap_analyzer.py --check
```

### Structure-based Analysis (Recommended)

```bash
# Basic analysis (text output)
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb af3/C45H4_coVL/Structure\ Prediction\ \(Boltz-2\)/rank_1.pdb

# With JSON output
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --json

# With Markdown output
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --markdown

# With Word output
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --word

# All formats (text + JSON + Markdown + Word)
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --all-formats \
    --output outputs/tap/
```

### Sequence-only Analysis

```bash
# Uses ANARCI for IMGT-based CDR length calculation
python3 tools/tap_profiler/tap_analyzer.py \
    --vh <VH_SEQUENCE> \
    --vl <VL_SEQUENCE> \
    --name "my_antibody"
```

### Batch Processing

```bash
# Analyze all PDB files in a directory
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb-dir af3/ \
    --output outputs/tap_batch/ \
    --all-formats
```

## TAP Thresholds (TAP2 2025)

| Metric | Green (Acceptable) | Amber (Caution) | Red (Danger) |
|--------|-------------------|-----------------|--------------|
| Ltot | 43-54 | 37-42 or 55-65 | <37 or >65 |
| PSH | 111-168 | 96-111 or 168-212 | <96 or >212 |
| PPC | 0-1.33 | 1.34-4.20 | >4.20 |
| PNC | 0-1.98 | 1.99-4.43 | >4.43 |
| SFvCSP | >-6.00 | -30.60 to -6.00 | <-30.60 |

## Important: PSH/PPC/PNC Calculation Differences

**Our implementation uses patch counting, NOT surface area calculation:**

- **Official TAP**: Uses MSMS to calculate molecular surface area (Å²)
- **Our tool**: Counts contiguous patches of 3+ residues with same property

This means:
- PSH/PPC/PNC values are **NOT directly comparable** to official TAP values
- Our values represent **patch counts** (small integers), not surface area
- For accurate PSH/PPC/PNC, use the official TAP webserver or FlashTAP

**Ltot (CDR length) is accurate** and uses the same IMGT numbering as official TAP.

## Risk Assessment

- **LOW RISK**: All metrics in green range
- **LOW-MEDIUM RISK**: One amber flag
- **MEDIUM RISK**: Two or more amber flags
- **MEDIUM-HIGH RISK**: One red flag
- **HIGH RISK**: Two or more red flags

## Amino Acid Classifications

### Hydrophobic (for PSH)
A, V, L, I, M, F, W, Y

### Positively Charged (for PPC)
K (+1), R (+1), H (+0.5)

### Negatively Charged (for PNC)
D (-1), E (-1)

## Output Example

```
================================================================================
  TAP Analysis: C45H4_coVL_rank_1
================================================================================

  Name: C45H4_coVL_rank_1
  Numbering: ANARCI IMGT
  Surface: MSMS
  VH length: 120 residues
  VL length: 107 residues

================================================================================
  TAP METRICS
================================================================================

  Metric          Value        Flag           Threshold (Green)
  --------------- ------------ -------------- ------------------------------
  Ltot            54           ✓ GREEN        43-54
  PSH             145.67       ✓ GREEN        111-168
  PPC             2.34         ⚠ AMBER        0-1.33
  PNC             1.89         ✓ GREEN        0-1.98
  SFvCSP          -8.45        ⚠ AMBER        >-6.00

================================================================================
  CDR LENGTHS (IMGT)
================================================================================

  CDR         VH       VL
  ---------- -------- --------
  CDR1        8        12
  CDR2        8        3
  CDR3        13       11
  Total       29       26

================================================================================
  RISK ASSESSMENT
================================================================================

  Overall: LOW-MEDIUM RISK
  Amber flags: PPC, SFvCSP
```

## Integration with Pipeline

This tool is independent of the main humanization pipeline. It can be used:

1. **Before humanization**: Assess developability risks in donor sequences
2. **After humanization**: Validate that optimized variants have acceptable TAP profiles
3. **For external structures**: Analyze any antibody PDB/CIF file
4. **Comparative analysis**: Compare TAP profiles across variants

## Notes

- Structure-based analysis provides accurate PSH/PPC/PNC/SFvCSP values
- Sequence-only analysis can only calculate Ltot; other metrics require structure
- TAP thresholds are based on 851 post Phase-I clinical-stage therapeutics (2025 update)
- SFvCSP measures charge asymmetry between VH and VL chains
- More negative SFvCSP values indicate greater charge asymmetry (higher risk)

## References

1. Raybould MIJ, et al. Five computational developability guidelines for therapeutic antibody profiling. PNAS. 2019.
2. Raybould MIJ, et al. Contextualising the developability risk of antibodies with lambda light chains using enhanced therapeutic antibody profiling. Communications Biology. 2024.
