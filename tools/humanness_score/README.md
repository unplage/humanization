# Humanness Score Tool

Evaluate the humanization degree of antibody sequences by comparing against human germline V/J genes.

## Features

- **Dual Numbering**: Both Kabat and IMGT numbering schemes
- **Germline Comparison**: Compare against all human V genes
- **FR/CDR Identity**: Calculate framework and CDR region identities
- **Humanization Classification**: Academic estimate based on sequence identity
- **Batch Processing**: Analyze multiple sequences from FASTA files
- **Word Reports**: Generate professional Word documents

## Requirements

```bash
pip install biopython anarci python-docx
```

## Usage

### Single Sequence Evaluation

```bash
# From FASTA file
python3 tools/humanness_score/evaluate_humanness.py \
    --input data/examples/amg110.fasta

# Direct sequences
python3 tools/humanness_score/evaluate_humanness.py \
    --vh EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS \
    --vl DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
```

### Batch Processing (Word Reports)

```bash
python3 tools/humanness_score/generate_all_reports.py \
    data/examples/antibody_panel.fasta
```

## Output

### Terminal Output

```
======================================================================
  H chain (120 residues) - Kabat numbering
======================================================================

  Top 10 human germline matches (Kabat):
  Rank   Gene                 FR id      CDR id     Overall
  ------ -------------------- ---------- ---------- ----------
  1      IGHV1-69*01          0.8500     0.4000     0.7500
  2      IGHV1-69*02          0.8400     0.4000     0.7400
  ...

  Best match (Kabat): IGHV1-69*01
    FR identity:  0.8500 (85.0%)
    CDR identity: 0.4000 (40.0%)

  Humanization estimate: Academic: likely Humanized-like
  (Note: Official USAN/INN naming requires knowledge of manufacturing process)
```

### Word Reports

- Individual report per antibody
- Combined summary report with all antibodies

## Numbering Schemes

### Kabat Numbering
- FR1: 1-30, CDR1: 31-35, FR2: 36-49, CDR2: 50-65, FR3: 66-94, CDR3: 95-102, FR4: 103-113
- Used internally by the pipeline

### IMGT Numbering
- FR1: 1-26, CDR1: 27-38, FR2: 39-55, CDR2: 56-65, FR3: 66-104, CDR3: 105-117, FR4: 118-128
- Used for WHO/INN/USAN drug naming

## Humanization Classification

**Note**: This is an **academic estimate** based on sequence identity. Official USAN/INN naming is based on the technology used to create the antibody.

| FR Identity | Classification | Typical Origin |
|-------------|----------------|----------------|
| ≥ 95% | Human-like | Transgenic mice, phage display |
| ≥ 85% | Humanized-like | CDR grafting |
| ≥ 70% | Chimeric-like | Mouse V + human C |
| < 70% | Murine-like | Full mouse antibody |

## Evaluation Criteria

- **Excluded Regions**: CDR3 (donor) and FR4 (human J gene)
- **FR Identity**: Identity in framework regions (FR1+FR2+FR3)
- **CDR Identity**: Identity in CDR1+CDR2 regions
- **Overall Identity**: Weighted average

## Example Output Files

```
outputs/reports/
├── Antibody1_report.docx      # Individual report
├── Antibody2_report.docx
└── combined_report.docx       # Summary report
```

## Integration with Pipeline

This tool is independent of the main humanization pipeline. Use it to:

1. **Pre-humanization**: Assess baseline humanization of donor sequences
2. **Post-humanization**: Validate humanization achieved by the pipeline
3. **Competitive analysis**: Compare with approved therapeutics
