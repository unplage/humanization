# Humanness Score Tool

Evaluate the humanization degree of antibody sequences by comparing against human germline V genes (CDR3/FR4 are excluded from headline scores).

## Features

- **Dual Numbering**: Both Kabat and IMGT numbering schemes (IMGT is mandatory)
- **Germline Comparison**: Compare against all human V genes
- **Druggability Tie-break**: When several germlines give identical (FR, CDR)
  identity, the tool deterministically prefers established druggable frameworks
  (VH: IGHV1-69 > IGHV3-23 > IGHV1-46; VK: IGKV1-39 > IGKV3-20 > IGKV4-1),
  then germline usage frequency, then gene id (see `DRUGGABLE_*_PRIORITY` in
  `evaluate_humanness.py`)
- **FR/CDR Identity**: Calculate framework and CDR region identities
- **Humanization Classification**: Academic estimate based on sequence identity
- **Batch Processing**: Analyze multiple sequences from FASTA files
- **Word Reports**: Generate professional Word documents

## Requirements

```bash
pip install biopython anarci python-docx
```

**ANARCI is required** for IMGT numbering (WHO/INN/USAN standard). The tool will error out if ANARCI is not installed.

## Usage

### Single Sequence Evaluation

```bash
# From FASTA file
python3 tools/humanness_score/evaluate_humanness.py \
    --input data/examples/mouse_4d5_fab.fasta

# Direct sequences
python3 tools/humanness_score/evaluate_humanness.py \
    --vh EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS \
    --vl DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
```

### Batch Processing (Word Reports)

```bash
# Positional FASTA (default output: outputs/reports)
python3 tools/humanness_score/generate_all_reports.py \
    data/examples/antibody_panel.fasta

# Explicit output directory
python3 tools/humanness_score/generate_all_reports.py \
    data/examples/antibody_panel.fasta \
    --outdir outputs/humanness
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
    Compared:     98/98 query positions (coverage 100.0%); query-only 0, germline-only 0

  Detailed region breakdown (Kabat):
  Region     Identity     Match    Total
  ---------- ------------ -------- --------
  FR1        88.0%     22       25
  FR2        86.7%     13       15
  ...

  Humanization estimate: Academic: likely Humanized-like
  (Academic estimate only; official USAN/INN naming requires manufacturing-process knowledge)

======================================================================
  H chain (120 residues) - IMGT numbering
======================================================================
  (Per WHO/INN/USAN standards for drug naming)

  Top 10 human germline matches (IMGT):
  Rank   Gene                 FR id      CDR id     Overall
  ------ -------------------- ---------- ---------- ----------
  1      IGHV1-69*01          0.8421     0.4444     0.7451
  ...

  Best match (IMGT): IGHV1-69*01
    FR identity:  0.8421 (84.2%)
    CDR identity: 0.4444 (44.4%)
    Overall:      0.7451 (74.5%)
    Compared:     111/112 query positions (coverage 99.1%); query-only 1, germline-only 0

  IMGT region breakdown:
  Region     Identity     Match    Total
  ---------- ------------ -------- --------
  FR1        88.5%     23       26
  CDR1       42.9%     3        7
  ...

  Humanization estimate: Academic: likely Humanized-like
  (Academic estimate only; official USAN/INN naming requires manufacturing-process knowledge)
```

### Word Reports

- Individual report per antibody (Kabat + IMGT sections)
- Combined summary report with all antibodies
- Report filenames are sanitized; duplicate FASTA names get unique suffixes
- Invalid FASTA records are skipped and reported; the command exits non-zero if any failed

## Numbering Schemes

### Kabat Numbering
- VH: FR1 1-30, CDR1 31-35, FR2 36-49, CDR2 50-65, FR3 66-94, CDR3 95-102, FR4 103-113
- VL: FR1 1-23, CDR1 24-34, FR2 35-49, CDR2 50-56, FR3 57-88, CDR3 89-97, FR4 98-107
- Used internally by the pipeline

### IMGT Numbering (Mandatory)
- FR1: 1-26, CDR1: 27-38, FR2: 39-55, CDR2: 56-65, FR3: 66-104, CDR3: 105-117, FR4: 118-128
  (identical for VH and VL, Lefranc 2003)
- **Used for WHO/INN/USAN drug naming** — this is the standard for regulatory submissions
- Requires ANARCI (installed via `pip install anarci`) and the bundled IMGT germline database
  (`data/germline/abnumber_human_imgt.json`)

## Humanization Classification

**Note**: This is an **academic estimate** based on sequence identity. Official USAN/INN naming is based on the technology used to create the antibody.

| FR Identity | Classification | Typical Origin |
|-------------|----------------|----------------|
| ≥ 95% | Human-like | Transgenic mice, phage display |
| ≥ 85% | Humanized-like | CDR grafting |
| ≥ 70% | Chimeric-like | Mouse V + human C |
| < 70% | Murine-like | Full mouse antibody |

## Evaluation Criteria

- **Excluded Regions**: CDR3 (donor-grafted) and FR4 (human J gene, highly conserved across species) — both excluded from the headline scores so results stay comparable with published V-region-only humanness metrics
- **FR Identity**: Identity in framework regions (FR1+FR2+FR3)
- **CDR Identity**: Identity in CDR1+CDR2 regions
- **Overall Identity**: Pooled identity over FR1+FR2+FR3+CDR1+CDR2
- **Match/Count**: Reported region values are exact position counts (`match/total`), not rounded back-computed percentages
- **Coverage**: Headline identity uses positions present in both query and germline; coverage/query-only/germline-only are reported separately so indels are visible without changing the historical metric
- **Note**: The terminal "Sequence identity summary" table is a raw per-region diagnostic that *includes* CDR3 and FR4 and is aligned by position label, so its numbers will not reconcile with the headline scores

## Error Handling

The tool **requires** ANARCI and the IMGT germline database for IMGT evaluation. If either is missing, the tool will error out with instructions:

```
ERROR: ANARCI not installed. Required for IMGT numbering (WHO/INN/USAN standard).
  Install with: pip install anarci
```

or

```
ERROR: IMGT germline database unavailable. Required for IMGT evaluation.
  Expected: data/germline/abnumber_human_imgt.json
```

There is **no `--no-imgt` option** — IMGT evaluation is mandatory for WHO/INN/USAN compliance.

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