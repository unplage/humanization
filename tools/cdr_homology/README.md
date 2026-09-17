# CDR Homology Analysis Tool

Compare antibody CDR sequences against a clinical antibody database (TheraSAbDab) to find homologous therapeutics.

## Features

- **CDR Extraction**: Uses ANARCI IMGT numbering to extract CDR1/CDR2/CDR3
- **Sequence Alignment**: BLOSUM62-based global alignment for accurate identity calculation
- **Database Search**: Compare against 1100+ clinical antibodies from TheraSAbDab
- **Multiple Reports**: Markdown, Word, and CSV output formats
- **Batch Processing**: Analyze multiple antibodies from FASTA files

## Requirements

```bash
pip install pandas biopython anarci python-docx openpyxl
```

## Workflow

### Step 0: Build a paired FASTA from pipeline output (optional)

If you have a humanization `variants.fasta`, generate a paired VH/VL FASTA
(one pair per ladder index V0, V1, V2a, ...):

```bash
python3 tools/cdr_homology/create_paired_fasta.py \
    --input outputs/amg110_humanization_step3/step3/variants.fasta \
    --output outputs/amg110_humanization_step3/step3/analysis_input/variants_paired.fasta \
    --name AMG110
```

This accepts the pipeline header convention `>AMG110.H|H_V0|desc`, plain
`>name_VH` / `>name_VL` FASTAs, and single-chain (VHH) inputs. Output records
are named `>AMG110_V0_VH` / `>AMG110_V0_VL`, which the query tool pairs
automatically.

### Step 1: Preprocess Database

First, preprocess the TheraSAbDab Excel database to extract CDRs:

```bash
python3 tools/cdr_homology/preprocess.py \
    --database data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx \
    --output data/therasabdab_cdrs.json
```

This generates:
- `data/therasabdab_cdrs.json` - CDR sequences for all clinical antibodies
- `data/therasabdab_stats.json` - Database statistics

### Step 2: Query

Run CDR homology search:

```bash
# Single query (VH + VL)
python3 tools/cdr_homology/query.py \
    --vh QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWIG... \
    --vl DIQMTQSPSSLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLI... \
    --top 5 \
    --output outputs/cdr_homology/

# From FASTA file
python3 tools/cdr_homology/query.py \
    --input query.fasta \
    --top 5 \
    --output outputs/cdr_homology/
```

## Output Files

| File | Description |
|------|-------------|
| `cdr_homology_report.md` | Markdown summary report |
| `cdr_homology_report.docx` | Word report with tables |
| `cdr_homology_results.json` | Raw JSON data |
| `detail/*.csv` | Per-antibody detailed results |

## FASTA Format

The tool expects FASTA files with VH/VL pairs. Use suffixes to identify chains:

```fasta
>Antibody1_VH
QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWIG...
>Antibody1_VL
DIQMTQSPSSLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLI...
```

Or use header keywords:
- VH/Heavy → heavy chain
- VL/Light/Kappa/Lambda → light chain

## Report Contents

### Summary Table
| # | Antibody | VH_CDR3 best | VH_CDR3 id | VL_CDR3 best | VL_CDR3 id | Mean id |
|---|----------|-------------|-----------|-------------|-----------|---------|

### Per-Antibody Details
- CDR sequences (IMGT numbering)
- Top N matches per CDR
- Alignment visualization

### Appendix
- Database statistics (format, genetics, targets)
- CDR length distribution
- Top therapeutic targets

## Example

```bash
# Analyze AMG110
python3 tools/cdr_homology/query.py \
    --input data/examples/amg110.fasta \
    --top 5 \
    --output outputs/cdr_homology_amg110/
```

Output shows most similar clinical antibodies for each CDR region, helping assess developability and regulatory landscape.
