# CDR Homology Analysis Tool

Compare antibody CDR sequences against a clinical antibody database
(TheraSAbDab) to find homologous therapeutics.

## Features

- **CDR extraction**: ANARCI IMGT numbering, CDR1 27-38 / CDR2 56-65 / CDR3 105-117
  (identical to `scripts/humanize/numbering.py::NATIVE_IMGT_CDRS`)
- **Length-aware alignment**: BLOSUM62 with *global* identity for CDR1/CDR2 and
  *local* (Smith-Waterman) core identity for CDR3, plus coverage and a
  length-normalised score
- **Statistical context**: background median identity and top-hit percentile
  per CDR, so a hit can be judged against the database rather than in isolation
- **VHH aware**: single-domain queries are detected from the FR2 hallmark
  (IMGT 42/49/50/52) and labelled instead of being silently reported as VH
- **Provenance**: every report records the timestamp, source-file SHA-256,
  database SHA-256, tool/library versions and every parameter used
- **Multiple reports**: Markdown, Word, JSON and per-antibody CSV

## Requirements

System packages (ANARCI shells out to HMMER):

```bash
# Debian/Ubuntu
sudo apt-get install hmmer
# conda
conda install -c bioconda hmmer
```

Python packages:

```bash
pip install pandas openpyxl biopython anarci python-docx
```

`anarci` alone is not sufficient: `hmmsearch` and `hmmbuild` must be on `PATH`
or numbering fails for every sequence. Verify with:

```bash
hmmsearch -h | head -1
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

Accepts the pipeline header convention `>AMG110.H|H_V0|desc`, plain
`>name_VH` / `>name_VL` FASTAs, and single-chain (VHH) inputs. Output records
are named `>AMG110_V0_VH` / `>AMG110_V0_VL`, which the query tool pairs
automatically.

The grouping key is **antibody name + ladder index**, so several antibodies in
one file (e.g. `>AMG110.H|H_V0` and `>OTHER.H|H_V0`) stay separate. A
duplicate chain for the same group is reported as a warning (exit code 2)
instead of being silently dropped.

### Step 1: Preprocess the database

```bash
python3 tools/cdr_homology/preprocess.py \
    --database data/benchmarks/TheraSAbDab_SeqStruc_OnlineDownload.xlsx \
    --output data/therasabdab_cdrs.json
```

This generates:

- `data/therasabdab_cdrs.json` - CDR sequences for all clinical antibodies
  (primary and bispecific chains)
- `data/therasabdab_stats.json` - database statistics **and the provenance of
  the source spreadsheet** (file name, SHA-256, row count, timestamp,
  ANARCI/Biopython versions)

Rows whose sequence ANARCI cannot number are reported on stderr and counted in
the statistics rather than being dropped silently.

### Step 2: Query

```bash
# Single query (VH + VL)
python3 tools/cdr_homology/query.py \
    --vh QVQLVQSGAEVKKPGASVKVSCKAS... \
    --vl DIQMTQSPSSLSASVGDRVTITC... \
    --name my_antibody \
    --top 5 \
    --output outputs/cdr_homology/

# From FASTA file
python3 tools/cdr_homology/query.py \
    --input query.fasta \
    --top 5 \
    --output outputs/cdr_homology/
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--top N` | 5 | Hits per CDR; must be >= 1 |
| `--align` | `auto` | `auto` = local for CDR3, global for CDR1/CDR2; or force `global` / `local` |
| `--min-score` | 0.0 | Minimum BLOSUM62 score for a hit to be reported |
| `--min-identity` | 0.0 | Minimum identity (0-1) for a hit to be reported |
| `--markdown` | off* | Write `cdr_homology_report.md` |
| `--json` | off* | Write `cdr_homology_results.json` |
| `--csv` | off* | Write `detail/*.csv` |
| `--docx` | off* | Write `cdr_homology_report.docx` |
| `--all-formats` | - | Write every format |

\* When no format flag is given, all formats are written.

## Output Files

| File | Description |
|------|-------------|
| `cdr_homology_report.md` | Markdown report (provenance, method, summary, details, appendix) |
| `cdr_homology_report.docx` | Word report mirroring the Markdown sections |
| `cdr_homology_results.json` | `{metadata, queries}`; metadata carries full provenance |
| `detail/*.csv` | Per-antibody results, collision-free file names |

## How identity is computed

| CDR | Alignment | `identity` |
|-----|-----------|------------|
| CDR1, CDR2 | global (Needleman-Wunsch) | identical residues / full alignment length |
| CDR3 | local (Smith-Waterman) | identical residues / homologous core length |

CDR3 loops range from 3 to more than 30 residues. A global denominator divides
the identity by the length difference, which systematically under-reports
homology: a query CDR3 of 8 residues found intact inside a 13-residue subject
CDR3 scores 8/13 = 62% globally but 100% over its core. The tool therefore
reports, for every hit:

- `identity` - the headline value defined above
- `identity_global` - the full-length value, for comparison
- `identity_over_query` - identical residues / query length ("how much of my
  loop is shared")
- `coverage` - fraction of the query CDR in the alignment
- `score` / `score_norm` - raw BLOSUM62 score and score / query self-score

The report also gives the **background median identity** (the query CDR against
every database CDR) and the **top-hit percentile**, so a 70% CDR3 hit can be
read against the database baseline instead of in a vacuum.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Fatal (database missing, no queries, unreadable input) |
| 2 | Partial: at least one chain could not be numbered or a CDR was empty. Reports are still written and the affected regions are marked N/A / listed under **Errors** |

A missing `python-docx` never aborts the run: the JSON, CSV and Markdown are
written first, the Word report is skipped with a warning.

## FASTA Format

The tool expects FASTA files with VH/VL pairs. Use suffixes to identify chains:

```fasta
>Antibody1_VH
QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWIG...
>Antibody1_VL
DIQMTQSPSSLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLI...
```

Or use header keywords: `VH`/`Heavy`/`IGH` → heavy; `VL`/`Light`/`Kappa`/
`Lambda`/`IGK`/`IGL` → light. When the header says nothing, ANARCI's own HMM
chain call (H / K / L) is used; sequence length is only a last resort and is
flagged in the JSON (`chain_source`).

## Report Contents

### Summary
| # | Antibody | Type | VH_CDR3 best | VH_CDR3 id | VL_CDR3 best | VL_CDR3 id | CDR3 mean id |

`CDR3 mean id` aggregates only the CDR3 hits. CDR1 and CDR2 are highly
conserved across clinical antibodies (a large fraction of the database is
already >= 80% identical to any given CDR1), so averaging them in would dilute
the only discriminative signal.

### Per-antibody details
- Chain type (Fab / VHH / `vhh_suspect`)
- CDR sequences (IMGT numbering) with coverage and background median
- Top-N matches per CDR with chain origin (primary vs bispecific)

### Appendix
- A. Antibody format distribution
- B. Genetics / humanization strategy
- C. Light chain type
- D. Highest clinical trial phase
- E. CDR length distribution (all chains, primary and bispecific)
- F. Top 30 therapeutic targets
- G. Top 30 drugs by entry count

## Example

```bash
# Analyse 4D5 (trastuzumab parent) from the bundled example
python3 tools/cdr_homology/create_paired_fasta.py \
    --input data/examples/mouse_4d5_fab.fasta \
    --output /tmp/4d5_paired.fasta --name 4D5

python3 tools/cdr_homology/query.py \
    --input /tmp/4d5_paired.fasta \
    --top 5 \
    --output outputs/cdr_homology_4d5/
```

Output shows the most similar clinical antibodies for each CDR region, helping
assess developability and the regulatory landscape.

## Caveat

Sequence homology is not functional homology. A high CDR3 identity means the
loop is shared with a clinical antibody; it says nothing about epitope,
affinity or mechanism. Use the result as a regulatory/developability signal,
not as evidence of functional equivalence.
