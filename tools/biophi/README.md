# BioPhi Humanness Analysis Tool

Standalone wrapper around [BioPhi](https://github.com/Merck/BioPhi) for independent humanness cross-checks and comprehensive variant evaluation. Two entry points:

| Script | Role |
|--------|------|
| `biophi_analyze.py` | CLI wrapper around the `biophi` binary (Sapiens mean, OASis identity, variant evaluation with cross-validation) |
| `biophi_comprehensive.py` | **Full BioPhi capability reproduction via the Python API** — every analysis defaults to the fullest output (Sapiens matrix/embeddings/best/masked, OASis peptides/curves/frequencies/XLSX, all four humanization modes, OASis before/after, cross-validation) |

- **Sapiens** — human antibody language model; mean per-residue probability (higher = more human-like), full 20-aa probability matrix, per-position scores, best-score/masked sequences, residue + sequence embeddings, and Sapiens-guided humanized FASTA.
- **OASis** — fraction of 9-mer peptides observed in human antibody repertoires (requires the local OASis peptide database), plus peptide-level tables, positional humanness, OASis curves, germline matches and residue-frequency tables.
- **Humanization** — Sapiens, CDR grafting, Manual/Designer, and coordinated full-antibody (VH+VL) modes with alignments, mutation lists and OASis before/after.
- **Evaluate** — comprehensive evaluation of pipeline variant panels with cross-validation against pipeline metrics.

This tool is deliberately **independent of `scripts/humanize/`**. It runs BioPhi from a dedicated conda environment, so the portable pipeline keeps its stdlib-only core and is not pinned to BioPhi's old ANARCI/abnumber stack.

> The existing `scripts/humanize/humanness.py` adapter is not used or modified by this tool. Do not edit the pipeline when changing this tool.

## Environment

BioPhi needs Python 3.9 and an older ANARCI/abnumber than the main pipeline. Create a dedicated environment:

```bash
conda create -n biophi python=3.9 pip -y
conda activate biophi
conda install -c bioconda -c conda-forge \
    "anarci=2020.04.23" "abnumber=0.3.2" "hmmer>=3.1" -y

pip install "transformers<5" "sapiens==1.1.0" openpyxl
pip install "git+https://github.com/Merck/BioPhi.git@v1.0.11"
```

The wrapper calls `conda run -n biophi biophi ...`, so the environment does not need to be active when running the tool.

System `hmmsearch`/`hmmbuild`/`hmmscan` are also present in `/usr/bin`. Verify:

```bash
python3 tools/biophi/biophi_analyze.py check
```

## OASis database

Place the **unzipped** database here:

```text
data/biophi/OASis_9mers_v1.db
```

Download source (CC-BY-4.0):

```bash
wget https://zenodo.org/record/5164685/files/OASis_9mers_v1.db.gz
md5sum OASis_9mers_v1.db.gz
# expected: 28fa18135b3ef1710970f72b5953cbc5
gunzip OASis_9mers_v1.db.gz
```

The database is ~22 GB uncompressed. `data/biophi/` is gitignored.

## Sapiens weights (offline)

Sapiens loads two HuggingFace model repositories plus a tokenizer. On this server direct `huggingface.co` access is blocked, so the weights are already cached under `~/.cache/huggingface/hub/` using the `hf-mirror.com` endpoint:

```bash
export HF_ENDPOINT=https://hf-mirror.com
```

Manual download (if the cache has to be rebuilt on another machine):

```text
# VH model
https://huggingface.co/prihodad/biophi-sapiens1-vh/resolve/main/config.json
https://huggingface.co/prihodad/biophi-sapiens1-vh/resolve/main/model.safetensors

# VL model
https://huggingface.co/prihodad/biophi-sapiens1-vl/resolve/main/config.json
https://huggingface.co/prihodad/biophi-sapiens1-vl/resolve/main/model.safetensors

# Tokenizer
https://huggingface.co/prihodad/biophi-sapiens1-tokenizer/resolve/main/vocab.json
https://huggingface.co/prihodad/biophi-sapiens1-tokenizer/resolve/main/merges.txt
https://huggingface.co/prihodad/biophi-sapiens1-tokenizer/resolve/main/tokenizer_config.json
https://huggingface.co/prihodad/biophi-sapiens1-tokenizer/resolve/main/special_tokens_map.json
```

Mirror URLs are identical with `huggingface.co` replaced by `hf-mirror.com`.

To populate the Hub cache from local directories:

```bash
HF_ENDPOINT=https://hf-mirror.com \
  ~/miniconda3/envs/biophi/bin/huggingface-cli download \
  prihodad/biophi-sapiens1-vh
HF_ENDPOINT=https://hf-mirror.com \
  ~/miniconda3/envs/biophi/bin/huggingface-cli download \
  prihodad/biophi-sapiens1-vl
HF_ENDPOINT=https://hf-mirror.com \
  ~/miniconda3/envs/biophi/bin/huggingface-cli download \
  prihodad/biophi-sapiens1-tokenizer
```

Weights are small (~2.3 MB per model); the total cache is ~4.5 MB.

The wrapper defaults to **offline mode** (`HF_HUB_OFFLINE=1`, `TRANSFORMERS_OFFLINE=1`) so a missing cache fails clearly instead of hanging on a blocked network. Pass `--hf-endpoint https://hf-mirror.com` to allow a one-time download.

## Usage

### Check everything

```bash
python3 tools/biophi/biophi_analyze.py check \
    --oasis-db data/biophi/OASis_9mers_v1.db

# same check for the comprehensive tool (auto-enters the conda env)
python3 tools/biophi/biophi_comprehensive.py check
```

Outputs `problems`/`warnings` and exits non-zero only on hard problems (missing conda env, missing HMMER, gzipped DB).

### Sapiens mean score (no OASis DB required)

```bash
python3 tools/biophi/biophi_analyze.py sapiens \
    --input data/examples/mouse_4d5_fab.fasta \
    --outdir outputs/biophi
```

By default this now also computes **per-position scores + the full 20-aa matrix**
(`*_biophi_scores.csv`, works around the BioPhi CLI `--scores-only` bug by
running the Python API inside the conda env). Add `--humanize` to also write
`*_biophi_humanized.fa`; `--no-per-position` skips the matrix.

### OASis identity (needs the DB)

```bash
python3 tools/biophi/biophi_analyze.py oasis \
    --input data/examples/mouse_4d5_fab.fasta \
    --oasis-db data/biophi/OASis_9mers_v1.db \
    --min-percent-subjects 10 \
    --outdir outputs/biophi
```

For a Fab, BioPhi pairs heavy and light chains only when both FASTA records share the same ID with an optional `_VH`/`_VL` or `_HC`/`_LC` suffix. A single heavy chain is treated as a VHH-like chain.

### Evaluate pipeline variants (comprehensive)

Comprehensive BioPhi evaluation of pipeline-generated variant panels. Requires the pipeline's `variants.fasta`, `variants_numbering.json`, and optionally `humanization_result.json` for cross-validation.

```bash
python3 tools/biophi/biophi_analyze.py evaluate-variants \
    --input outputs/step2/variants.fasta \
    --metadata outputs/step2/variants_numbering.json \
    --pipeline-result outputs/step2/humanization_result.json \
    --outdir outputs/biophi_eval \
    --oasis-db data/biophi/OASis_9mers_v1.db \
    --verbose --timeout 7200
```

**Options:**

| Flag | Description |
|------|-------------|
| `--input` | Variants FASTA from pipeline (`variants.fasta`) |
| `--metadata` | Variants numbering JSON (`variants_numbering.json`) |
| `--pipeline-result` | Pipeline result JSON (`humanization_result.json`) for cross-validation |
| `--outdir` | Output directory (default: `outputs/biophi_eval`) |
| `--run-sapiens` / `--no-sapiens` | Run/skip Sapiens scoring (default: on) |
| `--run-oasis` / `--no-oasis` | Run OASis 9-mer identity (requires `--oasis-db`, default: on) |
| `--run-humanize` / `--no-humanize` | Generate Sapiens-humanized FASTA (default: **on**) |
| `--heatmap` / `--no-heatmap` | Sapiens per-position heatmap JSON mapped to Kabat (default: **on**) |
| `--iterations` / `--humanize-cdrs` | Sapiens humanization options (default: 1 iteration, CDRs fixed) |
| `--min-percent-subjects` | OASis threshold percent (default: 10) |
| `--verbose` | Stream BioPhi subprocess output in real-time |
| `--timeout` | Subprocess timeout in seconds (default: 3600; recommend 7200+) |
| `--scheme` / `--cdr-definition` | Numbering scheme & CDR definition (kabat/chothia/imgt/aho) |

**Typical run (13 variants):**

```bash
python3 tools/biophi/biophi_analyze.py evaluate-variants \
    --input outputs/step2/variants.fasta \
    --metadata outputs/step2/variants_numbering.json \
    --pipeline-result outputs/step2/humanization_result.json \
    --outdir outputs/biophi_eval \
    --oasis-db data/biophi/OASis_9mers_v1.db \
    --verbose --timeout 7200 \
    2>&1 | tee biophi_eval.log
```

> **Note:** OASis is enabled by default (adds ~5-15 min). Use `--no-oasis` to skip. First run is slow (model loading + inference). Use `--verbose` to monitor progress. Increase `--timeout` for large panels.

## Full capability reproduction: `biophi_comprehensive.py`

Reproduces **every BioPhi analysis** through the Python API. The script
auto-re-execs into the `biophi` conda env, so plain `python3` works from the
repo root. **All features default to the fullest output**; disable with the
`--no-*` flags.

| Subcommand | What it runs (defaults) |
|------------|-------------------------|
| `check` | env versions, Sapiens cache, HMMER, OASis DB |
| `sapiens` | mean + per-position + full 20-aa matrix + best-score + masked + residue/sequence embeddings + humanized FASTA + heatmap |
| `humanness` | OASis identity/percentile/germline content, V/J germlines, peptide table, positional humanness, OASis curve, top frequencies, chain-type & germline-family residue frequencies, paired-antibody metrics, full XLSX workbook |
| `humanize` | all four modes: `sapiens`, `cdr_grafting`, `manual` (Designer), `full_antibody` (VH+VL) — with alignments, mutation lists, top proposals, germline candidates, OASis before/after, per-mode FASTAs |
| `cdr-grafting` | CDR grafting only (germline auto/select, vernier back-mutation, optional Sapiens pass) |
| `full-antibody` | coordinated VH+VL humanization only |
| `evaluate-all` | **everything**: sapiens + humanness + all humanization modes + heatmap + embeddings + scores/peptides CSV + OASis XLSX + humanized FASTAs + pipeline cross-validation |

```bash
# Full analysis of a FASTA (Fab pair or single VHH)
python3 tools/biophi/biophi_comprehensive.py evaluate-all \
    --input data/examples/mouse_4d5_fab.fasta \
    --outdir outputs/biophi_full

# Pipeline variant panel with numbering + cross-validation
python3 tools/biophi/biophi_comprehensive.py evaluate-all \
    --input outputs/step2/variants.fasta \
    --metadata outputs/step2/variants_numbering.json \
    --pipeline-result outputs/step2/humanization_result.json \
    --outdir outputs/biophi_full \
    --verbose

# Individual capabilities
python3 tools/biophi/biophi_comprehensive.py sapiens --input ab.fasta
python3 tools/biophi/biophi_comprehensive.py humanness --input ab.fasta \
    --oasis-db data/biophi/OASis_9mers_v1.db
python3 tools/biophi/biophi_comprehensive.py humanize --input ab.fasta
```

Common flags: `--scheme/--cdr-definition`, `--oasis-db`, `--min-fraction-subjects`
(0.1 = 10%), `--iterations`, `--humanize-cdrs`, `--modes` (comma list or `all`),
`--heavy-v-germline/--light-v-germline`, `--sapiens-iterations`,
`--no-embeddings`, `--no-humanization`, `--no-oasis`, `--limit N`, `--verbose`,
`--hf-endpoint`.

## Outputs

### Sapiens / OASis (single commands)

For an input `mouse_4d5_fab.fasta`:

| File | Contents |
|------|----------|
| `mouse_4d5_fab_biophi_sapiens.json` | provenance + per-record Sapiens score + per-position |
| `mouse_4d5_fab_biophi_sapiens.md` | human-readable Sapiens table |
| `mouse_4d5_fab_biophi_sapiens.csv` | `id,length,chain,sapiens_score,...` |
| `mouse_4d5_fab_biophi_scores.csv` | full position × 20-aa probability matrix |
| `mouse_4d5_fab_biophi_humanized.fa` | Sapiens-humanized FASTA (with `--humanize`) |
| `mouse_4d5_fab_biophi_oasis.xlsx` | raw BioPhi OASis workbook |
| `mouse_4d5_fab_biophi_oasis.json` | provenance + parsed Overview sheet |
| `mouse_4d5_fab_biophi_oasis.md` | human-readable OASis summary |
| `mouse_4d5_fab_biophi_oasis.csv` | parsed Overview sheet |

### Evaluate Variants (`biophi_analyze.py`)

For pipeline variants (`variants.fasta`):

| File | Contents |
|------|----------|
| `{stem}_biophi_eval.json` | **Complete machine-readable report**: all variants with Sapiens mean/per-position, OASis, pipeline metadata, cross-validation |
| `{stem}_biophi_eval.md` | Human-readable Markdown report with variant summary, per-variant details, cross-validation tables |
| `{stem}_biophi_eval.csv` | Variant-level summary CSV (incl. pipeline back-mutation counts/CVI) |
| `{stem}_biophi_scores.csv` | full Sapiens probability matrix |
| `{stem}_biophi_heatmap.json` | Sapiens per-position heatmap mapped to Kabat positions (default on) |
| `{stem}_biophi_humanized.fa` | Sapiens-humanized FASTA (default on) |
| `{stem}_biophi_oasis.xlsx` | Raw BioPhi OASis workbook (when `--run-oasis`) |

### Evaluate All (`biophi_comprehensive.py evaluate-all`)

| File | Contents |
|------|----------|
| `{stem}_biophi_comprehensive.json` | everything: sapiens (mean/matrix/best/masked), OASis per-chain + per-antibody, all humanization modes, cross-validation |
| `{stem}_biophi_comprehensive.md` | full Markdown report |
| `{stem}_biophi_comprehensive.csv` | variant summary CSV |
| `{stem}_biophi_scores.csv` | full Sapiens probability matrix |
| `{stem}_biophi_peptides.csv` | OASis peptide-level table (all chains) |
| `{stem}_biophi_oasis.xlsx` | OASis workbook (Overview/VH/VL/curves/per-antibody sheets) |
| `{stem}_biophi_heatmap.json` | Kabat-mapped per-position Sapiens scores |
| `{stem}_biophi_embeddings.json` | residue + sequence embeddings (all layers) |
| `{stem}_biophi_humanized_sapiens.fa` | Sapiens-humanized sequences |
| `{stem}_biophi_humanized_cdr_grafting.fa` | CDR-grafted sequences |
| `{stem}_biophi_humanized_manual.fa` | Designer (unchanged) sequences |
| `{stem}_biophi_humanized_full_antibody.fa` | paired VH+VL humanized sequences |

### Key report sections (Markdown)

1. **Variant Summary** — Sapiens mean, OASis identity/percentile, germline content, pipeline back-mutations, humanization mutation counts
2. **Per-Variant Details** — Sapiens hotspots, OASis peptides/germlines, humanization before/after OASis, pipeline metrics
3. **Cross-Validation: BioPhi vs Pipeline** — side-by-side metrics + Spearman rank correlations + rank disagreements

## Interpretation

- **Sapiens score**: mean probability of the input residues under the human antibody language model. Higher = more human-like. It is **not** a germline identity percentage and should be read as a relative cross-check between variants, not as an absolute pass/fail threshold.
- **Sapiens per-position**: residue-level human-likeness probabilities. Low scores pinpoint "local non-human" hotspots; map to Kabat positions via `variants_numbering.json` for targeted back-mutation decisions.
  - **Note**: The BioPhi CLI `--scores-only` has a known bug (empty position/residue fields). This tool works around it by calling the Sapiens Python API directly (`sapiens.predict.predict_scores`), which produces correct per-position scores.
- **Sapiens-guided humanization**: BioPhi's own humanization proposals (`--humanize`), usable as a third candidate route (`V_biophi`) alongside germline-based and structure-based designs.
- **OASis Identity**: fraction of overlapping 9-mers found in OAS repertoires of at least `--min-percent-subjects` human subjects. Higher = more human.
- **OASis Percentile**: position of that identity against background antibodies in the database.
- **Cross-validation**: compares BioPhi metrics (Sapiens/OASis) against pipeline's internal humanness/CVI/composite scores to flag discrepancies.

## Notes and limitations

- OASis queries a ~22 GB SQLite database. Large panels can be slow; run per-antibody or use `--summary`.
- VHH is processed as a heavy chain. Sapiens/OASis human-chain statistics are a reasonable approximation but are not VHH-specific; label VHH results as reference-only.
- First run is slow (model loading + inference). Subsequent runs reuse cached models and are much faster.
- Sapiens per-position output uses the Python API (`sapiens.predict.predict_scores`) to work around a BioPhi CLI bug (`--scores-only` produces empty position/residue fields). The CLI `--mean-score-only` and `--fasta-only` work correctly; the helper runs inside the conda env so the wrapper itself works under system `python3`.
- `biophi_comprehensive.py` needs the `biophi` conda env: it re-execs itself with the env interpreter automatically (or run it as `conda run -n biophi python tools/biophi/biophi_comprehensive.py ...`).
- Use `--verbose` to stream subprocess output in real-time (recommended for long runs).
- OASis DB is CC-BY-4.0. Do not commit it to Git.
- BioPhi is MIT-licensed. This wrapper is independent of the pipeline and only invokes the `biophi` CLI / BioPhi Python API.
