# BioPhi Humanness Analysis Tool

Standalone wrapper around [BioPhi](https://github.com/Merck/BioPhi) for two
independent humanness cross-checks:

- **Sapiens** — human antibody language model; mean per-residue probability
  (higher = more human-like), plus optional Sapiens-guided humanized FASTA.
- **OASis** — fraction of 9-mer peptides observed in human antibody
  repertoires (requires the local OASis peptide database).

This tool is deliberately **independent of `scripts/humanize/`**. It runs
BioPhi from a dedicated conda environment, so the portable pipeline keeps its
stdlib-only core and is not pinned to BioPhi's old ANARCI/abnumber stack.

> The existing `scripts/humanize/humanness.py` adapter is not used or modified
> by this tool. Do not edit the pipeline when changing this tool.

## Environment

BioPhi needs Python 3.9 and an older ANARCI/abnumber than the main pipeline.
Create a dedicated environment:

```bash
conda create -n biophi python=3.9 pip -y
conda activate biophi
conda install -c bioconda -c conda-forge \
    "anarci=2020.04.23" "abnumber=0.3.2" "hmmer>=3.1" -y

pip install "transformers<5" "sapiens==1.1.0" openpyxl
pip install "git+https://github.com/Merck/BioPhi.git@v1.0.11"
```

The wrapper calls `conda run -n biophi biophi ...`, so the environment does
not need to be active when running the tool.

System `hmmsearch`/`hmmbuild`/`hmmscan` are also present in `/usr/bin`.
Verify:

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

Sapiens loads two HuggingFace model repositories plus a tokenizer. On this
server direct `huggingface.co` access is blocked, so the weights are already
cached under `~/.cache/huggingface/hub/` using the `hf-mirror.com` endpoint:

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

The wrapper defaults to **offline mode** (`HF_HUB_OFFLINE=1`,
`TRANSFORMERS_OFFLINE=1`) so a missing cache fails clearly instead of hanging
on a blocked network. Pass `--hf-endpoint https://hf-mirror.com` to allow a
one-time download.

## Usage

### Check everything

```bash
python3 tools/biophi/biophi_analyze.py check \
    --oasis-db data/biophi/OASis_9mers_v1.db
```

Outputs `problems`/`warnings` and exits non-zero only on hard problems
(missing conda env, missing HMMER, gzipped DB).

### Sapiens mean score (no OASis DB required)

```bash
python3 tools/biophi/biophi_analyze.py sapiens \
    --input data/examples/mouse_4d5_fab.fasta \
    --outdir outputs/biophi
```

Add `--humanize` to also write `*_biophi_humanized.fa`.

### OASis identity (needs the DB)

```bash
python3 tools/biophi/biophi_analyze.py oasis \
    --input data/examples/mouse_4d5_fab.fasta \
    --oasis-db data/biophi/OASis_9mers_v1.db \
    --min-percent-subjects 10 \
    --outdir outputs/biophi
```

For a Fab, BioPhi pairs heavy and light chains only when both FASTA records
share the same ID with an optional `_VH`/`_VL` or `_HC`/`_LC` suffix. A single
heavy chain is treated as a VHH-like chain.

## Outputs

For an input `mouse_4d5_fab.fasta`:

| File | Contents |
|------|----------|
| `mouse_4d5_fab_biophi_sapiens.json` | provenance + per-record Sapiens score |
| `mouse_4d5_fab_biophi_sapiens.md` | human-readable Sapiens table |
| `mouse_4d5_fab_biophi_sapiens.csv` | `id,length,chain,sapiens_score` |
| `mouse_4d5_fab_biophi_oasis.xlsx` | raw BioPhi OASis workbook |
| `mouse_4d5_fab_biophi_oasis.json` | provenance + parsed Overview sheet |
| `mouse_4d5_fab_biophi_oasis.md` | human-readable OASis summary |
| `mouse_4d5_fab_biophi_oasis.csv` | parsed Overview sheet |

## Interpretation

- **Sapiens score**: mean probability of the input residues under the human
  antibody language model. Higher = more human-like. It is **not** a germline
  identity percentage and should be read as a relative cross-check between
  variants, not as an absolute pass/fail threshold.
- **OASis Identity**: fraction of overlapping 9-mers found in OAS repertoires
  of at least `--min-percent-subjects` human subjects. Higher = more human.
- **OASis Percentile**: position of that identity against background
  antibodies in the database.

## Notes and limitations

- OASis queries a ~22 GB SQLite database. Large panels can be slow; run
  per-antibody or use `--summary`.
- VHH is processed as a heavy chain. Sapiens/OASis human-chain statistics are
  a reasonable approximation but are not VHH-specific; label VHH results as
  reference-only.
- OASis DB is CC-BY-4.0. Do not commit it to Git.
- BioPhi is MIT-licensed. This wrapper is independent of the pipeline and only
  invokes the `biophi` CLI.
