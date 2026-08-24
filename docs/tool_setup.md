# Tool setup for server deployment (工具部署)

The pipeline is portable (stdlib-only) on the laptop. On the server, install
the optional tools in this order. `python3 scripts/humanize/cli.py setup-check`
reports what is available.

## 1. Core + exact numbering (recommended on every server)

```bash
conda create -n humanize python=3.11 -y
conda activate humanize
conda install -c bioconda -c conda-forge hmmer==3.3.2 -y   # ANARCI dependency
pip install anarci abnumber biopython
```

ANARCI provides exact IMGT/Chothia/Kabat numbering (HMMER-based) used to
cross-validate the portable engine's output (see validation.md).

## 2. IgBLAST + NCBI germline database (chosen germline source)

```bash
python3 scripts/humanize/cli.py setup-germline --dir data/germline
```

Downloads `ncbi-igblast-<ver>-x64-linux.tar.gz` from
`ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/LATEST/` and extracts
`human_gl_V.fasta` (IGHV), `human_gl_L.fasta` (IGKV+IGLV),
`human_gl_J.fasta` (IGHJ+IGKJ+IGLJ) plus the `igblastn` binary and
`internal_data`. The pipeline uses the FASTA files as its germline source
(identical origin to the bundled set). IgBLAST itself can additionally be used
for mutation analysis on the same data.

Note: NCBI FTP can be slow; the downloader retries. If it still fails, the
bundled `data/germline/human_germline_kabat.json` (373 V + 28 J human genes,
Kabat space, derived from the same IMGT data) is used automatically.

## 3. AlphaFold3 (structure scoring)

AF3 runs as a local binary or API. Only the PDB output format is needed
(`--output_format=pdb`).

```bash
# local (GPU server with AF3 installed):
python3 scripts/humanize/cli.py run --input seq.fasta \
  --af3-mode local --af3-binary /path/run_alphafold.py \
  --af3-workdir /path/models --af3-db /path/pdb_databases

# API (e.g. internal AF3 service):
python3 scripts/humanize/cli.py run --input seq.fasta \
  --af3-mode api --af3-api https://af3.internal/v1/predict \
  # set env AF3_TOKEN=...
```

What the pipeline does with AF3 models (per chain and per variant):

1. predict Fv (VH+VL or VHH) and, if `--antigen <seq>` is given, the complex;
2. compute per-position hints: buriedness proxy (heavy-atom contacts),
   CDR contacts (<4.5 A), antigen contacts;
3. feed the hints back into back-mutation scoring (structural score
   refinement);
4. per-variant: CDR-loop CA RMSD vs the donor model (same numbering),
   VH/VL interface contact count retention.

Recommended seeds: 2-3 per job; report the best-pLDDT model.

## 4. ProteinMPNN (alternative framework design)

ProteinMPNN re-designs the framework on the fixed backbone with CDRs,
interface, vernier and VHH hallmark locked:

```bash
python3 scripts/humanize/cli.py run --input seq.fasta \
  --mpnn-mode local --mpnn-script /path/protein_mpnn.py
```

The adapter generates the `--fixed_residues` list automatically, parses the
`*fa` outputs, filters designs by human germline identity of the designed
framework, and flags developability risks. Top designs are emitted as
additional variants. Without ProteinMPNN, a top-germline consensus framework
is produced as a sanity check instead.

## 5. Optional cross-checks

- BioPhi (https://github.com/oxpig/biophi or the `biophi` package):
  independent ML-based humanization scoring; use as a second opinion on
  germline choice and per-position human-likeness.
- IgBLAST `-humanize` mode: built-in "closest human germline" search for
  cross-validation of the germline selection.

## Quick start on a fresh server

```bash
bash scripts/install.sh --full
conda activate humanize
python3 tests/test_pipeline.py
python3 scripts/humanize/cli.py run --input data/examples/mouse_4d5_fab.fasta --outdir demo
```
