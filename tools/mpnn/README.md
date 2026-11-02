# ProteinMPNN Framework Design Tool

Adapter that lets the humanization pipeline use an upstream
[ProteinMPNN](https://github.com/dauparas/ProteinMPNN) checkout for
structure-guided framework design (Step 4 developability optimization).

The pipeline's MPNN integration still speaks the **pre-2023** ProteinMPNN CLI
(`--chains_to_design`, `--fixed_residues`). Upstream renamed the entry point to
`protein_mpnn_run.py` and replaced those flags with `--pdb_path_chains` +
`--fixed_positions_jsonl`. `proteinmpnn_adapter.py` translates between the two,
so the pipeline keeps working with a current checkout.

## Files

| File | Purpose |
|------|---------|
| `proteinmpnn_adapter.py` | **Use this.** Legacy CLI -> modern `protein_mpnn_run.py` translator. |
| `protein_mpnn.py` | Deprecated alias; forwards its args to the adapter. |

## Requirements

ProteinMPNN runs fine on **CPU-only** PyTorch. A dedicated env keeps the heavy
`torch` stack out of the portable pipeline env:

```bash
conda create -n proteinMPNN python=3.11 -y
conda activate proteinMPNN
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install pytorch-lightning omegaconf pandas numpy biopython
```

Then obtain the model repo and weights (weights are already bundled in the
release/main tarball under `vanilla_model_weights/`):

```bash
git clone https://github.com/dauparas/ProteinMPNN
# or download ProteinMPNN-main.zip and unzip
```

## Configuration

The adapter finds the repo via, in order:

1. `--proteinmpnn-dir <path>`
2. `PROTEINMPNN_DIR` environment variable
3. `~/ProteinMPNN`, `~/ProteinMPNN-main`, `~/ProteinMPNN/vanilla_proteinmpnn`

`python scripts/humanize/cli.py setup-check` reports whether it is detected.

## Usage with the pipeline

```bash
export PROTEINMPNN_DIR=/path/to/ProteinMPNN

python scripts/humanize/cli.py run --input ab.fasta --outdir outputs \
  --mpnn-mode local \
  --mpnn-script tools/mpnn/proteinmpnn_adapter.py \
  --mpnn-python /path/to/miniconda3/envs/proteinMPNN/bin/python \
  --donor-structure donor.pdb
```

- `--mpnn-python` selects the interpreter that has `torch`/ProteinMPNN
  (defaults to the current interpreter).
- The PDB chain to design is auto-detected from the structure. Override with
  `--mpnn-chains A` (or `--mpnn-chains "A B"`) if needed.

MPNN only runs when the V2 variant carries a high-risk developability motif and
a structure is available (`--donor-structure` or AF3). Without a structure the
pipeline falls back to a germline-consensus design.

## Standalone invocation

The adapter accepts the same legacy flags the pipeline passes, so you can test
it directly:

```bash
python tools/mpnn/proteinmpnn_adapter.py \
  --pdb_path structure.pdb --out_folder /tmp/mpnn_out \
  --num_seq_per_target 8 --batch_size 4 \
  --chains_to_design A --fixed_residues 1,2,3,5,10 \
  --seed 0 --suppress_print 1 \
  --proteinmpnn-dir /path/to/ProteinMPNN
```

## Behaviour notes

- **Output location**: upstream writes `<out>/seqs/<name>.fa`; the adapter also
  copies the FASTA to `<out>/<name>.fa`, which is where the pipeline globs.
- **Native record is dropped**: upstream writes the WT sequence first. The
  adapter strips it so the pipeline does not count it as a design
  (`--keep-native` disables this).
- **Seeds**: upstream treats `--seed 0` as "pick a random seed"; the adapter
  maps `0 -> 37` so runs are reproducible.
- **Out-of-range fixed positions** are clipped per chain with a warning
  (ProteinMPNN indexes them into a chain-length array and would otherwise crash
  with an opaque `IndexError`).
- **mmCIF input**: `.cif`/`.mmcif` is converted to a temporary PDB via
  biopython before design.
- **Missing residues**: if a chain has gaps in its numbering, ProteinMPNN
  inserts GAP slots, which can shift fixed positions that lie after a gap. The
  adapter prints a warning; `--strict-gaps` turns it into a hard error.
- Fixed positions are **1-based positions in the parsed chain sequence**, the
  same space the pipeline's `NumberedChain.index` uses.
