# IgFold Structure Prediction Tool

## Overview

This tool uses [IgFold](https://github.com/Graylab/IgFold) to predict antibody structures from FASTA sequences. It's designed to work with the humanization pipeline for RMSD evaluation.

## Installation

### 1. Create conda environment

```bash
conda create -n igfold python=3.11 -y
```

### 2. Install IgFold

```bash
conda run -n igfold pip install igfold
```

### 3. Download AntiBERTy weights (first run only)

The first time you run IgFold, it will download ~100MB of model weights from Hugging Face. This requires internet access.

If you need to run offline, you can manually download the weights:

```bash
# On a machine with internet access:
pip install huggingface_hub
python -c "from huggingface_hub import snapshot_download; snapshot_download('jeffruffolo/AntiBERTy', local_dir='./antiberty_weights')"

# Copy the weights to the target machine and set:
export ANTIbertY_WEIGHTS_DIR=/path/to/antiberty_weights
```

## Usage

### Single antibody prediction

```bash
conda run -n igfold python tools/igfold_predict.py \
  --input af3/amg110_variants/amg110_V1.fasta \
  --output af3/amg110_variants/
```

### Batch prediction (all variants)

```bash
conda run -n igfold python tools/igfold_predict.py \
  --input-dir af3/amg110_variants/ \
  --output af3/amg110_variants/
```

### With refinement and renumbering

```bash
conda run -n igfold python tools/igfold_predict.py \
  --input-dir af3/amg110_variants/ \
  --output af3/amg110_variants/ \
  --refine openmm \
  --renumber
```

### Using GPU

```bash
conda run -n igfold python tools/igfold_predict.py \
  --input-dir af3/amg110_variants/ \
  --output af3/amg110_variants/ \
  --device cuda:0
```

## Options

| Option | Description | Default |
|--------|-------------|---------|
| `--input` | Input FASTA file (single prediction) | - |
| `--input-dir` | Input directory with FASTA files (batch) | - |
| `--output` | Output directory for PDB files | Required |
| `--refine` | Refinement method: `none`, `openmm`, `pyrosetta` | `none` |
| `--renumber` | Renumber with Chothia scheme | False |
| `--num-models` | Number of ensemble models (1-4) | 4 |
| `--device` | Device: `cpu`, `cuda:0`, `cuda:1`, `mps` | `cpu` |
| `--python-api` | Use Python API (experimental) | False |

## Input FASTA Format

The script accepts FASTA files with VH and VL sequences:

```
>H
EVQLLQSGAEVKKPGESLKISCKGSGYSFTNYWLGWVKQMPGKGLEWIGDIFPGSGNIHYNEKFKGQATLSADKSISTAYLQWSSLKASDTAMYYCARLRNWDEPMDYWGQGTTVTVSS
>L
DLVMTQSPDSLAVSLGERATINCKSSQSLLNSGNQKNYLTWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQNDYSYPLTFGQGTKLEIK
```

Or with descriptive headers:

```
>V1_H Heavy chain
EVQLLQSGAEVKKPGESLKISCKGSGYSFTNYWLGWVKQMPGKGLEWIGDIFPGSGNIHYNEKFKGQATLSADKSISTAYLQWSSLKASDTAMYYCARLRNWDEPMDYWGQGTTVTVSS
>V1_L Light chain
DLVMTQSPDSLAVSLGERATINCKSSQSLLNSGNQKNYLTWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQNDYSYPLTFGQGTKLEIK
```

## Output

PDB files are saved to the specified output directory with the same name as the input FASTA file.

Example:
```
af3/amg110_variants/
├── amg110_V0.fasta  (input)
├── amg110_V0.pdb    (output)
├── amg110_V1.fasta
├── amg110_V1.pdb
└── ...
```

## RMSD Evaluation

After predicting structures, run RMSD evaluation with the humanization pipeline:

```bash
python3 scripts/humanize/cli.py rmsd \
  --input inputs/amg110_correct.fasta \
  --donor "af3/amg110_structure/Structure Prediction (Boltz-2)/rank_1.pdb" \
  --variants af3/amg110_variants/amg110_V0.pdb \
            af3/amg110_variants/amg110_V1.pdb \
            af3/amg110_variants/amg110_V2.pdb \
  --chain H
```

## Performance

| Metric | Value |
|--------|-------|
| Speed (CPU) | ~15-25 seconds per structure |
| Speed (GPU) | ~5-10 seconds per structure |
| CDR-H3 RMSD | ~4.25 Å (vs experimental) |
| Framework RMSD | 0.57-0.80 Å |

## Troubleshooting

### Network unreachable error

If you see "Network is unreachable" error, it means IgFold cannot download the AntiBERTy weights. Solutions:

1. Run on a machine with internet access
2. Download weights manually (see Installation section)
3. Set `ANTIbertY_WEIGHTS_DIR` environment variable to local weights

### CUDA out of memory

If you run out of GPU memory, use CPU mode:

```bash
--device cpu
```

### Slow prediction

For faster prediction:
1. Use `--num-models 1` (faster, slightly less accurate)
2. Use GPU if available (`--device cuda:0`)
3. Use `--refine none` (skip refinement)

## References

- IgFold paper: [Nature Communications (2023)](https://www.nature.com/articles/s41467-023-38063-x)
- GitHub: https://github.com/Graylab/IgFold
