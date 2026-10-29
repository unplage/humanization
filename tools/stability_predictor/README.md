# Protein Stability Predictor (ΔΔG) — ThermoMPNN

Antibody stability prediction tool for humanization pipeline. Predicts mutational effects (ΔΔG) on protein stability using **ThermoMPNN** (Kuhlman Lab, MIT license). Designed for evaluating back-mutation candidates in antibody humanization.

## Features

- **Mutational Effect (ΔΔG)**: Predicts ΔΔG for single-point mutations using ThermoMPNN (site-saturation mutagenesis)
- **Back-mutation Focus**: Evaluates stability impact of framework back-mutations (T1/T2 tier positions)
- **Kabat/IMGT Region Annotation**: Each mutation annotated with Kabat and IMGT region (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4)
- **Batch Analysis**: Analyze multiple variant structures (V0-V3) in a directory
- **Pipeline Integration**: Works with humanization pipeline AF3 structure outputs
- **Multiple Output Formats**: Text, JSON, Markdown, Word

## Requirements

- **ThermoMPNN** repo cloned locally (see Installation)
- **Conda environment** with ThermoMPNN dependencies (Python 3.10+)
  - PyTorch, pytorch-lightning, omegaconf, biopython, pandas, numpy, torchmetrics
- **Model checkpoint** (`.ckpt` or `.pt`) in `ThermoMPNN/models/` or `ThermoMPNN/vanilla_model_weights/`

## Installation

### 1. Clone ThermoMPNN

```bash
git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git ~/ThermoMPNN
```

### 2. Create conda environment

```bash
# Using the provided environment.yaml (recommended)
cd ~/ThermoMPNN
mamba env create -f environment.yaml
# or: conda env create -f environment.yaml

# Activate
conda activate thermoMPNN
```

### 3. Verify installation

```bash
# Should find the checkpoint
ls ~/ThermoMPNN/models/*.ckpt ~/ThermoMPNN/vanilla_model_weights/*.pt 2>/dev/null | head -5

# Test Python environment
~/miniconda3/envs/thermoMPNN/bin/python -c "import torch; import pytorch_lightning; import omegaconf; print('OK')"
```

## Quick Start

### Single Structure Analysis (with Kabat positions)

```bash
# Analyze back-mutation positions for heavy chain
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb \
    --chain H \
    --positions 37,44,45,47,48,49,67,69,71,73,76,78,80,82,87,94 \
    --thermompnn-dir ~/ThermoMPNN \
    --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/
```

### Batch Analysis of Variant Directory

```bash
# Analyze all AF3-predicted variant structures
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ \
    --chain H \
    --thermompnn-dir ~/ThermoMPNN \
    --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/ \
    --json
```

### Using `conda run` directly

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb --chain H \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output outputs/stability/
```

### Via Humanization Pipeline CLI

```bash
# After running humanization with --af3-mode
python3 scripts/humanize/cli.py stability \
    --pdb-dir outputs/step3/af3/ \
    --chain H \
    --positions 37,44,45,47,48,49,67,69,71,73,76,78,80,82,87,94 \
    --thermompnn-dir ~/ThermoMPNN \
    --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/ \
    --json
```

## Position Specification

- **Kabat positions**: Provide as comma-separated list (e.g., `37,44,45` or `H37,H44,H45`)
- **Chain prefix**: Optional; defaults to `H` (heavy chain). Use `L` for light chain (e.g., `L24,L27`)
- **Typical back-mutation positions** (VH framework):
  - FR1: H37 (VHH hallmark), H44, H45, H47
  - FR2: H48, H49, H67, H69, H71, H73, H76, H78
  - FR3: H80, H82, H87, H94

## Output Formats

### Text Output (default)

```
================================================================================
  ThermoMPNN Stability Analysis: V2_variant.pdb (chain H)
================================================================================

  Model: ThermoMPNN
  Position: H37 (PDB pos 37)
  Mutation: V -> F
  Kabat Region: FR2
  IMGT Region: FR2

  Mutational Effect (ΔΔG): -0.85 kcal/mol
  Effect: STABILIZING

================================================================================
  SUMMARY
================================================================================
  Total mutations analyzed:  16
  Stabilizing (ΔΔG < 0):     12
  Destabilizing (ΔΔG > 0):   4
  Significant (|ΔΔG| ≥ 1.0): 3
  Mean ΔΔG:                  -0.42 kcal/mol
  Min ΔΔG:                   -1.23 kcal/mol
  Max ΔΔG:                   0.87 kcal/mol

  Thresholds: |ΔΔG| ≥ 1.0 kcal/mol = significant; ΔΔG < 0 = stabilizing
  Model: ThermoMPNN (Kuhlman Lab, MIT license)
```

### JSON Output

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb --chain H --positions 37,44,45 \
    --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --json
```

```json
{
  "summary": {
    "n_mutations": 3,
    "n_successful": 3,
    "n_failed": 0,
    "ddg_mean": -0.42,
    "ddg_min": -1.23,
    "ddg_max": 0.15,
    "n_stabilizing": 2,
    "n_destabilizing": 1,
    "n_significant": 1
  },
  "results": [
    {
      "pdb_file": "structure.pdb",
      "chain_id": "H",
      "position_pdb": 37,
      "position_kabat": "H37",
      "wildtype": "V",
      "mutation": "F",
      "ddg": -1.23,
      "region_kabat": "FR2",
      "region_imgt": "FR2",
      "model": "ThermoMPNN",
      "error": null,
      "is_stabilizing": true,
      "is_significant": true
    }
  ]
}
```

### Save Reports to Files

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir variants/ --chain H \
    --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/ \
    --all-formats
```

Creates per-structure files:
- `outputs/stability/<stem>_stability.json`
- `outputs/stability/<stem>_stability.md`
- `outputs/stability/<stem>_stability.docx`
- `outputs/stability/stability_batch_results.json` (combined)

## Interpreting Results

| ΔΔG Value | Interpretation | Action |
|-----------|----------------|--------|
| ΔΔG < -1.0 | Strongly stabilizing | Favorable back-mutation |
| -1.0 ≤ ΔΔG < 0 | Mildly stabilizing | Acceptable back-mutation |
| 0 ≤ ΔΔG ≤ 1.0 | Mildly destabilizing | Caution; monitor |
| ΔΔG > 1.0 | Significantly destabilizing | Avoid or compensate |

**Note**: ThermoMPNN predicts *relative* mutational effects (ΔΔG), not absolute folding stability (ΔG). For absolute stability, consider complementary tools.

## Integration Examples

### With PTM Exposure Analysis

```bash
# Analyze both stability and PTM risks
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ --chain H \
    --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/

python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb-dir outputs/step3/af3/ \
    --output outputs/ptm/
```

### With TAP Developability

```bash
# Complete developability assessment
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb --chain H \
    --thermompnn-dir ~/ThermoMPNN --python ~/miniconda3/envs/thermoMPNN/bin/python \
    --output outputs/stability/

python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --output outputs/tap/
```

## Troubleshooting

### "ThermoMPNN not properly installed"

- Verify `--thermompnn-dir` points to the cloned repo
- Check `models/*.ckpt` or `vanilla_model_weights/*.pt` exists
- Verify `local.yaml` exists in repo root

### "Python env check failed"

- Ensure `--python` points to correct conda environment
- For `conda run`: use full command `conda run -n thermoMPNN python`
- Verify environment has: `torch`, `pytorch_lightning`, `omegaconf`, `biopython`

### "ThermoMPNN failed" / "No output CSV"

- Input PDB must have valid backbone atoms for the specified chain
- Chain ID must match PDB (case-sensitive: `H`, `L`, `A`, `B`, etc.)
- Check ThermoMPNN stderr with `--verbose`

### Long runtime

- ThermoMPNN runs site-saturation mutagenesis (L × 20 forward passes)
- For ~120 residue chain: ~2400 mutations, ~3-10 min on CPU
- Use `--positions` to limit to specific Kabat positions (much faster)

## Scientific Background

### ΔΔG (Mutational Effect)

The change in folding free energy upon mutation:

- **ΔΔG < 0**: Stabilizing mutation (good for humanization)
- **ΔΔG > 0**: Destabilizing mutation (may affect folding/expression)
- **Threshold**: |ΔΔG| ≥ 1.0 kcal/mol considered significant

ThermoMPNN uses a graph neural network (GNN) with transfer learning from ProteinMPNN, trained on mega-scale stability datasets (MegaScale, FireProtDB). Benchmarks: Spearman ρ ≈ 0.72 on MegaScale test set.

### References

1. Dieckhaus et al. (2024) "Transfer learning to leverage larger datasets for improved prediction of protein stability changes" PNAS 121(6): e2314853121
2. ThermoMPNN: https://github.com/Kuhlman-Lab/ThermoMPNN
3. MegaScale dataset: Tsuboyama et al. (2023)

## License

ThermoMPNN code and model weights: **MIT License** (Kuhlman Lab)

This wrapper: MIT License (part of humanization pipeline)

**Note**: Unlike ESM3ΔG (Cambrian Non-Commercial License, prohibits drug discovery), ThermoMPNN is MIT-licensed and suitable for therapeutic antibody development.