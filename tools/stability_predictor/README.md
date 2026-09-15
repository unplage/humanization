# Protein Stability Predictor (ΔG/ΔΔG)

Antibody stability prediction tool for humanization pipeline. Predicts absolute folding stability (ΔG) and mutational effects (ΔΔG) using free, open-source AI models.

## Features

- **Absolute Stability (ΔG)**: Predicts protein folding free energy
- **Mutational Effect (ΔΔG)**: Evaluates impact of back-mutations on stability
- **Batch Analysis**: Analyze multiple variants (V0-V3) in parallel
- **Pipeline Integration**: Works with humanization pipeline output
- **Multiple Models**: ESM3ΔG (recommended) and SaProtΔG

## Installation

```bash
# Install ESM3ΔG (recommended, MIT license)
pip install git+https://github.com/yehlincho/absolute-stability-predictor.git

# Or install SaProtΔG (alternative)
pip install git+https://github.com/yehlincho/absolute-stability-predictor.git
```

**Requirements:**
- Python 3.8+
- PyTorch (installed automatically)
- ~2GB disk space for model weights

## Quick Start

### Single Structure Analysis

```bash
# Analyze absolute stability ΔG
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb

# Specify chain (default: A)
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb antibody.pdb \
    --chain H
```

### Mutational Effect (ΔΔG)

```bash
# Compare wildtype vs mutant
python3 tools/stability_predictor/stability_analyzer.py \
    --wt wildtype.pdb \
    --mutant mutant.pdb
```

### Batch Analysis

```bash
# Analyze all PDB files in a directory
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir variants/ \
    --output results/
```

### Integration with Humanization Pipeline

```bash
# Analyze all variants from Step 3 output
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ \
    --output outputs/stability/
```

## Usage with Pipeline

### Step 1: Run Humanization Pipeline

```bash
# Generate variants with AF3 structure prediction
python3 scripts/humanize/cli.py run \
    --input inputs/antibody.fasta \
    --outdir outputs/step3 \
    --donor-structure donor.pdb \
    --af3-mode local
```

### Step 2: Analyze Stability

```bash
# Analyze all predicted structures
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ \
    --output outputs/stability/
```

### Step 3: Interpret Results

- **ΔG < 0**: Protein is predicted to be stable
- **ΔG > 0**: Protein may be unstable
- **ΔΔG < 0**: Mutation is stabilizing (good for humanization)
- **ΔΔG > 0**: Mutation is destabilizing (may affect binding)

## Output Formats

### Text Output (default)

```
================================================================================
  Stability Analysis: V2_variant.pdb
================================================================================

  Model: ESM3dG
  Chain: A

  Absolute Stability (ΔG): -12.34 kcal/mol
  Status: STABLE
```

### JSON Output

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb \
    --json
```

```json
{
  "pdb_file": "structure.pdb",
  "chain_id": "A",
  "dg": -12.34,
  "model": "ESM3dG",
  "error": null,
  "is_stable": true
}
```

### Save Results

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir variants/ \
    --output results/ \
    --json
```

Creates `results/stability_batch_results.json` with summary statistics.

## Model Comparison

| Model | Speed | Accuracy | License | Recommended |
|-------|-------|----------|---------|-------------|
| ESM3ΔG | Fast | High | MIT | Yes |
| SaProtΔG | Medium | High | MIT | Alternative |

## Integration Examples

### With PTM Exposure Analysis

```bash
# Analyze both stability and PTM risks
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ \
    --output outputs/stability/

python3 tools/ptm_exposure/analyze_ptm.py \
    --pdb-dir outputs/step3/af3/ \
    --output outputs/ptm/
```

### With TAP Developability

```bash
# Complete developability assessment
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb structure.pdb \
    --output outputs/stability/

python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb \
    --output outputs/tap/
```

## Troubleshooting

### "ModuleNotFoundError: No module named 'ESM3dG'"

```bash
pip install git+https://github.com/yehlincho/absolute-stability-predictor.git
```

### "CUDA out of memory"

ESM3ΔG uses GPU acceleration. If GPU memory is insufficient:

```bash
# Force CPU mode
export CUDA_VISIBLE_DEVICES=""
python3 tools/stability_predictor/stability_analyzer.py --pdb structure.pdb
```

### "Prediction failed"

1. Check PDB file format is valid
2. Ensure structure has valid backbone atoms
3. Try with `--chain A` to specify chain

## Scientific Background

### ΔG (Absolute Stability)

The folding free energy ΔG predicts whether a protein will fold spontaneously:

- **ΔG < 0**: Thermodynamically favorable (stable)
- **ΔG > 0**: Thermodynamically unfavorable (unstable)
- Typical antibody ΔG: -5 to -15 kcal/mol

### ΔΔG (Mutational Effect)

The change in stability upon mutation:

- **ΔΔG < 0**: Stabilizing mutation (good)
- **ΔΔG > 0**: Destabilizing mutation (may affect function)
- Threshold: |ΔΔG| > 1 kcal/mol is significant

## References

1. Cho et al. (2026) "Accurate protein stability prediction for small domains using mega-scale experiments" bioRxiv
2. ESM3: EvolutionaryScale's protein language model
3. SaProt: Westlake's structure-aware protein model

## License

MIT License - same as ESM3ΔG and SaProtΔG models.
