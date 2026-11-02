# Boltz-2 structure prediction

Predict antibody/protein structures with Boltz-2 and write PDB files with
per-residue pLDDT in the B-factor column. Designed for the humanization
pipeline's structure validation steps.

## Requirements

- conda environment `boltz` (Python >= 3.11, PyTorch >= 2.2):

  ```bash
  conda create -n boltz python=3.11
  conda run -n boltz pip install "torch==2.5.1+cu121" -f https://download.pytorch.org/whl/cu121
  conda run -n boltz pip install boltz==2.2.1
  ```

- Model weights in `~/.boltz/` (or `$BOLTZ_CACHE`):
  `boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols/` (or `mols.tar`).
- MSA server access for the default online mode; use `--offline` for
  single-sequence prediction.

Check the environment before running:

```bash
conda run -n boltz python -c "import boltz, torch; print(boltz.__version__, torch.cuda.is_available())"
```

## Usage

Run from the repository root with the `boltz` conda environment:

```bash
# Fab/Fv from FASTA (headers H and L)
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input antibody.fasta --output predictions/

# VHH/nanobody (single chain)
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input vhh.fasta --output predictions/ --single-chain

# Batch all FASTA/YAML files in a directory
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input-dir variants/ --output predictions/

# Offline mode (no MSA server)
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input fv.yaml --output predictions/ --offline

# Longer timeout for large complexes (default 600 s per structure)
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input fv.yaml --output predictions/ --timeout 1800
```

### Kernels

cuequivariance kernels are **disabled by default** because
`cuequivariance_torch` is an optional dependency. The default invocation
therefore works on a plain `boltz==2.2.1` install. If you have installed
`cuequivariance_torch`, enable the faster kernel path explicitly:

```bash
conda run -n boltz python tools/boltz/boltz_predict.py \
  --input fv.yaml --output predictions/ --kernels
```

`--no-kernels` is still accepted and is equivalent to the default.

## Input

- **FASTA**: chain IDs are inferred from headers. `>H` / `>L` / `>A` are used
  directly; headers such as `>V1_H` or `>...H` map to heavy chain, and
  `>...L` to light chain. A single sequence with `--single-chain` is treated
  as a VHH.
- **YAML**: native Boltz input. Protein entries must provide `id` and
  `sequence`; `msa: empty` is recommended for antibody single-sequence runs.

## Output

Boltz writes to:

```
<output>/boltz_results_<stem>/predictions/<stem>/<stem>_model_0.pdb
```

The wrapper prints the selected PDB path on success. The PDB contains pLDDT
in the B-factor column. Confidence metrics are written alongside it in
`confidence_<stem>_model_0.json` (including `confidence_score`,
`complex_plddt`, `ptm`, `iptm`).

## Notes

- The wrapper invokes the `boltz` CLI from `PATH`; run it inside the `boltz`
  environment (or via `conda run -n boltz`).
- On failure, the full Boltz stderr is printed to aid diagnosis.
- The wrapper does not run the downstream humanization pipeline. Use the
  generated PDB with `scripts/humanize/cli.py rmsd` or pass it as
  `--donor-structure`/structure input to the relevant pipeline stage.
