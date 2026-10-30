# Protein Stability Predictor (ΔΔG) — ThermoMPNN

Antibody stability prediction tool for the humanization pipeline. It calls
**ThermoMPNN** (Kuhlman Lab, MIT license) to predict folding-stability effects
of point mutations, then maps ThermoMPNN's internal sequence positions to
**true Kabat/IMGT positions**. It is designed for deciding which framework
back-mutations (T1/T2/T3) to keep.

Run all commands from the **repository root** so that the bundled numbering
modules and relative paths resolve correctly.

---

## What it does

1. Runs ThermoMPNN on a PDB/CIF structure.
2. Maps ThermoMPNN's **0-based chain-sequence index** to the real Kabat/IMGT
   position and region using the pipeline's numbering engine.
3. In back-mutation mode, evaluates exactly the human↔donor substitutions from
   the humanization pipeline and reports a direction-aware
   `ddG_backmutation`, verdict, and recommendation.

Two modes:

| Mode | When to use | How |
|---|---|---|
| `backmutation` (recommended) | Decide which pipeline T1/T2/T3 back-mutations to keep | `--pipeline-json` or `--backmutations-csv` |
| `ssm` | Scan all substitutions at selected positions, or full SSM | `--positions H67,H69` or no position filter |

---

## Requirements

- **ThermoMPNN repository** cloned locally.
- **Python environment** with ThermoMPNN dependencies (Python 3.10+):
  PyTorch, pytorch-lightning, omegaconf, biopython, pandas, numpy, torchmetrics.
- **Stability checkpoint** in `ThermoMPNN/models/`
  (`thermoMPNN_default.pt`, `*.pt`, or `*.ckpt`).

> `vanilla_model_weights/*.pt` are **ProteinMPNN** weights, **not** a stability
> model. The wrapper refuses to use them.

---

## Installation

```bash
git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git ~/ThermoMPNN

cd ~/ThermoMPNN
mamba env create -f environment.yaml   # or: conda env create -f environment.yaml
conda activate thermoMPNN

# Verify checkpoint and Python dependencies
ls ~/ThermoMPNN/models/*.ckpt ~/ThermoMPNN/models/*.pt
python -c "import torch, pytorch_lightning, omegaconf; print('OK')"
```

If you already have an environment with these dependencies, you can use it
directly and skip creating a new one.

---

## Before you start

### 1. Check the installation

```bash
python3 scripts/humanize/cli.py setup-check
```

Expected line:

```text
ThermoMPNN repo     : found at /path/to/ThermoMPNN (ckpt: thermoMPNN_default.pt)
```

### 2. Prepare the input structure

A **V0 structure** is strongly recommended:

- V0 = human framework + donor CDRs, with human residues at back-mutation
  positions.
- ThermoMPNN then measures **human → donor** directly.
- `ddG_backmutation = ddG`; no sign inversion is needed.

If you do not have V0, V2/V3 can be used, but the sign is inferred by
inverting the donor→human prediction. Prefer V0 for decisions.

Typical inputs:

```text
outputs/.../igfold_variants/<name>_V0.pdb
outputs/.../step3/humanization_result.json
# or
outputs/.../step3/backmutations_<name>.H.csv
```

### 3. Identify the chain IDs

```bash
# Inspect chain IDs if unsure
grep '^ATOM' V0.pdb | cut -c22 | sort -u
```

Use `--chain all` for Fab (H and L) or when chain IDs are `A`/`B`. The tool can
map `H`/`L` to `A`/`B` by checking heavy/light numbering.

### 4. Point `--python` at the ThermoMPNN environment

```bash
# Direct interpreter
--python ~/miniconda3/envs/thermoMPNN/bin/python

# Or conda run (quote the whole command)
--python "conda run -n thermoMPNN python"

# Any existing environment with the dependencies also works
--python "conda run -n my_torchem_env python"
```

`~/...` is expanded automatically, and `conda run ...` is shell-split correctly.

---

## Step-by-step: recommended workflow

### Step 1 — Evaluate pipeline back-mutations on V0

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb outputs/step3/igfold_variants/<name>_V0.pdb \
    --chain all \
    --pipeline-json outputs/step3/humanization_result.json \
    --tiers all \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output outputs/stability_<name>_V0 \
    --all-formats
```

What happens:

- Loads all back-mutation candidates (T1/T2/T3) from the pipeline JSON.
- Maps each candidate to the structure residue and decides the direction:
  - structure residue == human → evaluate human→donor (`ddG_bm = ddG`)
  - structure residue == donor → evaluate donor→human (`ddG_bm = -ddG`)
- Runs ThermoMPNN only for those exact substitutions (fast targeted mode).
- Prints a table and writes JSON/Markdown/Word reports.

### Step 2 — Read the terminal table

```text
Kabat    Tier  H>D    WT  Mut        dG    dG_bm Region    Verdict                        Recommendation
H48      T1    M>I    M   I       -1.02    -1.02 FR2       strongly supports back-mutation keep
H67      T1    V>A    V   A        0.16     0.16 FR3       neutral                        keep; structural evidence decides
```

Column meanings:

| Column | Meaning |
|---|---|
| `Kabat` | True Kabat position |
| `Tier` | Pipeline tier (T1/T2/T3) |
| `H>D` | Human → donor substitution |
| `WT` | Residue in the input structure |
| `Mut` | Residue evaluated in this run |
| `dG` | Raw ThermoMPNN prediction |
| `dG_bm` | Direction-aware effect of **applying the back-mutation** |
| `Region` | Kabat region (FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4) |
| `Verdict` | Automatic interpretation |
| `Recommendation` | Suggested action |

### Step 3 — Interpret `ddG_backmutation`

| `ddG_backmutation` (kcal/mol) | Verdict | Action |
|---|---|---|
| `<= -1.0` | Strongly supports | Favor keeping the back-mutation |
| `-1.0` to `-0.2` | Supports | Favor keeping |
| `-0.2` to `0.2` | Neutral | Keep; let structural/contact evidence decide |
| `0.2` to `1.0` | Cautions | Consider dropping |
| `>= 1.0` | Argues against | Drop unless experimentally required |

### Step 4 — Open the reports

```bash
# Markdown is the easiest to read
cat outputs/stability_<name>_V0/<name>_V0_H_stability.md

# Machine-readable results
cat outputs/stability_<name>_V0/<name>_V0_H_stability.json
```

### Step 5 — Optional: audit the exact mutations evaluated

```bash
cat outputs/stability_<name>_V0/thermompnn_raw/<name>_V0/H/targeted_mutations.json
cat outputs/stability_<name>_V0/thermompnn_raw/<name>_V0/H/ThermoMPNN_inference_<name>_V0.csv
```

This is the best way to verify that a Kabat label maps to the expected
`seq_pos` and amino acid.

---

## Run modes

### A. Back-mutations from the pipeline JSON (recommended)

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb <V0.pdb> --chain all \
    --pipeline-json <step3>/humanization_result.json \
    --tiers all \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output <outdir> --all-formats
```

- `--tiers all` is the default and evaluates every tier (T1/T2/T3). Restrict with
  `--tiers T1,T2` when only the must-keep and recommended positions matter.
- `--pipeline-chain H` (or `L`, or a full chain name) restricts a multi-chain
  JSON to one chain. Matching is by chain `name` or `chain_type`.
- Do **not** combine `--pipeline-json` with `--positions`; the tool will reject
  this because the two modes answer different questions.

### B. Back-mutations from a CSV

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb <V0.pdb> --chain H \
    --backmutations-csv <step3>/backmutations_<name>.H.csv \
    --tiers all \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output <outdir> --json
```

Use either `--pipeline-json` or `--backmutations-csv`, not both.

### C. SSM at true Kabat positions

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb <structure.pdb> --chain H \
    --positions H67,H69,H44 \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output <outdir> --json
```

- `--positions` takes **true Kabat positions**, not ThermoMPNN 0-based indices.
- Insertion codes are supported: `H35A`, `H52A`.
- Unprefixed tokens get the `--chain-type` prefix: `37,44` → `H37,H44` for a
  heavy chain.
- Each position evaluates 19 non-self amino acids.
- Without `--positions`, full site-saturation mutagenesis runs at every residue
  (slow).

### D. Batch analysis

```bash
python3 tools/stability_predictor/stability_analyzer.py \
    --pdb-dir outputs/step3/af3/ --chain all \
    --pipeline-json outputs/step3/humanization_result.json \
    --tiers all \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output outputs/stability_batch --json
```

- Recursively finds `.pdb`, `.cif`, and `.mmcif`.
- Each structure/chain gets its own report and raw CSV directory.
- Writes `combined_stability_results.json` when JSON output is enabled.

### E. Via the humanize CLI

The CLI accepts the same options:

```bash
python3 scripts/humanize/cli.py stability \
    --pdb <V0.pdb> --chain all \
    --pipeline-json <step3>/humanization_result.json \
    --tiers all \
    --thermompnn-dir ~/ThermoMPNN \
    --python "conda run -n thermoMPNN python" \
    --output <outdir> --all-formats
```

---

## Parameter reference

| Parameter | Description |
|---|---|
| `--pdb FILE` | Single PDB/CIF to analyze |
| `--pdb-dir DIR` | Directory for batch analysis |
| `--chain H` | PDB chain ID (default `H`); accepts comma lists or `all` |
| `--chain-type auto/H/L` | Force heavy/light numbering when chain IDs are `A`/`B` |
| `--positions H67,H35A` | SSM at true Kabat positions |
| `--pipeline-json FILE` | Targeted back-mutation mode from `humanization_result.json` |
| `--pipeline-chain H` | Restrict the JSON to one chain (`name` or `chain_type`) |
| `--backmutations-csv FILE` | Targeted back-mutation mode from a CSV |
| `--tiers all` | Tiers to evaluate (default `all` = T1+T2+T3; use `T1,T2` to restrict) |
| `--thermompnn-dir DIR` | Path to the ThermoMPNN repository |
| `--python ...` | ThermoMPNN Python interpreter or `conda run ...` command |
| `--timeout 1800` | Inference timeout in seconds (default 1800) |
| `--output DIR` | Output directory (default `outputs/stability`) |
| `--json` | Write JSON report |
| `--markdown` | Write Markdown report |
| `--word` | Write Word report (requires `python-docx`) |
| `--all-formats` | Write JSON + Markdown + Word |
| `--verbose` | Print the exact ThermoMPNN command and stderr |

Mutually exclusive combinations:

- `--pdb` vs `--pdb-dir`
- `--pipeline-json` vs `--backmutations-csv`
- `--positions` vs `--pipeline-json`/`--backmutations-csv`

---

## Position and chain handling

### Why ThermoMPNN positions need mapping

ThermoMPNN's `analysis/custom_inference.py` reports `position` as a **0-based
index into its parsed, gap-aware chain sequence** — not PDB numbering and not
Kabat numbering. The wrapper:

1. Rebuilds the same gapped sequence ThermoMPNN builds (missing PDB numbers
   become `-`; MSE is treated as MET).
2. Maps each ThermoMPNN index to the ungapped residue.
3. Numbers the ungapped sequence with the pipeline's Kabat or IMGT engine.
4. Attaches `position_kabat`, `position_pdb`, `insertion_code`, and regions.

Consequences:

- `--positions H67` always means Kabat H67.
- Insertion codes such as `H35A` work.
- The same Kabat position can have different `seq_pos` values in V0 and V2 if
  the structures have different indel patterns; the mapping is per structure.
- If a residue cannot be mapped, it is reported with a warning and a synthetic
  label rather than silently mislabeled.

### Chain IDs and chain types

- `--chain` is the **PDB chain ID**, not the chain type.
- `--chain all` analyzes every chain.
- `--chain H` will also resolve to a chain named `A`/`B` whose sequence numbers
  as a heavy chain, provided `--chain-type` is `H` or `auto` (and similarly for
  `L`).
- `--chain-type auto` detects heavy/light by numbering and FR4 score.

---

## Reading the results

### Back-mutation mode

Each evaluated candidate reports:

- `position_kabat`: true Kabat label
- `position_seq`: ThermoMPNN 0-based sequence index
- `position_pdb`: PDB author residue number
- `wildtype`: residue in the input structure
- `mutation`: residue evaluated in this run
- `donor_aa` / `human_aa`: the two states from the pipeline
- `tier`: T1/T2/T3
- `direction`: `human_to_donor` or `donor_to_human`
- `ddG`: raw ThermoMPNN prediction
- `ddG_backmutation`: effect of applying the back-mutation
- `verdict` / `recommendation`: automatic interpretation
- `error` / `note`: reason when the candidate could not be evaluated

### V0 vs V2/V3 direction and sign

| Structure | WT state | Evaluated change | `ddG_backmutation` |
|---|---|---|---|
| V0 | human | human → donor | `ddG` (direct) |
| V1/V2/V3 | donor | donor → human | `-ddG` (approximation) |

Use a **V0** structure for decisions. On V2/V3 the sign inversion assumes the
mutation effect is approximately symmetric; it is a heuristic, not an exact
thermodynamic cycle. The same position can also give different magnitudes in
V0 and V2 because the structural context differs.

### Skipped candidates

Skipped candidates remain in the output with `ddG_backmutation: null` and a
reason. Common reasons:

- The position maps to a CDR (CDR back-mutations are not allowed).
- The structure residue matches neither the human nor the donor amino acid.
- The candidate is an FR indel (`human_aa` or `donor_aa` is `-`).
- The position is missing from the structure (numbering/structure mismatch).
- The target is not a standard amino acid.
- ThermoMPNN did not return a prediction for that mutation.

Always check the "Not evaluated" table before concluding that a candidate is
safe.

---

## Output files

Text output is always printed to stdout. With `--output DIR`, reports are
written as:

```text
DIR/
├── <stem>_stability.json
├── <stem>_stability.md
├── <stem>_stability.docx
├── combined_stability_results.json          # multi-chain/batch + JSON
├── targeted_mutations.json                  # single-chain raw audit
├── ThermoMPNN_inference_<stem>.csv           # single-chain raw audit
└── thermompnn_raw/
    └── <pdb_stem>/
        └── <chain>/
            ├── targeted_mutations.json
            └── ThermoMPNN_inference_<pdb_stem>.csv
```

- With multiple chains, the report stem becomes `<pdb_stem>_<chain>`.
- Multi-chain/batch raw files are kept in separate per-chain directories so
  they never overwrite each other.
- `targeted_mutations.json` lists the exact `[WT, seq_pos, Mut]` entries sent to
  ThermoMPNN.
- The raw CSV is the unmodified ThermoMPNN output.
- `combined_stability_results.json` is written when there are multiple reports
  and JSON/all-formats output is enabled.

### JSON example (back-mutation mode)

```json
{
  "pdb_file": "V0.pdb",
  "chain": "H",
  "chain_type": "H",
  "mode": "backmutation",
  "context": "human-framework context (V0-like)",
  "backmutations": [
    {
      "position_kabat": "H67",
      "position_seq": 67,
      "position_pdb": 68,
      "wildtype": "V",
      "mutation": "A",
      "donor_aa": "A",
      "human_aa": "V",
      "tier": "T1",
      "direction": "human_to_donor",
      "ddG": 0.16,
      "ddG_backmutation": 0.16,
      "region_kabat": "FR3",
      "region_imgt": "FR3",
      "verdict": "neutral",
      "recommendation": "keep; structural evidence decides",
      "error": null
    }
  ],
  "summary": {
    "candidates": 12,
    "evaluated": 12,
    "skipped": 0,
    "supports": 4,
    "neutral": 3,
    "cautions": 5,
    "opposes": 0,
    "mean_ddg_backmutation": -0.05
  },
  "warnings": []
}
```

### SSM mode

`--positions` or full SSM produces one row per mutation with the same
`position_kabat` / `position_seq` / `region_kabat` / `ddG` annotations. The text
and Markdown tables do not show back-mutation fields; in JSON those fields are
present as `null`/defaults (`ddG_backmutation`, `tier`, `donor_aa`, `human_aa`,
`direction`).

---

## Notes and limitations

1. **Single-chain calculation**: each chain is scored as a monomer. For Fab,
   the VH/VL interface is absent, so interface positions (e.g. H44/H45/H47,
   L35/L37/L38/L44/L45/L47) may be mis-scored. Analyze both chains
   (`--chain all`) and treat interface positions with caution; use structural
   contact/buriedness information from the pipeline as the primary evidence.
2. **Folding stability only**: ThermoMPNN does **not** predict affinity,
   expression, aggregation, immunogenicity, or VH/VL pairing.
3. **Relative effects only**: ThermoMPNN predicts ΔΔG, not absolute ΔG.
4. **CDR candidates are skipped**: the pipeline does not generate CDR/FR4
   back-mutations, and the tool additionally skips any candidate that maps to a
   CDR.
5. **Missing residues / indels are skipped**: check the "Not evaluated" table.
6. **Use V0 for decisions**: V2/V3 sign inversion is an approximation.
7. **Runtime**: targeted mode is fast (seconds to minutes); full SSM is roughly
   `L × 20` forward passes and can take much longer. Increase `--timeout` if
   needed; use a GPU-enabled environment when possible.
8. **Word output**: requires `python-docx`; if missing, the tool warns and skips
   the Word report without failing the analysis.

---

## Troubleshooting

### "ThermoMPNN not properly installed"

Check:

```bash
ls ~/ThermoMPNN/analysis/custom_inference.py
ls ~/ThermoMPNN/local.yaml
ls ~/ThermoMPNN/models/*.pt ~/ThermoMPNN/models/*.ckpt
```

- `--thermompnn-dir` must point at the cloned repo.
- A checkpoint must exist in `<repo>/models/`.
- `vanilla_model_weights/` is not accepted as a stability checkpoint.

### "Python env check failed"

- `--python` must point at the ThermoMPNN environment.
- `~/...` is expanded automatically.
- `conda run -n thermoMPNN python` must be quoted as one argument.
- Verify `torch`, `pytorch_lightning`, `omegaconf`, and `biopython`.

### "ThermoMPNN failed" / "No output CSV found"

- The PDB/CIF must contain valid backbone atoms for the requested chain.
- Chain IDs are case-sensitive; use `--chain all` if unsure.
- Re-run with `--verbose` to see the exact command and stderr.

### Wrong position labels

- Confirm the structure contains the requested chain.
- Confirm `--chain-type` is correct for `A`/`B` structures.
- Audit `targeted_mutations.json` and the raw CSV.
- Remember that `--positions` takes true Kabat positions, not 0-based indices.

### Timeout

```bash
--timeout 3600
```

Or use targeted mode / `--positions` to reduce the number of mutations.

---

## Validation

Run the stability test suite:

```bash
python3 tests/test_stability.py
```

The tests cover position mapping (including the 4D5 and AMG110 off-by-one
cases), gapped structures, targeted direction/sign, skipped candidates, report
formatting, invocation helpers, and an end-to-end path with a fake ThermoMPNN
environment.

---

## Scientific background

### ΔΔG (mutational effect)

- **ΔΔG < 0**: stabilizing (favors the mutation)
- **ΔΔG > 0**: destabilizing
- `|ΔΔG| >= 1.0 kcal/mol` is conventionally treated as significant

ThermoMPNN uses a graph neural network with transfer learning from ProteinMPNN,
trained on MegaScale/FireProtDB. Reported benchmark: Spearman ρ ≈ 0.72 on
MegaScale.

### References

1. Dieckhaus et al. (2024) "Transfer learning to leverage larger datasets for
   improved prediction of protein stability changes" PNAS 121(6): e2314853121
2. ThermoMPNN: https://github.com/Kuhlman-Lab/ThermoMPNN
3. MegaScale: Tsuboyama et al. (2023)

---

## License

ThermoMPNN code and model weights: **MIT License** (Kuhlman Lab). This wrapper:
MIT License (part of the humanization pipeline).

Unlike ESM3ΔG (Cambrian Non-Commercial License), ThermoMPNN is MIT-licensed and
suitable for therapeutic antibody development.
