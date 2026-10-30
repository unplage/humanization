# Immunogenicity Analysis Tool

Post-humanization immunogenicity analysis for antibody variants. Predicts MHC-II T-cell epitope burden (the core ADA driver), comparing humanized variants to the donor.

## Two-Step Analysis Workflow

**Critical Finding**: Sequence-based predictions alone are insufficient. Buriedness changes significantly between donor and variant structures due to backmutations.

### Step 1: Sequence-Based Analysis (HLAIIPred/Heuristic)

- **Input**: Antibody sequence (FASTA)
- **Method**: Predicts MHC-II binding probability from sequence
- **Output**: Per-peptide binding scores (0-1), epitope maps
- **Limitation**: Does not consider 3D structure → may overestimate risk for buried positions

### Step 2: Structure-Based Confirmation (FreeSASA)

- **Input**: Variant PDB structure (AF3/IgFold/Boltz/Protenix)
- **Method**: Computes relative SASA (relSASA) for each residue
- **Output**: Peptide-level adjustment factor; the adjusted mean **rewrites the final score**
- **Key insight**: Buried residues (relSASA < 0.15) are less accessible to MHC-II antigen processing

### Risk Adjustment Logic

The factor is a core-weighted average over all residues of each 15-mer
peptide (core 1.0, flanks 0.5). Back-mutation tiers are never used as
structural weights.

| relSASA | Classification | Risk Factor | Rationale |
|---------|----------------|-------------|-----------|
| < 0.05 | Deeply buried | 0.05 | Minimal accessibility |
| 0.05-0.15 | Partially buried | 0.05-0.30 | Low accessibility |
| 0.15-0.25 | Uncertain zone | 0.30-0.70 | Linear interpolation |
| 0.25-0.40 | Partially exposed | 0.70-0.90 | Moderate accessibility |
| > 0.40 | Fully exposed | 1.00 | Full risk |

Missing structure data contributes a neutral factor of 1.0, never a silent
0.5. The final score is `100 * adjusted_variant_mean / adjusted_reference_mean`
and is not capped at 100.

## Backend

| Backend | Method | License | Notes |
|---------|--------|---------|-------|
| **HLAIIPred** (default) | Transformer cross-attention, pan-allele | Apache-2.0 | Commercial use OK |
| **Heuristic** (fallback) | Germline rarity × exposure proxy | MIT | Always available, non-ML |

## Features

- **ML MHC-II prediction**: HLAIIPred Transformer model (Pfizer, Apache-2.0)
- **Pure-Python heuristic fallback**: always works, no external deps
- **Per-position epitope map**: which residues fall in strong/weak predicted T-cell epitopes
- **Donor comparison**: epitope load reduction vs mouse donor and V0 graft
- **Back-mutation impact**: Δepitope for each reverted position
- **Multi-format output**: text + JSON + Markdown + Word + CSV heatmap
- **Structural risk adjustment**: Two-step workflow with variant-specific buriedness analysis

## Installation

### Quick Install (Recommended)

```bash
# 1. Create conda environment (Python 3.11)
conda create -n hlapred python=3.11 -y
conda activate hlapred

# 2. Install PyTorch CPU
pip install torch --index-url https://download.pytorch.org/whl/cpu

# 3. Install HLAIIPred dependencies
pip install scipy numpy pandas tqdm biopython pyyaml

# 4. Clone HLAIIPred (Apache-2.0 license) INTO this tool directory and
#    install it from there. The clone must stay at this path: the install
#    below is editable, so the package body is read from this directory at
#    import time. Deleting or moving it breaks `import hlapred`.
HLAIIPRED_DIR=tools/immunogenicity/HLAIIPred   # from the repository root
git clone https://github.com/pfizer-opensource/HLAIIPred.git "$HLAIIPRED_DIR"
pip install -e "$HLAIIPRED_DIR"

# 5. Copy model files to project directory
mkdir -p tools/immunogenicity/models
cp "$HLAIIPRED_DIR"/models/epT_0.pt tools/immunogenicity/models/
cp "$HLAIIPRED_DIR"/models/epT_1.pt tools/immunogenicity/models/

# 6. Install structural risk assessment dependencies
pip install freesasa python-docx

# 7. Verify installation (run from the repository root)
python3 tools/immunogenicity/immunogenicity_analyzer.py --check
#   Expected:
#     hlaiipred       available
#     heuristic       available
#     freesasa        available
```

> **Do not delete `tools/immunogenicity/HLAIIPred/`.** An editable install
> (`pip install -e .`) keeps the package body in the clone directory and only
> writes a path pointer into `site-packages`. Removing the directory leaves the
> pointer dangling and `import hlapred` fails, which silently downgrades the
> tool to the heuristic backend. The clone is third-party code and is
> git-ignored; it is not meant to be committed or deleted.

If the editable install is ever broken, re-point it at the existing clone
instead of re-downloading anything:

```bash
pip install -e tools/immunogenicity/HLAIIPred --force-reinstall --no-deps
```

### Dependencies

| Dependency | Version | Purpose | Required? |
|------------|---------|---------|-----------|
| Python | ≥ 3.9 | Runtime | Yes |
| conda | ≥ 4.10 | Environment management | Recommended |
| PyTorch | ≥ 2.0 | HLAIIPred inference | Yes (if using HLAIIPred) |
| scipy | ≥ 1.7 | HLAIIPred dependency | Yes |
| numpy | ≥ 1.21 | Numerical computation | Yes |
| pandas | ≥ 1.3 | Data processing | Yes |
| tqdm | ≥ 4.60 | Progress bars | No |
| biopython | ≥ 1.79 | Sequence processing | Yes |
| pyyaml | ≥ 5.4 | Config file parsing | Yes |
| freesasa | ≥ 2.0 | SASA computation (structural risk) | No (degraded) |
| python-docx | ≥ 0.8 | Word report generation | No (degraded) |

### Verification

```bash
# Activate environment
conda activate hlapred

# Check backend availability
python3 tools/immunogenicity/immunogenicity_analyzer.py --check
# Output should show:
#   HLAIIPred: ✓ (Apache-2.0)
#   Heuristic: ✓ (MIT)
#   FreeSASA: ✓

# Run quick test
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --input data/examples/mouse_4d5_fab.fasta \
    --backend hlaiipred \
    --all-formats
```

### Windows Notes

```bash
# Windows users may need to install Visual C++ Redistributable
# Download: https://learn.microsoft.com/en-us/cpp/windows/latest-supported-vc-redist

# If freesasa installation fails, you can skip it (functionality degrades but still usable)
pip install python-docx
```

## Quick Start

### Check backend availability

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py --check
```

### Analyze pipeline output (auto-discovers variants.fasta)

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --outdir outputs/amg110_step3/step3 \
    --all-formats
```

### Analyze with structural risk adjustment

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --outdir outputs/amg110_step3/step3 \
    --structural-risk \
    --structure-dir af3/amg110_variants \
    --backmutation-dir outputs/amg110_step3/step3 \
    --all-formats
```

### Analyze specific FASTA

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --input variants.fasta \
    --all-formats
```

### Heuristic fallback (no external deps)

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --input variants.fasta \
    --backend heuristic \
    --all-formats
```

### With donor reference

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --input variants.fasta \
    --donor donor.fasta \
    --all-formats
```

## Key CLI options

| Option | Default | Meaning |
|--------|---------|---------|
| `--backend` | `auto` | `hlaiipred`, `heuristic`, or `auto` (falls back with a warning) |
| `--allele-panel` | `9common` | `9common` (official HLAIIPred example) or `extended` (24 alleles) |
| `--alleles` | panel | Explicit comma-separated alleles; max 14 per HLAIIPred call (batched automatically) |
| `--pep-len` | `15` | MHC-II peptide window length, 9-30 |
| `--strong-threshold` | `0.5` | Peptide window is strong at/above this score |
| `--weak-threshold` | `0.3` | Peptide window is weak in `[weak, strong)` |
| `--donor` | none | Always scored from its own sequence as the reference |
| `--reference-variant` | auto | Exact variant name to use as reference when no donor is given |
| `--structural-risk` | off | Enable structure-based rewriting of final scores |
| `--structure-dir` | required with structural | Directory containing per-variant PDB/CIF files |
| `--backmutation-dir` | none | Pipeline step3 dir with back-mutation CSVs and `variants_numbering.json` |
| `--numbering-json` | auto-discovered | Explicit per-variant Kabat/IMGT labels |

`extended` adds common DQ/DP alleles. HLAIIPred expects DP/DQ heterodimers as
`DPA1*01:03-DPB1*02:01` / `DQA1*01:02-DQB1*06:02`; the tool batches alleles
and pads missing slots with `0`, matching the official API.

## Position numbering (FR4 and insertion letters)

`structural_risk.py` matches back-mutation positions from the pipeline CSVs
by their exact label. The label space follows the pipeline run: **Kabat**
(`H5`, `H82A`, `H103`, `L100B`, ...) for a default run, or **IMGT**
(`H5`, `H111A`, `H118`, ...) for `--numbering imgt`. To do this exactly it
needs the per-residue labels of each variant. The pipeline writes them
automatically:

```
outputs/<run>/step3/variants_numbering.json   # Kabat run (flat)
{
  "H_V0": ["H1", "H2", ..., "H113"],
  "H_V2a": ["H1", "H2", ..., "H113"],
  ...
}

outputs/<run>/step3/variants_numbering.json   # IMGT run (wrapped)
{
  "numbering_scheme": "imgt",
  "variants": {"H_V0": ["H1", ..., "H128"], ...}
}
```

`structural_risk.py` auto-discovers this file in `--backmutation-dir`; you can
also pass it explicitly with `--numbering-json`.

Without it, the tool falls back to matching by **sequence index**, which
diverges from the label wherever insertions exist (Kabat FR3 82A-C / CDR
27A-/100A-K, IMGT 111A-112F), so FR4 and insertion-letter back-mutations
would be missed. A warning is printed in verbose mode.

Note: `positions_of_interest` in the sequence report uses **1-based sequence
indices**, not Kabat/IMGT labels (the MHC backend scores residue offsets).
Use the structural-risk report for label-aware positions.

**FR4 coverage across the tools:**

| Tool | FR4 handled? |
|------|--------------|
| `immunogenicity_analyzer.py` (sequence) | Yes — MHC-II peptides are generated over the whole chain, FR4 included |
| `structural_risk.py` (structure) | Yes — with `variants_numbering.json`; FR4 uses the same core-weighted buriedness as other positions |
| `../ptm_exposure/analyze_ptm.py` | Yes — PTM motifs scanned over the whole chain (e.g. `H103` W oxidation) |
| `../igfold/igfold_predict.py` | Predicts the full chain (FR4 in the model); does no region analysis itself |

## Usage with Humanization Pipeline

### Step 1: Run humanization

```bash
python3 scripts/humanize/cli.py run \
    --input data/examples/mouse_4d5_fab.fasta \
    --outdir outputs/amg110_step3 \
    --germline-strategy adimab_frequency \
    --donor-structure af3/amg110_structure/rank_1.pdb
```

### Step 2: Analyze immunogenicity (sequence-based)

```bash
# Auto-discovers variants.fasta in step3 output
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --outdir outputs/amg110_step3/step3 \
    --donor data/examples/mouse_4d5_fab.fasta \
    --all-formats
```

### Step 3: Structural risk adjustment (structure-based)

Recommended: run analysis and structural adjustment in one pass. The adjusted
mean directly rewrites `immunogenicity_score` in every report.

```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --outdir outputs/amg110_step3/step3 \
    --donor data/examples/mouse_4d5_fab.fasta \
    --structural-risk \
    --structure-dir af3/amg110_variants \
    --backmutation-dir outputs/amg110_step3/step3 \
    --all-formats
```

The standalone structural step is still available and writes
`*_with_structural_risk.json`:

```bash
python3 tools/immunogenicity/structural_risk.py \
    --immunogenicity-json outputs/immunogenicity/variants_immunogenicity.json \
    --structure-dir af3/amg110_variants \
    --backmutation-dir outputs/amg110_step3/step3 \
    --output-dir outputs/immunogenicity \
    --verbose
```

### Step 4: Interpret results

- **Score**: percent of reference mean; `100%` = same as reference; values **above 100%** mean higher predicted presentation. With `--structural-risk`, both variant and reference means are structurally adjusted.
- **Strong / Weak**: number of **peptide windows** at/above the strong or weak threshold (not residue counts)
- **Density**: `(strong + weak) / total peptide windows`
- **ΔEPI**: epitope score change at back-mutated positions (negative = reduced)
- **StructAdj**: mean buriedness-derived factor (0.1 = buried, 1.0 = exposed); `*` marks a structurally rewritten score
- **Reference**: with `--donor`, the donor is always scored from its own sequence. Without donor, an exact `V0` / `H_V0` / `L_V0` / `combo_V0` variant is used. Use `--reference-variant` to override.

## Output Formats

### Text (stdout)

```
================================================================================
  Immunogenicity Analysis: amg110_variants
================================================================================
  Input:    outputs/amg110_step3/step3/variants.fasta
  Backend:  hlaiipred (apache-2.0)
  Alleles:  DRB1*01:01, DRB1*03:01, DRB1*04:01 ...
  Pep len:  15

  Chain H (7 variants)
  ------------------------------------------------------------------------------
  Variant              Score  Strong    Weak  Density   Mean EPI
  ------------------------------------------------------------------------------
  H_V0                100.0%      12      8   0.245     0.421
  H_V1                 85.3%       9      6   0.198     0.359
```

### JSON

```json
{
  "input_fasta": "variants.fasta",
  "mhc_backend": "hlaiipred",
  "alleles": ["DRB1*01:01", ...],
  "chains": {
    "H": {
      "variants": {
        "H_V0": {
          "immunogenicity_score": 100.0,
          "strong_threshold": 0.5,
          "weak_threshold": 0.3,
          "n_strong": 12,
          "n_weak": 8,
          "n_total_windows": 81,
          "epitope_density": 0.247,
          "mean_epitope_load": 0.421,
          "raw_mean_epitope_load": 0.532,
          "structural_adjustment": 0.79,
          "structural_score_applied": true,
          "adjusted_peptides": [...],
          "positions_of_interest": [...]
        }
      }
    }
  }
}
```

### Markdown + Word + CSV

Generated with `--all-formats`:
- `{stem}_immunogenicity.json`
- `{stem}_immunogenicity.md`
- `{stem}_immunogenicity.docx`
- `{stem}_epitope_heatmap.csv`
- `structural_risk_assessment.md` (from structural_risk.py)
- `structural_risk_assessment.docx` (from structural_risk.py)

## Backend Details

### HLAIIPred (default)

- **Method**: Transformer cross-attention, pan-allele
- **License**: Apache-2.0 (commercial use OK)
- **Reference**: Wang et al., *Commun Biol* (2025)
- **Performance**: 16% increase in prediction of presented peptides vs second-best on unseen alleles; 3% AUC improvement for immunogenic vs non-immunogenic antibodies
- **Download**: ~410 MB (PyTorch CPU + deps); model weights ~9 MB (included in repo)

### Heuristic (fallback)

- **Method**: Germline rarity × surface exposure proxy
- **License**: MIT
- **Note**: Explicitly a non-ML approximation. Always available. Use only when HLAIIPred is unavailable.

### FreeSASA (structural analysis)

- **Method**: Lee-Richards algorithm for SASA computation (FreeSASA default)
- **Reference**: Beer et al., *J Chem Inf Model* (2016)
- **Note**: Used for structural risk adjustment in Step 2

## Scientific Background

### Why MHC-II epitopes matter

The dominant mechanism of anti-drug antibody (ADA) formation is:
1. Antibody is internalized by antigen-presenting cells
2. Proteolyzed into peptide fragments
3. Peptides presented on MHC class II molecules
4. CD4+ T cells recognize MHC-II–peptide complexes
5. T-cell help drives B-cell activation and ADA production

**MHC-II binding is the rate-limiting step.** Predicting which antibody regions form strong MHC-II epitopes allows:
- Ranking humanized variants by ADA risk
- Identifying positions where back-mutations retain/remove epitopes
- Comparing to donor (mouse) for immunogenicity reduction

### Why structural confirmation is critical

**Key Finding**: Buriedness changes significantly between donor and variant structures:

| Position | Donor (V0) | V2 (T1+T2) | V3 (T1+T2+T3) | Risk Change |
|----------|-----------|------------|---------------|-------------|
| H5 | Exposed (0.433) | **Buried** (0.107) | **Buried** (0.113) | Risk DECREASES |
| H12 | Buried (0.145) | **Exposed** (0.654) | **Exposed** (0.552) | Risk INCREASES |
| H16 | Uncertain (0.242) | **Exposed** (0.714) | **Exposed** (0.710) | Risk INCREASES |
| H18 | Uncertain (0.239) | **Exposed** (0.693) | **Exposed** (0.676) | Risk INCREASES |

**Implications**:
- Using donor structure alone UNDERESTIMATES risk for positions that become exposed in variants
- Backmutations can cause local structural changes that alter buriedness
- Variant-specific structures are essential for accurate risk assessment

### Immunogenicity score

The per-variant score compares epitope load to the reference:
- **< 50%**: substantially reduced epitope load (strong improvement)
- **50–80%**: moderately reduced
- **80–100%**: similar to reference (minimal improvement)
- **> 100%**: increased epitope load (regression)

There is no upper cap; a variant with twice the reference load reports 200%.

### Structural adjustment

The structural adjustment factor (0–1) is applied to each peptide's raw score:
- **0.0–0.2**: buried, minimal risk
- **0.2–0.5**: partially buried, reduced risk
- **0.5–0.8**: partially exposed, moderate risk
- **0.8–1.0**: fully exposed, full risk

`mean_epitope_load`, `n_strong`, `n_weak`, `epitope_density`, and
`immunogenicity_score` are all recomputed from the adjusted peptides.

## Key Findings from AMG110 Analysis

### H Chain

| Peptide | Raw Score | Adjustment | Adjusted Score | Risk Level |
|---------|-----------|------------|----------------|------------|
| EVQLLEQSGAEVVKP | 0.998 | ×0.50 | **0.499** | MEDIUM |
| VQLLEQSGAEVVKPG | 0.995 | ×0.60 | **0.597** | MEDIUM-HIGH |
| STAYMQLSSLKASDT | 0.983 | ×0.00 | **0.000** | FALSE POSITIVE |

### L Chain

| Peptide | Raw Score | Adjustment | Adjusted Score | Risk Level |
|---------|-----------|------------|----------------|------------|
| LLIYWASTRESGVPD | 0.997 | ×1.00 | **0.997** | HIGH |
| LIYWASTRESGVPDR | 0.997 | ×1.00 | **0.997** | HIGH |

### Key Conclusions

1. **H5 (L→V, T1)**: Becomes buried in variants → NOT a real risk
2. **H12 (V→K, T2)**: Becomes exposed in variants → IS a real risk
3. **STAYMQLSSLKASDT**: No backmutation positions → FALSE POSITIVE
4. **L63 (T→S, T2)**: Surface-exposed + CDR contact → HIGHEST PRIORITY

## References

1. Wang et al. HLAIIPred: Cross-attention for HLA-II peptide presentation. *Commun Biol* (2025)
2. Raybould et al. TAP: Five developability guidelines. *PNAS* (2019)
3. Tien et al. A dataset and consensus approach for identifying the relative solvent accessibility of amino acids in proteins. *PLOS ONE* (2013)
4. Beer et al. FreeSASA: An implementation of a solvent accessible surface area calculation. *J Chem Inf Model* (2016)
