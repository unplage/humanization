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

- **Input**: Variant PDB structure (AF3/IgFold)
- **Method**: Computes relative SASA (relSASA) for each residue
- **Output**: Risk adjustment factors based on buriedness
- **Key insight**: Buried residues (relSASA < 0.20) are less accessible to MHC-II antigen processing

### Risk Adjustment Logic

| relSASA | Classification | Risk Factor | Rationale |
|---------|----------------|-------------|-----------|
| < 0.10 | Deeply buried | 0.05 | Minimal accessibility |
| 0.10-0.15 | Partially buried | 0.05-0.30 | Low accessibility |
| 0.15-0.25 | Uncertain zone | 0.30-0.70 | Linear interpolation |
| 0.25-0.40 | Partially exposed | 0.70-0.90 | Moderate accessibility |
| > 0.40 | Fully exposed | 1.00 | Full risk |

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

# 4. Clone and install HLAIIPred (Apache-2.0 license)
git clone https://github.com/pfizer-opensource/HLAIIPred.git /tmp/HLAIIPred
cd /tmp/HLAIIPred && pip install -e .

# 5. Copy model files to project directory
mkdir -p tools/immunogenicity/models
cp /tmp/HLAIIPred/models/epT_0.pt tools/immunogenicity/models/
cp /tmp/HLAIIPred/models/epT_1.pt tools/immunogenicity/models/

# 6. Install structural risk assessment dependencies
pip install freesasa python-docx

# 7. Cleanup
rm -rf /tmp/HLAIIPred

# 8. Verify installation
python3 tools/immunogenicity/immunogenicity_analyzer.py --check
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
# Step 1: Run sequence-based analysis
python3 tools/immunogenicity/immunogenicity_analyzer.py \
    --outdir outputs/amg110_step3/step3 \
    --all-formats

# Step 2: Run structural risk adjustment
python3 tools/immunogenicity/structural_risk.py \
    --immunogenicity-json outputs/immunogenicity/variants_immunogenicity.json \
    --structure-dir af3/amg110_variants \
    --backmutation-dir outputs/amg110_step3/step3 \
    --output-dir outputs/immunogenicity \
    --verbose
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

```bash
# Use variant-specific structures (AF3/IgFold)
python3 tools/immunogenicity/structural_risk.py \
    --immunogenicity-json outputs/immunogenicity/variants_immunogenicity.json \
    --structure-dir af3/amg110_variants \
    --backmutation-dir outputs/amg110_step3/step3 \
    --output-dir outputs/immunogenicity \
    --verbose
```

### Step 4: Interpret results

- **Score**: variant immunogenicity relative to donor (100% = same as donor, <100% = reduced)
- **Strong**: number of residues in strong predicted epitopes (sigmoid > 0.5)
- **Density**: fraction of residues in any epitope window
- **ΔEPI**: epitope score change at back-mutated positions (negative = reduced)
- **Structural Adjustment**: risk factor based on buriedness (0.1 = buried, 1.0 = exposed)

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
          "epitope_density": 0.245,
          "structural_adjustment": 0.75,
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

- **Method**: Shrake-Rupley algorithm for SASA computation
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

The per-variant score (0–100%) compares epitope load to the donor:
- **< 50%**: substantially reduced epitope load (strong improvement)
- **50–80%**: moderately reduced
- **80–100%**: similar to donor (minimal improvement)
- **> 100%**: increased epitope load (regression)

### Structural adjustment

The structural adjustment factor (0–1) modifies the immunogenicity score:
- **0.0–0.2**: position is buried, minimal risk
- **0.2–0.5**: position is partially buried, reduced risk
- **0.5–0.8**: position is partially exposed, moderate risk
- **0.8–1.0**: position is fully exposed, full risk

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
