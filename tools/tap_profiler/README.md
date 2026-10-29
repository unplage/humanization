# TAP (Therapeutic Antibody Profiler) Analysis Tool

Implements the five published TAP developability guidelines on antibody
structures or sequences:

| Metric | Meaning |
|--------|---------|
| **Ltot** | Total IMGT CDR length of the Fv (six CDR loops) |
| **PSH** | Patches of surface hydrophobicity over the CDR vicinity |
| **PPC** | Patches of positive charge over the CDR vicinity |
| **PNC** | Patches of negative charge over the CDR vicinity |
| **SFvCSP** | Structural Fv charge symmetry parameter (Q_VH × Q_VL) |

References:

1. Raybould MIJ et al. *Five computational developability guidelines for
   therapeutic antibody profiling.* PNAS 2019;116(10):4025–4030.
2. Raybould MIJ et al. *Contextualising the developability risk of antibodies
   with lambda light chains using enhanced therapeutic antibody profiling.*
   Commun Biol 2024;7:62 (thresholds, `--thresholds 2024`, the default).

The metric definitions in `tap_core.py` follow the published Methods
verbatim: pairwise `H(R1)·H(R2)/r12²` sums over surface residue pairs within
7.5 Å in the CDR vicinity, charges zeroed in salt bridges (≤ 3.2 Å), and
SFvCSP as the **product** of the net VH and VL surface charges.

## What changed in the rewrite

The previous implementation produced numbers that could not be compared with
any published TAP value, so every structure was mis-flagged. Fixed:

| # | Legacy bug | Effect | Fix |
|---|-----------|--------|-----|
| 1 | PSH/PPC/PNC counted *sequence-order* runs of ≥3 exposed residues (after deleting buried residues) | Values were integer patch counts (0–5) compared against published area/score bands (96–212) → **PSH red on every Fv** | Published pairwise `x1·x2/r²` sums over the 3D CDR vicinity |
| 2 | SFvCSP computed `Q_VH − Q_VL` (difference) with product-based thresholds | Wrong sign/magnitude for any charged Fv | `Q_VH × Q_VL` (product), per PNAS 2019 |
| 3 | Ltot summed CDRs of **all** chains (antigen, CH1, CL included) | e.g. Fab+antigen: Ltot 90 instead of 47 | Only ANARCI-classified Fv chains are analysed |
| 4 | VH/VL assigned by chain letter or chain length (`A` or len>115 → VH) | Antigen chains entered SFvCSP/PSH | ANARCI chain-type detection, `--vh-chain/--vl-chain` override |
| 5 | MSMS `.area` parser looked for `ATOM`-prefixed lines (MSMS writes ` REMARK`) | With MSMS installed all SASA silently became 0 | FreeSASA/JOY psa backends; MSMS support dropped |
| 6 | Patches merged across chain boundaries and buried gaps | Inflated/merged patches | 3D proximity defines patches; no sequence runs |
| 7 | Thresholds: Ltot red >65 (published >63), PSH red >212 (published >201.6), amber/green gaps defaulting to red | Wrong flags near boundaries | Exact published bands, amber→green→red precedence |
| 8 | README example output did not even match the code (`PSH 3.00 → GREEN`) | Misleading docs | README regenerated from actual runs |

Validation: on the Exscientia reference structure the rewrite reproduces the
independent MIT implementation of the same method **to full precision**
(Ltot 64, PSH 199.813, PPC 0.323, PNC 0.466, SFvCSP −2.880 — see
`tests/test_tap_profiler.py::test_reference_parity`).

## Requirements

```bash
pip install anarci        # IMGT numbering (required for structure/sequence modes)
pip install freesasa      # surface accessibility (required unless JOY psa is used)
pip install python-docx   # optional, Word reports
```

Optional, for exact parity with the official TAP surface reference:

```bash
# JOY psa (Ala-X-Ala side-chain relative ASA); academic JOY distribution
# contains the `psa` binary - point at it with --psa-binary or put it on PATH
```

## Usage

### Check available tools

```bash
python3 tools/tap_profiler/tap_analyzer.py --check
```

### Structure analysis (recommended)

```bash
# Default: ANARCI IMGT numbering, FreeSASA surface, Fv-only scope, 2024 thresholds
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb "af3/c45_af3/Structure Prediction (Boltz-2)/rank_1.pdb"

# All report formats into one directory
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb --all-formats --output outputs/tap/

# Exact TAP surface reference (JOY psa) instead of the FreeSASA fallback
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb structure.pdb --sasa psa --psa-binary /path/to/psa

# Structure is already IMGT-numbered (skips the ANARCI lookup)
python3 tools/tap_profiler/tap_analyzer.py --pdb structure.pdb --numbering native

# Keep antigen/constant chains so they occlude the Fv surface
python3 tools/tap_profiler/tap_analyzer.py --pdb structure.pdb --scope complex

# Explicit chain assignment / older threshold set
python3 tools/tap_profiler/tap_analyzer.py --pdb structure.pdb \
    --vh-chain H --vl-chain L --thresholds 2019
```

### Sequence-only analysis (Ltot only)

```bash
python3 tools/tap_profiler/tap_analyzer.py \
    --vh <VH_SEQUENCE> --vl <VL_SEQUENCE> --name my_antibody
```

PSH/PPC/PNC/SFvCSP require a structure; sequence mode reports them as `n/a`.

### Batch processing

```bash
python3 tools/tap_profiler/tap_analyzer.py \
    --pdb-dir af3/ --output outputs/tap_batch/ --all-formats
```

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `--sasa {auto,psa,freesasa}` | `auto` | `auto` uses JOY psa when found, else FreeSASA |
| `--psa-binary PATH` | — | JOY psa executable (exact Ala-X-Ala reference) |
| `--numbering {anarci,native}` | `anarci` | `native` = structure already IMGT-numbered |
| `--scope {fv,complex}` | `fv` | `fv` = Fv chains only (official TAP protocol) |
| `--thresholds {2024,2019,exscientia,legacy}` | `2024` | Flagging bands |
| `--vh-chain / --vl-chain` | auto | Explicit chain ids (ANARCI detection otherwise) |
| `--json / --markdown / --word / --all-formats` | off | Report formats |
| `--output DIR` | `outputs/tap_profiler` | Report directory |

## Thresholds

`--thresholds 2024` (default) — Raybould 2024, ABodyBuilder2 models of 664
post-Phase-I clinical-stage therapeutics (CSTs). Amber/green boundaries;
anything outside is red.

| Metric | Green | Amber | Red |
|--------|-------|-------|-----|
| Ltot | 43–54 | 37–42, 55–63 | <37, >63 |
| PSH | 110.11–168.06 | 95.58–110.11, 168.06–201.59 | <95.58, >201.59 |
| PPC | 0–1.32 | 1.32–4.22 | >4.22 |
| PNC | 0–2.00 | 2.00–4.42 | >4.42 |
| SFvCSP | ≥ −6.00 | −30.60 to −6.00 | < −30.60 |

Other sets: `2019` (ABodyBuilder1 bands tabulated in the 2024 paper),
`exscientia` (bands hard-coded in Exscientia/ab-characterisation),
`legacy` (the pre-rewrite tool's values — kept only so old reports stay
interpretable; **do not use for new analyses**).

### Risk summary

| Condition | Level |
|-----------|-------|
| all green | LOW RISK |
| 1 amber | LOW-MEDIUM RISK |
| 2+ amber | MEDIUM RISK |
| 1 red | MEDIUM-HIGH RISK |
| 2+ red | HIGH RISK |

## Surface engine

TAP defines a residue as exposed when its **side-chain** relative accessible
surface area (Ala-X-Ala reference) is ≥ 7.5 %, as reported by the JOY `psa`
program used by the official TAP pipeline.

* `--sasa psa` (or `auto` with `psa` on PATH): exact TAP definition.
* FreeSASA fallback (default when psa is absent): whole-residue rSASA vs
  Tien et al. 2013 maxima with the same 7.5 % cut-off. Measured agreement
  with the psa surface/buried decision: **92.9 %** of residues on a
  Fab+antigen test structure; PSH moves by a few units, flags rarely change.
  The report always states which engine produced the numbers.

## Scope

* `--scope fv` (default): only the ANARCI-classified Fv chains are written
  for SASA and annotated — matches the official TAP protocol (Fv models).
  Constant regions (CH1/CL) and antigen chains are dropped.
* `--scope complex`: all chains are kept, so antigen contacts bury Fv
  residues; CDR-vicinity scores can then include neighbouring antigen
  surface residues. Useful for checking a bound-Fv context.

## Methodology (implemented exactly as published)

```
surface   : relASA >= 7.5 % (psa: side chain / Ala-X-Ala; FreeSASA: residue / Tien)
neighbour : closest heavy-atom distance <= 7.5 A
CDR       : IMGT 27-38, 56-65, 105-117 (insertions included)
anchors   : IMGT 25, 26, 39, 40, 54, 55, 66, 67, 103, 104, 118, 119
vicinity  : exposed CDR+anchor seeds, plus exposed residues with a heavy
            atom within 4 A of a seed
salt bridge: K/R N vs D/E O <= 3.2 A, both exposed -> charge 0,
            hydrophobicity := glycine
H(R)      : Kyte-Doolittle normalised to [1, 2]
Q(R)      : D/E -1, K/R +1, H +0.1, else 0

PSH  = sum over ordered vicinity pairs (r12 <= 7.5 A) of H(R1)*H(R2)/r12^2
PPC  = same with |Q(R1)|*|Q(R2)|/r12^2, both charges > 0
PNC  = same, both charges < 0
SFvCSP = (sum of Q over exposed VH residues) * (sum over exposed VL)
Ltot = number of IMGT CDR residues in the Fv
```

## Tests

```bash
python3 tests/test_tap_profiler.py
# optional exact-parity check against the Exscientia reference structure:
TAP_PARITY_PDB=/path/test_complex_reference.pdb \
TAP_PARITY_PSA=/path/to/psa python3 tests/test_tap_profiler.py
```

Covers threshold logic, pairwise metric arithmetic on hand-computed cases,
salt-bridge geometry, structure-vs-sequence Ltot consistency, antigen
exclusion, determinism, all report formats and the CLI end-to-end.

## Relation to other open-source tools

| Tool | Licence | Status | Use it for |
|------|---------|--------|-----------|
| **This tool** | repo licence | maintained | TAP metrics inside the humanization pipeline |
| [Exscientia/ab-characterisation](https://github.com/Exscientia/ab-characterisation) | MIT | 2024 | Independent TAP reference implementation (co-authored by the TAP authors); full pipeline also needs Rosetta/ChimeraX |
| [Ginkgo abdev-benchmark](https://github.com/ginkgobioworks/abdev-benchmark) | MIT | 2025 | ML model benchmark on HIC/Tm2/Titer/PR_CHO/AC-SINS — consumes TAP features, does not replace them |
| TDC `Develop(name='TAP')` | Apache-2.0 | 2022+ | Published TAP values for ~100 named therapeutics (validation data) |
| OPIG TAP web server | academic only | — | Official tool; not licensed for commercial pipelines |

Methodology (formulas, thresholds) is from the open papers above; only
software licences restrict redistribution. `JOY psa` and `ABodyBuilder` are
academic-licensed — this tool therefore treats them as optional and keeps
the default path on permissive dependencies: ANARCI (BSD-3-Clause), FreeSASA
(MIT), Biopython (BSD), python-docx (MIT).
