# Structure-aware pI / Fv Charge Analysis

`pi_analyzer.py` estimates the isoelectric point (pI), net charge, charge
asymmetry, charge patches and developability risks of an antibody **Fv** or
**VHH**. It uses a PDB/mmCIF structure when available and falls back to a
sequence-only calculation when no structure exists.

The structure path reuses the TAP/ANARCI/FreeSASA annotations in
`tools/tap_profiler/tap_core.py`; the sequence path is pure Python stdlib and
runs anywhere.

## Why structure-aware pI?

Sequence pI assumes every ionizable group has its model pKa. In a folded Fv:

- buried D/E/K/R groups are less solvated and shift pKa;
- salt bridges and H-bonds stabilise the charged form;
- disulfide Cys is no longer ionizable;
- N/C termini may be free, capped or blocked.

These effects can move an Fv pI by several tenths of a pH unit and change the
net charge at the formulation pH. The structure-aware value is therefore more
useful for formulation, IEX and charge-patch engineering, but it remains a
computational estimate - confirm by cIEF/IEX.

## pKa sets

`--pka-set` accepts the canonical names or aliases:

| Name | Aliases | N-term | C-term | D | E | C | Y | H | K | R |
|------|---------|--------|--------|---|---|---|---|---|---|---|
| `Bjellqvist` | `exPASy`, `expasy` | 7.50 | 3.55 | 4.05 | 4.45 | 9.00 | 10.00 | 5.98 | 10.00 | 12.00 |
| `IPC_protein` | `ipc` | 9.094 | 2.869 | 3.872 | 4.412 | 7.555 | 10.85 | 5.637 | 9.052 | 11.84 |
| `EMBOSS` | - | 8.60 | 3.60 | 3.90 | 4.10 | 8.50 | 10.10 | 6.50 | 10.80 | 12.50 |
| `Sillero` | - | 8.20 | 3.20 | 4.00 | 4.50 | 9.00 | 10.00 | 6.40 | 10.40 | 12.00 |
| `Nozaki_Tanford` | `nozaki`, `nozaki-tanford` | 7.50 | 3.80 | 4.00 | 4.40 | 9.50 | 9.60 | 6.30 | 10.40 | 12.00 |

Bjellqvist also uses residue-specific N-terminal overrides
(A 7.59, M 7.00, S 6.93, P 8.36, T 6.82, V 7.44, E 7.70) and C-terminal
overrides (D 4.55, E 4.75).

Default: `IPC_protein`. Use `--all-pka-sets` to report every set and see how
sensitive the result is to the pKa model. For comparison with classic ExPASy
numbers, use `--pka-set Bjellqvist`.

## Structure model

### `--structure-model heuristic` (fallback)

For every ionizable side chain:

1. Compute relative SASA with FreeSASA (or JOY `psa` if available).
2. Apply a burial penalty: acidic groups raise pKa when buried; basic groups
   lower pKa when buried.
3. Apply salt-bridge and H-bond corrections (acidic pKa down, basic pKa up
   when the charged form is stabilised), screened by Debye-Huckel ionic
   strength (`--ionic-strength`, default 0.15 M).
4. Clamp the total shift to +/-2 pH units.

Disulfide Cys is detected geometrically (SG-SG <= 2.5 A) or, if no SG atoms
are present, from the conserved Kabat Cys positions (H 22/92, L 22/23/88).
Disulfide Cys is treated as non-ionizable.

Termini are reported as:

- **free**: both N- and C-termini charged (default for VHH);
- **capped**: C-terminus amidated/blocked (default for Fab);
- **blocked**: both termini blocked.

`--termini auto` selects capped for Fab and free for VHH. The selected model is
used for pI, net charge, patches and recommendations; free/capped values are
still reported so the construct-dependent difference is visible. If a chain is
truncated before IMGT 1/127, that end is treated as blocked and a warning is
emitted.

### `--structure-model propka` (default)

Runs PROPKA3 if `propka3`/`propka` is on `PATH`, parses the
`SUMMARY OF THIS PREDICTION` section, uses per-residue pKa values and the
predicted N/C-terminal pKa values when present.
If PROPKA is missing the run **falls back to the heuristic model** with a
warning (use `heuristic` explicitly to suppress the warning, or install
`propka3`).

### `--structure-model none`

Uses the model pKa values without microenvironment corrections. Useful as a
control to see how much burial/salt-bridge/H-bond effects change the result.

### Structure quality

Before interpreting a structure, the tool reports:

- mean/min/max CA B-factor (a pLDDT proxy for AF2/AF3 models);
- missing side-chain heavy atoms, especially missing ionizable atoms;
- CA-CA chain breaks (> 4.5 A);
- truncated N/C termini.

Low pLDDT/B-factor models, missing ionizable atoms and chain breaks produce
warnings because they directly affect SASA and pI.

### Salt bridges and pKa uncertainty

Two salt-bridge definitions are reported deliberately:

- **TAP**: 3.2 A, both partners surface-exposed, one partner per residue;
- **pKa microenvironment**: 4.0 A, burial allowed, multiple partners.

They answer different questions and can differ substantially (for example,
A reference Fab has 1 TAP salt bridge but 7 microenvironment pairs). The pKa model uses the
broader definition; the TAP metric remains comparable with published TAP
output. Every run also reports the pI range across all five bundled pKa sets
and a net-charge-vs-pH curve (`--ph-curve` prints it; JSON always contains it).

## Sequence fallback

If no structure is supplied:

```bash
python3 tools/pi_analysis/pi_analyzer.py --fasta inputs/antibody.fasta
python3 tools/pi_analysis/pi_analyzer.py --vh <VH_SEQ> --vl <VL_SEQ>
python3 tools/pi_analysis/pi_analyzer.py --vh <VHH_SEQ>
```

Sequence mode reports pI, capped pI, net charge and per-chain charge. It does
not report TAP, SASA, charge patches or structure-aware pKa shifts.

FASTA handling:

- one record -> VHH;
- two records -> VH/VL; headers containing VH/VL, H/L, heavy/light, kappa or
  lambda are used to assign chains even when VL appears first;
- ambiguous two-record files are treated as VH then VL.

## Requirements

Structure mode:

```bash
pip install anarci        # IMGT numbering
pip install freesasa      # SASA
pip install propka        # default structure pKa model (falls back to heuristic)
```

Sequence mode: no third-party dependencies.

JOY `psa` is optional and gives the exact TAP Ala-X-Ala surface reference.
Without it, FreeSASA uses whole-residue Tien rSASA; approximately 7 % of
surface/buried decisions can differ from official TAP.

## Usage

```bash
# Check availability
python3 tools/pi_analysis/pi_analyzer.py --check

# Structure-aware analysis, all report formats
python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --all-formats

# Compare every pKa set
python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --all-pka-sets

# PROPKA pKa values
python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --structure-model propka

# Explicit chains and formulation pH
python3 tools/pi_analysis/pi_analyzer.py --pdb structure.pdb \
    --vh-chain A --vl-chain B --formulation-ph 5.5

# Batch
python3 tools/pi_analysis/pi_analyzer.py --pdb-dir af3/ --all-formats

# Sequence fallback
python3 tools/pi_analysis/pi_analyzer.py --fasta inputs/antibody.fasta --all-formats
```

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `--pdb PATH` | - | PDB/mmCIF structure |
| `--pdb-dir DIR` | - | Recursive structure batch |
| `--fasta PATH` | - | One VHH or two VH/VL FASTA records |
| `--vh`, `--vl` | - | Direct sequences (omit `--vl` for VHH) |
| `--pka-set SET` | `IPC_protein` | pKa model |
| `--all-pka-sets` | off | Report all five pKa sets |
| `--ph` | `7.4` | Physiological pH for net charge |
| `--formulation-ph` | `6.0` | Formulation pH for charge/patch reporting |
| `--disulfide {auto,exclude,free}` | `auto` | Cys handling |
| `--structure-model {propka,heuristic,none}` | `propka` | Structure pKa model (propka falls back to heuristic if PROPKA3 is missing) |
| `--termini {auto,free,capped,blocked}` | `auto` | Terminus model (Fab capped, VHH free) |
| `--ionic-strength M` | `0.15` | Debye-Huckel screening ionic strength |
| `--ph-grid START:STOP:STEP` | `4:10:0.5` | pH grid for the charge curve |
| `--ph-curve` | off | Print charge-vs-pH in the text report |
| `--vh-chain`, `--vl-chain` | auto | Explicit structure chain IDs (`--vh-chain` alone = VHH) |
| `--sasa {auto,psa,freesasa}` | `auto` | Surface engine |
| `--psa-binary PATH` | - | JOY `psa` executable |
| `--tap-thresholds {2024,2019,exscientia,legacy}` | `2024` | TAP flag bands |
| `--max-patches N` | `12` | Patches shown in text/Markdown; 0 = all |
| `--ionizable` | off | Include per-residue pKa table |
| `--json`, `--markdown`, `--csv`, `--text`, `--all-formats` | off | Reports |
| `--output DIR` | `outputs/pi_analysis` | Report directory |

## Outputs

`--all-formats` writes:

| File | Contents |
|------|----------|
| `<stem>_pi.txt` | Human-readable report |
| `<stem>_pi.json` | Full machine-readable result |
| `<stem>_pi.md` | Markdown report |
| `<stem>_pi.csv` | One-row summary |
| `<stem>_pi_patches.csv` | One row per charge patch |

Batch runs also write `combined_pi_results.json`.

The JSON contains:

- `sequence_pI` and `structure_pI` per pKa set (`pI_selected`, `pI_free`,
  `pI_capped`, `pI_blocked`, `delta_vs_sequence`);
- `sensitivity`: selected pI for every bundled pKa set and the resulting range;
- `charge_curve`: net charge sampled across pH;
- `structure_quality`: B-factor/pLDDT proxy, missing side chains, chain breaks,
  truncated termini;
- per-chain pI (`free`/`capped`/`selected`) and net charge;
- `fv`: selected/free/capped pI, pI range, net charge, total and surface
  `Fv_chml`/`SFvCSP`, total and surface per-chain charges;
- `tap`: Ltot, PSH, PPC, PNC, SFvCSP, flags, risk summary, and both salt-bridge
  counts (`n_salt_bridges_tap` and `n_salt_bridges_pka`);
- `patches`: 3D charge patches with residue rows, local risk basis and
  mutation candidates with sequence-level what-if `delta_pI_estimate` /
  `delta_charge_formulation_estimate`;
- `ionizable_residues`: structure-aware pKa, exposure, salt bridge distance,
  H-bond distance, disulfide and shift reason;
- `recommendations` and `warnings`.

## Interpretation

### Fv pI windows

| Fv pI | Interpretation | Practical direction |
|-------|----------------|---------------------|
| `< 6.5` | acidic Fv | Raise toward 7.5-9.0 by neutralising exposed framework D/E; check IEX/HCP separation and Fc charge match |
| `6.5-7.5` | near-neutral | Avoid formulation within ~0.5 pH units; check self-interaction and viscosity |
| `7.5-9.0` | favourable weakly basic | Maintain during affinity/developability optimisation |
| `> 9.0` | strongly basic | Lower toward 8.0-8.5 by neutralising exposed framework K/R; check heparin/polyreactivity |

The 7.5-9.0 window follows Gupta et al. 2022 and Thorsteinson et al. 2021.
It is a developability heuristic, not a hard requirement.

### Charge asymmetry

- `Fv_chml = Q_VH - Q_VL` at pH 7.4. Large absolute values (`|Fv_chml| > 4`)
  indicate VH/VL charge imbalance and are associated with attractive
  self-interactions and higher viscosity.
- `SFvCSP = Q_VH * Q_VL` is the TAP charge-symmetry parameter. The TAP 2024
  green region is `>= -6.0`; values below `-6.0` are amber.
- Structure mode reports **total** and **surface** versions separately:
  total uses all modelled ionizable groups, surface uses TAP's exposed-residue
  definition. Sequence mode reports total charge (surface values are `null`).
  Use the surface `SFvCSP` for TAP comparability and the total values for
  biochemical charge balance.
- Balance the more charged chain using **exposed framework** positions; keep
  CDR, anchors, VHH hallmarks and FR4 untouched.

### TAP flags

| Metric | Meaning |
|--------|---------|
| `Ltot` | Total IMGT CDR length |
| `PSH` | Surface hydrophobic patch score in the CDR vicinity |
| `PPC` | Positive-charge patch score in the CDR vicinity |
| `PNC` | Negative-charge patch score in the CDR vicinity |
| `SFvCSP` | VH x VL surface-charge product |

TAP flags are structure-only. See `tools/tap_profiler/README.md` for the
published thresholds and the exact surface-reference caveat.

### Charge patches

Patches are built at both the formulation pH and pH 7.4:

- exposed charged residues with rSASA >= 7.5 %;
- clustered by 3D proximity (<= 8 A);
- charge computed at the requested pH (salt-bridged groups zeroed);
- score is the pairwise `|q1*q2|/r^2` sum;
- `score_tap` is the TAP-style ordered integer-charge score (`PPC`/`PNC`);
- local risk is recalibrated for continuous partial charges: positive patches
  use 0.50/2.00 bands, negative patches 0.60/2.50, with a 1.25x CDR-vicinity
  penalty; `risk_basis` states which rule fired;
- first/last residues of each chain are flagged and excluded from mutation
  candidates because they change the terminal pKa and can form
  pyroglutamate/processing liabilities.

Each patch reports its residues, CDR/framework composition, VH/VL interface
flag and candidate substitutions.

## Engineering recommendations

### By pI

- **Low pI (`< 6.5`)**: neutralise exposed framework D/E to Q/N/S/T/A; check
  whether the Fc isotype matches the Fv charge (IgG1 Fc is positive,
  IgG2/IgG4 Fc is negative); consider anion-exchange/HCP separation early.
- **Near-neutral (`6.5-7.5`)**: avoid formulating within 0.5 pH units of the
  pI; screen colloidal stability (AC-SINS, BVP, DLS) across pH 5.0-7.5.
- **Favourable (`7.5-9.0`)**: keep the window during optimisation; still check
  patches, PTMs and hydrophobicity.
- **High pI (`> 9.0`)**: neutralise exposed framework K/R; test heparin
  binding and polyreactivity; avoid charge reversal unless necessary.

### By patch

- **Positive patch**: prefer framework K/R -> Q/N/S/T/A.
- **Negative patch**: prefer framework D/E -> Q/N/S/T/A.
- **Mixed patch**: only act when one sign dominates (`|total charge| > 0.5`);
  otherwise the patch is approximately balanced.
- Never mutate CDR residues, anchors, VHH hallmark positions (IMGT 42/49/50/52;
  Kabat 37/44/45/47) or FR4.
- VH/VL interface positions are flagged **high risk**: mutation can perturb
  chain pairing. Validate pairing and affinity.
- Charge reversal (K/R -> D/E or D/E -> K/R) is reported separately and should
  be treated as a last resort.

### Follow-up workflow

After selecting candidate positions:

1. **TAP** - confirm PPC/PNC/SFvCSP and overall risk:
   `python3 tools/tap_profiler/tap_analyzer.py --pdb <pdb> --all-formats`
2. **PTM exposure** - check deamidation/oxidation/glycosylation motifs:
   `python3 tools/ptm_exposure/analyze_ptm.py --pdb <pdb> --all-formats`
3. **ThermoMPNN** - filter destabilising substitutions:
   `python3 tools/stability_predictor/stability_analyzer.py --pdb <pdb> --chain all --positions <labels> --thermompnn-dir <repo>`
4. **BioPhi/Sapiens** - check humanness after framework changes:
   `python3 tools/biophi/biophi_analyze.py sapiens --input <fasta> --outdir <dir>`
5. **ProteinMPNN** - generate/score alternative framework sequences:
   `python3 tools/mpnn/proteinmpnn_adapter.py --pdb_path <pdb> --out_folder <dir>`
6. **Structure validation** - predict variants (AF3/Boltz/IgFold) and compare
   CDR RMSD, pLDDT and clashes before ordering.

## Example: reference Fab

On the AF3 reference Fab structure (`af3/antibody.pdb`, VH=A, VL=B):

```text
Fv pI (selected, capped):    5.40
Fv pI range (5 pKa sets):    5.40-5.71
Fv pI (free termini):        5.00
Net charge @ pH 7.4:        -2.38
Fv_chml total / surface:    -8.17 / -9.00
SFvCSP total / surface:    -15.26 / -18.00
TAP: Ltot 49 green, PSH 126.3 green, PPC 0.00 green,
     PNC 0.24 green, SFvCSP -18.0 amber
Salt bridges: TAP 1 / pKa microenvironment 7 (12 residues)
```

Interpretation: acidic Fv with strong VH/VL charge asymmetry. The VH is
strongly negative (`-5.27`) and the VL is positive (`+2.89`). The main levers
are neutralising exposed VH framework acids and checking whether VL positive
patches are needed for antigen binding. The capped-Fab pI is ~5.4, so a pH 6.0
formulation is only ~0.6 units away; the free/capped difference (0.4 units)
matters for construct design. Confirm colloidal stability and IEX/HCP
behaviour experimentally.

## Validation and limitations

- Bjellqvist per-chain pI matches Biopython `IsoelectricPoint` to < 0.02 pH
  units for normal antibody domains. Biopython clamps its search to
  `[4.05, 12.0]`, so extremely acidic/basic peptides can differ.
- The heuristic structure model is a screening tool, not a replacement for
  PROPKA or experimental pI. Electrostatic coupling beyond the screened
  Debye-Huckel term, specific-ion effects, excipients and conformational
  ensembles are not modelled.
- FreeSASA vs JOY `psa` changes ~7 % of surface/buried decisions; TAP parity
  requires `--sasa psa`.
- Structure mode analyses the Fv/VHH only. It does not model the Fc, glycan,
  antigen or formulation buffer.
- Sequence mode cannot detect burial, salt bridges, disulfides by geometry or
  charge patches.
- Always confirm pI by cIEF and charge variants by IEX; confirm self-interaction
  by AC-SINS/BVP/SEC and viscosity at the intended concentration.

## References

1. Bjellqvist B et al. *The focusing positions of polypeptides in immobilized
   pH gradients can be predicted from their amino acid sequences.*
   Electrophoresis 1993;14:1023-1031.
2. Kozlowski LP. *IPC - Isoelectric Point Calculator.* Biol Direct 2016;11:55.
3. Kozlowski LP. *IPC 2.0: prediction of isoelectric point and pKa
   dissociation constants.* Nucleic Acids Res 2021;49(W1):W285-W292.
4. Thorsteinson N et al. *Structure-based charge calculations for predicting
   antibody developability.* mAbs 2021;13(1):1863582.
5. Gupta P et al. *Antibody developability profiles and Fv pI.*
   mAbs 2022;14(1):2020203.
6. Raybould MIJ et al. *Five computational developability guidelines for
   therapeutic antibody profiling.* PNAS 2019;116(10):4025-4030.
7. Raybould MIJ et al. *Contextualising the developability risk of antibodies
   with lambda light chains using enhanced therapeutic antibody profiling.*
   Commun Biol 2024;7:62.
8. Olsson MHM et al. *PROPKA3: consistent treatment of internal and surface
   residues in empirical pKa predictions.* J Chem Theory Comput
   2011;7:525-537.

## Files

| File | Purpose |
|------|---------|
| `pi_core.py` | pI/pKa/patch/TAP logic |
| `pi_analyzer.py` | CLI and report renderers |
| `README.md` | This document |
| `../../tests/test_pi_analysis.py` | Regression tests |
