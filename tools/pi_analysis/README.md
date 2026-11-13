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

The N- and C-terminal pKa are **residue-specific in every bundled set**: the
side chain of the terminal residue shifts the pKa of the adjacent group (an
N-terminal Asp is less basic, an N-terminal Pro more basic, a C-terminal
Asp/Glu raises the pKa of the alpha-carboxyl because the resulting dianion
is destabilised). The offsets are the residue-specific model pKa values of
Grimsley et al. (2009) for the N-terminus and of Bjellqvist (1993) for the
C-terminus, expressed relative to the generic value of each set, and they
are applied to **all** sets:

| terminal residue | N-term offset | C-term offset |
|---|---|---|
| A | +0.09 | - |
| C | 0.00 | - |
| D | -0.50 | +1.00 |
| E | +0.20 | +1.20 |
| M | -0.50 | - |
| P | +0.86 | - |
| S | -0.57 | - |
| T | -0.68 | - |
| V | -0.06 | - |

Because the offset is a property of the residue and not of the reference
set, the terminal behaviour of a given antibody no longer depends on the
arbitrary choice of pKa set. For `Bjellqvist` this reproduces the published
table exactly (A 7.59, M 7.00, S 6.93, P 8.36, T 6.82, V 7.44, E 7.70,
D 7.00 for the N-terminus; D 4.55, E 4.75 for the C-terminus). The only
visible consequence is that an **N-terminal Asp** now uses the Grimsley
model pKa of 7.00 where ExPASy/Biopython use the generic 7.50, which moves
the pI of such a chain by about 0.02 pH units.

Default: `IPC_protein`. Use `--all-pka-sets` to report every set and see how
sensitive the result is to the pKa model. For comparison with classic ExPASy
numbers, use `--pka-set Bjellqvist`.

## Structure model

### `--structure-model heuristic` (fallback)

For every ionizable side chain:

1. Compute relative SASA with FreeSASA (or JOY `psa` if available).
2. Apply a burial penalty: acidic groups raise pKa when buried (the charged
   anion is destabilised in a hydrophobic environment), basic groups lower
   pKa when buried (the charged cation is destabilised). Histidine follows
   the basic convention.
3. Apply salt-bridge and H-bond corrections (acidic pKa down, basic pKa up
   when the charged form is stabilised), screened by Debye-Huckel ionic
   strength (`--ionic-strength`, default 0.15 M).
4. Clamp the total shift to +/-2 pH units.

Disulfide Cys is detected geometrically (SG-SG <= 2.5 A) or, if no SG atoms
are present, from the conserved Kabat Cys positions (H 22/92, L 22/23/88).
A conserved Cys is only assigned to the disulfide when **both** partners are
present: a lone conserved Cys is a free thiol and stays ionizable.
Disulfide Cys is treated as non-ionizable.

Termini are reported as:

- **free**: both N- and C-termini charged (default for VHH);
- **capped**: C-terminus amidated/blocked (default for Fab);
- **blocked**: both termini blocked.

`--termini auto` selects capped for Fab and free for VHH. The selected model is
used for pI, net charge, patches and recommendations; free/capped values are
still reported so the construct-dependent difference is visible.

The completeness of the modelled C-terminus is checked per chain type, because
the last FR4-IMGT position differs: a VH/VHH runs to IMGT 128 while a VL ends
at IMGT 127 (as ANARCI numbers them). A chain truncated before its own
boundary is treated as blocked at that end and a warning is emitted, so a
missing FR4 residue never becomes a spurious charged terminus.

### `--structure-model propka` (default)

Runs PROPKA3 if `propka3`/`propka` is on `PATH`, parses the
`SUMMARY OF THIS PREDICTION` section, and uses the per-residue and
N/C-terminal pKa values. PROPKA is run **on the trimmed Fv** (the same
residues the SASA-based model sees), not on the whole input file, so that the
two structure models describe the same molecular system and the
structure-vs-sequence delta stays interpretable.

PROPKA reports both the final pKa and the *model* pKa it started from. The
difference between the two is the purely environmental contribution
(desolvation, hydrogen bonds, salt bridges, charge interactions), which is
independent of the free-solution reference. The tool therefore transfers that
environmental contribution onto each bundled pKa set:

```
pKa_effective = pKa_set(model) + [pKa_PROPKA - pKa_model_PROPKA]
```

This keeps the pKa-set sensitivity analysis meaningful in PROPKA mode. (Using
the raw PROPKA value for every set would collapse all five sets onto one
number and report a spurious zero-width uncertainty.)

If PROPKA is missing the run **falls back to the heuristic model** with a
warning (use `heuristic` explicitly to suppress the warning, or install
`propka3`).

### `--structure-model none`

Uses the model pKa values without microenvironment corrections. Useful as a
control to see how much burial/salt-bridge/H-bond effects change the result.

### Structure quality

Before interpreting a structure, the tool reports:

- mean/min/max CA B-factor, together with the **provenance** of the file;
- missing side-chain heavy atoms, especially missing ionizable atoms;
- CA-CA chain breaks (> 4.5 A);
- truncated N/C termini.

The B-factor column means different things in different files, so the tool
establishes the **provenance** of the file first and interprets the numbers
accordingly:

| provenance | B-factor column | low value means | warning |
|---|---|---|---|
| experimental | crystallographic B-factor (A^2) | well ordered | mean > 60 A^2: partly disordered |
| predicted | pLDDT (0-100) | low confidence | mean < 70: low confidence |
| unknown | unknown | - | reported, not interpreted |

Interpreting a crystal structure's low mean B-factor as "low confidence"
(inverting its meaning) would have flagged every well-refined structure as
unreliable, so the provenance is established rather than guessed from the
numeric range.

**How the provenance is established.** The file header is consulted first
(`EXPDTA` / `REMARK 2 RESOLUTION` / `_exptl.method` -> experimental;
AlphaFold/AF3 markers -> predicted). A header is **not** required, though:
AF2, AF3 and ESMFold all write bare `ATOM` records with no `HEADER`/`REMARK`
block at all, so a header-only test would silently disable the pLDDT reading
for exactly the files the tool is most often pointed at. When the header is
inconclusive the file content is inspected instead:

- a standard crystallographic entry always carries `CRYST1` and usually
  `SSBOND` and non-water `HETATM` records - any of those means
  `experimental`;
- a file with **none** of those records whose B-factors all lie inside the
  pLDDT range [0, 100] is treated as `predicted`.

The content call is a documented inference, not a certainty, so the report
says so explicitly and records it in `structure_quality.provenance_basis`
(`header`, `content` or `none`). Use `--provenance predicted` /
`--provenance experimental` to override the inference in either direction.

The completeness of the modelled C-terminus is checked per chain type: a
VH/VHH runs to IMGT 128 while a VL ends at IMGT 127, so a chain truncated
before its own boundary is treated as blocked at that end. A single heavy
chain is reported as a VHH (free C-terminus) with an explicit warning,
because an isolated Fab VH or a VH fusion construct would instead have a
blocked C-terminus, and the two readings differ by up to ~2 pH units of Fv
pI for a near-neutral Fv.

### Salt bridges and pKa uncertainty

Two salt-bridge definitions are reported deliberately:

- **TAP**: 3.2 A, both partners surface-exposed, one partner per residue;
- **pKa microenvironment**: 4.0 A, burial allowed, multiple partners,
  and histidine included as a donor (a protonated imidazole does form
  charge-assisted hydrogen bonds to carboxylates).

They answer different questions and can differ substantially (for example,
a reference Fab has 1 TAP salt bridge but 7 microenvironment pairs). The pKa model uses the
broader definition; the TAP metric remains comparable with published TAP
output. Every run also reports the pI range across all five bundled pKa sets
and a net-charge-vs-pH curve (`--ph-curve` prints it; JSON always contains it).

Hydrogen bonds are detected with explicit donor/acceptor polarity: a
carboxylate O or thiolate S is an acceptor, an ammonium/guanidinium/imidazole
N is a donor, and a contact only counts when the two partners have
complementary roles and the antecedent-donor...acceptor angle exceeds 120
degrees. A D...D or K...K contact is therefore no longer reported as a
hydrogen bond.

## Sequence fallback

If no structure is supplied:

```bash
python3 tools/pi_analysis/pi_analyzer.py --fasta inputs/antibody.fasta
python3 tools/pi_analysis/pi_analyzer.py --vh <VH_SEQ> --vl <VL_SEQ>
python3 tools/pi_analysis/pi_analyzer.py --vh <VHH_SEQ>
```

Sequence mode reports pI, capped pI, net charge and per-chain charge. It does
not report TAP, SASA, charge patches or structure-aware pKa shifts.

### What sequence mode can and cannot tell you

A sequence determines **which** ionizable groups exist and their intrinsic
pKa. It does not determine **where** they sit in the folded protein, and it
is exactly the second part that moves the pI: a buried carboxylate gains
~1.4 pH units, a salt-bridged one loses ~0.8.

| Output | Sequence-only | Why |
|--------|---------------|-----|
| Fv/VHH pI, pKa-set sensitivity range | yes | free-solution titration |
| net charge vs pH, per-chain pI and charge | yes | free-solution titration |
| total-charge asymmetry (`Fv_chml`, `SFvCSP_total`) | yes | charge balance only |
| `Ltot` (CDR length) | yes | exact from IMGT numbering |
| TAP `PSH` / `PPC` / `PNC` / `SFvCSP` | **no** | defined on solvent-exposed charges and 3D distances |
| structure-aware pI (burial / salt bridge / H-bond shifts) | **no** | needs the microenvironment |
| 3D charge patches and the mutation candidates from them | **no** | patches are 3D clusters of exposed residues |
| VH/VL interface, antigen contacts | **no** | needs the 3D complex |

Measured on the 4D5 Fab (same terminus model, all five pKa sets), the
structure-aware pI differs from the sequence pI by **-0.2 to -0.7 pH units**,
with no consistent sign, and the structure mode additionally raised a real
`PSH = red` hydrophobic-patch flag that sequence mode cannot raise at all.

Because an absent flag reads like a clean flag, sequence mode states the
omission explicitly:

- the JSON carries `analysis_type` (`"sequence_only"` / `"structure"`) and
  `unavailable_metrics`, a list of `{key, label, reason}` for every output
  that needs a structure (with the physical reason, not just "n/a");
- the text and Markdown reports open with a **SCOPE OF THIS ANALYSIS** block
  listing them, and close with the warning that *an absent flag is not a
  clean flag*;
- the recommendation list carries the same statement as a `scope` entry, so
  a machine-readable consumer sees it too.

Structure mode reports `analysis_type: "structure"` and an empty
`unavailable_metrics` list.

**Practical route from sequence to structure:** a single-sequence predicted
model (ESMFold / AF2 / AF3) is enough to obtain the missing metrics. The
tool interprets the B-factor column of a predicted model as pLDDT and
reports the low-confidence warning, so a predicted model can be fed in
directly. AF2/AF3/ESMFold files carry no header at all, so the provenance
is inferred from the file content (no crystal records + pLDDT-range
B-factors) and the inference is stated in the report; `--provenance
predicted` forces it.

FASTA handling:

- one record -> VHH;
- two records -> VH/VL; headers containing VH/VL, H/L, heavy/light, kappa or
  lambda are used to assign chains even when VL appears first;
- when the header does **not** state the chain type, the sequence itself
  decides: ANARCI's HMM assignment is used when available, otherwise the
  repo's portable numbering engine, whose conserved FR3 Cys anchor
  separates the two (a heavy V domain carries Cys92, a light one Cys88).
  A two-record file whose light chain comes first is therefore assigned
  correctly even with uninformative headers;
- when neither the header nor the sequence can classify a chain, the file
  order is used (VH then VL).

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
| `--temperature K` | `298.15` | Temperature for the Debye-Huckel screening (kappa ~ 1/sqrt(T)) |
| `--ph-grid START:STOP:STEP` | `4:10:0.5` | pH grid for the charge curve |
| `--ph-curve` | off | Print charge-vs-pH in the text report |
| `--vh-chain`, `--vl-chain` | auto | Explicit structure chain IDs (`--vh-chain` alone = VHH) |
| `--sasa {auto,psa,freesasa}` | `auto` | Surface engine |
| `--psa-binary PATH` | - | JOY `psa` executable |
| `--tap-thresholds {2024,2019,exscientia,legacy}` | `2024` | TAP flag bands |
| `--patch-calibration FILE` | - | JSON file with local patch risk bands (replaces the documented defaults) |
| `--patch-cutoff A` | `8.0` | 3D proximity cutoff of the charge-patch graph |
| `--propka-inactive-threshold PKA` | `20.0` | PROPKA pKa at or above which a group is non-ionizable |
| `--provenance {auto,predicted,experimental}` | `auto` | How to read the B-factor column (auto: header, then file content) |
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

Batch runs also write `combined_pi_results.json`. It carries the same
contract as a single report - `analysis_type`, `unavailable_metrics`,
`structure_model`, `termini_mode`, `chains`, `charge_curve`,
`patch_calibration` and the full `patches` list - so a consumer of the
batch summary can tell which pKa model produced a number and whether a
metric was computed at all. (An omitted field reads like a clean flag.)

The JSON contains:

- `analysis_type`: `"sequence_only"` or `"structure"`, so a consumer can
  branch on what was actually computed;
- `unavailable_metrics`: the structure-dependent outputs that this run could
  not produce (`{key, label, reason}` entries; empty in structure mode);
- `sequence_pI` and `structure_pI` per pKa set (`pI_selected`, `pI_free`,
  `pI_capped`, `pI_blocked`, `delta_vs_sequence`);
- `sensitivity`: selected pI for every bundled pKa set and the resulting range;
- `charge_curve`: net charge sampled across pH;
- `structure_quality`: B-factor statistics **with the file provenance**
  (`provenance`, `provenance_basis`, `bfactor_meaning`), missing side
  chains, chain breaks, truncated termini;
- per-chain pI (`free`/`capped`/`selected`) and net charge;
- `fv`: selected/free/capped pI, pI range, net charge, total and surface
  `Fv_chml`/`SFvCSP` (plus `SFvCSP_basis`, which records whether the `SFvCSP`
  alias holds the surface or the total value), total and surface per-chain
  charges;
- `tap`: Ltot, PSH, PPC, PNC, SFvCSP, flags, risk summary, and both salt-bridge
  counts (`n_salt_bridges_tap` and `n_salt_bridges_pka`);
- `patches`: 3D charge patches with residue rows, local risk basis and
  mutation candidates with sequence-level what-if `delta_pI_estimate` /
  `delta_charge_formulation_estimate` (patch IDs are unique across the pH
  values);
- `patch_calibration`: the risk bands actually applied and where they came from;
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
- The published thresholds are defined on the **surface** descriptors, so the
  tool only compares them against `SFvCSP_surface` / `Fv_chml_surface`. The
  `SFvCSP` alias follows the value that exists in the current mode and
  `SFvCSP_basis` records which one it is (`surface` or `total`), so a
  total-charge product is never silently judged against a surface threshold.
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

The thresholds, the 7.5-9.0 pI window and the `|Fv_chml| <= 4` limit were all
derived on mammalian IgG/Fab. For a VHH they are applied by analogy only: a
camelid VHH has a longer CDR3, a different solvent-exposure profile and no VL
partner, so the flag values are relative guidance rather than calibrated
cut-offs. The tool states this explicitly in the VHH report.

### Charge patches

Patches are built at both the formulation pH and pH 7.4:

- exposed charged residues with rSASA >= 7.5 %;
- clustered by 3D proximity (<= 8 A; the patch graph uses its own 8 A
  neighbour map, while the TAP metrics keep their published 7.5 A map);
- charge computed at the requested pH (salt-bridged groups zeroed);
- score is the pairwise `|q1*q2|/r^2` sum over like-charged pairs;
- `score_tap` is the TAP-style ordered integer-charge score (`PPC`/`PNC`);
- local risk uses bands on the fractional-charge score: positive patches
  0.50/2.00, negative patches 0.60/2.50, with a 1.25x CDR-vicinity penalty;
  `risk_basis` states which rule fired;
- first/last residues of each chain are flagged and excluded from mutation
  candidates because they change the terminal pKa and can form
  pyroglutamate/processing liabilities.

**The band values are a documented prior, not a validated threshold set.** The
local score is an unordered fractional-charge sum, whereas the published TAP
PPC/PNC bands are ordered integer-charge sums over the whole CDR vicinity, so
the two scales are not interchangeable. A site with its own experimental data
can replace the bands with `--patch-calibration FILE`:

```json
{
  "positive": {"elevated": 0.5, "high": 2.0},
  "negative": {"elevated": 0.6, "high": 2.5},
  "cdr_factor": 1.25
}
```

The bands actually used are reported in the JSON output under
`patch_calibration`.

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
  Kabat 37/44/45/47) or FR4. The IMGT anchors (26/39, 55/66, 104/118) are
  **excluded** from the candidate list, not merely flagged: they are the
  conserved paratope framework that delimits each CDR loop.
- VH/VL interface positions are flagged **high risk**: mutation can perturb
  chain pairing. Validate pairing and affinity.
- Charge reversal (K/R -> D/E or D/E -> K/R) is reported separately and should
  be treated as a last resort.

### Scope of the numbers

- The pI reported is the **Fv** (variable domain) pI. cIEF and IEX measure the
  intact molecule, whose charge also contains the constant regions (an IgG1 Fc
  contributes roughly +7 at pH 7.4). The Fv pI is a variable-domain quantity
  and is not directly comparable with an experimental IgG pI.
- The pI window, the TAP thresholds and the charge-asymmetry limit come from
  mammalian IgG/Fab work. For a VHH they are applied by analogy only, and the
  report states this.
- The pI is a **variable-domain** quantity and the termini model is part of
  it: `free` (both termini ionizable), `capped` (C-terminus amidated or
  fused, the Fab default) and `blocked` (both termini blocked) are all
  reported. For a near-neutral Fv the free/capped difference can reach ~2
  pH units, so the construct's actual termini state has to be stated.
- The B-factor semantics are established from the file, not guessed from the
  numeric range. When the file carries no usable marker and no crystal
  records, the interpretation is an **inference** from the content and is
  labelled as such (`provenance_basis: content`); `--provenance` overrides it.
- The salt-bridge/H-bond microenvironment model is a screening heuristic.
  In `propka` mode the environmental contribution is transferred onto every
  bundled pKa set, so the reported pKa-set spread is the spread of the
  *reference* pKa values, not a spread of PROPKA runs.

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
