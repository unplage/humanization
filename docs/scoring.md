# Back-mutation scoring rationale (评分依据)

## 0. Germline selection strategy (auto routing)

After large-scale backtest optimization (mouse→humanized pairs, 9 strategies;
driver: `tests/backtest_large.py`, which requires the non-bundled
`data/benchmarks/humanized_pairs.csv`):

| Strategy | Precision | Recall | F1 |
|----------|-----------|--------|-----|
| **adimab_frequency** | **0.448** | **0.326** | **0.377** |
| current | 0.363 | 0.311 | 0.335 |
| fr_best | 0.335 | 0.325 | 0.330 |

**Default auto routing:** VH → `adimab_frequency`, VL → `current`

## 1. Score components

Every framework candidate (donor != human germline) receives three components
and one blended composite:

```
composite = 100 * (0.55 * structural + 0.30 * immunogenicity + 0.15 * chemical)
```

The composite is **bounded to [0, 100]**. A positive `chemical` term is
capped at 1.0 and added, as in the formula above. A negative term (a
reversion that re-introduces a donor liability, e.g. N-glycan) instead
rescales the non-chemical base by `(1 - |chemical|/1.10)`: a stronger liability
always ranks lower, but the score can no longer become negative. The
practical range with the current weights is ~0-94 (`structural` tops out at
1.0, `immunogenicity` spans 0.3-0.8, `chemical` is bounded by ±1.10 after
exposure scaling). The signed chemical score itself remains in the report and
CSV, so the penalty is never hidden.

### 1.1 Structural score (0-1) — weights in config.py

| Feature | Weight | Source |
|---------|--------|--------|
| interface_core | 1.00 | Chothia 1985 packing set |
| vernier | 0.85 | Foote & Winter 1992 |
| canonical | 0.80 | Chothia & Lesk 1987 |
| interface_extended | 0.65 | Vargas-Madrazo 2003 |
| buried | 0.70 | AF3 structure (relSASA < 0.20) |
| loop_shell | 0.35 | measured second-shell loop proximity (D3, side chain 4.5-8 A; ranking only) |
| cdr_contact_sc | 0.85 | AF3 structure: **side-chain** heavy atom < 4.5 A to CDR |
| cdr_contact_bb | 0.30 | AF3 structure: backbone-only contact (β-sandwich geometry) |
| antigen_contact | 0.60 | AF3 complex structure |
| vhvl_contact | 0.70 | **measured** VH/VL packing contact (Fab): side chain within 4.5 Å of the partner chain; functional evidence for interface residues that touch no CDR |
| vhh_hallmark | 1.00 | keep-donor, never revert |
| disulfide_cys | 1.00 | keep-donor, never revert |

The structural score combines the present features with a **noisy-OR**
(`1 - Π(1 - w_i)`), so independent pieces of evidence reinforce each other
instead of collapsing to the single strongest feature. It is then made
**structure-adaptive** (per antibody) in four ways:

1. **Fractional contacts.** A CDR contact is weighted by how it is formed:
   side-chain-mediated (`cdr_contact_sc` = 0.85) is functional — a side-chain
   swap can break it — while a backbone-only contact (`cdr_contact_bb` = 0.30)
   is fixed β-sandwich geometry that no side-chain substitution changes.
2. **Exposure modulation** (`config.STRUCTURAL_EXPOSURE_FACTOR`, applied only
   to *non-functional* positions): deeply buried positions are boosted (1.15)
   and exposed ones down-weighted (0.70), continuously in relSASA when
   available, else from the binary buried hint. The same literature position
   therefore scores differently across antibodies whose local packing differs.
   Positions with a functional contact are left untouched.
3. **pLDDT factor**: continuous confidence weighting
   (`0.2 + 0.8·clip((pLDDT-50)/40, 0, 1)`, unknown pLDDT = 1.0), so
   low-confidence structural evidence is down-weighted smoothly.
4. **Per-antibody loop-class correction** (`loopclass.py`, see
   `docs/theory.md` §2.6): the feature *set* is corrected before scoring, not
   just reweighted. A `canonical` claim whose predicted loop is > 8 A away is
   removed; a `vernier` claim > 12 A from every loop is removed; a buried
   position whose **side chain** lies 4.5-8 A from a loop, absent from every
   other structural feature, is added as `loop_shell` (0.35, ranking only - a
   side-chain contact is already `cdr_contact_sc`, and an 8 A proximity is a
   neighbourhood, not a packing contact). When the canonical claim survives,
   its weight is scaled by the modelled confidence of the loop it supports
   (mean pLDDT: 1.0 at >= 80, 0.60 at <= 50).

Without structure all four reduce to identity, so the portable path is
unchanged. The correction is skipped for a position whose model is missing a
whole loop (or part of it), and for pLDDT < 50 positions; the Kabat-only
membership table is not applied in IMGT mode. In a multi-model consensus the
measured shell relation must hold in >= 60% of the models.

### 1.2 Immunogenicity benefit (0-1)

```
benefit = 0.3 + 0.5 * exposure * (1 - conservation)            # default proxy
benefit = 0.3 + 0.5 * (0.4 * exposure*(1-conservation) + 0.6 * epitope)   # if MHC-II supplied
```

- exposure: **measured relSASA clamped to 0.15-0.85** (AF3/Fab structure, the
  same continuous value used by the structural and chemical terms); the binary
  buried hint (0.85 exposed / 0.15 buried) is kept only when relSASA is
  absent; 0.5 without any structure.
- conservation: fraction of the **entire human germline repertoire** carrying
  the DONOR residue at that position, **weighted by each germline's
  therapeutic-antibody usage frequency** (`germline_frequency.get_frequency`).
  The reference panel is the full human V set (not the top-20 homologs), so a
  residue shared only with rare germlines counts as less "human". A rare donor
  residue (low conservation) = more foreign = higher benefit of reverting.
- Positions with no structural feature are capped at 0.50 (surface
  humanization is still worth something, but low priority; capped to
  preserve conservation-based ranking without flattening the gradient).
- **Optional MHC-II override** (`--immunogenicity-json` or `--netmhciipan`):
  when a per-position epitope score is available it dominates the proxy
  (weight 0.6). A donor residue inside a strong predicted T-cell epitope gets a
  high benefit for reverting to the human germline (lower ADA risk). The score
  is stored per candidate (`immunogenicity_score`) and shown as an `immuno`
  column in the report. Falls back silently to the proxy when the tool/map is
  absent.

### 1.3 Chemical / developability

The chemical term is a **symmetric liability delta** computed on the sequence
pattern of the donor vs the human (graft) state at the candidate position:

```
chem = (liability(human residue at pos) - liability(donor residue at pos))
       * CHEMICAL_EXPOSURE_WEIGHT[exposure_class(pos)]
```

- `chem > 0` → reverting to the donor **removes** a liability the graft would
  carry (reward);
- `chem < 0` → reverting **introduces / retains** a donor liability (penalty).
- **Exposure scaling** (`config.CHEMICAL_EXPOSURE_WEIGHT`): developability
  liabilities are surface chemistry, so the delta magnitude is scaled by the
  substituted position's solvent exposure — exposed 1.00, intermediate 0.60,
  buried 0.25. A buried NG/DG/M is far less reactive than an exposed one, so a
  buried liability no longer scores as if it were surface-exposed. Without
  structure the class is `unknown` (weight 1.0), keeping the sequence-only
  value unchanged.

Liability motifs and weights (`config.LIABILITY_MOTIFS`):

| Motif | Weight | Notes |
|-------|--------|-------|
| N-glycan (NxS/T) | 0.80 | N-X-S/T (X≠P); glycosylation site |
| deamidation (NG) | 0.55 | Asn-Gly hotspot (Lu 2018) |
| deamidation (NS) | 0.50 | Asn-Ser |
| deamidation (NH) | 0.40 | Asn-His |
| deamidation (ND) | 0.35 | Asn-Asp |
| isomerization (DG) | 0.50 | Asp-Gly (FRIDA 2024) |
| isomerization (DS/DT) | 0.40 | Asp-Ser/Thr |
| isomerization (DH) | 0.35 | Asp-His |
| acid hydrolysis (DD) | 0.55 | Asp-Asp; dual risk |
| acid hydrolysis (D-X) | 0.45 | D + small residue (A/V/L/I/P) |
| oxidation (M/W) | 0.30 | Met/Trp oxidation |
| base hydrolysis (K-X) | 0.25 | Lys + D/E (rare) |
| met-lyscleavage (MK) | 0.25 | Metalloprotease cleavage |

This supersedes the old asymmetric `removes_*` / `introduces_nglycan` keys,
which only handled N-glycan and had inconsistent signs. Chemical weights are
from the literature (Lu et al. 2018 deamidation; FRIDA 2024 isomerization;
Booster et al. 2022 stability) plus expert consensus.

**Overlap merge** (`backmut._liability_score`): overlapping motif matches
are merged into one interval scored at the highest weight of its members —
a region matching both N-glycan (0.80) and deamidation NG (0.55) counts
once at 0.80 (one physical region, one risk). Adjacent but non-overlapping
motifs stay separate, and isolated motifs keep their original magnitudes.

## 2. Tier assignment (final recommendation)

| Tier | Rule | Meaning |
|------|------|---------|
| KEEP_DONOR | VHH hallmark OR framework Cys (disulfide partner) | must not revert |
| T1 | interface_core (and not clearly exposed); vernier ∩ buried; vernier ∩ canonical | must revert |
| T2 | vernier; canonical; cdr_contact; antigen_contact; buried; interface_extended | strongly recommended |
| T3 | no structural role; exposed or unknown | optional, immunogenicity only |
| (no entry) | donor == human | not a candidate |

## 3. Interpretation rules (production)

- **V2 (T1+T2) is the default production candidate.** Historical cases
  (trastuzumab, palivizumab, bevacizumab lineage) kept reverions of exactly
  this class.
- A T3 reversion is only justified by clinical immunogenicity data or a
  regulatory ask for maximal human content; it carries affinity risk.
- With AF3 data: treat `cdr_contact` + `buried` as nearly mandatory;
  `antigen_contact` reversions are discouraged unless the contact is
  backbone-only and side-chain buried.
- The composite score is a ranking aid, not a hard filter: positions in
  T1/T2 with composite >= 60 are "high confidence"; 40-60 "moderate";
  < 40 "borderline" (typically T3).

## 3.5 Structure-driven tier adjustment (Step 3, AF3/结构数据)

Once per-position structure hints are available, tiers are adjusted as follows
(implemented in `backmut.py::analyze_backmutations`):

| Evidence | Effect | Rationale |
|----------|--------|-----------|
| `cdr_contact` / `antigen_contact` / `cdr_partners` | keep/raise to T1/T2 | direct functional support of the paratope |
| `partner_contact_sc` (measured VH/VL contact) | keep/raise to T1/T2 | interface packing is functional selection; the partner chain is in the donor structure and is now included in the contact pool instead of discarded |
| `buried = True`, relSASA < 0.10 (deep core) | keep T2/T1 (no automatic demotion) | deep packing is structurally constrained; even an isosteric swap can matter |
| `buried = True`, relSASA 0.10-0.20 AND disruptive swap (e.g. Lys→Ala) | keep T2 | shallow burial alone is weak, but a large volume/charge/H-bond change still creates a cavity / breaks packing |
| `buried = True`, relSASA 0.10-0.20 AND near-isosteric swap | demote T2 → **T3** | shallow burial is not functional evidence; the human residue is preferable when the swap cannot disrupt packing |
| `buried = False` AND no contact AND pLDDT ≥ 50 | demote T1/T2 → **T3** | surface-exposed, non-contacting residue: no structural reason to stay donor, and it carries real immunogenicity risk |
| near-isosteric swap at a deep-core site (relSASA < 0.10) with ≥3 side-chain framework neighbours and ≥20 A³ volume strain | treated as **disruptive** | cavity creation in a tightly packed core is not visible to the sequence-only 40 A³ rule |
| `plddt < 50` | no demotion | low-confidence structure is not trusted |
| `interface_core` | exempt from T3 demotion | literature-only label; the measured `partner_contact_sc` now supplies the same evidence when the structure permits |

FR4 structural reversions now require **joint evidence** (`FR4_REVERSION_THRESHOLD` 0.55): a side-chain CDR3 contact (0.45) alone is insufficient; it must be corroborated by burial, antigen contact or interface packing (0.45 + 0.15 = 0.60 passes). FR4 is J-encoded and human by construction, so a single weak signal must not override the domain rule.

Literature-only claims (vernier/canonical/interface tables) with no direct measured corroboration in the structure are discounted by `LITERATURE_ONLY_FACTOR = 0.6` in the composite; tiers are unchanged but ranking no longer over-trusts the static 1980s tables.

This is symmetric: the earlier rule demoted *buried* positions but kept
*exposed* literature positions, which could silently drop structurally
required back-mutations (e.g. the trastuzumab H78 reversion).

### 3.5.1 Gold-standard empirical demotion (金标准实测降级)

内置"实测无效应位点"表（`config.py::EMPIRICAL_NO_EFFECT`），来自回测金标准：

| 位点 | 证据 | 无结构时的处理 |
|------|------|----------------|
| VL L87 | 曲妥珠单抗保留 germline Y87 且 KD ~0.1 nM 完好（docs/backtest_report.md） | T1/T2 → **T3**（score 封顶 40） |

- 仅当无结构证据时生效：若 AF3 判定该位点埋藏或接触 CDR（`buried=True` /
  `cdr_contact=True`），降级被覆盖，维持原分级
- 新增此类条目需附回测/实验证据

### 3.5.2 FR4 structural reversion (T_FR4)

FR4 is human (J gene) by construction, but when a J-region residue contacts
CDR3 or the antigen it may need to stay donor. This is evaluated against the
**J gene** map (`backmut.py`, requires `j_gene` **and** structure data) and the
resulting `T_FR4` candidates are included in the V2 variant.

Trigger: a donor/J mismatch at a J-region position whose structure-derived
score meets `FR4_REVERSION_THRESHOLD` (0.4) — `cdr3_contact` 0.45,
`antigen_contact` 0.35, `buried` 0.15, `interface_core` 0.05. An exposed
mismatch with no CDR3/antigen contact is kept human (reverting it would only
add a non-human surface residue).

J gene choice itself is already identity-maximizing (highest FR4 identity, then
frequency), so mismatches are minimized before reversion is considered. Every
report always lists the FR4/J mismatches and their disposition ("FR4 / J-region
Analysis"), so a no-reversion outcome is explicit rather than silent.

## 4. Known limitations

- Without AF3, buriedness/contact hints are absent: tiering relies on the
  literature position sets only (still the industry baseline).
- Structure-adaptive weighting modulates the literature weights by measured
  exposure and contact type; it is **not** an energy calculation (no
  Rosetta/FoldX ΔΔG). Exposure is applied at the substituted residue, which is
  an approximation for liability motifs spanning the adjacent residue.
- Removing an N-glycan (or introducing one) is still reported at full motif
  weight through the liability scan even when buried; only the composite's
  chemical term is exposure-scaled.
- The Kabat insertion lettering of the portable engine is approximate for
  exotic loop lengths; run ANARCI mode on the server before synthesis
  (see validation.md).
- Germline identity uses the bundled human IMGT-derived set; the NCBI
  IgBLAST FASTA set (identical origin) is used when downloaded.
- When no human germline reaches the 60% FR-identity threshold, selection
  automatically degrades (FR ≥ 50% → 40% → …) with a warning; the chosen
  gene is the highest-FR feasible one, but the many back-mutations expected
  (FR < 60%) mean affinity retention is at risk (BI 2024).
