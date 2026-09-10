# Back-mutation scoring rationale (评分依据)

## 0. Germline selection strategy (auto routing)

After large-scale backtest optimization (31 mouse→humanized pairs, 9 strategies):

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

The `chemical` term is capped at 1.0 on the positive side, but **negative
penalties pass through unclamped** — in particular `introduces_nglycan`
(-0.80) must lower the composite whenever a reversion would create a new
N-glycan motif (regression-tested in tests/test_pipeline.py). Note the
practical range of the composite is therefore ~9-92 rather than the full
0-100: `structural` tops out at 1.0, `immunogenicity` (benefit) spans
0.3-0.725, and `chemical` is bounded by ±0.80 with the current weights.

### 1.1 Structural score (0-1) — weights in config.py

| Feature | Weight | Source |
|---------|--------|--------|
| interface_core | 1.00 | Chothia 1985 packing set |
| vernier | 0.85 | Foote & Winter 1992 |
| canonical | 0.80 | Chothia & Lesk 1987 |
| interface_extended | 0.65 | Vargas-Madrazo 2003 |
| buried | 0.70 | AF3 structure (relSASA < 0.20) |
| cdr_contact | 0.85 | AF3 structure (heavy atom < 4.5 A) |
| antigen_contact | 0.60 | AF3 complex structure |
| vhh_hallmark | 1.00 | keep-donor, never revert |
| disulfide_cys | 1.00 | keep-donor, never revert |

The structural score combines the present features with a **noisy-OR**
(`1 - Π(1 - w_i)`), so independent pieces of evidence reinforce each other
instead of collapsing to the single strongest feature. A confirmed buried
position raises it to at least 0.70; a confirmed exposed position scales it by
0.85. The result is then multiplied by a **continuous pLDDT factor**
(`0.2 + 0.8·clip((pLDDT-50)/40, 0, 1)`, unknown pLDDT = 1.0), so low-confidence
structural evidence is down-weighted smoothly rather than by a hard cutoff.

### 1.2 Immunogenicity benefit (0-1)

```
benefit = 0.3 + 0.5 * exposure * (1 - conservation)
```

- exposure: 0.85 exposed / 0.15 buried (AF3); 0.5 default without structure.
- conservation: fraction of the **entire human germline repertoire** carrying
  the DONOR residue at that position, **weighted by each germline's
  therapeutic-antibody usage frequency** (`germline_frequency.get_frequency`).
  The reference panel is the full human V set (not the top-20 homologs), so a
  residue shared only with rare germlines counts as less "human". A rare donor
  residue (low conservation) = more foreign = higher benefit of reverting.
- Positions with no structural feature are capped at 0.50 (surface
  humanization is still worth something, but low priority; capped to
  preserve conservation-based ranking without flattening the gradient).

### 1.3 Chemical / developability

The chemical term is a **symmetric liability delta** computed on the sequence
pattern of the donor vs the human (graft) state at the candidate position:

```
chem = liability(human residue at pos) - liability(donor residue at pos)
```

- `chem > 0` → reverting to the donor **removes** a liability the graft would
  carry (reward);
- `chem < 0` → reverting **introduces / retains** a donor liability (penalty).

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
| `buried = True` | **keep** T1/T2 (never demote) | core/packing residues can still support the fold; a buried reversion is immunologically silent |
| `buried = False` AND no contact AND pLDDT ≥ 50 | demote T1/T2 → **T3** | surface-exposed, non-contacting residue: no structural reason to stay donor, and it carries real immunogenicity risk |
| `plddt < 50` | no demotion | low-confidence structure is not trusted |
| `interface_core` | exempt from demotion | the contact detector can miss packing partners |

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
**J gene** map (`backmut.py`, requires `j_gene`) and the resulting `T_FR4`
candidates are included in the V2 variant.

## 4. Known limitations

- Without AF3, buriedness/contact hints are absent: tiering relies on the
  literature position sets only (still the industry baseline).
- The Kabat insertion lettering of the portable engine is approximate for
  exotic loop lengths; run ANARCI mode on the server before synthesis
  (see validation.md).
- Germline identity uses the bundled human IMGT-derived set; the NCBI
  IgBLAST FASTA set (identical origin) is used when downloaded.
- When no human germline reaches the 60% FR-identity threshold, selection
  automatically degrades (FR ≥ 50% → 40% → …) with a warning; the chosen
  gene is the highest-FR feasible one, but the many back-mutations expected
  (FR < 60%) mean affinity retention is at risk (BI 2024).
