# Back-mutation scoring rationale (评分依据)

## 1. Score components

Every framework candidate (donor != human germline) receives three components
and one blended composite (0-100):

```
composite = 100 * (0.55 * structural + 0.30 * immunogenicity + 0.15 * chemical)
```

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

The structural score is the max over present features; a confirmed buried
position raises it, a confirmed exposed position caps it at 85%.

### 1.2 Immunogenicity benefit (0-1)

```
benefit = 0.3 + 0.5 * exposure * (1 - conservation)
```

- exposure: 0.85 exposed / 0.15 buried (AF3); 0.5 default without structure.
- conservation: fraction of the top-20 human germlines (by FR identity)
  carrying the DONOR residue at that position. A rare donor residue (low
  conservation) = more foreign = higher benefit of reverting. The panel is
  the top-N homologous germlines, not the per-strategy winners, so the
  estimate is unbiased across the human repertoire.
- Positions with no structural feature are capped at 0.50 (surface
  humanization is still worth something, but low priority; capped to
  preserve conservation-based ranking without flattening the gradient).

### 1.3 Chemical / developability (0-1)

Bonus for reverting away from a donor residue that creates a risk motif in a
3-residue window (anchored position, not full-window scan):

| Motif | Weight | Notes |
|-------|--------|-------|
| N-glycan (NxS/T) | 0.80 | N-X-S/T (X≠P); glycosylation site |
| deamidation (NG) | 0.55 | Asn-Gly hotspot (Lu 2018) |
| deamidation (NS) | 0.50 | Asn-Ser |
| deamidation (NH) | 0.40 | Asn-His |
| deamidation (ND) | 0.35 | Asn-Asp |
| isomerization (DG) | 0.50 | Asp-Gly (FRIDA 2024) |
| isomerization (DS) | 0.40 | Asp-Ser |
| isomerization (DT) | 0.40 | Asp-Thr |
| isomerization (DH) | 0.35 | Asp-His |
| acid hydrolysis (D-X) | 0.45 | D + small residue (A/V/L/I/P) |
| acid hydrolysis (DD) | 0.55 | Asp-Asp; dual risk (isomerization + hydrolysis) |
| oxidation (M/W/C) | 0.30 | Met/Trp/Cys oxidation |
| base hydrolysis (K-X) | 0.25 | Lys + D/E (rare; alkaline conditions) |
| met-lyscleavage (MK) | 0.25 | Metalloprotease cleavage |
| **introduces_nglycan** | **-0.80** | reversion would CREATE an N-glycan (penalty) |

Chemical weights are from the literature: Lu et al. 2018 (deamidation), FRIDA
2024 (isomerization), Boosted et al. 2022 (general stability), and expert
consensus. Positions where the reversion would break an existing risk motif
get the corresponding positive weight; positions where the reversion would
INTRODUCE a glycan get a -0.80 penalty.

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

## 3.5 Gold-standard empirical demotion (金标准实测降级)

内置"实测无效应位点"表（`config.py::EMPIRICAL_NO_EFFECT`），来自回测金标准：

| 位点 | 证据 | 无结构时的处理 |
|------|------|----------------|
| VL L87 | 曲妥珠单抗保留 germline Y87 且 KD ~0.1 nM 完好（docs/backtest_report.md） | T1/T2 → **T3**（score 封顶 40） |

- 仅当无结构证据时生效：若 AF3 判定该位点埋藏或接触 CDR（`buried=True` /
  `cdr_contact=True`），降级被覆盖，维持原分级
- 结构模式（AF3）开启时，该表自动失效（结构数据优先于历史先验）
- 新增此类条目需附回测/实验证据

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
