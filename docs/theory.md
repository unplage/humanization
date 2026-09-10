# Humanization theory: back-mutation design (回复突变理论基础)

## 1. Why back-mutations are needed

CDR grafting transfers only the loops; the donor framework's structural
"pillars" that support CDR conformation (vernier residues), the VH/VL packing
interface, and buried core residues are NOT transplanted. If a human framework
residue at such a position differs from the donor's, the loop conformation —
and therefore affinity — can change. Back-mutations restore the donor residue
at exactly those positions. The art of humanization is to revert ONLY what the
structure needs: every reverted position reduces human-likeness and increases
immunogenicity risk.

## 2. Position classes (Kabat numbering)

### 2.1 Vernier zone — Foote & Winter, JMB 224:487 (1992)

Framework residues that adjust CDR loop conformations ("vernier", like a
calibration scale). They are ~1-2 residues away from the CDRs, mostly buried,
and stabilize loop packing.

| Chain | Positions (Kabat) |
|-------|-------------------|
| VH    | 2, 27-30, 47-49, 67, 69, 71, 73, 78, 93-94, 103 |
| VL    | 2, 4, 35-38, 43-46, 48-49, 58, 62-63, 66-69, 71, 87-88, 98 |

> **Update (2026-08-29):** L37 and L45 have been added to VL vernier zone
> based on large-scale backtest analysis showing these positions are
> consistently missed by the pipeline (FN analysis across 31 pairs).

### 2.2 VH/VL interface — Chothia et al., JMB 186:651 (1985); Vargas-Madrazo & Paz-Garcia, JMB 330:783 (2003)

Core packing residues between VH and VL (mostly buried, hydrophobic); they
maintain the relative orientation of the two domains, i.e. the geometry of the
combining site. Wrong residues here can shift the paratope even if CDRs are
unchanged.

| Chain | Core (Kabat) | Extended |
|-------|--------------|----------|
| VH    | 37, 39, 45, 47, 91, 93, 95, 103 | 32, 34, 50, 58, 104 |
| VL    | 34, 36, 38, 44, 46, 87, 89, 91, 96, 98 | 32, 50, 53 |

(Interface positions inside CDRs — e.g. H95, L89/91/96, H100K — are grafted
with the CDRs and never back-mutated; H103/L98 belong to the human J region.)

### 2.3 Canonical-structure residues — Chothia & Lesk, JMB 196:901 (1987); Al-Lazikani et al., JMB 273:927 (1997)

Residues that determine the canonical conformation class of each CDR loop.
Framework support residues (excluding loop residues themselves):

| Loop | Framework support residues (Kabat) |
|------|------------------------------------|
| H1   | 24, 27, 29, 33, 34, 94 |
| H2   | 47, 48, 49, 57, 58, 71, 78 |
| L1   | 2, 25, 33, 71 |
| L2   | 48, 49, 51, 53, 55, 58, 64 |
| L3   | 34, 36, 46, 49, 90, 94, 97, 98 |

Note: Chothia CDR definitions graft the H1 stem (Kabat 26-35), so several
"canonical" H1 residues are carried with the loop in that scheme.

### 2.4 VHH hallmark — Muyldermans, J Biotechnol 74:277 (2001)

Camelid VHH keep a hydrophilic FR2 patch (no VL partner; the hydrophobic VH3
patch would be exposed). Positions (Kabat): **37 (F/Y), 44 (E/Q), 45 (R),
47 (G/F/S)** vs human VH3 (V/G/L/W). **Never revert these** — they are
essential for single-domain stability and solubility. (In Chothia numbering
these are 42/49/50/52.)

### 2.5 Structure-derived classes (AF3 mode)

- **Buried** (relative SASA < 20%): core packing; reverting a buried
  mismatch is usually safe structurally but low immunogenic value.
  Buriedness is computed by FreeSASA (accurate SASA calculation) when
  available, with a contact-counting heuristic fallback.
- **CDR contact** (heavy atom < 4.5 A to any CDR atom): directly stabilizes
  loop conformation → high reversion priority.
- **Antigen contact** (heavy atom < 4.5 A to antigen, complex model):
  can affect binding directly; keep donor unless proven neutral.
- **Disulfide Cys** (e.g. VHH CDR3 disulfide partners in FR): keep donor.

## 3. Decision logic

Per framework position where donor != human germline:

1. If VHH hallmark or disulfide-Cys → **KEEP_DONOR** (never revert).
2. Interface-core (and buried, or with other structural features) → **T1**.
3. Vernier + buried, or vernier + canonical → **T1**; vernier alone → **T2**.
4. Canonical, CDR-contact, antigen-contact, buried, interface-extended → **T2**.
5. No structural role → **T3** (optional, immunogenicity).
6. Fully human-identical → no candidate.

**Step 3 structure adjustment** (only when AF3/experimental hints exist):
CDR/antigen contact keeps or raises the tier; a *buried* position is retained
(packing/core); a *surface-exposed* position with no CDR/antigen contact and a
reliable pLDDT (≥ 50) is demoted T1/T2 → T3. `interface_core` is exempt from
demotion. This is the symmetric form of the previous (buried-only) demotion and
avoids dropping structurally required back-mutations.

Recommended production variant: **V2 (T1+T2)**. V1 when framework identity is
high and risk tolerance low; V3 only when additional humanization is demanded
(e.g. repeated anti-drug-antibody responses), accepting affinity risk.

## 4. Germline selection (multi-strategy)

The pipeline offers 9 selection strategies (`--germline-strategy`), all
computed on the same per-position comparison (FR / CDR / CVI identity):

| Strategy | Rule | Rationale |
|----------|------|-----------|
| `fr_best` | highest FR identity | framework stability priority |
| `cdr_best` | highest CDR1+2 identity | CDR conservation, fewer reversions |
| `composite` | 0.7·FR + 0.3·CDR | balanced |
| `cvi_best` | highest CVI identity | BI 2024: CVI predicts expression + affinity retention |
| `min_backmutations` | fewest T1+T2 (V2) | minimal change (actual back-mutation analysis, not a raw diff count) |
| `adimab_frequency` | Adimab-recommended + usage frequency | therapeutic-antibody track record |
| `pioneer_frequency` | Pioneer-library gene + frequency | 600+ clinical-stage antibodies |
| `composite_3axis` | 0.5·CVI + 0.3·frequency + 0.2·FR | combined optimum |
| `auto` (default) | VH=adimab_frequency, VL=current | backtest-optimized routing (see docs/scoring.md) |

Selection notes:
- **FR threshold**: default 60% (graft viability). When no human gene reaches
  it (highly non-human sequences), the threshold auto-degrades
  (50% → 40% → …) with a warning, and the highest-FR feasible gene is chosen.
- **J gene**: highest FR4 identity, matching CDR3/J junction length.
- The legacy single-strategy "hybrid" (top-30% FR, then max CDR1+2) is kept
  as the `current` strategy for backward compatibility.

## 4.5 CDR3 loop residues at Kabat H93/H94 (严格 Kabat 语义)

严格 Kabat 编号中 CDR3 = H95-102，而 CDR3 环的前两个残基占据 H93/H94
（形式上属 FR3；IMGT 105-106）。**移植时整个环（含 93/94）必须来自供体**；
编号表中 H93/H94 仍标注为 FR3。回复突变分析将 H93/H94 视为环残基排除
（不是框架候选），见 docs/backtest_report.md 的 P0 修复记录。

## 5. Validation case

Mouse anti-HER2 4D5 (PDB 1N8Z/1FVE family) → trastuzumab (Carter et al.,
PNAS 89:4285, 1992). The pipeline reproduces the structurally critical
reversions (VH 48/67/69/71/73; VL 4/66/87) in T1/T2, with the same germline
families used historically. Run `tests/test_pipeline.py` for regression checks.

## 6. Numbering edge cases (known & handled)

The anchor-based Kabat engine handles non-canonical inputs with explicit
warnings; ANARCI cross-validation is recommended before synthesis:

- **Cys22 anchor shift** (e.g. an extra N-terminal residue moves Cys to
  position 23): FR1 is capped to the standard 30 residues and CDR1 absorbs
  the extra residue, so H31 is never double-assigned and no residue is
  silently dropped on grafting.
- **V-region-only germline sequences** (no J region, e.g. bundled JSON V
  genes): numbered without a Trp103/Phe98 anchor; the tail is treated as CDR3
  and FR4 stays empty. J genes are loaded from the curated position map.
- **Germline gap positions** (e.g. VH3-family H49 empty for VHH): the donor
  residue is kept on grafting so the domain stays complete; dropping it would
  shift every downstream position and corrupt CDR2/CDR3.

## 7. FR insertions / deletions (donor framework length differences)

`fr_indel.detect_fr_indels` finds donor FR insertions/deletions per region and
returns, for each, a ranked list of candidate insertion positions with
confidence scores. In interactive mode the user confirms the insertion point
(`select_insertion_interactive`); the choice flows through
`update_indel_selection` into the correspondence and is honoured everywhere.

The donor<->germline correspondence (`fr_indel.build_fr_correspondence`) is the
single source of truth for the alignment: the confirmed insertion positions are
removed from the ordered donor FR list and the remainder is paired 1:1 with the
ordered germline FR list. Example::

    donor    H5 H6(E, insertion) H6A H7
    germline     H5 H6             H7
    -> germline H6 maps to donor H6A; H6 is the insertion

Consequences (all fixed by this layer):

- grafting copies the germline residue to the *aligned* donor label, so the
  inserted residue never overwrites it and no residue is duplicated;
- back-mutation candidates are walked in germline numbering and resolved to the
  aligned donor residue, so downstream positions are not shifted;
- the inserted residue itself is scored as an indel candidate (not a
  substitution), and is included in V2; V0/V1 exclude it (pure graft);
- `human_likeness_percent` compares aligned columns, so insertions are not
  counted as framework mismatches.

Donor *deletions* (germline residues absent from the donor) are also detected;
the graft fills them with the human residue and they are reported.

## 8. AF3 structure validation of variants (CDR-RMSD)

When AF3 is available (`--af3-mode ... --af3-rmsd`), the pipeline validates the
designed variants structurally:

1. Predict the donor Fv (reference) — already done for Step 3 hints.
2. AF3-predict each variant. For a **Fab**, the VH and VL variants sharing a
   suffix (H_V2 / L_V2) are predicted together as one variant Fv; if only one
   side has that variant, the donor partner chain is used. VHH/single-chain
   inputs are predicted as monomers. The antigen is added when provided.
3. Match each output chain to its input sequence (`match_pdb_chain`), map PDB
   residue numbers to Kabat positions, and superpose the variant on the donor
   **using high-confidence framework CA atoms** (Kabsch,
   `structure.kabsch_superpose`; only FR residues with pLDDT ≥ 70 when AF3
   confidence is available).
4. Report CDR CA-RMSD per loop (CDR1/CDR2/CDR3) and overall, plus the framework
   RMSD and an inter-residue heavy-atom **clash count** (< 2.0 Å,
   `structure.count_clashes`), in the markdown report and
   `humanization_result.json` (`structure_validation`).

Superposing on the framework is essential: independently predicted models have
arbitrary global orientation, so a raw RMSD would be meaningless. After the
framework fit, the CDR-RMSD isolates the loop-conformation change induced by
the humanized framework — the quantity that correlates with affinity loss.
This is the lightweight (no new dependency) score from
`docs/improvement_roadmap.md` §2.2 D.

Notes:

- Step 3 tiers need **only the donor** structure. The CDR-RMSD step is opt-in
  (`--af3-rmsd`) and never runs otherwise; a donor-only AF3 structure is enough
  for Step 3.
- The donor structure is always the reference. The variant structure **does not
  need to contain a donor chain** — only the chain being evaluated is used for
  the framework fit and CDR measurement. In the automatic path the variant is
  modelled with the donor partner chain purely to isolate the single-chain
  effect.
- If the variant structures already exist, skip AF3 re-prediction and use the
  standalone command (no pipeline run):

  ```bash
  python3 scripts/humanize/cli.py rmsd --input ab.fasta --donor donor.pdb \
    --variants v0.pdb v1.pdb v2.pdb --chain H [--scheme kabat] [--out rmsd.json]
  ```

## 9. Variant enumeration: structure-guided panel

Besides the fixed V0–V3 ladder, `variants.structure_guided_panel` emits a
graded **V_opt panel** (opt-in via `--design-panel`, `--design-panel-steps N`).
Starting from the V2 set it greedily adds the candidate that recovers the most
currently-uncovered framework→CDR contacts (donor `cdr_partners`), breaking
ties by composite score and falling back to composite ranking without
structure. Each step yields one additional back-mutation, giving a screening
gradient of increasing contact recovery / mutation count for experimental
titration.
