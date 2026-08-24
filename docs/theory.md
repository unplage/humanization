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
5. No structural role: exposed mismatch → **T3** (optional, immunogenicity);
   buried without structure data → **T3**.
6. Fully human-identical → no candidate.

Recommended production variant: **V2 (T1+T2)**. V1 when framework identity is
high and risk tolerance low; V3 only when additional humanization is demanded
(e.g. repeated anti-drug-antibody responses), accepting affinity risk.

## 4. Germline selection (why "hybrid")

- Primary: **framework identity** (FR1+FR2+FR3) — a poor framework match means
  many back-mutations.
- Tie-break within top-FR genes: **highest CDR1+2 identity** — preserves more
  CDR residues as-is and minimizes structural perturbation (Olimpieri et al.,
  Bioinformatics 31:434, 2015).
- J gene: highest FR4 identity, matching CDR3/J junction length.

## 5. Validation case

Mouse anti-HER2 4D5 (PDB 1N8Z/1FVE family) → trastuzumab (Carter et al.,
PNAS 89:4285, 1992). The pipeline reproduces the structurally critical
reversions (VH 48/67/69/71/73; VL 4/66/87) in T1/T2, with the same germline
families used historically. Run `tests/test_pipeline.py` for regression checks.
