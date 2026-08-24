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
- conservation: fraction of the top-10 human germlines carrying the DONOR
  residue at that position. A rare donor residue (low conservation) = more
  foreign = higher benefit of reverting.
- Positions with no structural feature are capped at 0.45 (surface
  humanization is still worth something, but low priority).

### 1.3 Chemical / developability (0-1)

Bonus for reverting away from a donor residue that creates a risk motif in a
3-residue window:

- N-glycan (NxS/T): +0.8
- deamidation (NG / NS): +0.5
- isomerization (DG): +0.4
- oxidation-prone (M/W): +0.3

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

## 4. Known limitations

- Without AF3, buriedness/contact hints are absent: tiering relies on the
  literature position sets only (still the industry baseline).
- The Kabat insertion lettering of the portable engine is approximate for
  exotic loop lengths; run ANARCI mode on the server before synthesis
  (see validation.md).
- Germline identity uses the bundled human IMGT-derived set; the NCBI
  IgBLAST FASTA set (identical origin) is used when downloaded.
