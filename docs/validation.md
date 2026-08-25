# Validation plan (验证方案)

How to gain confidence in a humanization result before committing to synthesis.

## 1. Numbering cross-validation (server, ANARCI)

The portable engine uses anchor-based Kabat numbering. ANARCI (exact, HMMER)
must confirm the region map for the input and the final variants:

```python
from anarci import anarci
from humanize.numbering import number_heavy

seq = "EVQLQQSGPELVKPGASVKMSCKASGYTFTDYYMYWVKQSHGKSLEWIGYINP..."
# our numbering
our = number_heavy(seq)
# ANARCI (IMGT + chothia schemes)
result = anarci([("x", seq)], scheme="chothia", output=False)
# compare CDR boundaries and framework position mapping
```

Checks:
- CDR1/2/3 start/end residues identical for both chains.
- No framework position is classified differently (vernier/interface lookups
  must hit the same residues).
- Insertion-letter placements in long loops (CDR-H1 >5, CDR-H2 >16,
  CDR-H3 >8, CDR-L1 >11) must match.
- Known engine edge cases to double-check in ANARCI:
  - Cys22 anchor shifted (N-terminal extra residue) → FR1 capped to 30;
    confirm the CDR1 boundary and that no residue was dropped.
  - V-region-only germline / VHH with a germline gap (e.g. VH3 H49):
    confirm the donor residue is kept and the grafted sequence length
    equals the donor's.

Expected agreement: >99% of positions for standard V domains; any divergence
must be resolved before synthesis (adjust numbering or ANARCI-mode config).

## 2. Germline selection sanity

- The chosen V gene should also be among IgBLAST `-humanize` or BioPhi
  top hits for the same chain (identical IMGT data; expect same family).
- Framework identity values in `humanization_result.json` should be
  reproducible across runs (deterministic pipeline).

## 3. Back-mutation benchmark (regression)

`tests/test_pipeline.py` includes the 4D5 → trastuzumab case:
- VH germline should be an IGHV1-family gene (FR id >= 0.7).
- T1/T2 sets must include the historically critical reversions
  (H48/H67/H69/H71/H73 for VH; L4/L66/L87 for VL).
- VHH (1BZQ): hallmark positions must be KEEP_DONOR; CDR3 (incl. any Cys)
  must be grafted intact.

## 4. Structure-based acceptance (AF3 mode, server)

Per variant:
1. CDR loop CA-RMSD vs donor model: H1/H2/L1/L2 < 1.0 A,
   H3/L3 < 1.5 A (AF3 accuracy limits; higher = redesign needed).
2. VH/VL interface: buried interface area within 15% of donor Fv;
   identical interface-core residue side-chain conformations.
3. Antigen complex (if `--antigen`): paratope contact set of the variant
   must cover >= 90% of donor contacts; no new steric clashes
   (AF3 clash score / pLDDT per residue).
4. pLDDT of grafted CDRs >= pLDDT of donor CDRs - 5.

If any check fails: revert additional positions (from T3/contact lists),
re-run AF3, iterate (max 2 rounds before consulting experimental data).

## 5. Cross-validation with BioPhi (optional)

If available, compare per-position human-likeness scores; positions where
BioPhi disagrees with our germline choice should be reviewed manually.

## 6. Final gate

A humanized design is "production-ready" only after passing checks 1-3
(mandatory) and, on the server, check 4. Experimental validation then follows
docs/experimental_SOP.md.
