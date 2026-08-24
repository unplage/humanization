# Experimental validation SOP (实验验证全流程)

Goal: confirm that humanized variants retain affinity, specificity, epitope,
and developability of the parental antibody. Timeline ~6-10 weeks.

## 0. Input from the pipeline

Recommended lead: **V2 (T1+T2)**. Backup: V1 (conservative), V3 (aggressive),
and (if generated) the ProteinMPNN/consensus design. The parental
(murine/camelid) antibody is the positive control throughout.

## 1. Gene synthesis & cloning (Week 0-1)

- Codon-optimize VH/VL or VHH for the expression host (e.g. Expi293; use
  standard codons; avoid repeats >8 nt; keep 5' Kozak).
- Formats per target:
  - Fab: VH-CH1 + kappa/lambda CL in separate or bicistronic vectors.
  - VHH: VHH-hFc (human IgG1 Fc, LALA or N297A to reduce Fc effector
    function if not needed) and/or VHH-HiBiT for kinetics.
  - Full IgG: VH-IgG1 + VL-CL.
- Synthesis QC: Sanger confirm full ORF; plasmid miniprep, sequence 3 clones.
- Suggested clone panel: V2 (lead) + V0 (pure graft, activity probe) + V1 +
  parental. V3 only if immunogenicity is a known issue.

## 2. Expression & purification (Week 2-4)

- Expi293F (or CHO-S), 1 mL test scale (96-well deep block) for all variants
  in parallel; then 50-100 mL scale-up for the lead.
- Transfect 1:1 H:L DNA ratio (Fab/IgG) or single chain (VHH).
- Harvest day 5-6; Protein A (VHH-Fc/IgG) or KappaSelect/λ-select (Fab);
  SEC polishing step (Superdex 200 Increase).
- QC panel per variant:
  - SEC: % monomer (target >95%).
  - SDS-PAGE (reduced/non-reduced): intact chain masses, no clipping.
  - Mass spec (LC-MS intact or IdeS digest): verify sequence mass, no
    unexpected PTMs.
  - Endotoxin < 0.5 EU/mg if in vivo work planned.

## 3. Functional validation (Week 4-6)

### 3.1 Affinity (kinetics) — mandatory

- BLI (Octet/ForteBio) or SPR (Biacore 8K/S200): capture anti-human Fc
  (IgG/VHH-Fc) or anti-Fab (Fab), titrate antigen 5-7 concentrations
  (2-3x dilution), fit 1:1 (with mass-transport check) or bivalent model.
- **Acceptance:** KD within 2-3x of parental (ideally <2x). For V2
  deviations >5x trigger structural re-analysis (AF3 complex + contact
  comparison) and possibly V1 instead.
- Report kon/koff: koff changes dominate; slowed koff usually OK, faster koff
  is a warning.

### 3.2 Activity / potency — depends on MOA

- Reporter assay, target-cell binding (FACS with titrated Ab), competition
  ELISA vs parental, or antigen-binding ELISA. Acceptance: EC50 within 3x.
- Epitope binning: cross-block parental vs variant vs competitor mAbs
  (BLI or FACS); must remain in the same bin.

### 3.3 Specificity (if required)

- SPR off-rate panel or membrane FACS against counter-targets; human
  proteome panel optional (e.g. Retrogenix) if de-risking immunogenicity.

## 4. Developability panel (Week 6-7)

- DSF (thermal stability, Tm1/Tm2; target Tm1 > 65C for VHH-Fc/IgG).
- Accelerated stability: 40C/2 weeks, SEC + CE-SDS (aggregation/fragmentation).
- HIC retention (surface hydrophobicity; compare with parental).
- CIC (charge variants; IEX) — optional.
- VHH-specific: CDR3 disulfide integrity (mass spec under reducing conditions;
  check the grafted CDR3 Cys pairing), solubility (PEG precipitation or 100
  mg/mL concentration test).
- Any new risk motifs flagged by the pipeline (N-glycan/deamidation in FR
  after reversion) must be checked in the sequence before synthesis and, if
  present, the position excluded from the variant.

## 5. Immunogenicity risk assessment (Week 6-8, optional but recommended)

- MHC-II binding prediction (netMHCIIpan 4.0) for peptides spanning the
  remaining non-human residues; compare V2 vs V0 vs parental.
- In vitro T-cell assay (ELISpot / CD4 proliferation) on donor PBMC panel
  (n >= 20 donors) if regulatory submission is planned.
- Expectation: V2's T-cell epitope score close to fully-human controls.

## 6. Go/No-Go criteria

| Criterion | Go | Conditional | No-Go |
|-----------|-----|-------------|-------|
| Affinity (KD vs parental) | <= 2x | 2-5x (consider V1) | > 5x (structural fix required) |
| Potency (EC50 vs parental) | <= 3x | 3-5x | > 5x |
| SEC monomer | > 95% | 90-95% | < 90% |
| Tm1 (DSF) | >= parental - 3C | -3 to -6C | < -6C |
| Epitope bin | same | - | changed |
| Developability flags | none new | 1 flag (fixable) | >= 2 flags |

## 7. If affinity is lost (troubleshooting loop)

1. Run the pipeline again with `--antigen <seq>` and `--donor-structure`
   (or AF3 donor model): compare CDR loop RMSD + paratope contacts between
   V2 and donor.
2. Identify lost contacts: revert additional contact-adjacent positions
   (from the T3 list that contact the CDR backbone, not just solvent).
3. If VH/VL orientation changed: revert interface-core positions from the
   T3 list (L89/L91/L96-class positions that were borderline).
4. Optional: ProteinMPNN redesign of FR2/FR3 patches with CDR+interface
   fixed, rescore with AF3, iterate.
5. Re-test lead at 2-4 variants per round.

## 8. Documentation & release

- Archive: input FASTA, report md/json/csv, clone sequences, synthesis QC,
  assay data, decision log. Naming: `<Ab>_hV2_r1` (round 1).
- Transfer to GMP vector / stable cell line only after criteria met.
