---
name: antibody-humanization
description: Antibody humanization pipeline (Fab and VHH). Use when the user provides a non-human antibody amino acid sequence (VH/VL, Fab, or single-domain VHH) and wants CDR grafting onto human germline frameworks, back-mutation (回复突变) candidate design with scoring, variant generation, structure validation via AlphaFold3/ProteinMPNN, or the full experimental SOP. Triggers on keywords: humanization, 人源化, back-mutation, 回复突变, CDR grafting, germline, Fab, VHH, nanobody, 纳米抗体.
---

# Antibody Humanization (Fab / VHH)

Production-grade CDR-grafting + back-mutation design pipeline. Runs fully on the
user's laptop (stdlib only) and enables ProteinMPNN/AF3 structure mode on the
server.

## When to use

- User supplies a non-human antibody sequence (mouse/rat/rabbit/camelid etc.)
  as FASTA (VH+VL for Fab, or a single VHH chain) and asks for humanization.
- User asks for back-mutation (回复突变) candidates with scoring, or for the
  complete humanization + experimental validation workflow.

## How to run

```bash
# laptop (portable mode, no installation):
python3 scripts/humanize/cli.py run --input my_antibody.fasta --outdir outputs
# server (full mode with structure validation):
python3 scripts/humanize/cli.py run --input my_antibody.fasta --outdir outputs \
  --af3-mode local --af3-binary /path/run_alphafold.py \
  --mpnn-mode local --mpnn-script /path/protein_mpnn.py
# tooling status:
python3 scripts/humanize/cli.py setup-check
```

Input FASTA convention: one record per V domain (VH + VL in one file, or a
single VHH). Chain types are auto-detected; camelid VHH is detected by the
FR2 hallmark (Kabat 37/44/45/47) plus single-chain status.

## Pipeline stages (what the skill does for the user)

1. **Numbering** — anchor-based Kabat engine (zero-dependency) in
   `scripts/humanize/numbering.py`. Exact mode (ANARCI) recommended on server
   for cross-validation (`humanize setup-check`; see docs/validation.md).
2. **Germline selection** — NCBI IgBLAST human germline (server download via
   `humanize setup-germline`) or the bundled human IMGT-derived set
   (373 V + 28 J genes, `data/germline/human_germline_kabat.json`).
   Strategy: rank by framework identity, then pick the top-FR gene with the
   highest CDR1+2 identity (hybrid approach, Olimpieri et al. 2015).
3. **CDR grafting** — per scheme (kabat/chothia/abm/imgt), all four are
   reported; the variant ladder uses the scheme chosen with `--scheme`.
   VHH: the FR2 hallmark residues 37/44/45/47 are ALWAYS kept from the donor.
4. **Back-mutation analysis** — the core deliverable. For every framework
   position where donor != human germline:
   - classify by structural role: vernier zone (Foote & Winter 1992),
     VH/VL interface (Chothia 1985), canonical residues (Chothia & Lesk 1987),
     buried/CDR-contact/antigen-contact (from AF3 when available),
     VHH hallmark, disulfide Cys
   - score: structural + immunogenicity + developability components
     (weights in `scripts/humanize/config.py`, rationale in docs/scoring.md)
   - tier: T1 must-revert / T2 recommended / T3 optional / KEEP_HUMAN /
     KEEP_DONOR (VHH hallmark, Cys — never revert)
5. **Variant ladder** — V0 pure graft, V1 (T1), V2 (T1+T2), V3 (T1+T2+exposed
   T3). All written to `outputs/variants.fasta`.
6. **Structure validation (server, optional)** — AF3 prediction of donor,
   variants and (optionally) complex with antigen; computes CDR loop RMSD,
   interface contacts, buriedness → refines back-mutation scores.
   ProteinMPNN: framework re-design with fixed CDRs/interface as an
   alternative design source, filtered by human-likeness.
7. **Experimental SOP** — docs/experimental_SOP.md: gene synthesis, expression,
   affinity (BLI/SPR), developability panel, Go/No-Go criteria.

## Outputs

`outputs/humanization_report.md` (full report), `humanization_result.json`
(machine-readable), `backmutations_*.csv` (per-position scoring),
`variants.fasta`. Key tables in the report: germline choice + alternatives,
back-mutation candidates with tier/score/features, variant sequences, and
the tier summary. **The recommended production choice is V2** (T1+T2).

## Rules that must never be violated

- Never revert VHH hallmark positions (37/44/45/47) or CDR-disulfide Cys.
- Never design back-mutations inside CDRs (they are donor by definition) or
  in FR4 (J region is human by construction).
- For immunogenicity claims, only surface-exposed mismatches count; buried
  mismatches are structurally tolerated and reverting them is low-value.
- Always check `outputs/humanization_report.md` warnings (e.g. unusual loop
  lengths, truncated J region, format suspicion) before recommending variants.

## References

- Foote & Winter, JMB 224:487 (1992) — vernier zone
- Chothia et al., JMB 186:651 (1985) — VH/VL interface
- Chothia & Lesk, JMB 196:901 (1987); Al-Lazikani et al., JMB 273:927 (1997) — canonical residues
- Muyldermans, J Biotechnol 74:277 (2001) — VHH hallmark
- Carter et al., PNAS 89:4285 (1992) — 4D5 → trastuzumab (validation case)
- Olimpieri et al., Bioinformatics 31:434 (2015) — germline selection strategy
