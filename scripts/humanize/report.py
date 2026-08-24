"""Report generation: Markdown summary + per-position CSV + JSON dump."""

from __future__ import annotations

import csv
import json
import os
from typing import Dict, List, Optional

from .backmut import BackMutationResult
from .config import TIER_LABELS
from .pipeline import ChainReport, RunResult


def _num(v, nd=2):
    if v is None:
        return "-"
    return f"{v:.{nd}f}"


def write_csv(path: str, backmut: BackMutationResult) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "position", "donor_aa", "human_aa", "tier", "composite",
            "structural", "benefit", "chemical", "buried", "cdr_contact",
            "antigen_contact", "features", "rationale",
        ])
        for c in backmut.candidates:
            w.writerow([
                c.position, c.donor_aa, c.human_aa, c.tier, c.composite,
                c.structural_score, c.benefit_score, c.chemical_score,
                c.buried, c.cdr_contact, c.antigen_contact,
                "+".join(c.features), "; ".join(c.rationale),
            ])


def write_json(path: str, result: RunResult) -> None:
    payload = {
        "version": "0.1.0",
        "format": result.format,
        "warnings": result.warnings,
        "chains": [],
    }
    for rep in result.chains:
        chain = rep.input_chain
        chain_payload = {
            "name": chain.name,
            "chain_type": chain.chain_type,
            "is_vhh": chain.is_vhh,
            "sequence": chain.sequence,
            "germline": {
                "v_gene": rep.germline.v_gene.gene_id if rep.germline.v_gene else None,
                "j_gene": rep.germline.j_gene.gene_id if rep.germline.j_gene else None,
                "fr_identity": rep.germline.scores.get("fr_identity"),
                "cdr_identity": rep.germline.scores.get("cdr_identity"),
                "alternatives": [
                    (g.gene_id, s) for g, s in rep.germline.alternatives[:5]
                ],
            },
            "human_likeness": rep.human_likeness,
            "backmutations": [
                {
                    "position": c.position,
                    "donor_aa": c.donor_aa,
                    "human_aa": c.human_aa,
                    "tier": c.tier,
                    "composite": c.composite,
                    "structural": c.structural_score,
                    "benefit": c.benefit_score,
                    "chemical": c.chemical_score,
                    "features": c.features,
                    "rationale": c.rationale,
                }
                for c in rep.backmut.candidates
            ],
            "variants": [
                {
                    "name": v.name,
                    "description": v.description,
                    "sequence": v.sequence,
                    "backmutations": v.backmutations,
                }
                for v in rep.variants
            ],
        }
        payload["chains"].append(chain_payload)
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)


def write_markdown(path: str, result: RunResult) -> None:
    L = []
    L.append("# Antibody Humanization Report")
    L.append("")
    L.append(f"- Format: **{result.format.upper()}**")
    L.append(f"- CDR graft scheme (variants): **{result.config.cdr_scheme if result.config else 'kabat'}**")
    if result.warnings:
        L.append("- Warnings:")
        for w in result.warnings:
            L.append(f"  - {w}")
    L.append("")
    for rep in result.chains:
        chain = rep.input_chain
        L.append(f"## Chain: {chain.name} ({chain.chain_type}{' | VHH' if chain.is_vhh else ''})")
        L.append("")
        if rep.germline.v_gene:
            v, j = rep.germline.v_gene, rep.germline.j_gene
            s = rep.germline.scores
            L.append(f"**Germline choice:** {v.gene_id} (FR identity {_num(s.get('fr_identity'), 3)}, "
                     f"CDR identity {_num(s.get('cdr_identity'), 3)}) + {j.gene_id if j else '?'}")
            L.append("")
            L.append("| rank | gene | FR id | CDR id |")
            L.append("|------|------|-------|--------|")
            for g, sc in rep.germline.alternatives[:5]:
                L.append(f"| | {g.gene_id} | {_num(sc.get('fr_identity'), 3)} | {_num(sc.get('cdr_identity'), 3)} |")
            L.append("")
        if rep.human_likeness:
            L.append(f"**Human-likeness (FR germline identity, pure graft):** "
                     + ", ".join(f"{k} {v:.1f}%" for k, v in rep.human_likeness.items()))
            L.append("")

        L.append("### Back-mutation candidates (framework positions)")
        L.append("")
        L.append("| pos | donor | human | tier | score | features |")
        L.append("|-----|-------|-------|------|-------|----------|")
        for c in sorted(rep.backmut.candidates, key=lambda c: -c.composite):
            tier_txt = c.tier + (" " + TIER_LABELS[c.tier][:24] if c.tier in TIER_LABELS else "")
            L.append(f"| {c.position} | {c.donor_aa} | {c.human_aa} | {tier_txt} | {c.composite} | "
                     f"{'+'.join(c.features) or '-'} |")
        L.append("")
        counts = {}
        for c in rep.backmut.candidates:
            counts[c.tier] = counts.get(c.tier, 0) + 1
        L.append("**Tier summary:** " + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())) + "")
        L.append("")

        L.append("### Recommended variants")
        L.append("")
        L.append("| variant | description | # back-mutations | sequence |")
        L.append("|---------|-------------|------------------|----------|")
        for v in rep.variants:
            L.append(f"| {v.name} | {v.description} | {len(v.backmutations)} | `{v.sequence}` |")
        L.append("")
        L.append("### Variant differences (vs pure graft V0)")
        L.append("")
        L.append("| variant | positions reverted |")
        L.append("|---------|-------------------|")
        v0_map = rep.variants[0].graft.origin
        for v in rep.variants[1:]:
            diff = [p for p in v.graft.origin if v.graft.origin[p] != v0_map.get(p)]
            L.append(f"| {v.name} | {', '.join(diff) or '-'} |")
        L.append("")
    L.append("## Notes")
    L.append("")
    L.append("- Position numbering: Kabat (portable engine); ANARCI (exact) "
             "recommended on the server for verification.")
    L.append("- `KEEP_DONOR` positions (VHH hallmark / Cys) must NOT be reverted.")
    L.append("- Structure-based hints (buriedness, contacts) appear only when "
             "AF3 mode is enabled.")
    L.append("- See docs/experimental_SOP.md for the experimental validation path.")
    L.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


def write_all(outdir: str, result: RunResult) -> Dict[str, str]:
    paths = {}
    md = os.path.join(outdir, "humanization_report.md")
    js = os.path.join(outdir, "humanization_result.json")
    write_markdown(md, result)
    write_json(js, result)
    paths["markdown"] = md
    paths["json"] = js
    for rep in result.chains:
        csvp = os.path.join(outdir, f"backmutations_{rep.input_chain.name}.csv")
        write_csv(csvp, rep.backmut)
        paths[f"csv_{rep.input_chain.name}"] = csvp
    # FASTA of all variants
    fasta = os.path.join(outdir, "variants.fasta")
    with open(fasta, "w") as fh:
        for rep in result.chains:
            for v in rep.variants:
                fh.write(f">{rep.input_chain.name}|{v.name}|{v.description}\n{v.sequence}\n")
    paths["fasta"] = fasta
    return paths
