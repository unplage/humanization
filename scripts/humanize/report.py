"""Report generation: Markdown summary + per-position CSV + JSON dump."""

from __future__ import annotations

import csv
import json
import os
from typing import Dict, List, Optional

from .backmut import BackMutationResult
from .config import TIER_LABELS
from .minimal import MatrixEntry, MinimalReversion
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
            "antigen_contact", "empirical_ddG", "empirical_n", "features", "rationale",
        ])
        for c in backmut.candidates:
            w.writerow([
                c.position, c.donor_aa, c.human_aa, c.tier, c.composite,
                c.structural_score, c.benefit_score, c.chemical_score,
                c.buried, c.cdr_contact, c.antigen_contact,
                c.empirical_ddG, c.empirical_n,
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
            "cvi_homology": rep.cvi_homology,
            "germline": {
                "v_gene": rep.germline.v_gene.gene_id if rep.germline.v_gene else None,
                "j_gene": rep.germline.j_gene.gene_id if rep.germline.j_gene else None,
                "fr_identity": rep.germline.scores.get("fr_identity"),
                "cdr_identity": rep.germline.scores.get("cdr_identity"),
                "alternatives": [
                    (g.gene_id, s) for g, s in rep.germline.alternatives[:5]
                ],
            },
            "minimal_reversion": (
                {
                    "positions": rep.minimal_reversion.positions,
                    "method": rep.minimal_reversion.method,
                    "covered_contacts": rep.minimal_reversion.covered_contacts,
                    "total_contacts": rep.minimal_reversion.total_contacts,
                    "note": rep.minimal_reversion.note,
                } if rep.minimal_reversion else None
            ),
            "framework_matrix": [
                {
                    "germline": e.germline.gene_id,
                    "cvi_homology": e.cvi,
                    "n_backmutations": len(e.backmut.revert_positions(("T1", "T2"))),
                }
                for e in rep.matrix
            ],
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
                    "empirical_ddG": c.empirical_ddG,
                    "empirical_n": c.empirical_n,
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
            "developability_optimization": (
                {
                    "risks": [
                        {
                            "seq_pos": r.seq_pos,
                            "kabat_pos": r.kabat_pos,
                            "motif": r.motif,
                            "motif_name": r.motif_name,
                            "risk": r.risk,
                            "issue": r.issue,
                            "context": r.context,
                            "skip_reason": r.skip_reason,
                            "rel_sasa": r.rel_sasa,
                        }
                        for r in (rep.developability_optimization.risks if hasattr(rep.developability_optimization, 'risks') else [])
                    ],
                    "skipped_risks": [
                        {
                            "seq_pos": r.seq_pos,
                            "kabat_pos": r.kabat_pos,
                            "motif": r.motif,
                            "motif_name": r.motif_name,
                            "risk": r.risk,
                            "issue": r.issue,
                            "context": r.context,
                            "skip_reason": r.skip_reason,
                            "rel_sasa": r.rel_sasa,
                        }
                        for r in (rep.developability_optimization.skipped_risks if hasattr(rep.developability_optimization, 'skipped_risks') else [])
                    ],
                    "optimized_designs": rep.developability_optimization.optimized_designs if hasattr(rep.developability_optimization, 'optimized_designs') else [],
                    "warnings": rep.developability_optimization.warnings if hasattr(rep.developability_optimization, 'warnings') else [],
                    "summary": rep.developability_optimization.summary if hasattr(rep.developability_optimization, 'summary') else "",
                } if rep.developability_optimization else None
            ),
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
        L.append(f"**CVI homology (canonical+vernier+interface vs donor, BI 2024):** "
                 f"{rep.cvi_homology:.3f}")
        if rep.minimal_reversion is not None:
            mr = rep.minimal_reversion
            L.append("")
            L.append(f"**Minimal reversion set** ({mr.method}): "
                     f"{', '.join(mr.positions) or '(none beyond graft)'} "
                     f"| contacts preserved {mr.covered_contacts}/{mr.total_contacts} "
                     f"| {mr.note}")
        if rep.matrix:
            L.append("")
            L.append("**Framework matrix (alternative germlines, V2-class):**")
            L.append("")
            L.append("| germline | CVI homology | # back-mutations |")
            L.append("|----------|--------------|-----------------|")
            for e in rep.matrix:
                L.append(f"| {e.germline.gene_id} | {e.cvi:.3f} | "
                         f"{len(e.backmut.revert_positions(('T1', 'T2')))} |")
            L.append("")
        L.append("")

        L.append("### Back-mutation candidates (framework positions)")
        L.append("")
        L.append("| pos | donor | human | tier | score | empirical ddG | features |")
        L.append("|-----|-------|-------|------|-------|---------------|----------|")
        for c in sorted(rep.backmut.candidates, key=lambda c: -c.composite):
            tier_txt = c.tier + (" " + TIER_LABELS[c.tier][:24] if c.tier in TIER_LABELS else "")
            emp = f"{c.empirical_ddG:+.2f} (n={c.empirical_n})" if c.empirical_ddG is not None else "-"
            L.append(f"| {c.position} | {c.donor_aa} | {c.human_aa} | {tier_txt} | {c.composite} | {emp} | "
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
        if rep.humanness:
            L.append("### BioPhi/Sapiens humanness cross-check")
            L.append("")
            L.append("| variant | Sapiens mean | OASis identity |")
            L.append("|---------|--------------|----------------|")
            for k in sorted(rep.humanness):
                h = rep.humanness[k]
                sap = h.get("sapiens_mean")
                oas = h.get("oasis_identity")
                L.append(f"| {k} | {sap if sap is not None else '-'} | "
                         f"{oas if oas is not None else '-'} |")
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
        
        # ---- Developability optimization report (only show if structure data exists) ----
        has_structure = rep.structure_hints and rep.structure_hints.data.get("buried")
        if has_structure and rep.developability_optimization and (rep.developability_optimization.risks or rep.developability_optimization.skipped_risks):
            dev_opt = rep.developability_optimization
            L.append("### Developability optimization (Step 4)")
            L.append("")
            L.append(f"**Summary:** {dev_opt.summary}")
            L.append("")
            
            # Optimizable positions (surface-exposed)
            if dev_opt.risks:
                L.append("**Optimizable positions (surface-exposed):**")
                L.append("")
                L.append("| Position | Motif | Risk | Issue | relSASA | Context |")
                L.append("|----------|-------|------|-------|---------|---------|")
                for r in dev_opt.risks:
                    sasa = f"{r.rel_sasa:.2f}" if r.rel_sasa is not None else "unknown"
                    L.append(f"| {r.kabat_pos or r.seq_pos} | {r.motif} | {r.risk} | {r.issue} | {sasa} | {r.context} |")
                L.append("")
            
            # Skipped positions (structural constraints)
            if dev_opt.skipped_risks:
                L.append("**Skipped positions (structural constraints):**")
                L.append("")
                L.append("| Position | Motif | Risk | Issue | Reason |")
                L.append("|----------|-------|------|-------|--------|")
                for r in dev_opt.skipped_risks:
                    L.append(f"| {r.kabat_pos or r.seq_pos} | {r.motif} | {r.risk} | {r.issue} | {r.skip_reason} |")
                L.append("")
            
            if dev_opt.optimized_designs:
                L.append("**Optimized designs (risk-free):**")
                L.append("")
                L.append("ProteinMPNN has generated alternative sequences that eliminate the high-risk motifs while preserving the framework structure.")
                L.append("")
            else:
                L.append("**Note:** No risk-free designs generated. Consider:")
                L.append("- Running ProteinMPNN with `--mpnn-mode local` for actual optimization")
                L.append("- Manually reviewing the high-risk positions for experimental validation")
                L.append("")
            
            if dev_opt.warnings:
                L.append("**Warnings:**")
                for w in dev_opt.warnings:
                    L.append(f"- {w}")
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
    
    # 增强版报告 (借鉴 WeMol 格式)
    try:
        from .report_enhanced import generate_enhanced_report, write_enhanced_docx
        enhanced_md = os.path.join(outdir, "enhanced_report.md")
        enhanced_content = generate_enhanced_report(result, outdir)
        with open(enhanced_md, "w") as fh:
            fh.write(enhanced_content)
        paths["enhanced_markdown"] = enhanced_md
        
        # 生成 Word 格式增强报告
        enhanced_docx = write_enhanced_docx(result, outdir)
        if enhanced_docx:
            paths["enhanced_docx"] = enhanced_docx
    except Exception as e:
        paths["_errors"] = paths.get("_errors", [])
        paths["_errors"].append(f"enhanced report: {e}")

    # Word report (python-docx; skipped when unavailable)
    try:
        from .report_docx import write_docx_report
        docx_path = write_docx_report(result, outdir)
        if docx_path:
            paths["docx"] = docx_path
    except Exception as e:
        paths.setdefault("_errors", []).append(f"docx: {e}")

    # FASTA of all variants
    fasta = os.path.join(outdir, "variants.fasta")
    with open(fasta, "w") as fh:
        for rep in result.chains:
            for v in rep.variants:
                fh.write(f">{rep.input_chain.name}|{v.name}|{v.description}\n{v.sequence}\n")
    paths["fasta"] = fasta
    return paths
