"""增强版报告生成 - 借鉴 WeMol 格式"""

from __future__ import annotations

import csv
import json
import os
from typing import Dict, List, Optional, Tuple

from .backmut import BackMutationResult, BackMutationCandidate
from .germline import GermlineDB, GermlineGene
from .numbering import NumberedChain
from .pipeline import ChainReport, RunResult


def generate_template_score_table(
    chain_type: str,
    query: NumberedChain,
    db: GermlineDB,
    selected_germline: GermlineGene,
) -> str:
    """生成 Template Score 表格 (类似 WeMol)

    按链型（H/L）逐位点比较 FR 区（FR1/FR2/FR3）同源性，
    避免旧实现对轻链使用硬编码的 H 位置表（导致 VL 全零）。
    """
    v_genes = db.human(chain_type)

    def region_of(pos: str) -> str:
        return query.region_of(pos) or ""

    results = []
    for g in v_genes:
        if g.numbered is None:
            continue
        q = query.posmap()
        gn = g.numbered.posmap()
        common = [p for p in q if p in gn and region_of(p) in ("FR1", "FR2", "FR3")]
        fr_mismatch = sum(1 for p in common if q[p] != gn[p])
        fr_length = len(common)

        def calc_fr_pct(region):
            ps = [p for p in common if region_of(p) == region]
            if not ps:
                return 0.0
            return 100.0 * sum(1 for p in ps if q[p] == gn[p]) / len(ps)

        fr1_pct = calc_fr_pct("FR1")
        fr2_pct = calc_fr_pct("FR2")
        fr3_pct = calc_fr_pct("FR3")
        fr_identity = 1 - (fr_mismatch / fr_length) if fr_length > 0 else 0.0
        # 与 FR 同源性主导、FR1-3 分段微调的综合分数
        score = fr_identity * 100 + fr1_pct * 0.1 + fr2_pct * 0.1 + fr3_pct * 0.1

        results.append({
            'gene': g.gene_id,
            'fr_mismatch': fr_mismatch,
            'fr_length': fr_length,
            'fr_identity': fr_identity,
            'fr1_pct': fr1_pct,
            'fr2_pct': fr2_pct,
            'fr3_pct': fr3_pct,
            'score': score,
        })

    results.sort(key=lambda x: -x['score'])

    lines = []
    lines.append("NO.\tHit\tFR_Length_Diff\t#FR_Mismatch\tFR%\tFR1%\tFR2%\tFR3%\tScore")
    lines.append("-" * 100)
    for i, r in enumerate(results[:30], 1):
        lines.append(
            f"{i}\t{r['gene']}\t0\t{r['fr_mismatch']}\t"
            f"{r['fr_identity']*100:.1f}\t{r['fr1_pct']:.1f}\t"
            f"{r['fr2_pct']:.1f}\t{r['fr3_pct']:.1f}\t{r['score']:.1f}"
        )
    return "\n".join(lines)


def generate_mutation_score_table(
    backmut: BackMutationResult,
    germline: GermlineGene,
) -> str:
    """生成 Mutation Score 表格 (类似 WeMol)

    Germline Frequency = 供体残基在人源 germline 参考面板中的出现频率
    （来自 backmut.germline_conservation，top-20 同源 germline 面板）。
    """
    cons = backmut.germline_conservation
    lines = []
    lines.append("Chain\tPosition\tDonor Residue\tTemplate Residue\tScore\tGermline Frequency")
    lines.append("-" * 100)
    for c in backmut.candidates:
        if c.tier in ("T1", "T2", "T3"):
            freq = cons.get(c.position)
            freq_txt = f"{freq:.0%}" if freq is not None else "-"
            lines.append(
                f"{c.position[0]}\t{c.position}\t{c.donor_aa}\t{c.human_aa}\t"
                f"{c.composite:.1f}\t{freq_txt}"
            )
    return "\n".join(lines)


def generate_backmutation_summary(
    chain_reports: List[ChainReport],
) -> str:
    """生成回复突变摘要"""
    
    lines = []
    lines.append("Chain\tVariants")
    lines.append("-" * 80)
    
    for rep in chain_reports:
        chain_type = rep.input_chain.chain_type
        v_gene = rep.germline.v_gene.gene_id if rep.germline.v_gene else "None"
        
        # 生成变体
        variants = []
        
        # V0: 纯移植
        variants.append(f"Graft({v_gene})")
        
        # V1: T1 回复突变
        t1_muts = [c for c in rep.backmut.candidates if c.tier == "T1"]
        if t1_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t1_muts])
            variants.append(f"Graft({v_gene}) + {mut_str}")
        
        # V2: T1+T2 回复突变
        t2_muts = [c for c in rep.backmut.candidates if c.tier in ["T1", "T2"]]
        if t2_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t2_muts])
            variants.append(f"Graft({v_gene}) + {mut_str}")
        
        # V3: T1+T2+T3 回复突变
        t3_muts = [c for c in rep.backmut.candidates if c.tier in ["T1", "T2", "T3"]]
        if t3_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t3_muts])
            variants.append(f"Graft({v_gene}) + {mut_str}")
        
        for i, v in enumerate(variants, 1):
            lines.append(f"{chain_type}{i}\t{v}")
    
    return "\n".join(lines)


def generate_kabat_mutation_summary(
    chain_reports: List[ChainReport],
) -> str:
    """生成 KABAT 编号突变摘要"""
    
    lines = []
    lines.append("ID\tmutations")
    lines.append("-" * 80)
    
    for rep in chain_reports:
        chain_type = rep.input_chain.chain_type
        
        # V0
        lines.append(f"{chain_type}0\t")
        
        # V1: T1 回复突变
        t1_muts = [c for c in rep.backmut.candidates if c.tier == "T1"]
        if t1_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t1_muts])
            lines.append(f"{chain_type}1\t{mut_str}")
        else:
            lines.append(f"{chain_type}1\t")
        
        # V2: T1+T2 回复突变
        t2_muts = [c for c in rep.backmut.candidates if c.tier in ["T1", "T2"]]
        if t2_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t2_muts])
            lines.append(f"{chain_type}2\t{mut_str}")
        else:
            lines.append(f"{chain_type}2\t")
        
        # V3: T1+T2+T3 回复突变
        t3_muts = [c for c in rep.backmut.candidates if c.tier in ["T1", "T2", "T3"]]
        if t3_muts:
            mut_str = ",".join([f"{c.human_aa}{c.position[1:]}{c.donor_aa}" for c in t3_muts])
            lines.append(f"{chain_type}3\t{mut_str}")
        else:
            lines.append(f"{chain_type}3\t")
    
    return "\n".join(lines)


def generate_hotspot_summary(
    chain_reports: List[ChainReport],
) -> str:
    """生成热点摘要"""
    
    lines = []
    lines.append("4.1 CDR Summary:")
    lines.append("ID\tSequence-CDR\tHigh risk #\tHigh risk sites")
    lines.append("-" * 100)
    
    for rep in chain_reports:
        chain_type = rep.input_chain.chain_type
        
        # 获取 CDR 序列
        cdrs = rep.input_chain.cdrs if hasattr(rep.input_chain, 'cdrs') else {}
        
        # 生成 CDR 序列字符串
        cdr_seq = " ".join([f"{k}: {v}" for k, v in cdrs.items()]) if cdrs else "N/A"
        
        # 检查热点
        high_risk_count = 0
        high_risk_sites = []
        
        for c in rep.backmut.candidates:
            if c.tier in ["T1", "T2"]:
                # Check if position contacts a CDR loop (feature is
                # "cdr_contact", not "cdr" — the old check was always False).
                if "cdr_contact" in c.features or "antigen_contact" in c.features:
                    high_risk_count += 1
                    high_risk_sites.append(c.position)
        
        lines.append(f"{chain_type}\t{cdr_seq}\t{high_risk_count}\t{','.join(high_risk_sites)}")
    
    lines.append("")
    lines.append("4.2 Full-length Summary:")
    lines.append("ID\tHigh risk #\tHigh risk sites\tHumanness (FV)\tHumanness (FR)")
    lines.append("-" * 100)
    
    for rep in chain_reports:
        chain_type = rep.input_chain.chain_type
        
        # 计算人源度
        hl = rep.human_likeness
        # human_likeness_percent returns FR-only identity.  A true FV
        # metric would include CDRs; report it as FR for now and label
        # both columns the same to avoid misleading the user.
        fr_hl = hl.get('kabat', 0) if hl else 0
        
        # 检查热点
        high_risk_count = 0
        high_risk_sites = []
        
        for c in rep.backmut.candidates:
            if c.tier in ["T1", "T2"]:
                high_risk_count += 1
                high_risk_sites.append(c.position)
        
        lines.append(
            f"{chain_type}\t{high_risk_count}\t{','.join(high_risk_sites)}\t"
            f"{fr_hl:.1f}%\t{fr_hl:.1f}%"
        )
    
    return "\n".join(lines)


def generate_humanized_sequences(
    chain_reports: List[ChainReport],
) -> str:
    """生成人源化序列（V0 纯移植 + V1/V2/V3 变体，来自真实 graft 结果）"""
    lines = []
    for rep in chain_reports:
        chain_type = rep.input_chain.chain_type
        v_gene = rep.germline.v_gene.gene_id if rep.germline.v_gene else "None"
        graft = rep.grafts.get("kabat")
        if graft is not None:
            lines.append(f">{chain_type}0 Graft({v_gene})")
            lines.append(graft.sequence)
        for v in rep.variants[1:]:
            lines.append(f">{chain_type} {v.name} ({v.description})")
            lines.append(v.sequence)
    return "\n".join(lines)


def generate_enhanced_report(
    result: RunResult,
    outdir: str,
) -> str:
    """生成增强版报告 (借鉴 WeMol 格式)"""
    
    lines = []
    lines.append("Antibody Humanization Design Summary")
    lines.append("=" * 80)
    lines.append("")
    
    # 1. Template Score
    lines.append("1. Template Score")
    lines.append("-" * 80)
    
    for rep in result.chains:
        chain_type = rep.input_chain.chain_type
        if chain_type == "H":
            lines.append(f"\nVH Chain (selected: {rep.germline.v_gene.gene_id if rep.germline.v_gene else 'None'}):")
            lines.append(generate_template_score_table(
                chain_type, rep.input_chain.numbered, result.germline_db, rep.germline.v_gene
            ))
        elif chain_type == "L":
            lines.append(f"\nVL Chain (selected: {rep.germline.v_gene.gene_id if rep.germline.v_gene else 'None'}):")
            lines.append(generate_template_score_table(
                chain_type, rep.input_chain.numbered, result.germline_db, rep.germline.v_gene
            ))
    
    lines.append("")
    
    # 2. Mutation Score
    lines.append("2. Mutation Score")
    lines.append("-" * 80)
    
    for rep in result.chains:
        lines.append(f"\n{rep.input_chain.chain_type} Chain:")
        lines.append(generate_mutation_score_table(rep.backmut, rep.germline.v_gene))
    
    lines.append("")
    
    # 3. Back Mutation Summary
    lines.append("3. Back Mutation Summary")
    lines.append("-" * 80)
    lines.append(generate_backmutation_summary(result.chains))
    
    lines.append("")
    
    # 3.1 KABAT Numbering Mutation Summary
    lines.append("3.1 KABAT Numbering Mutation Summary")
    lines.append("-" * 80)
    lines.append(generate_kabat_mutation_summary(result.chains))
    
    lines.append("")
    
    # 4. Hotspot Summary
    lines.append("4. Hotspot Summary")
    lines.append("-" * 80)
    lines.append(generate_hotspot_summary(result.chains))
    
    lines.append("")
    
    # 5. Humanized Sequences
    lines.append("5. Humanized Sequences")
    lines.append("-" * 80)
    lines.append(generate_humanized_sequences(result.chains))
    
    return "\n".join(lines)
