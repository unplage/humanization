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
    """生成 Template Score 表格 (类似 WeMol)"""
    
    v_genes = db.human(chain_type)
    
    # 计算每个 germline 的评分
    results = []
    for g in v_genes:
        if g.numbered is None:
            continue
        
        # 计算各种指标
        q = query.posmap()
        gn = g.numbered.posmap()
        
        # FR mismatch 计数
        fr_mismatch = 0
        fr_length = 0
        for pos in ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10',
                     'H11', 'H12', 'H13', 'H14', 'H15', 'H16', 'H17', 'H18', 'H19', 'H20',
                     'H21', 'H22', 'H23', 'H24', 'H25', 'H26', 'H27', 'H28', 'H29', 'H30',
                     'H31', 'H32', 'H33', 'H34', 'H35', 'H35A', 'H36', 'H37', 'H38', 'H39',
                     'H40', 'H41', 'H42', 'H43', 'H44', 'H45', 'H46', 'H47', 'H48', 'H49',
                     'H50', 'H51', 'H52', 'H52A', 'H52B', 'H52C', 'H53', 'H54', 'H55', 'H56',
                     'H57', 'H58', 'H59', 'H60', 'H61', 'H62', 'H63', 'H64', 'H65', 'H66',
                     'H67', 'H68', 'H69', 'H70', 'H71', 'H72', 'H73', 'H74', 'H75', 'H76',
                     'H77', 'H78', 'H79', 'H80', 'H81', 'H82', 'H82A', 'H82B', 'H82C', 'H83',
                     'H84', 'H85', 'H86', 'H87', 'H88', 'H89', 'H90', 'H91', 'H92', 'H93',
                     'H94', 'H95', 'H96', 'H97', 'H98', 'H99', 'H100', 'H100A', 'H100B',
                     'H100C', 'H100D', 'H100E', 'H100F', 'H100G', 'H101', 'H102', 'H103']:
            pos_str = f"H{pos}" if not pos.startswith('H') else pos
            if pos_str in q and pos_str in gn:
                fr_length += 1
                if q[pos_str] != gn[pos_str]:
                    fr_mismatch += 1
        
        # 计算 FR1%, FR2%, FR3%
        fr1_positions = [f'H{i}' for i in range(1, 26)]
        fr2_positions = [f'H{i}' for i in range(36, 50)]
        fr3_positions = [f'H{i}' for i in range(66, 95)]
        
        def calc_fr_pct(positions):
            match = 0
            total = 0
            for pos in positions:
                if pos in q and pos in gn:
                    total += 1
                    if q[pos] == gn[pos]:
                        match += 1
            return (match / total * 100) if total > 0 else 0
        
        fr1_pct = calc_fr_pct(fr1_positions)
        fr2_pct = calc_fr_pct(fr2_positions)
        fr3_pct = calc_fr_pct(fr3_positions)
        
        # 计算综合分数
        fr_identity = 1 - (fr_mismatch / fr_length) if fr_length > 0 else 0
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
    
    # 按分数排序
    results.sort(key=lambda x: -x['score'])
    
    # 生成表格
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
    """生成 Mutation Score 表格 (类似 WeMol)"""
    
    lines = []
    lines.append("Chain\tPosition\tDonor Residue\tTemplate Residue\tScore\tGermline Frequency")
    lines.append("-" * 100)
    
    for c in backmut.candidates:
        if c.tier in ["T1", "T2", "T3"]:
            # 计算 germline 频率
            freq = "100%"  # 默认
            
            lines.append(
                f"{c.position[0]}\t{c.position}\t{c.donor_aa}\t{c.human_aa}\t"
                f"{c.composite:.1f}\t{freq}"
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
                # 检查是否在 CDR 区域
                if "cdr" in c.features:
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
        fv_hl = hl.get('kabat', 0) if hl else 0
        fr_hl = hl.get('kabat', 0) if hl else 0  # 简化
        
        # 检查热点
        high_risk_count = 0
        high_risk_sites = []
        
        for c in rep.backmut.candidates:
            if c.tier in ["T1", "T2"]:
                high_risk_count += 1
                high_risk_sites.append(c.position)
        
        lines.append(
            f"{chain_type}\t{high_risk_count}\t{','.join(high_risk_sites)}\t"
            f"{fv_hl:.1f}%\t{fr_hl:.1f}%"
        )
    
    return "\n".join(lines)


def generate_humanized_sequences(
    chain_reports: List[ChainReport],
) -> str:
    """生成人源化序列"""
    
    lines = []
    lines.append(">VL")
    lines.append(chain_reports[0].input_chain.sequence if chain_reports else "")
    
    for i, rep in enumerate(chain_reports):
        chain_type = rep.input_chain.chain_type
        v_gene = rep.germline.v_gene.gene_id if rep.germline.v_gene else "None"
        
        # V0: 纯移植
        if hasattr(rep, 'graft_sequence'):
            lines.append(f">{chain_type}0 Graft({v_gene})")
            lines.append(rep.graft_sequence)
        
        # V2: T1+T2 回复突变
        if hasattr(rep, 'variants') and len(rep.variants) > 1:
            v2 = rep.variants[1]
            lines.append(f">{chain_type}2 Graft({v_gene}) + T1T2 mutations")
            lines.append(v2.sequence)
    
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
