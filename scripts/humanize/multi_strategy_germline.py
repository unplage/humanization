"""多策略 Germline 选择模块

支持多种 germline 选择策略：
1. FR 优先：FR 同源性最高
2. CDR 优先：CDR 同源性最高
3. 综合评分：0.7*FR + 0.3*CDR
4. CVI 优先：考虑 canonical + vernier + interface
5. 回复突变最少：回复突变数最少
6. Adimab 频率：基于 Adimab 治疗性抗体的使用频率
7. Pioneer 频率：基于 Pioneer 库的使用频率 (600+ 临床阶段抗体)
8. 复合评分：0.5*CVI + 0.3*频率 + 0.2*FR
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .germline import (
    GermlineDB,
    GermlineGene,
    compare_to_germline,
    score_j_match,
)
from .germline_frequency import (
    get_frequency,
    is_recommended,
)
from .numbering import NumberedChain


@dataclass
class GermlineCandidate:
    """Germline 候选"""
    gene: GermlineGene
    fr_identity: float
    cdr_identity: float
    composite_score: float
    cvi_score: float  # canonical + vernier + interface 同源性
    n_backmutations_est: int  # 估计的回复突变数
    strategy: str  # 选择策略
    frequency_score: float = 0.0  # 使用频率评分
    composite_3axis: float = 0.0  # 3 轴综合评分 (CVI + 频率 + FR)


@dataclass
class MultiStrategyResult:
    """多策略选择结果"""
    query_chain: str
    candidates: Dict[str, List[GermlineCandidate]] = field(default_factory=dict)
    
    def get_best(self, strategy: str) -> Optional[GermlineCandidate]:
        """获取特定策略的最佳候选"""
        lst = self.candidates.get(strategy, [])
        return lst[0] if lst else None
    
    def get_top(self, strategy: str, n: int = 5) -> List[GermlineCandidate]:
        """获取特定策略的前N个候选"""
        return self.candidates.get(strategy, [])[:n]
    
    def summary(self) -> str:
        """生成摘要"""
        lines = [f"Germline 选择结果 ({self.query_chain}):"]
        lines.append("-" * 70)
        
        for strategy, candidate_list in self.candidates.items():
            lines.append(f"\n策略: {strategy}")
            for i, c in enumerate(candidate_list[:3], 1):
                lines.append(f"  #{i} 基因: {c.gene.gene_id}")
                lines.append(f"      FR 同源性: {c.fr_identity:.4f}")
                lines.append(f"      CDR 同源性: {c.cdr_identity:.4f}")
                lines.append(f"      综合评分: {c.composite_score:.4f}")
                lines.append(f"      CVI 同源性: {c.cvi_score:.4f}")
                lines.append(f"      估计回复突变数: {c.n_backmutations_est}")
        
        return "\n".join(lines)


# CVI (canonical + vernier + interface) 位置集
CVI_POSITIONS = {
    "H": {
        "canonical": {24, 26, 27, 29, 33, 34, 47, 48, 49, 57, 58, 71, 78, 94},
        "vernier": {2, 27, 28, 29, 30, 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103},
        "interface_core": {37, 39, 45, 47, 91, 93, 95, 103},
    },
    "L": {
        "canonical": {2, 25, 33, 34, 36, 37, 45, 46, 48, 49, 58, 64, 71, 90, 94, 97, 98},
        "vernier": {2, 4, 35, 36, 37, 38, 43, 44, 45, 46, 48, 49, 58, 62, 63, 66, 67, 68, 69, 71, 87, 88, 98},
        "interface_core": {34, 36, 37, 38, 44, 45, 46, 87, 89, 91, 96, 98},
    },
}


def _calculate_cvi_score(query: NumberedChain, gene: GermlineGene) -> float:
    """计算 CVI (canonical + vernier + interface) 同源性
    
    基于 BI 2024 研究：CVI 同源性与表达量和亲和力保留显著相关
    """
    if gene.numbered is None:
        return 0.0
    
    q = query.posmap()
    g = gene.numbered.posmap()
    ctype = query.chain_type
    
    # 合并所有 CVI 位置
    all_cvi_positions = set()
    for category in CVI_POSITIONS.get(ctype, {}).values():
        all_cvi_positions.update(category)
    
    n = 0
    same = 0
    for pos in all_cvi_positions:
        pos_str = f"{ctype}{pos}"
        if pos_str in q and pos_str in g:
            if q[pos_str] and g[pos_str]:
                n += 1
                if q[pos_str] == g[pos_str]:
                    same += 1
    
    return same / n if n > 0 else 0.0


def _estimate_n_backmutations(
    query: NumberedChain,
    gene: GermlineGene,
    is_vhh: bool = False,
) -> int:
    """估计回复突变数量（基于实际的 tier 分析）。

    使用与下游完全一致的 analyze_backmutations 统计 **V2 变体（T1+T2）**
    的回复突变数——即生产推荐变体的真实回复数。旧实现统计 T1+T2+T3，
    而 T3 数量 ≈ 框架差异数，导致 min_backmutations 策略与 fr_best 趋同。
    """
    if gene.numbered is None:
        return 999
    try:
        from .backmut import analyze_backmutations
        bm = analyze_backmutations(query, gene, is_vhh=is_vhh)
        return sum(1 for c in bm.candidates if c.tier in ("T1", "T2"))
    except Exception:
        return 999


def choose_germlines_multi_strategy(
    query: NumberedChain,
    db: GermlineDB,
    n_alternatives: int = 5,
    is_vhh: bool = False,
) -> MultiStrategyResult:
    """多策略 Germline 选择
    
    提供多种选择策略，让用户可以根据具体需求选择：
    1. fr_best: FR 同源性最高
    2. cdr_best: CDR 同源性最高
    3. composite: 综合评分最高 (0.7*FR + 0.3*CDR)
    4. cvi_best: CVI 同源性最高
    5. min_backmutations: 估计回复突变最少
    6. current: 当前系统策略 (top 30% FR 中 CDR 最高)
    7. adimab_frequency: 基于 Adimab 使用频率
    8. pioneer_frequency: 基于 Pioneer 库使用频率
    9. composite_3axis: 3 轴综合评分 (0.5*CVI + 0.3*频率 + 0.2*FR)
    """
    ctype = query.chain_type
    v_genes = db.human(ctype)
    
    # 计算所有候选的评分
    scored = []
    for g in v_genes:
        s = compare_to_germline(query, g)
        cvi = _calculate_cvi_score(query, g)
        n_backmut = _estimate_n_backmutations(query, g, is_vhh=is_vhh)
        composite = 0.7 * s['fr_identity'] + 0.3 * s['cdr_identity']
        
        scored.append({
            'gene': g,
            'fr': s['fr_identity'],
            'cdr': s['cdr_identity'],
            'composite': composite,
            'cvi': cvi,
            'n_backmut': n_backmut,
        })
    
    # FR 阈值：默认 60%（可移植性下限）。当没有候选达到该阈值时
    # （高度非人源序列，如部分鼠源 VH），逐级降级以保证流程可继续，
    # 并在结果中记录实际使用的阈值。优先选出 FR 最高的可行 germline。
    fr_threshold = 0.60
    pool = [x for x in scored if x['fr'] >= fr_threshold]
    for t in (0.50, 0.40, 0.30, 0.20, 0.0):
        if pool:
            break
        fr_threshold = t
        pool = [x for x in scored if x['fr'] >= t]
    scored = pool
    
    if not scored:
        # 仍然没有满足条件的（理论上不会发生，因为 t=0.0 时全部通过）
        return MultiStrategyResult(query_chain=ctype)
    if fr_threshold < 0.60:
        import warnings as _w
        _w.warn(
            f"[{ctype}] 无人源 germline 达到 FR 60% 阈值；"
            f"自动降级至 FR >= {fr_threshold:.0%} 继续选择。"
            f"该序列人源化难度较高，回复突变数会显著增加。"
        )
    
    result = MultiStrategyResult(query_chain=ctype)
    
    # 先计算所有候选的频率
    for item in scored:
        frequency = get_frequency(ctype, item['gene'].gene_id)
        item['frequency'] = frequency
        item['composite_3axis'] = 0.5 * item['cvi'] + 0.3 * frequency + 0.2 * item['fr']
    
    def _make_candidate(item, strategy):
        return GermlineCandidate(
            gene=item['gene'],
            fr_identity=item['fr'],
            cdr_identity=item['cdr'],
            composite_score=item['composite'],
            cvi_score=item['cvi'],
            n_backmutations_est=item['n_backmut'],
            strategy=strategy,
            frequency_score=item['frequency'],
            composite_3axis=item['composite_3axis'],
        )
    
    def _top_n(items, n=5):
        return [_make_candidate(x, '') for x in items[:n]]
    
    # 策略1: FR 同源性最高（FR 相同时用使用频率作为 tiebreaker）
    sorted_fr = sorted(scored, key=lambda x: (-x['fr'], -x['frequency']))
    result.candidates['fr_best'] = _top_n(sorted_fr)
    for c in result.candidates['fr_best']:
        c.strategy = 'fr_best'
    
    # 策略2: CDR 同源性最高（CDR 相同时用使用频率作为 tiebreaker）
    sorted_cdr = sorted(scored, key=lambda x: (-x['cdr'], -x['frequency']))
    result.candidates['cdr_best'] = _top_n(sorted_cdr)
    for c in result.candidates['cdr_best']:
        c.strategy = 'cdr_best'
    
    # 策略3: 综合评分最高（composite 相同时用使用频率作为 tiebreaker）
    sorted_composite = sorted(scored, key=lambda x: (-x['composite'], -x['frequency']))
    result.candidates['composite'] = _top_n(sorted_composite)
    for c in result.candidates['composite']:
        c.strategy = 'composite'
    
    # 策略4: CVI 同源性最高（CVI 相同时用使用频率作为 tiebreaker）
    sorted_cvi = sorted(scored, key=lambda x: (-x['cvi'], -x['frequency']))
    result.candidates['cvi_best'] = _top_n(sorted_cvi)
    for c in result.candidates['cvi_best']:
        c.strategy = 'cvi_best'
    
    # 策略5: 估计回复突变最少（回突变相同时用使用频率作为 tiebreaker）
    sorted_min = sorted(scored, key=lambda x: (x['n_backmut'], -x['frequency']))
    result.candidates['min_backmutations'] = _top_n(sorted_min)
    for c in result.candidates['min_backmutations']:
        c.strategy = 'min_backmutations'
    
    # 策略6: 当前系统策略 (top 30% FR 中 CDR 最高)
    sorted_by_fr = sorted(scored, key=lambda x: -x['fr'])
    top_30 = sorted_by_fr[:max(1, int(len(sorted_by_fr) * 0.3))]
    sorted_current = sorted(top_30, key=lambda x: (-x['cdr'], -x['frequency']))
    result.candidates['current'] = _top_n(sorted_current)
    for c in result.candidates['current']:
        c.strategy = 'current'
    
    # 策略7: Adimab 频率优先 (基于使用频率加权)
    for item in scored:
        is_adimab_rec = is_recommended(ctype, item['gene'].gene_id, "adimab")
        is_pioneer_rec = is_recommended(ctype, item['gene'].gene_id, "pioneer")
        item['adimab_score'] = item['frequency'] * (1.2 if is_adimab_rec else 1.0)
        item['pioneer_score'] = item['frequency'] * (1.2 if is_pioneer_rec else 1.0)
    
    sorted_adimab = sorted(scored, key=lambda x: -x['adimab_score'])
    result.candidates['adimab_frequency'] = _top_n(sorted_adimab)
    for c in result.candidates['adimab_frequency']:
        c.strategy = 'adimab_frequency'
    
    # 策略8: Pioneer 频率优先
    sorted_pioneer = sorted(scored, key=lambda x: -x['pioneer_score'])
    result.candidates['pioneer_frequency'] = _top_n(sorted_pioneer)
    for c in result.candidates['pioneer_frequency']:
        c.strategy = 'pioneer_frequency'
    
    # 策略9: 3 轴综合评分 (CVI + 频率 + FR)
    sorted_3axis = sorted(scored, key=lambda x: -x['composite_3axis'])
    result.candidates['composite_3axis'] = _top_n(sorted_3axis)
    for c in result.candidates['composite_3axis']:
        c.strategy = 'composite_3axis'
    
    return result


def format_multi_strategy_report(result: MultiStrategyResult) -> str:
    """格式化多策略选择报告"""
    lines = []
    lines.append("=" * 70)
    lines.append(f"多策略 Germline 选择报告 ({result.query_chain} 链)")
    lines.append("=" * 70)
    
    # 表头
    lines.append(f"\n{'策略':<20} {'排名':<4} {'Germline':<18} {'FR':<10} {'CDR':<10} {'综合':<10} {'CVI':<10} {'频率':<10} {'3轴':<10} {'回复突变'}")
    lines.append("-" * 110)
    
    for strategy, candidate_list in result.candidates.items():
        for i, c in enumerate(candidate_list[:5], 1):
            lines.append(
                f"{strategy:<20} {i:<4} {c.gene.gene_id:<18} "
                f"{c.fr_identity:<10.4f} {c.cdr_identity:<10.4f} "
                f"{c.composite_score:<10.4f} {c.cvi_score:<10.4f} "
                f"{c.frequency_score:<10.4f} {c.composite_3axis:<10.4f} "
                f"{c.n_backmutations_est}"
            )
    
    lines.append("\n" + "=" * 70)
    lines.append("策略说明:")
    lines.append("  fr_best: FR 同源性最高 - 框架区域最接近人源")
    lines.append("  cdr_best: CDR 同源性最高 - CDR 区域最接近人源")
    lines.append("  composite: 综合评分最高 (0.7*FR + 0.3*CDR)")
    lines.append("  cvi_best: CVI 同源性最高 - canonical+vernier+interface 最优")
    lines.append("  min_backmutations: 估计回复突变最少")
    lines.append("  current: 当前系统策略 (top 30% FR 中 CDR 最高)")
    lines.append("  adimab_frequency: Adimab 推荐 germline 优先 + 频率加权")
    lines.append("  pioneer_frequency: Pioneer 库 germline 优先 + 频率加权")
    lines.append("  composite_3axis: 3 轴综合评分 (0.5*CVI + 0.3*频率 + 0.2*FR)")
    lines.append("=" * 70)
    
    lines.append("\n基于 BI 2024 研究建议:")
    lines.append("  - CVI 同源性与表达量和亲和力保留显著相关")
    lines.append("  - VL 框架同源性 < 60% 会导致 98% 亲和力丢失")
    lines.append("  - 最高同源性的框架不一定是最优选择")
    lines.append("  - 推荐使用 cvi_best 或 fr_best 策略")
    lines.append("\n基于行业数据:")
    lines.append("  - IGHV3-23 是治疗性抗体中最常用的 VH germline (~18%)")
    lines.append("  - IGHV1-69 是第二常用的 VH germline (~12%)")
    lines.append("  - IGKV1-39 是治疗性抗体中最常用的 VK germline (~15%)")
    
    return "\n".join(lines)
