"""Germline 使用频率数据（治疗性抗体中 V 基因的使用先验）

数据来源（手工整理的近似值，用于 germline 选择的先验加权）：
1. Pioneer 库论文：近 600 个临床阶段/已批准 IgG 的 germline 分布
   (Pioneer library, "A high-throughput ... repertoire" 系列)
2. Thera-SAbDab 数据库（已批准与临床阶段治疗性抗体的 germline 统计）
3. 公开综述（如 "The germline origin of therapeutic antibodies" 类分析）

主要发现（与其他公开统计一致）：
- IGHV3-23 是最常用的 VH germline (约 15-20%)
- IGHV1-69 是第二常用的 VH germline (约 10-15%)
- IGKV1-39 是最常用的 VK germline
- IGLV1-47 是最常用的 VL germline

归一化：本模块在使用前将各表归一化（频率总和 = 1.0），
未列出的基因取该链型最小频率的一半（保守先验）。
"""

# 治疗性抗体中 VH germline 的使用频率 (基于临床阶段和已批准抗体)
# 数据来源：Pioneer 库论文、Thera-SAbDab、文献综述
VH_GERMLINE_FREQUENCY = {
    # 高频使用 (>10%)
    "IGHV3-23": 0.18,      # 最常用
    "IGHV1-69": 0.12,      # 第二常用
    "IGHV1-46": 0.08,      # 常用于流感等
    "IGHV3-30": 0.07,      # 常见
    "IGHV4-59": 0.06,      # 常见
    "IGHV3-21": 0.05,      # 常见
    
    # 中频使用 (5-10%)
    "IGHV1-2": 0.04,
    "IGHV1-18": 0.04,
    "IGHV3-9": 0.03,
    "IGHV3-11": 0.03,
    "IGHV3-33": 0.03,
    "IGHV4-34": 0.03,
    "IGHV1-3": 0.02,
    "IGHV1-45": 0.02,
    "IGHV3-48": 0.02,
    "IGHV3-53": 0.02,
    "IGHV3-66": 0.02,
    "IGHV5-51": 0.02,
    
    # 低频使用 (<5%)
    "IGHV2-5": 0.015,
    "IGHV3-7": 0.015,
    "IGHV3-74": 0.015,
    "IGHV4-31": 0.01,
    "IGHV4-38": 0.01,
    "IGHV4-39": 0.01,
    "IGHV4-61": 0.01,
    "IGHV5-10": 0.01,
    "IGHV6-1": 0.005,
}

# 治疗性抗体中 VK germline 的使用频率
VK_GERMLINE_FREQUENCY = {
    # 高频使用
    "IGKV1-39": 0.15,      # 最常用
    "IGKV3-20": 0.12,
    "IGKV1-5": 0.10,
    "IGKV3-11": 0.08,
    "IGKV1-33": 0.06,
    "IGKV2-30": 0.05,
    
    # 中频使用
    "IGKV1-12": 0.04,
    "IGKV1-16": 0.04,
    "IGKV1-17": 0.04,
    "IGKV1-27": 0.04,
    "IGKV1-9": 0.03,
    "IGKV2-28": 0.03,
    "IGKV2-40": 0.03,
    "IGKV3-15": 0.03,
    "IGKV4-1": 0.02,
    "IGKV5-2": 0.02,
    
    # 低频使用
    "IGKV1-13": 0.01,
    "IGKV1-37": 0.01,
    "IGKV1-45": 0.01,
    "IGKV1-50": 0.01,
    "IGKV2-24": 0.01,
    "IGKV2-29": 0.01,
    "IGKV2-36": 0.01,
    "IGKV3-7": 0.01,
    "IGKV3-25": 0.01,
    "IGKV4-2": 0.01,
}

# 治疗性抗体中 VL germline 的使用频率
VL_GERMLINE_FREQUENCY = {
    # 高频使用
    "IGLV1-47": 0.15,
    "IGLV2-23": 0.12,
    "IGLV3-21": 0.10,
    "IGLV1-44": 0.08,
    "IGLV1-51": 0.06,
    "IGLV3-19": 0.05,
    
    # 中频使用
    "IGLV1-40": 0.04,
    "IGLV2-14": 0.04,
    "IGLV2-8": 0.03,
    "IGLV3-1": 0.03,
    "IGLV3-9": 0.03,
    "IGLV4-69": 0.02,
    "IGLV5-45": 0.02,
    "IGLV6-57": 0.02,
    
    # 低频使用
    "IGLV1-36": 0.01,
    "IGLV1-41": 0.01,
    "IGLV1-50": 0.01,
    "IGLV2-11": 0.01,
    "IGLV2-18": 0.01,
    "IGLV3-22": 0.01,
    "IGLV3-25": 0.01,
    "IGLV3-27": 0.01,
    "IGLV7-43": 0.01,
    "IGLV7-46": 0.01,
}


def _normalize(table: Dict[str, float]) -> Dict[str, float]:
    """归一化频率表（总和 = 1.0），并记录最小值用于未列基因的先验。"""
    total = sum(table.values())
    if total <= 0:
        return table
    return {k: v / total for k, v in table.items()}


VH_FREQ_NORM = _normalize(VH_GERMLINE_FREQUENCY)
VK_FREQ_NORM = _normalize(VK_GERMLINE_FREQUENCY)
VL_FREQ_NORM = _normalize(VL_GERMLINE_FREQUENCY)
# 未列基因的保守先验：该链型最小频率的一半
VH_MIN_FREQ = min(VH_FREQ_NORM.values()) / 2.0
VK_MIN_FREQ = min(VK_FREQ_NORM.values()) / 2.0
VL_MIN_FREQ = min(VL_FREQ_NORM.values()) / 2.0


def get_vh_frequency(gene_id: str) -> float:
    """获取 VH germline 的归一化使用频率（总和 = 1.0）"""
    gene_name = gene_id.split("*")[0] if "*" in gene_id else gene_id
    return VH_FREQ_NORM.get(gene_name, VH_MIN_FREQ)


def get_vk_frequency(gene_id: str) -> float:
    """获取 VK germline 的归一化使用频率（总和 = 1.0）"""
    gene_name = gene_id.split("*")[0] if "*" in gene_id else gene_id
    return VK_FREQ_NORM.get(gene_name, VK_MIN_FREQ)


def get_vl_frequency(gene_id: str) -> float:
    """获取 VL germline 的归一化使用频率（总和 = 1.0）"""
    gene_name = gene_id.split("*")[0] if "*" in gene_id else gene_id
    return VL_FREQ_NORM.get(gene_name, VL_MIN_FREQ)


def get_frequency(chain_type: str, gene_id: str) -> float:
    """获取 germline 的归一化使用频率（0-1，各链型总和 = 1）"""
    if chain_type == "H":
        return get_vh_frequency(gene_id)
    elif chain_type == "K" or gene_id.startswith("IGKV"):
        return get_vk_frequency(gene_id)
    else:
        return get_vl_frequency(gene_id)


# Adimab 推荐的 germline（基于其酵母展示库设计）
ADIMAB_RECOMMENDED_VH = {
    "IGHV1-69",  # 高频使用，适合噬菌体展示
    "IGHV3-23",  # 最常用，平衡的属性
}

ADIMAB_RECOMMENDED_VK = {
    "IGKV1-39",  # 高频使用
}

ADIMAB_RECOMMENDED_VL = {
    "IGLV3-1",   # Lambda 轻链
}

# Pioneer 库推荐的 germline（基于近 600 个临床阶段抗体分析）
PIONEER_RECOMMENDED_VH = {
    "IGHV1-69",  # 高频使用，适合展示
    "IGHV3-23",  # 最常用
}

PIONEER_RECOMMENDED_VK = {
    "IGKV1-39",  # 高频使用
}

PIONEER_RECOMMENDED_VL = {
    "IGLV3-1",   # Lambda 轻链
}


def is_recommended(chain_type: str, gene_id: str, source: str = "adimab") -> bool:
    """检查是否为推荐的 germline"""
    gene_name = gene_id.split("*")[0] if "*" in gene_id else gene_id
    
    if source == "adimab":
        if chain_type == "H":
            return gene_name in ADIMAB_RECOMMENDED_VH
        elif gene_name.startswith("IGKV"):
            return gene_name in ADIMAB_RECOMMENDED_VK
        else:
            return gene_name in ADIMAB_RECOMMENDED_VL
    elif source == "pioneer":
        if chain_type == "H":
            return gene_name in PIONEER_RECOMMENDED_VH
        elif gene_name.startswith("IGKV"):
            return gene_name in PIONEER_RECOMMENDED_VK
        else:
            return gene_name in PIONEER_RECOMMENDED_VL
    
    return False
