# Pipeline 差距分析：评分系统与回复突变设计的不足

> 基于对 `backmut.py`、`config.py`、`minimal.py`、`graft.py`、`structure.py`、`pipeline.py`、`variants.py`、`germline.py` 等核心模块的代码审查，从科学和事实角度进行的系统性分析。

## 一、回复突变评分系统的不足

### 1.1 结构评分用 `max()` 而非累加

`backmut.py:200` 中 structural score 取所有匹配特征的最大权重，而非求和：

```python
struct_scores = [w[f] for f in features if f in w]
structural = max(struct_scores) if struct_scores else 0.0
```

**问题**：
- 一个 `interface_core`（权重 1.0）直接封顶，无论该位置是否同时是 canonical、vernier、CDR-contact
- 一个多角色位置（如 H39：interface_core + canonical + vernier）与单一 interface_core 位置得分相同
- 实际上，同时承担多重结构角色的位置比单一角色位置更重要，当前算法无法区分

### 1.2 权重是静态的，不随抗体结构变化

所有权重（interface_core=1.0, vernier=0.85, canonical=0.80...）是固定常数（`config.py:55-66`）。

**问题**：
- 不同抗体的 VH/VL interface 几何差异很大，某些 interface_core 位置可能实际不埋藏
- CDR canonical class 不同（如 H1 type 1 vs type 3），支撑残基的重要性不同
- 缺少基于实际 3D 结构的动态权重调整（如 Rosetta interface energy 或 FoldX ΔΔG）

### 1.3 化学评分仅做序列 motif 扫描，不考虑暴露

`backmut.py:227-282` 的化学风险评分在 3-residue 窗口内检测 deamidation/isomerization/oxidation motif。

**问题**：
- 一个埋藏的 NG motif（不可及溶剂）几乎不发生脱酰胺，但仍被计入高分
- 一个表面暴露的 Met 氧化风险远高于埋藏 Met，但权重相同（0.30）
- 缺少基于结构的溶剂可及性加权（SASA-weighted chemical risk）

### 1.4 免疫原性 benefit 评分过于简化

当前公式（`backmut.py:220`）：

```
benefit = 0.3 + 0.5 * exposure * (1 - conservation)
```

**问题**：
- **无 T 细胞表位预测**：免疫原性的核心驱动因素是 MHC-II 呈递 + T 细胞激活，当前仅用"表面暴露 + 稀有性"代理，完全忽略了 MHC-II binding affinity 预测（如 NetMHCIIpan）、T 细胞表位免疫原性评分（如 IEDB immunogenicity score）
- **conservation 计算有偏**：用 top-20 germline 的携带率，但不同 V 基因家族的多态性差异很大
- **exposure 粗糙**：只有 0.85/0.15/0.5 三档，无连续 SASA 值

## 二、位置分类体系的不足

### 2.1 静态位置集来自 1985-2003 年文献

| 位置集 | 来源年份 | 问题 |
|--------|---------|------|
| vernier zone | Foote & Winter 1992 | 基于有限的抗体结构集（~10个），可能遗漏位置 |
| interface core | Chothia 1985 | 基于更早期的结构，VH/VL packing 理解不完整 |
| canonical | Chothia & Lesk 1987 | CDR canonical 分类已扩展（Martin & Thornton 1996, North et al. 2011），但位置集未更新 |

**已知遗漏**：`config.py` 中 VL vernier zone 缺少 L37 和 L45（`theory.md:27-31` 已注明），这两个位置在 Foote & Winter 原文中有列出。

### 2.2 接口位置集可能不适用于所有抗体

Chothia interface core 是"平均化"的 packing 集，但：
- 某些抗体的 VH/VL interface 有独特几何（如 llama VHH 的 extended interface）
- 某些"interface core"位置在特定抗体中实际暴露或远离 interface
- 缺少基于实际结构的 per-antibody interface 分析

## 三、Germline 选择的不足

### 3.1 无结构兼容性评估

当前 germline 选择（`germline.py:338-389`）仅基于序列 identity（FR/CDR/CVI），不考虑：
- 选定 germline 的 FR 与 donor CDR 的结构兼容性（某些 germline FR 可能导致 CDR 扭转）
- CDR3-J junction 的长度和序列偏好
- 某些 V 基因家族与特定 CDR canonical class 的适配性

### 3.2 CVI 同源性是必要但不充分条件

BI 2024 论文表明 CVI 与表达量/亲和力相关，但：
- CVI 只覆盖 ~40% 的 FR 位置，剩余 60% 的非 CVI 位置也可能影响结构
- CVI 相同的两个 germline，非 CVI 位置可能差异很大
- 缺少综合考虑整个 FR 的结构稳定性评估

## 四、变体组装的不足

### 4.1 线性 V0→V1→V2→V3 梯度缺乏组合优化

当前变体是简单的 tier 叠加（`variants.py:49-74`）：
- V1 = T1 only
- V2 = T1 + T2
- V3 = T1 + T2 + top T3

**问题**：
- 某些 T2 位置可能相互冲突（如同时回突变两个相邻 interface 位置可能破坏 packing）
- 某些 T3 位置在特定组合下可能比 T2 位置更重要
- 缺少基于结构能量的组合优化（如 Rosetta ΔΔG scanning）
- 缺少考虑位置间协同/拮抗效应

### 4.2 T3 选择仅按 composite 降序 + 暴露筛选

`variants.py:59-65` 中 V3 的 T3 选择逻辑：
```python
t3_exposed = [c.position for c in sorted(
    (c for c in backmut.candidates if c.tier == "T3" and (c.buried is False or c.buried is None)),
    key=lambda c: (-c.composite, c.position),
)][:extra_t3_max]
```

- 无结构数据时，`buried is None` 使所有 T3 都进入候选池
- 仅按 composite 排序，不考虑 T3 位置间的独立性

## 五、验证体系的不足

### 5.1 金标准回测用例太少

- `EMPIRICAL_NO_EFFECT` 仅含 VL L87 一个位点
- `backtest.py` 仅验证 3 个案例（4D5→trastuzumab, A4.6.1→bevacizumab, VHH gold standard）
- `backtest_scale.py` 用 25 个获批抗体，但缺乏系统性的对照实验
- 缺少与已发表的人源化抗体临床数据的系统对比

### 5.2 无负面案例验证

- 未测试"人源化失败"的案例（如某些抗体人源化后亲和力严重丢失）
- 未验证 pipeline 是否能正确识别"不应回突变"的位置
- 缺少假阳性/假阴性率的量化评估

## 六、缺失的现代方法

### 6.1 无机器学习/深度学习评分

当前完全基于规则 + 文献权重，未使用：
- 基于结构的 ML 亲和力预测（如 IgG lexer, antibody-mutation-effect predictors）
- 基于序列的免疫原性预测模型
- 大规模抗体-抗原复合物数据训练的 scoring function

### 6.2 无分子动力学/柔性考虑

- 所有分析基于静态结构（或无结构）
- CDR loop 的动态构象变化未被考虑
- 某些"稳定"的 back-mutation 可能影响 loop flexibility

### 6.3 可开发性预测仅限于序列 motif

`developability.py` 和化学评分仅检测序列 pattern，缺少：
- 基于结构的聚集倾向预测（如 COLIBRI, Aggrescan3D）
- 表面电荷分布分析（net charge, charge patches）
- 热稳定性预测（如 ThermoNet, ProTherm）
- 粘度预测（与 surface hydrophobicity 相关）

---

## 七、Step 3（AF3 结构模式）的改善与遗留

### 7.1 Step 3 能显著改善的

| 原始问题 | Step 3 如何改善 | 改善程度 |
|---------|----------------|---------|
| exposure 只有三档 0.85/0.15/0.5 | FreeSASA 提供连续 relSASA 值，且 0.15-0.25 标记为 uncertain（`structure.py:304-307`） | **显著** |
| 埋藏/暴露误判导致 tier 错误 | `buried=True/False` 精确覆盖静态位置集，interface_core + exposed → T2 而非 T1 | **显著** |
| 无 CDR 接触数据 | `cdr_contact` 检测（< 4.5 Å）→ 高优先级 T2 回突变 | **显著** |
| 无抗原接触数据 | `antigen_contact` 检测 → 评估回突变是否直接影响结合 | **有** |
| pLDDT 低置信度区域 | pLDDT < 50 → structural × 0.3, < 70 → × 0.7，避免不可靠结构误导 | **有** |
| 单模型噪声 | 多模型共识（rank_1-5, min_consensus=3）减少假阳性 | **有** |
| Vmin 退化为 T1-only | 有结构时用 greedy set-cover 算法精确选择最少回突变覆盖所有 CDR 接触 | **显著** |
| V3 无暴露筛选 | 有结构时 `buried=False` 才进入 T3 候选 | **有** |

### 7.2 Step 3 仍然解决不了的

| 问题 | 为什么 AF3 帮不上 | 影响 |
|------|----------------|------|
| structural score 用 max() 而非累加 | max() 逻辑在有结构数据时不变，多角色位置仍无法叠加 | **高** |
| 免疫原性无 T 细胞表位预测 | AF3 只预测结构，不预测 MHC-II 呈递或 T 细胞激活；benefit 公式不变 | **高** |
| 化学风险不考虑 SASA | 化学评分仍做纯序列 motif 扫描，不乘以 relSASA | **中** |
| 位置集是静态文献值 | AF3 结构只覆盖 buried/cdr_contact/antigen_contact，不更新文献集本身 | **中** |
| germline 选择无结构兼容性 | germline 选择在结构预测之前完成，不考虑 FR-CDR 结构适配 | **中** |
| 变体组装是线性叠加 | V0→V1→V2→V3 仍是 tier 叠加，不考虑位置间协同/拮抗 | **中** |
| 金标准回测用例少 | AF3 结构数据不增加验证用例 | **中** |
| 无 ML scoring | AF3 是结构预测模型，不是亲和力/免疫原性评分模型 | **低**（短期） |
| 无 MD/柔性考虑 | AF3 输出静态结构，CDR loop 构象动态不变 | **低** |

### 7.3 结论

**Step 3 本质上解决的是"输入数据质量"问题**——让 buried/contact/exposed 的判断从静态文献猜测变为结构实测。这对 tier 分配和 benefit 评分的准确性有显著帮助。

**但它不解决"评分函数本身"的问题**——max() 逻辑、免疫原性公式、化学评分方式、变体组装策略等算法层面的缺陷，在有结构数据时依然存在。

可以这样理解：
- **Step 3 无结构** → 数据差 + 算法差 → 双重劣势
- **Step 3 有结构** → 数据好 + 算法差 → 仍有改进空间
- **理想状态** → 数据好 + 算法好 → 需要同时改进评分函数

---

## 八、优先改进方向

| 优先级 | 改进 | 科学依据 | 实现难度 |
|--------|------|---------|---------|
| **高** | 用实际结构的 SASA 替代静态 exposure 三档 | 消除暴露/埋藏误判 | 中（需 FreeSASA） |
| **高** | 累加结构评分替代 max() | 反映多重角色的叠加重要性 | 低 |
| **高** | 补全 VL vernier zone（L37, L45） | 修复已知遗漏 | 低 |
| **中** | 引入 MHC-II T 细胞表位预测 | 免疫原性的核心机制 | 中（需 NetMHCIIpan） |
| **中** | 化学风险乘以 SASA 权重 | 埋藏 motif 不应计同等风险 | 中 |
| **中** | 位置间协同效应评估 | 回突变组合可能冲突 | 高（需结构能量计算） |
| **低** | ML scoring function | 数据驱动替代规则 | 高（需训练数据） |
| **低** | 负面案例回测 | 量化假阳性/假阴性率 | 中 |
