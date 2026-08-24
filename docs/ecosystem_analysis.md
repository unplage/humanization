# Antibody humanization ecosystem analysis（同类项目全景分析）

全网（GitHub 为主）抗体人源化/设计工具调研：本项目的定位、差距与可借鉴/可集成
的项目。调研日期：2026-08。

---

## 1. 同类项目清单与定位

### 1.1 直接同类（人源化流程/平台）

| 项目 | 组织 | ⭐ | 定位 | 核心技术 |
|------|------|-----|------|----------|
| **BioPhi** | Merck（原 oxpig） | 263 | 开源抗体设计平台（人源化+人源度评估） | Sapiens 抗体语言模型（逐位点人源化概率）+ OASis 人源度（OAS 9-mer 数据库）+ 移植/回复突变建议 |
| **HuDiff** | TencentAI4S | 94 | 抗体+纳米抗体人源化扩散模型 | HuDiff-Ab / HuDiff-Nb（自回归扩散 + 结合位点保留），提供 HumAb25 基准 |
| **AnthroAb** | nagarh | 1 | 人源化语言模型 | LM 生成式人源化 |
| **Llamanade** | Lefrunila | 2 | 开源纳米抗体人源化管线（归档恢复） | VHH 特异人源化 |
| **antibody-humanization-tool / AntibodyHumanization** | 个人 | 1 | ML 人源化 | — |
| **GraftJ / autograft** | 个人 | 0-3 | CDR 环移植工具 | 结构环移植 |
| **immunocheck** | salinas2000 | 0 | VHH 免疫原性/可开发性 API | MHC-II 表位 |

### 1.2 支撑工具（编号/结构/数据）

| 项目 | ⭐ | 用途 | 与我们的关系 |
|------|-----|------|--------------|
| ANARCI / abnumber | — | 抗体编号（IMGT/Chothia/Kabat） | 已作为服务器端交叉验证 |
| **SAbDab**（牛津 OPIG） | — | 结构抗体数据库 | 回测金标准数据源（本报告 1N8Z/1CZ8 等均来自 PDB/SAbDab） |
| **SNAC-DB**（Sanofi） | 35 | 结构纳米抗体数据库 + 管线 | VHH 验证数据源 |
| **abdev-benchmark**（Ginkgo） | 34 | 抗体可开发性基准 | 可开发性评估对标 |
| DeepViscosity / antibody_developability | 28/9 | 高浓度粘度/可开发性预测 | 后续扩展 |

### 1.3 语言模型 / 生成模型（人源度与序列生成）

| 项目 | ⭐ | 说明 | 采用建议 |
|------|-----|------|----------|
| **Sapiens**（Merck） | 80 | 人类抗体 BERT 语言模型，BioPhi 人源化引擎 | **已集成适配器**（`--biophi-env`） |
| **IgLM**（Graylab） | 194 | 抗体序列生成 LM（infilling） | 备选生成路线 |
| **AbLang**（oxpig） | 170 | 抗体 LM（掩码重建） | 备选生成/修复 |
| **abmap / AntiBERTy** | 155/80 | 抗体表示学习 | 特征嵌入（远期） |
| **HuDiff**（腾讯） | 94 | 扩散人源化（含 VHH） | 生成式替代路线；基准数据已采用 |

### 1.4 paratope / 界面预测（支撑 V_SDR）

| 项目 | ⭐ | 说明 | 采用建议 |
|------|-----|------|----------|
| **Parapred** | 62 | 深度学习的 paratope 预测（仅序列） | 无抗原结构时的 V_SDR 备用 |
| PECAN / ParaSurf / Paragraph | 35/29/7 | 图神经网络 paratope | 同上（对标） |

---

## 2. 借鉴与已落地的改进

| 借鉴来源 | 内容 | 本项目落地 |
|----------|------|------------|
| **BioPhi**（Sapiens/OASis） | 数据驱动的人源度评估（非仅 germline 匹配） | `humanness.py` 适配器：`--biophi-env` 运行 `biophi sapiens --mean-score-only`，报告逐变体 Sapiens 均分 + OASis 身份（22GB DB 可选，服务器模式） |
| **HuDiff / HumAb25** | 25 个真实人源化药物亲本序列基准 | `data/benchmarks/humab25_parental_mouse.csv` + `tests/backtest_scale.py`：25 例规模化验证全过（编号/选架/变体生成 100%） |
| **HuDiff 评估口径** | FR 突变精确率（share_precision） | 已纳入回测指标体系（tests/backtest.py 的 precision/recall） |
| **BI 矩阵移植** | 框架组合多样性 | 已实现框架矩阵变体（minimal.py） |
| **CUMAb/AbEL** | 能量排序 | 路线图 P1（AF3 逐位点回复模拟/Rosetta ddG） |
| **abdev-benchmark** | 可开发性量化 | 路线图 P2（扩展 developability.py 为完整面板） |
| **SNAC-DB** | VHH 结构验证数据 | 路线图 P2（VHH 专用回测集） |
| **Parapred 类** | 无抗原时的 paratope 预测 | 路线图 P2（V_SDR 的序列级备用） |

---

## 3. 差距与竞争定位

### 3.1 相对 BioPhi（最强同类）

**BioPhi 优势**：
- Sapiens/OASis 的**数据驱动人源度**（OAS 千亿序列统计）显著强于我们的
  germline 匹配近似
- 完整 Web 平台 + 社区

**本项目优势**：
- **结构驱动**：AF3 接触对 + 最小回复集（set-cover）+ paratope 移植（V_SDR），
  BioPhi 的回复突变建议无此结构分辨率
- **实验闭环**：KD 数据 → 逐位点 ΔΔG 校准（learning.py），BioPhi 无此机制
- **VHH 全程支持**（hallmark 保护、VHH 参数）
- **可移植**：纯 stdlib 便携模式，BioPhi 需 22GB DB 或 conda 环境
- **框架矩阵 / germline hopping**（BI 2024 实证）

### 3.2 相对 HuDiff（生成式）

HuDiff 是"生成式替代路线"：不显式做回复突变，而是让扩散模型直接输出
人源化序列。优势：无需人工规则；劣势：黑盒、无逐位点依据、需 GPU 训练/部署。
**本项目定位**：可解释、结构可验证、实验可闭环的"白盒"路线；
HuDiff 输出可作为交叉验证/备选面板（路线图 P2 集成其 HumAb25 评测口径）。

---

## 4. 采用路线图（按价值排序）

| 优先级 | 采用项 | 状态 |
|--------|--------|------|
| P0 | HumAb25 规模化回测 | ✅ 已交付（25/25 通过） |
| P0 | BioPhi/Sapiens 人源度适配器 | ✅ 已交付（`--biophi-env`，服务器模式） |
| P1 | Sapiens **逐位点**概率矩阵回注评分（IMGT→Kabat 映射） | 待办（需 ANARCI 映射层） |
| P1 | promb（轻量 OASis，MSDLLCpapers）集成 | 待办（pip install promb） |
| P2 | Parapred 序列级 paratope（无抗原时启用 V_SDR） | 待办 |
| P2 | HuDiff 生成式备选面板（推理接口封装） | 待办 |
| P2 | SNAC-DB VHH 结构回测集 | 待办 |
| P2 | abdev-benchmark 指标面板（TM/粘度/HIC 预测） | 待办 |
| P3 | IgLM/AbLang infilling 做框架"再人源化" | 探索 |

---

## 5. 结论

1. **本项目在当前开源生态中的差异化优势**：结构接触驱动的回复突变设计 +
   实验数据闭环校准 + VHH 全支持 + 零依赖可移植 —— 在 GitHub 开源项目中
   （BioPhi 之外）没有同时具备这四者的。
2. **最大短板**：人源度评估停留在 germline 匹配层面；已通过 BioPhi/Sapiens
   适配器补齐（服务器模式），逐位点回注列为 P1。
3. **验证优势**：回测金标准（曲妥珠/贝伐珠/25 例 HumAb25）在同类开源项目中
   属于最严格的一档。
