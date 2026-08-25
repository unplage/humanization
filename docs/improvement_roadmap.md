# Humanization improvement roadmap & industry comparison
（人源化改进路线图与行业对比）

本文回答三个问题：
1. 当前版本（v0.1）在回复突变阶段还有哪些改进空间？
2. 如何以**最少的回复突变**保留**最大的亲和力**？（业界证据 + 具体算法）
3. 与行业方案（WeMol / BioPhi / CUMAb / 矩阵移植 / paratope 移植等）、专利与
   公开案例的对比。

> 已落地实现的部分（`scripts/humanize/minimal.py`，本版已交付）：
> CVI 同源性指标、最小回复突变集（set-cover）、框架矩阵变体、
> paratope（SDR）精确移植变体。其余列为路线图。

---

## 1. 当前版本（v0.1）的差距分析

### 1.1 已具备的能力

- 文献位置表驱动的分级（T1/T2/T3）：vernier / 界面 / canonical / VHH hallmark / Cys
- AF3 结构提示回注（buried / CDR 接触 / 抗原接触）
- 三维评分（结构 0.55 + 免疫原性 0.30 + 化学 0.15）
- V0–V3 梯度 + 本次新增：Vmin（最小回复集）、V_SDR（paratope 移植）、框架矩阵

### 1.2 主要差距（按影响排序）

| # | 差距 | 影响 | 业界对标 |
|---|------|------|----------|
| 1 | ~~单 germline 决策~~（已实现多策略）：9 种 germline 选择策略（FR/CDR/综合/CVI/最小回复突变/频率加权/3 轴）+ 框架矩阵变体 | 已通过多策略选择 + matrix 备选变体覆盖 | BI 矩阵移植、CUMAb 千框架能量排序 |
| 2 | **无能量/结构稳定性评分**：回复突变选择基于文献特征，不直接评估"该突变对 CDR 构象/界面/稳定性的贡献" | 可能多回/漏回突变；不能区分功能等价设计 | Rosetta ddG（CUMAb/AbEL）、MOE、AF3 能量代理 |
| 3 | **无免疫原性 ML 打分**：仅用 germline 匹配 + 暴露度代理 | 暴露的非人源残基不一定免疫原性最强；T 细胞表位（MHC-II）未建模 | BioPhi OAS 打分、SITA、netMHCIIpan、Abzena EpiScreen |
| 4 | **变体数量少（4-6 个）**：无法覆盖组合空间 | 亲和力丢失后需人工迭代 | 组合回复突变库 + 高通量筛选（BI 210 graft 矩阵、噬菌体/酵母展示） |
| 5 | **无 VH/VL 不对称处理**：VL 框架对亲和力的影响权重未单独体现 | VL 框架选择失误是亲和力丢失主因（BI：VL <60% 同源 → 98% 失活） | BI 的 VL 优先策略 |
| 6 | **接触集缺失**：之前没有 framework→CDR 接触对数据 | 无法做"恢复接触"的最小集计算（本版已补 cdr_partners） | — |
| 7 | **CDR3 内残基无结构级支持判断**：CDR3 全部照搬 | 部分 CDR3 尾段由 J 区提供更优 | J 交界精细处理 |
| 8 | **无 ML/生成式人源化** | 无法自动探索"非回复突变"式人源化路径 | HuAbDiffusion、语言模型人源化、ProteinMPNN 框架再设计（已接适配器） |

---

## 2. 如何用最少的回复突变保留最大的亲和力

### 2.1 业界证据链（为何"最少"可行）

1. **CVI 残基同源性是核心指标**（Boehringer Ingelheim，JBC 300:105555，2024）：
   对同一套鼠源 CDR 做 210 个 VH×VL 人类框架组合移植，发现：
   - 表达量与亲和力均与 **CVI（canonical+vernier+interface）残基同源性**显著相关；
   - **VL 框架同源性 < 60% 的移植 98% 丢失结合**；VH 同源性 30–70% 均可保留结合；
   - **最高同源性的移植（Vκ1:VH1）不是最终最优**——Vκ4:VH9、Vκ8:VH11 等
     中低同源组合在亲和力、功能、可开发性上全面更优。
   → 含义：回复突变的"目标"不是恢复所有差异残基，而是**保住 CVI 残基**；
     框架选择本身比回复突变数量更决定成败。

2. **能量排序优于同源性直觉**（CUMAb/AbEL，Weizmann，bioRxiv 2022；
   专利已申请）：在**数千个人类框架**上移植 CDR，用 Rosetta 原子模拟按能量与
   结构完整性排序。结论：
   - **非高同源框架常被选中**（与 BI 一致）；
   - 多个编码数十个互不相同突变的框架组合**功能等价** → 存在"等效解空间"，
     应输出多样性面板而非单点解。

3. **AI 精确 paratope 移植**（XtalPi Ailux，bioRxiv 2025）：不只移植整段 CDR，
   而是用 AI 复合物结构鉴定**关键 paratope 残基**，只移植这些残基 + 必要支柱，
   免疫原性显著降低且活性保持。
   → "最少回复突变"的极致形态：CDR 内部也做减法（本版 V_SDR 已实现骨架）。

4. **亲和力丢失的实测机理**：BI 的结构分析显示 HCDR3 RMSD 在 1.08–5.01 Å 之间
   波动，与同源性无关、与 VH/VL 框架配对相关 → CDR 构象主要由**框架支柱残基**
   决定，与"vernier 位点"理论完全一致，也说明回复突变应聚焦**接触 CDR 的
   框架残基**（本版 Vmin 的 set-cover 即按此设计）。

### 2.2 具体算法（已实现 / 路线图）

**A. 最小回复突变集（已实现，Vmin）**
```
问题：在 donor 结构中，候选回复位置 p 的 donor 残基与若干 CDR 残基有接触
      (<4.5 Å)。若 p 保留人源残基，这些接触丢失。
求解：集合覆盖。universe = {(p, cdr) 接触对 | p 是候选};
      p 覆盖其全部接触对；贪心选择覆盖未覆盖接触对最多的 p；
      平局按 composite 分数。
无结构时回退：Tier-1（文献支柱）即最小安全集。
```
- 输出：Vmin 变体 + "保留接触 x/y" 指标
- 扩展（路线图）：把 universe 从 framework→CDR 接触扩展到
  framework→antigen 接触（复合物结构）与 VH/VL 界面双侧接触（成对回复）。

**B. 框架矩阵 / germline hopping（已实现，报告输出）**
- 对 top-3 备选 germline 各生成一套 T1+T2 回复突变与 V2 类变体；
- 每套报告 CVI 同源性 + 回复突变数；
- 扩展（路线图）：按 BI 结论给 **VL 的 CDR 同源性更高权重**（VL 驱动亲和力），
  按 VH 的 CVI 同源性预测可开发性；输出 3×3 VH×VL 矩阵。

**C. paratope（SDR）精确移植（已实现，V_SDR，需抗原复合物）**
- 移植残基 = 抗原接触的 CDR 残基 + T1/T2 支柱；
- 其余 CDR 残基回人源 germline（CDR 内部"减法"）；
- 扩展（路线图）：若 donor 无复合物结构，用 AF3 复合物预测并校验
  移植前后接触保留率 ≥ 90%。

**D. 能量/结构评分（路线图，优先级最高）**
- 选项 1（轻量，无新依赖）：AF3 批量预测每个单点回复突变后的变体，
  计算 (a) CDR 环 CA-RMSD vs donor，(b) VH/VL 界面面积变化，
  (c) paratope 接触保留率 → 作为每位的"回复收益"分；
- 选项 2（重）：接入 Rosetta `cartesian_ddg` / Flex ddG 逐位打分
  （CUMAb 方法学），或 AbEL 的公开流程（https://github.com/Fleishman-Lab）；
- 决策规则：当 p 回复后 CDR RMSD 无改善且无接触恢复 → 降级为"保持人源"。

**E. 免疫原性 ML 打分（路线图）**
- 接入 BioPhi 的 OAS 人源化打分（per-position humanness）；
- 或 netMHCIIpan 4.0 对 9-mer/15-mer 肽段预测 MHC-II 结合（限服务器，
  与 germline 匹配结合成"暴露 × 表位"复合免疫原性分）；
- VHH 特例：FR2 hallmark 区被证实在 T 细胞激活中贡献表位
  （Giraudet et al., MAbs 2025），人源化时应评估"部分 germline 化"
  （如 37V、44G、45L、47W 逐个尝试的保守变体面板）。

**F. 组合变体设计（路线图）**
- 固定 T1，对 T2/T3 用贪心/遗传算法按结构指标搜索 20–50 个组合变体，
  输出"筛选面板"（酵母/噬菌体展示或 96 孔小试）；
- 与 BI 的"矩阵移植 + liability 修复"流程对齐。

---

## 3. 行业方案对比

| 方案 | 主体 | 核心方法 | 回复突变策略 | 对本文的启示 |
|------|------|----------|--------------|--------------|
| **WeMol**（药明康德，2023 上线） | 药明康德 | Web 抗体设计平台：IMGT/OAS 数据 + 结构建模 + 可开发性预测；人源化模块 = CDR 移植 + 回复突变 + 表位风险与 developability 报告 | 经典移植为主，平台整合结构/AI 评分；结合药明湿实验（表达/亲和力/可开发性）闭环反馈 | 平台化、报告化、与湿实验闭环是差距 |
| **BioPhi**（Prihoda et al., PLoS Comput Biol 2022） | OPIG 开源 | OAS 人类抗体库训练的逐位点人源化打分；germline 选择用 OAS 覆盖度 | 建议回复突变 = 高 humanness 优先 + 结构常识位点 | **逐位点 humanness 打分**可替代纯 germline 匹配，建议作为交叉验证集成 |
| **CUMAb / AbEL**（Weizmann，2022） | 学术 | 数千框架移植 + Rosetta 能量排序 | "能量 + 结构完整性"联合选择框架与突变，几乎不依赖同源性 | 能量评分框架；非高同源框架优先 |
| **矩阵移植**（Boehringer Ingelheim，JBC 2024） | 药企 | 210 组合（14 VL × 15 VH）湿实验矩阵 + PCA 选优 | 回复突变被"框架组合多样性"取代；CVI 同源性作 KPI | CVI 指标 + VH/VL 不对称 + 组合面板 |
| **Paratope 移植**（XtalPi Ailux，bioRxiv 2025） | CRO/AI 公司 | AI 复合物结构 → 只移植关键 paratope 残基 | CDR 内做减法，最小免疫原性 | 本版 V_SDR 对应此路线 |
| **HuAbDiffusion**（Brief Bioinform 2025） | 学术 | 离散语言扩散模型直接生成人源化序列 | 无需人工回复突变规则 | 生成式人源化可作为替代路线（需大量验证） |
| **SITA**（J Pharm Anal 2025） | 学术 | 位点特异性免疫原性预测（MHC-II 表位） | 无（评估工具） | 免疫原性评分集成 |
| **Abzena Composite Human Antibodies** | CRO | 多个人源 germline 片段拼接的"复合框架"，先验去除 T 细胞表位 | 框架层面避免回复突变需求 | "免回复突变"框架构建思路 |
| **IgBLAST -humanize / OAS 社区** | NCBI/社区 | 最接近人源 germline 搜索 | — | 交叉验证 |

> 注：WeMol 为商业平台，具体算法细节未公开；以上为公开资料（官网/会议）层面的描述。

---

## 4. 专利与文献要点（自由实施注意事项）

| 专利/文献 | 覆盖内容 | 对管线的影响 |
|-----------|----------|--------------|
| **Queen 等人源化抗体专利**（US 5,530,101 等，Protein Design Labs） | CDR 移植 + 框架残基替换标准（影响 CDR 构象的位点、界面位点等）——即回复突变方法的专利根基 | 方法论属经典实践；**管线自身不构成侵权**，但产品化路径应做 FTO |
| **Adair/Celltech 人源化专利**（US 5,585,089 等） | 变体框架替换策略 | 同上 |
| **SDR grafting / resurfacing 专利族**（Padlan 相关工作） | 只移植 SDR 或仅表面残基人源化 | V_SDR 变体的设计空间提示 |
| **CUMAb 方法专利**（Tennenhouse, Fleishman 等申请） | 能量排序的人源化流程 | 若采用 Rosetta 能量排序，注意与申请中权利要求的边界 |
| **VHH 人源化相关专利**（多家 nanobody 公司） | VHH hallmark 保留策略 | 本管线默认保护 hallmark，符合业界共识 |
| **BI anti-LAG3 案例专利**（US20170334995A1） | 矩阵移植具体序列 | 案例学习 |

> 建议：任何进入临床的候选都应在法务指导下做 FTO 检索；本文不构成法律意见。

---

## 5. 经典案例对照（验证我们设计逻辑的合理性）

| 案例 | 人源化方式 | 关键回复突变 | 与本管线对照 |
|------|------------|--------------|--------------|
| **曲妥珠单抗**（4D5 → trastuzumab，Carter 1992） | CDR 移植 + 回复突变 | VH 48/67/69/71/73；VL 4/66/87（本管线 T1/T2 全部命中，见 tests） | T1/T2 分级与金标准一致 |
| **帕利珠单抗**（MEDI-493） | 框架替换 + 回复突变 | 界面/vernier 位点 | 同上 |
| **贝伐珠单抗**（A4.6.1 人源化） | CDR 移植 | vernier 支柱 | 同上 |
| **anti-LAG3 矩阵移植**（BI 2024） | 210 组合矩阵 | 无传统回复突变，靠框架组合 | 本版框架矩阵 + CVI 指标对应 |
| **VHH germline 化**（Giraudet et al., MAbs 2025） | FR2 部分 germline 化降低 T 细胞表位 | FR2 表位区 | 提示 VHH 的 FR2 需"分步"人源化评估 |

---

## 6. 优先路线图（P0→P2）

| 优先级 | 改进 | 交付物 | 依赖 |
|--------|------|--------|------|
| **P0（本版已交付）** | CVI 同源性指标 | 报告 + JSON | 无 |
| **P0（本版已交付）** | 最小回复突变集（set-cover） | Vmin 变体 + 接触保留指标 | AF3（无结构时回退 T1） |
| **P0（本版已交付）** | 框架矩阵变体 | 备选 germline 的 V2 面板 + CVI 对比 | 无 |
| **P0（本版已交付）** | Paratope 精确移植 | V_SDR 变体 | 抗原复合物结构 |
| **P1** | 能量评分：AF3 逐位点回复模拟（CDR RMSD/界面/接触保留） | 每位"回复收益"分 + 决策规则 | AF3（服务器） |
| **P1** | BioPhi OAS humanness 集成 | 逐位点 humanness 交叉验证 | pip install biophi（GitHub） |
| **P1** | VH/VL 不对称权重 | germline 选择加权（VL 驱动亲和力、VH 驱动可开发性） | 无 |
| **P1** | VHH FR2 分步人源化面板 | VHH 专属变体（37/44/45/47 单点尝试） | 无 |
| **P2** | 组合变体筛选面板（20–50 个） | 筛选库 + 排序 | 计算资源 |
| **P2** | netMHCIIpan / SITA 集成 | T 细胞表位分 | 服务器 |
| **P2** | 生成式人源化（HuAbDiffusion/语言模型） | 替代路线候选 | 模型部署 |

---

## 7. 参考文献

1. Gupta P, et al. Matrixed CDR grafting: A neoclassical framework for antibody
   humanization and developability. *J Biol Chem* 300:105555 (2024)
2. Tennenhouse A, et al. Reliable energy-based antibody humanization and
   stabilization (CUMAb). *bioRxiv* 2022.08.14.503891 (2022)
3. Wang H, et al. AI-Guided Precision in Antibody Humanization. *bioRxiv*
   2025.03.09.641110 (2025, XtalPi/Ailux)
4. Prihoda D, et al. BioPhi: A platform for antibody design, humanization, and
   humanness evaluation. *PLoS Comput Biol* 18:e1009884 (2022)
5. Liu D, et al. HuAbDiffusion: discrete language diffusion model for antibody
   humanization. *Brief Bioinform* 26:bbaf658 (2025)
6. Cun Y, et al. SITA: Predicting site-specific immunogenicity for therapeutic
   antibodies. *J Pharm Anal* 15:101316 (2025)
7. Giraudet R, et al. Immunogenicity of single-chain antibodies: germlining of a
   VHH lowers T-cell activation from epitopes in FR2 and CDR regions. *MAbs*
   17:2571406 (2025)
8. Carter P, et al. Humanization of an anti-p185HER2 antibody. *PNAS* 89:4285 (1992)
9. Queen C, et al. US 5,530,101 (humanized immunoglobulins)
10. Boehringer Ingelheim. US20170334995A1 (anti-LAG3 antibodies)
