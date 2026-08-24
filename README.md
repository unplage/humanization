# Antibody Humanization Pipeline · 抗体人源化流程

**Production-grade CDR grafting + back-mutation (回复突变) design for Fab and VHH antibodies.**

一个可直接进入生产环节的抗体人源化系统：从非人源（鼠/兔/驼等）抗体序列出发，
完成人源 germline 框架匹配 → CDR 移植（4 套 CDR 定义）→ **回复突变候选设计与
量化评分** → 变体梯度设计（V0–V3）→ （服务器端）AlphaFold3 结构验证 +
ProteinMPNN 框架再设计 → 实验验证 SOP（基因合成 → 表达 → 亲和力 → 可开发性 →
Go/No-Go 判定）。

本机（便携模式）**零安装、零系统依赖**（纯 Python 标准库即可端到端运行）；
服务器端可启用 ANARCI 精确编号交叉验证、IgBLAST germline 数据库、
AlphaFold3 结构评分与 ProteinMPNN 设计。

---

## 目录

- [特性](#特性)
- [架构总览](#架构总览)
- [快速开始](#快速开始)
- [输入输出](#输入输出)
- [管线各阶段详解](#管线各阶段详解)
- [回复突变评分体系](#回复突变评分体系)
- [变体梯度](#变体梯度)
- [结构模式（服务器）](#结构模式服务器)
- [验证与金标准](#验证与金标准)
- [实验验证 SOP](#实验验证-sop)
- [服务器部署](#服务器部署)
- [FAQ](#faq)
- [引用文献](#引用文献)
- [许可证与致谢](#许可证与致谢)

---

## 特性

| 特性 | 说明 |
|------|------|
| Fab / VHH 双格式 | 自动检测链型（VH+VL 或单链 VHH）；VHH hallmark（Kabat 37/44/45/47）识别并全程保护 |
| 零依赖便携核心 | anchor-based Kabat 编号引擎为纯标准库实现；germline 数据内置（373 V + 28 J 人源基因） |
| 多套 CDR 定义 | Kabat / Chothia / AbM / IMGT 边界同时报告；变体梯度可任选一套 |
| 回复突变量化评分 | 结构（vernier/界面/canonical/接触/埋藏）+ 免疫原性 + 可开发性三维评分，输出 T1/T2/T3/KEEP_DONOR 分级与依据 |
| 变体梯度 | V0 纯移植 → V1(T1) → V2(T1+T2，推荐) → V3(+暴露 T3) |
| 结构模式 | AlphaFold3 预测 → 埋藏度/CDR 接触/抗原接触回注评分；变体 CDR 环 RMSD、界面与抗原接触保留率验收 |
| 可开发性检查 | N-糖基化 / 脱酰胺 / 异构化 / 氧化风险基序自动扫描 |
| 完整报告 | Markdown + JSON + 逐位点 CSV + variants FASTA |
| 实验 SOP | 基因合成 → Expi293 表达 → BLI/SPR → 可开发性面板 → Go/No-Go 标准（6–10 周） |

---

## 架构总览

```
humanization/
├── .opencode/skills/antibody-humanization/SKILL.md   # opencode skill 入口（工作流编排）
├── scripts/
│   ├── humanize/                 # 核心 Python 包
│   │   ├── numbering.py          # ★ 自研 anchor-based Kabat 编号引擎（零依赖）
│   │   ├── sequences.py          # FASTA 解析、链型分类、VHH 检测
│   │   ├── germline.py           # germline 加载/选择（NCBI FASTA 优先，内置 JSON 后备）
│   │   ├── graft.py              # CDR 移植（多方案）+ 变体序列组装
│   │   ├── backmut.py            # ★ 回复突变候选生成与评分（核心）
│   │   ├── variants.py           # V0–V3 变体梯度
│   │   ├── structure.py          # AlphaFold3 适配器 + 无依赖 PDB/CIF 解析 + 结构指标
│   │   ├── mpnn.py               # ProteinMPNN 适配器（固定 CDR/界面再设计框架）
│   │   ├── developability.py     # 序列级可开发性风险扫描
│   │   ├── report.py             # Markdown/JSON/CSV/FASTA 报告
│   │   ├── pipeline.py           # 端到端编排
│   │   ├── cli.py                # 命令行入口（run / setup-germline / setup-check）
│   │   └── config.py             # 领域知识表（位置集合）与评分权重
│   ├── build_germline_kabat.py   # 从 abnumber IMGT 数据重建内置 germline
│   ├── install.sh                # 服务器一键安装（--full 启用 ANARCI/IgBLAST）
│   └── requirements.txt
├── data/
│   ├── germline/                 # 内置人源 germline（Kabat 坐标 JSON + 来源说明）
│   └── examples/                 # 示例：鼠源 4D5 Fab + 驼源 1BZQ VHH
├── docs/
│   ├── theory.md                 # 回复突变理论基础（位置表 + 文献）
│   ├── scoring.md                # 评分公式与权重依据
│   ├── tool_setup.md             # 服务器工具部署
│   ├── experimental_SOP.md       # 实验验证全流程
│   └── validation.md             # 验证方案（编号/金标准/结构验收）
└── tests/test_pipeline.py        # 回归测试（含 4D5→曲妥珠单抗金标准）
```

---

## 快速开始

### 本机（零安装）

```bash
# 1. 人源化一个抗体（Fab：VH+VL 写入同一 FASTA；或单个 VHH）
python3 scripts/humanize/cli.py run \
  --input data/examples/mouse_4d5_fab.fasta \
  --outdir outputs

# 2. 查看报告（germline 选择、回复突变表、变体序列）
open outputs/humanization_report.md

# 3. 运行回归测试
python3 tests/test_pipeline.py
```

示例（VHH）：

```bash
python3 scripts/humanize/cli.py run \
  --input data/examples/cab_rn05_vhh.fasta \
  --outdir outputs_vhh
```

### 服务器（完整模式）

```bash
bash scripts/install.sh --full          # 安装 ANARCI / IgBLAST / NCBI germline
conda activate humanize
python3 scripts/humanize/cli.py run --input seq.fasta --outdir outputs \
  --af3-mode local --af3-binary /path/run_alphafold.py \
  --mpnn-mode local --mpnn-script /path/protein_mpnn.py
```

工具状态自检：`python3 scripts/humanize/cli.py setup-check`

---

## 输入输出

### 输入

- **FASTA 文件**，每条记录一个 V 结构域：
  - Fab：`>VH` + `>VL` 两条记录（顺序任意）
  - VHH：单条记录（自动按 hallmark 检测驼源）
- 可选：`--antigen SEQ`（AF3 复合物预测）、`--donor-structure PDB`（CDR RMSD 参考）

### 输出（`--outdir`）

| 文件 | 内容 |
|------|------|
| `humanization_report.md` | 完整报告：germline 选择与备选、逐位点回复突变表（分级/分数/特征/依据）、变体序列、tier 汇总 |
| `humanization_result.json` | 机器可读完整结果 |
| `backmutations_<chain>.csv` | 逐位点评分明细 |
| `variants.fasta` | V0–V3 全部变体序列 |

---

## 管线各阶段详解

### 阶段 1：编号（Numbering）

自研 **anchor-based Kabat 编号引擎**（`numbering.py`），纯标准库实现，基于结构上不变的锚残基：

```
VH:  Cys22 (FR1) · Trp36 (FR2) · Cys92 (FR3) · Trp103 (FR4)
VL:  Cys23 (FR1) · Trp35 (FR2) · Cys88 (FR3) · Phe98 (FR4)
```

- 框架区长度固定，可变环（CDR1/2/3）按标准插入字母规则（35A/B、52A-D、82A-C、100A-K、27A-E…）编号
- 已针对 VH3-23/VH1/mouse VH1/VHH/kappa/lambda 逐区验证（`tests/`）
- 服务器端用 **ANARCI**（HMMER）做精确模式交叉验证（`docs/validation.md`）

### 阶段 2：Germline 选择（IgBLAST 数据源）

- **数据源**：NCBI IgBLAST 人类 germline FASTA（`humanize setup-germline` 下载，管线优先读取）；未下载时自动回退到内置的 373 V + 28 J 人源基因（Kabat 坐标 JSON，与 IgBLAST 同源的 IMGT 数据）
- **选择策略**（业界"hybrid"方案，Olimpieri et al. 2015）：
  1. 按**框架区同源性**（FR1+FR2+FR3，CDR 除外）排序，要求 ≥ 60%
  2. 在 top 框架命中内，选 **CDR1+2 同源性最高**者 → 需要回复突变的位点最少、结构扰动最小
  3. J 基因按 FR4 同源性 + CDR3/J 交界长度匹配
- 报告输出前 5 名备选及各自 FR/CDR 同源性

### 阶段 3：CDR 移植（CDR Grafting）

- 逐位点组装：`FR1-FR3 ← 人源 germline V`；`CDR1-3 ← 供体`；`FR4 ← 人源 J`
- **四套 CDR 边界**同时报告（kabat / chothia / abm / imgt），变体梯度默认用 Kabat（`--scheme` 可换）
- **VHH 特例**：FR2 hallmark 位点（Kabat 37/44/45/47）与二硫键 Cys 永远保留供体，确保单域折叠与溶解性

### 阶段 4：回复突变分析与评分（核心）

见下节。

### 阶段 5：变体梯度设计

见「变体梯度」节。

### 阶段 6：结构模式（服务器，可选）

见「结构模式（服务器）」节。

---

## 回复突变评分体系

对每个框架位点（FR1–FR3）且供体 ≠ 人源 germline 的位置：

### 特征分类（来源：文献位置表 + 结构）

| 特征 | 定义 | 文献 |
|------|------|------|
| `vernier` | 支撑 CDR 构象的框架残基 | Foote & Winter, JMB 224:487 (1992) |
| `interface_core` | VH/VL 核心堆叠界面 | Chothia et al., JMB 186:651 (1985) |
| `interface_extended` | 扩展界面 | Vargas-Madrazo & Paz-Garcia, JMB 330:783 (2003) |
| `canonical` | CDR 规范结构支撑残基 | Chothia & Lesk, JMB 196:901 (1987)；Al-Lazikani et al., JMB 273:927 (1997) |
| `buried` | 埋藏于核心（relSASA < 20%，AF3） | — |
| `cdr_contact` | 与 CDR 重原子 < 4.5 Å（AF3） | — |
| `antigen_contact` | 与抗原重原子 < 4.5 Å（AF3 复合物） | — |
| `vhh_hallmark` | 37/44/45/47（VHH 专有，**禁止回复**） | Muyldermans, J Biotechnol 74:277 (2001) |
| `disulfide_cys` | 框架 Cys（CDR3 二硫键配对，**禁止回复**） | — |

### 评分公式

```
composite = 100 × (0.55 × structural + 0.30 × immunogenicity + 0.15 × chemical)

structural     = max(特征权重)      # 权重见 scripts/humanize/config.py
immunogenicity = 0.3 + 0.5 × exposure × (1 − conservation)
chemical       = 去除 NxS/T(+0.8) / NG,NS(+0.5) / DG(+0.4) / M,W(+0.3)
```

- `exposure`：AF3 判定的表面暴露度（无结构时默认 0.5）
- `conservation`：供体残基在 top-10 人源 germline 中的出现率（越稀有 → 回复收益越高）

### Tier 分级（最终建议）

| Tier | 规则 | 含义 |
|------|------|------|
| **KEEP_DONOR** | VHH hallmark 或二硫键 Cys | 禁止回复 |
| **T1** | interface_core（未明确暴露）；vernier ∩ buried；vernier ∩ canonical | 必须回复 |
| **T2** | vernier；canonical；cdr_contact；antigen_contact；buried；interface_extended | 强烈建议 |
| **T3** | 无结构角色（暴露或不明确） | 可选，仅免疫原性收益 |
| — | 供体 == 人源 | 非候选 |

完整权重与公式依据见 `docs/scoring.md`；位置表与文献解读见 `docs/theory.md`。

---

## 变体梯度

| 变体 | 内容 | 用途 |
|------|------|------|
| **V0** | 纯移植：全人源框架 + 供体 CDR | 活性探针（验证 CDR 本身是否决定结合） |
| **V1** | V0 + T1 回复突变 | 保守方案（框架同源性高时首选） |
| **V2** | V0 + T1 + T2 | **★ 推荐生产选择**（业界金标准案例的回复突变即属此类） |
| **V3** | V0 + T1 + T2 + 暴露 T3 | 最大化人源化（如抗药抗体问题），接受亲和力风险 |

---

## 结构模式（服务器）

### AlphaFold3

- 预测 donor Fv（可选加抗原复合物，`--antigen SEQ`）
- **回注评分**：buriedness、CDR 接触、抗原接触 → 精化阶段 4 的 structural score
- **变体验收**（`docs/validation.md` 第 4 节）：
  - CDR 环 CA-RMSD vs donor：H1/H2/L1/L2 < 1.0 Å，H3/L3 < 1.5 Å
  - VH/VL 界面埋藏面积偏差 < 15%
  - 抗原接触保留 ≥ 90% 供体接触，无新冲突

### ProteinMPNN

- 对移植变体（或 donor）骨架，**锁定 CDR + 界面 + vernier + VHH hallmark**，再设计其余框架
- 按设计框架的人源 germline 同源性过滤 + 可开发性风险筛查 → 输出 MPNN 设计变体（与回复突变路线并列的备选）
- 无 ProteinMPNN 时回退为 top-germline 共识框架（便携模式）

---

## 验证与金标准

`tests/test_pipeline.py`（60+ 断言，`python3 tests/test_pipeline.py` 或 pytest）：

- 编号引擎 vs 已知序列（IGHV3-23、IGKV1-39、mouse VH1、VHH、kappa）
- germline 选择、graft CDR 完整性、变体长度一致性
- **4D5 → 曲妥珠单抗基准**：历史上关键回复突变（VH 48/67/69/71/73；VL 4/66/87）必须落在 T1/T2
- **VHH 基准（1BZQ）**：hallmark 位点必须 KEEP_DONOR，且所有变体保留

完整验证方案（含 ANARCI 交叉验证与结构验收）见 `docs/validation.md`。

---

## 实验验证 SOP

`docs/experimental_SOP.md` 提供 6–10 周的生产路径：

1. **基因合成**（第 0–1 周）：密码子优化（Expi293）、Fab/VHH-Fc/IgG 格式、克隆面板 = V2 + V0 + V1 + 亲本
2. **表达纯化**（第 2–4 周）：Expi293F 小试→放大；Protein A / KappaSelect；SEC ≥ 95% 单体；LC-MS 质控
3. **功能验证**（第 4–6 周）：BLI/SPR 亲和力（**KD ≤ 2× 亲本**）、报告基因/竞争 ELISA 活性、表位分箱
4. **可开发性面板**（第 6–7 周）：DSF（Tm1）、40°C 加速稳定性、HIC、VHH 二硫键完整性
5. **免疫原性评估**（可选）：netMHCIIpan 预测、PBMC 体外 T 细胞实验
6. **Go/No-Go 判定表** + 亲和力丢失的排查回路（AF3 接触分析 → 定向补回突变 → 迭代）

---

## 服务器部署

详见 `docs/tool_setup.md`：

```bash
# 核心 + 精确编号（推荐每台服务器）
conda create -n humanize python=3.11 -y && conda activate humanize
conda install -c bioconda -c conda-forge hmmer==3.3.2 -y
pip install anarci abnumber biopython

# IgBLAST + NCBI germline（选择的数据源）
python3 scripts/humanize/cli.py setup-germline --dir data/germline

# AF3（GPU）与 ProteinMPNN 另行部署，通过 --af3-*/--mpnn-* 参数接入
```

---

## FAQ

**Q: 便携模式与服务器模式结果会不同吗？**
A: 阶段 1–5（编号/germline/移植/评分/变体）完全一致。服务器模式额外通过 AF3 提供 buriedness/CDR 接触等结构提示，精化 T1/T2 判定；无结构时按文献位置表分级（仍是业界基线）。

**Q: 为什么推荐 V2 而不是 V1/V3？**
A: 历史获批案例（曲妥珠单抗、帕利珠单抗等）的回复突变恰好对应 T1+T2 类；V1 可能活性保留差，V3 引入无谓的亲和力风险。

**Q: VHH 与 Fab 人源化有何不同？**
A: VHH 的 FR2 hallmark（37/44/45/47）必须保留驼源残基（无 VL 伙伴时的亲水补丁），人源框架选择通常用 IGHV3 家族；CDR3 常含二硫键 Cys，必须完整移植。

**Q: 如何换用 NCBI 官方 germline？**
A: `humanize setup-germline` 下载后，管线自动优先读取 `human_gl_*.fasta`；内置 JSON 为同源数据的 Kabat 坐标版本，互为后备。

**Q: 报告里的分数怎么看？**
A: composite ≥ 60 且 T1/T2 → 高置信回复；40–60 → 中等；< 40 通常为 T3。分数是排序辅助，不是硬过滤。

---

## 引用文献

- Foote J, Winter G. Antibody framework residues affecting the conformation of the hypervariable loops. *J Mol Biol* 224:487–499 (1992)
- Chothia C, et al. Domain association in immunoglobulin molecules. *J Mol Biol* 186:651–663 (1985)
- Chothia C, Lesk AM. Canonical structures for the hypervariable regions of immunoglobulins. *J Mol Biol* 196:901–917 (1987)
- Al-Lazikani B, et al. Standard conformations for the canonical structures of immunoglobulins. *J Mol Biol* 273:927–948 (1997)
- Muyldermans S. Single domain camel antibodies. *J Biotechnol* 74:277–302 (2001)
- Carter P, et al. Humanization of an anti-p185HER2 antibody for human cancer therapy. *PNAS* 89:4285–4289 (1992)
- Olimpieri PP, et al. Prediction of sites of humanization… *Bioinformatics* 31:434–435 (2015)
- Dunbar J, Deane CM. ANARCI: antigen receptor numbering and receptor classification. *Bioinformatics* 32:298–300 (2016)

---

## 许可证与致谢

- 代码：MIT License（见 LICENSE）
- 内置 germline 数据：源自 [abnumber](https://github.com/prihoda/AbNumber)（MIT）内嵌的人类 IMGT germline 集合，与 NCBI IgBLAST `human_gl_*.fasta` 同源；NCBI 下载版使用时请遵循 NCBI/IMGT 引用要求
- 结构/设计工具：AlphaFold3（Google DeepMind）、ProteinMPNN（Baker Lab）均为各自许可证
