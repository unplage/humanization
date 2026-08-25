# Antibody Humanization Pipeline · 抗体人源化流程

**生产级 CDR grafting + 回复突变（back-mutation）设计系统，支持 Fab 与 VHH。**

从非人源抗体序列出发，完成：人源 germline 匹配 → CDR 移植（4 套定义）→
**回复突变量化评分与分级（T1/T2/T3/KEEP_DONOR）** → 变体梯度设计 →
（服务器端）AlphaFold3 结构验证 + ProteinMPNN 框架再设计 → 实验数据闭环校准
→ 专业 Word/Markdown 报告 → 实验 SOP（Go/No-Go 判定）。

本机便携模式**零安装、零系统依赖**（纯 Python 标准库端到端运行）；
已验证于真实获批案例（曲妥珠单抗、贝伐珠单抗、25 例 HumAb25 规模化基准）。

---

## 快速开始

```bash
# 1. 人源化（Fab：VH+VL 同一 FASTA；VHH：单链）
python3 scripts/humanize/cli.py run \
  --input data/examples/mouse_4d5_fab.fasta --outdir outputs

# 2. 使用多策略 Germline 选择（推荐）
python3 scripts/humanize/cli.py run \
  --input data/examples/mouse_4d5_fab.fasta --outdir outputs \
  --germline-strategy auto  # 默认：VH=cvi_best, VL=cdr_best

# 3. 查看所有策略的选择结果
python3 scripts/custom_humanize.py

# 4. 查看报告（Markdown / Word / JSON）
open outputs/humanization_report.md
open outputs/humanization_report.docx        # 需 pip install python-docx

# 5. 自检与测试
python3 scripts/humanize/cli.py setup-check
python3 tests/test_pipeline.py
```

详细使用手册见 **[docs/USER_GUIDE.md](docs/USER_GUIDE.md)**。

---

## Germline 选择策略

支持 9 种选择策略，基于 BI 2024 矩阵移植研究和治疗性抗体数据优化：

| 策略 | 说明 | 适用场景 |
|------|------|----------|
| `fr_best` | FR 同源性最高 | 框架稳定性优先 |
| `cdr_best` | CDR 同源性最高 | CDR 保守性优先 |
| `composite` | 综合评分 (0.7*FR + 0.3*CDR) | 平衡考虑 |
| `cvi_best` | CVI 同源性最高 | **推荐：BI 2024 策略** |
| `min_backmutations` | 估计回复突变最少 | 最小改动 |
| `adimab_frequency` | Adimab 推荐 germline + 频率 | 治疗性抗体优化 |
| `pioneer_frequency` | Pioneer 库 germline + 频率 | 600+ 临床阶段抗体 |
| `composite_3axis` | 0.5*CVI + 0.3*频率 + 0.2*FR | **综合最优** |
| `auto` | 自动：VH=cvi_best, VL=cdr_best | **默认推荐** |

**CVI 同源性** = Canonical + Vernier + Interface 同源性，基于 BI 2024 研究，
与表达量和亲和力保留显著相关。

**使用频率数据** 基于已获批和临床阶段的治疗性抗体分析：
- IGHV3-23 是最常用的 VH germline (~18%)
- IGHV1-69 是第二常用的 VH germline (~12%)
- IGKV1-39 是最常用的 VK germline (~15%)

---

## 核心能力

| 能力 | 说明 |
|------|------|
| **Fab / VHH 双格式** | 自动链型检测；VHH hallmark（Kabat 37/44/45/47）识别并全程保护 |
| **零依赖便携核心** | anchor-based Kabat 编号引擎（纯标准库）；内置 373 V + 28 J 人源 germline |
| **多策略 Germline 选择** | 9 种策略（FR/CDR/综合/CVI/最小突变/Adimab频率/Pioneer频率/3轴/自动），基于 BI 2024 和治疗性抗体数据优化 |
| **多套 CDR 定义** | Kabat / Chothia / AbM / IMGT 四套边界同时报告 |
| **回复突变量化评分** | 结构（vernier/界面/canonical/接触/埋藏）+ 免疫原性 + 可开发性三维评分 → T1/T2/T3/KEEP_DONOR 分级与逐条依据 |
| **变体梯度** | V0 纯移植 → V1(T1) → **V2(T1+T2，推荐)** → V3(+暴露T3) → Vmin(最小回复集) → V_SDR(paratope 精确移植) |
| **结构模式（服务器）** | AF3 预测 → buried/接触/抗原接触回注评分；CDR 环 RMSD、界面与接触保留率验收 |
| **框架矩阵** | 备选 germline 面板 + CVI 同源性指标（Boehringer Ingelheim JBC 2024 实证） |
| **实验数据闭环** | `humanize learn`：KD 数据 → 逐位点 ΔΔG → 自动校准分级与评分 |
| **可开发性检查** | N-糖基化 / 脱酰胺 / 异构化 / 氧化风险自动扫描 |
| **专业报告** | Markdown + **Word（.docx）** + JSON + 逐位点 CSV + variants FASTA |
| **验证基准** | 4D5→曲妥珠单抗、A4.6.1→贝伐珠单抗回测 + **25 例 HumAb25 规模化验证** |

---

## 典型输出（一次 run 的产物）

| 文件 | 内容 |
|------|------|
| `humanization_report.md` | 完整报告：执行摘要、编号、germline 选择、回复突变表、变体序列 |
| `humanization_report.docx` | **Word 专业版**：封面 + 12 章 + 2 附录（含全部完整序列与逐位点对照表） |
| `enhanced_report.md` | **增强报告**（借鉴 WeMol 格式）：Template Score、Mutation Score、回复突变摘要、Hotspot 摘要、人源化序列 |
| `humanization_result.json` | 机器可读完整结果 |
| `backmutations_<链>.csv` | 逐位点评分明细（含实验 ddG 列） |
| `variants.fasta` | 全部变体序列（直接用于基因合成） |

**生产推荐：V2（T1+T2 回复突变）**——金标准案例（曲妥珠单抗等）的
历史回复突变恰好属于该分级；Vmin 提供"最少回复突变"的安全下限。

---

## 服务器完整模式（一键安装）

```bash
bash scripts/install.sh --full     # ANARCI + IgBLAST + NCBI germline
conda activate humanize

# 结构模式（AlphaFold3 + ProteinMPNN + 抗原复合物）
python3 scripts/humanize/cli.py run --input seq.fasta --outdir outputs \
  --af3-mode local --af3-binary /path/run_alphafold.py \
  --mpnn-mode local --mpnn-script /path/protein_mpnn.py \
  --antigen <抗原序列>            # 启用 V_SDR paratope 移植

# 人源度交叉验证（BioPhi/Sapiens，Merck 开源）
python3 scripts/humanize/cli.py run --input seq.fasta --biophi-env biophi

# 实验数据闭环
python3 scripts/humanize/cli.py learn --experiments experiments.json --out calibration.json
python3 scripts/humanize/cli.py run --input seq.fasta --calibration calibration.json
```

---

## 回复突变设计逻辑（核心）

对每个框架位点（供体 ≠ 人源 germline）：

1. **特征分类**：vernier 区（Foote & Winter 1992）、VH/VL 界面（Chothia 1985）、
   canonical 残基（Chothia & Lesk 1987）、AF3 接触/埋藏、VHH hallmark、二硫键 Cys
2. **三维评分**：`composite = 100 × (0.55×结构 + 0.30×免疫原性 + 0.15×化学)`
3. **分级**：T1 必须回复 / T2 建议 / T3 可选 / KEEP_DONOR（hallmark、Cys）禁止回复
4. **实验校准**：KD 数据 → 逐位点 ΔΔG → 自动提升/降级

```bash
# 生成一个实际的回复突变表（示例）
python3 scripts/humanize/cli.py run --input data/examples/mouse_4d5_fab.fasta --outdir demo
grep -A25 "Back-mutation candidates" demo/humanization_report.md
```

评分公式与权重见 `docs/scoring.md`；理论依据见 `docs/theory.md`。

---

## 文档索引

| 文档 | 用途 |
|------|------|
| **[docs/USER_GUIDE.md](docs/USER_GUIDE.md)** | ★ 完整使用指导（安装/输入/命令/示例/故障排查/生产工作流） |
| [docs/theory.md](docs/theory.md) | 回复突变理论：位置表 + 文献依据 |
| [docs/scoring.md](docs/scoring.md) | 评分公式、权重与解释规则 |
| [docs/backtest_report.md](docs/backtest_report.md) | 回测方法与结论（曲妥珠/贝伐珠金标准） |
| [docs/improvement_roadmap.md](docs/improvement_roadmap.md) | 差距分析 + 行业方案对比 + P0-P2 路线图 |
| [docs/learning_loop.md](docs/learning_loop.md) | 实验数据闭环校准工作流 |
| [docs/ecosystem_analysis.md](docs/ecosystem_analysis.md) | GitHub 同类项目全景与采用分析 |
| [docs/tool_setup.md](docs/tool_setup.md) | 服务器部署（AF3/MPNN/BioPhi/IgBLAST） |
| [docs/validation.md](docs/validation.md) | 验证方案（ANARCI 交叉验证/结构验收） |
| [docs/experimental_SOP.md](docs/experimental_SOP.md) | 实验 SOP：合成→表达→亲和力→可开发性→Go/No-Go |

---

## 架构

```
scripts/humanize/
├── numbering.py       # anchor-based Kabat 编号引擎（零依赖，核心）
├── germline.py        # germline 加载/选择（NCBI FASTA 优先，内置后备）
├── graft.py           # CDR 移植（4 方案）+ 变体组装
├── backmut.py         # 回复突变候选 + 三维评分 + 分级（核心）
├── minimal.py         # Vmin 最小回复集 / CVI / 框架矩阵 / V_SDR
├── variants.py        # V0-V3 变体梯度
├── learning.py        # 实验数据闭环（ΔΔG 校准）
├── humanness.py       # BioPhi/Sapiens 人源度适配器
├── structure.py       # AlphaFold3 适配器 + 无依赖 PDB/CIF 解析
├── mpnn.py            # ProteinMPNN 适配器
├── developability.py  # 可开发性风险扫描
├── report.py          # Markdown/JSON/CSV/FASTA 报告
├── report_docx.py     # ★ Word (.docx) 专业报告
├── pipeline.py        # 端到端编排
└── cli.py             # 命令行入口
```

---

## 验证与基准

```bash
python3 tests/test_pipeline.py      # 75+ 项（编号/graft/分级/VHH/闭环/适配器）
python3 tests/backtest.py           # 深度回测（4D5→曲妥珠、A4.6.1→贝伐珠）
python3 tests/backtest_scale.py     # HumAb25：25 个真实药物亲本全过
```

- 结构关键位点召回率 **1.0**（无关键漏报），过回复均为保守方向的低风险差异
- 修正了流传的鼠源 4D5 CDR 序列错误（真实 CDR = DTYIH / SRWGGDGFYAMDY，
  经 4 个 PDB 结构 + HumAb25 双重印证）

---

## 许可证与致谢

- 代码：MIT License
- 内置 germline 数据：源自 abnumber（MIT）内嵌的人类 IMGT germline 集合
- 基准数据：HumAb25（TencentAI4S/HuDiff，Zenodo DOI 10.5281/zenodo.16974296）
- 结构/设计工具：AlphaFold3（Google DeepMind）、ProteinMPNN（Baker Lab）、
  BioPhi/Sapiens（Merck）、ANARCI（OPIG）
