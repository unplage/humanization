# User Guide 使用指导（全面版）

本文是抗体人源化流程的完整使用手册：安装、输入规范、全部命令、端到端示例、
报告解读、进阶模式（AF3/ProteinMPNN/实验闭环/BioPhi）、故障排查与生产工作流。

---

## 1. 项目概览

```
输入（非人源抗体序列）
      │
      ▼
┌─────────────────────────────────────────────────────────┐
│ 1. 编号       anchor-based Kabat 编号（零依赖）           │
│ 2. Germline   NCBI IgBLAST 数据 / 内置 373V+28J，hybrid 选架 │
│ 3. CDR 移植   kabat/chothia/abm/imgt 四套边界             │
│ 4. 回复突变   vernier/界面/canonical/接触 → T1/T2/T3 评分  │
│ 5. 变体梯度   V0/V1/V2/V3 + Vmin + V_SDR + 框架矩阵       │
│ 6. 验证       AF3 结构提示、实验数据闭环校准、回测基准      │
└─────────────────────────────────────────────────────────┘
      │
      ▼
输出：humanization_report.md / .docx / .json / *.csv / variants.fasta
```

支持格式：**Fab**（VH+VL）、**VHH**（单链驼源纳米抗体，自动检测 hallmark 并保护）。

---

## 2. 安装

### 2.1 便携模式（本机，零依赖）

```bash
# 直接运行即可（Python 3.8+，仅标准库）
python3 scripts/humanize/cli.py run --input 你的序列.fasta --outdir outputs
```

内置 germline 数据（373 个人源 V 基因 + 28 个 J 基因）随仓库分发，
无需联网。可选安装 `pip install python-docx` 以生成 Word 报告。

### 2.2 服务器完整模式

```bash
bash scripts/install.sh --full          # ANARCI + IgBLAST + NCBI germline
conda activate humanize
python3 tests/test_pipeline.py          # 自检
```

按需追加：
- **AlphaFold3**（结构评分）：`--af3-mode local --af3-binary ...`
- **ProteinMPNN**（框架再设计）：`--mpnn-mode local --mpnn-script ...`
- **BioPhi**（人源度交叉验证）：
  ```bash
  conda create -n biophi python=3.9 && conda install -n biophi biophi -c bioconda -c conda-forge
  ```
- **OASis 数据库**（可选，22GB）：见 docs/tool_setup.md

---

## 3. 输入规范

### 3.1 FASTA 约定

```fasta
>4D5_VH
EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTKYNQKFKGKATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS
>4D5_VL
DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK
```

- Fab：VH 与 VL 各一条记录（顺序任意）；VHH：单条记录
- 序列应为**完整的 V 结构域**（约 110-130 aa）；含 CH1/CL 尾巴也可
  （引擎自动忽略 FR4 之后的序列并告警）
- 名称建议以 `_VH`/`_VL` 结尾（便于识别，非强制）

### 3.2 命令行完整参考

```bash
# 主命令
python3 scripts/humanize/cli.py run \
  --input <fasta> \                    # 必需：输入序列
  --outdir outputs \                   # 输出目录（默认 outputs）
  --format auto|fab|vhh \              # 格式强制（默认自动检测）
  --scheme kabat|chothia|abm|imgt \    # 变体梯度的 CDR 定义（默认 kabat）
  --germline-dir <dir> \               # NCBI germline FASTA 目录（可选）
  --germline-strategy <strategy> \     # germline 选择策略（见下表）
  --antigen <seq> \                    # 抗原序列（AF3 复合物 + V_SDR）
  --calibration calibration.json \     # 实验数据校准文件
  --biophi-env biophi \                # BioPhi conda 环境（人源度交叉验证）
  --oasis-db <path> \                  # OASis 数据库（可选）
  --af3-mode off|local|api \           # AlphaFold3 模式
  --af3-binary <path> \                # run_alphafold.py 路径
  --mpnn-mode off|local \              # ProteinMPNN 模式
  --mpnn-script <path> \               # protein_mpnn.py 路径

# 工具与数据
python3 scripts/humanize/cli.py setup-check          # 环境自检
python3 scripts/humanize/cli.py setup-germline       # 下载 NCBI germline

# 实验闭环
python3 scripts/humanize/cli.py learn \
  --experiments experiments.json --out calibration.json
```

### Germline 选择策略 (--germline-strategy)

| 策略 | 说明 | 推荐场景 |
|------|------|----------|
| `fr_best` | FR 同源性最高 | 框架稳定性优先 |
| `cdr_best` | CDR 同源性最高 | CDR 保守性优先 |
| `composite` | 0.7*FR + 0.3*CDR | 平衡考虑 |
| `cvi_best` | CVI 同源性最高 | **推荐：BI 2024 策略** |
| `min_backmutations` | 估计回复突变最少 | 最小改动 |
| `adimab_frequency` | Adimab 推荐 germline + 频率 | 治疗性抗体优化 |
| `pioneer_frequency` | Pioneer 库 germline + 频率 | 600+ 临床阶段抗体 |
| `composite_3axis` | 0.5*CVI + 0.3*频率 + 0.2*FR | **综合最优** |
| `auto` (默认) | VH=cvi_best, VL=cdr_best | **默认推荐** |

**使用频率数据**：基于已获批和临床阶段的治疗性抗体分析，IGHV3-23 (~18%)、IGHV1-69 (~12%)、IGKV1-39 (~15%) 是最常用的 germline。

---

## 4. 端到端示例（鼠源 4D5 → 人源化）

```bash
# 1. 运行
python3 scripts/humanize/cli.py run \
  --input data/examples/mouse_4d5_fab.fasta --outdir outputs

# 2. 查看控制台摘要
#   [humanize] format: FAB
#   4D5_VH (H)
#     germline: IGHV1-3*01 + IGHJ1*01 (FR id 0.73, CDR id 0.55)
#     back-mutations: 20 (T1:2, T2:3, T3:15)
#     H_V0: 0 / H_V1: 2 / H_V2: 5 / H_V3: 13
#   ...

# 3. 查看报告（三选一）
open outputs/humanization_report.md
open outputs/humanization_report.docx     # Word 版（需 python-docx）
# outputs/humanization_result.json        # 机器可读
```

### 4.1 输出文件清单

| 文件 | 内容 |
|------|------|
| `humanization_report.md` | 完整报告（文本） |
| `humanization_report.docx` | 完整报告（Word，专业排版：封面/摘要/表格/附录） |
| `enhanced_report.md` | 增强报告（借鉴 WeMol 格式：Template Score 表、Mutation Score、Back-Mutation 摘要、Hotspot 摘要、人源化序列） |
| `humanization_result.json` | 结构化数据（脚本/下游处理用） |
| `backmutations_<链>.csv` | 逐位点回复突变明细 |
| `variants.fasta` | 全部变体序列（基因合成直接可用） |

### 4.2 报告解读要点

1. **执行摘要**：germline 选择、候选数（按 tier）、推荐变体（V2）
2. **编号表**：每链 FR1-CDR1-...-FR4 的位置范围与序列（Kabat 编号）
3. **Germline 表**：主选 + 前 5 备选的 FR/CDR 同源性（**CDR 同源性高 = 回复突变少**）
4. **回复突变表**：每行的含义——
   - `tier`：T1 必须回复 / T2 建议 / T3 可选 / KEEP_DONOR 禁止回复
   - `score`：0-100 复合分（结构 0.55 + 免疫原性 0.30 + 化学 0.15）
   - `features`：vernier / interface_core / canonical / buried / cdr_contact 等
   - `ddG(emp)`：实验数据校准后的实测效应（kcal/mol，正值 = 保留 donor 有益）
5. **变体表**：V0（纯移植）→ V2（推荐）→ V3（激进），序列完整列出
6. **增强报告**（`enhanced_report.md`）：WeMol 风格——
   - Template Score：候选 germline 按 FR%/FR1%/FR2%/FR3% 排名
   - Mutation Score：回复突变评分（可与 WeMol 平台结果对比核对）
   - Back-Mutation Summary + KABAT 编号突变摘要
   - Hotspot Summary：CDR 与全长的人源度、风险位点
   - Humanized Sequences：人源化序列总览

---

## 5. 进阶模式

### 5.1 结构模式（AF3，服务器）

```bash
python3 scripts/humanize/cli.py run --input seq.fasta --outdir outputs \
  --af3-mode local --af3-binary /path/run_alphafold.py \
  --antigen <抗原序列>          # 需要 V_SDR 变体时必填
```

启用后：
- 每个框架位点获得 buriedness / CDR 接触 / 抗原接触提示 → 精化 tier 判定
- 产出 **Vmin**（最小回复集，set-cover 保证接触恢复）
- 提供抗原时产出 **V_SDR**（仅移植抗原接触的 CDR 残基 + 结构支柱）

### 5.2 ProteinMPNN 框架再设计

```bash
python3 scripts/humanize/cli.py run --input seq.fasta \
  --mpnn-mode local --mpnn-script /path/protein_mpnn.py
```

锁定 CDR+界面+vernier+hallmark 后重新设计框架；按人源同源性过滤输出。
无 MPNN 时自动回退为 top-germline 共识框架。

### 5.3 实验数据闭环（推荐生产使用）

```json
// experiments.json
[{"name": "Ab1",
  "parent_vh": "...", "parent_vl": "...", "parent_kd": 0.15,
  "variants": [
    {"name": "Ab1_hV0", "vh": "...", "vl": "...", "kd": 0.9},
    {"name": "Ab1_hV2", "vh": "...", "vl": "...", "kd": 0.17}]}]
```

```bash
python3 scripts/humanize/cli.py learn --experiments experiments.json --out calibration.json
python3 scripts/humanize/cli.py run --input new_ab.fasta --calibration calibration.json
```

效果：实测 ΔΔG 自动修正分级（如 L87 类位点实测无效应 → 自动降级 T3）。

### 5.4 BioPhi/Sapiens 人源度交叉验证

```bash
python3 scripts/humanize/cli.py run --input seq.fasta --biophi-env biophi
```

报告新增逐变体 Sapiens 均分与 OASis 身份，与 germline 同源性指标互相印证。

---

## 6. 故障排查

| 现象 | 原因 | 处理 |
|------|------|------|
| `could not number as VH or VL` | 序列残缺/非 V 区/特殊环长 | 检查序列完整性；用 ANARCI 模式交叉验证 |
| 警告 `自动降级至 FR >= xx% 继续选择` | 亲本与人源框架同源性过低（<60%） | 属正常警示：人源化难度高、回复突变会很多；考虑 SDR/paratope 路线 |
| `FR identity 0.5x`（VHH 常见） | VHH 与人类框架天然差距大 | 正常；hallmark 位点已自动保护 |
| 报告出现 `unusual CDR2/FR3` 警告 | 长环/插入 | 服务器 ANARCI 精确编号复核 |
| 报告出现 `Cys22 anchor shifted` | N 端多残基使 Cys 不在 22 位 | 引擎已自动 cap FR1=30 并吸收多余残基；建议 ANARCI 复核 |
| 变体数 < 5 | Vmin 与 V2 集相同时不重复输出 | 正常 |
| `trailing residue(s) ignored` | 序列含 CH1/CL 或标签 | 正常（报告已忽略） |
| docx 未生成 | 未安装 python-docx | `pip install python-docx` |

---

## 7. 生产工作流（推荐）

```
第 0 轮   humanize run → V0/V1/V2 面板 → 合成+表达+BLI
第 1 轮   humanize learn（V0/V2 KD 数据）→ calibration.json
         humanize run --calibration → 修订版 V2（实测感知）
第 2 轮   单点变体面板（修订 V2 的边界位点）→ 精细校准
第 3 轮   组合面板（T3 子集）→ 组织级 calibration.json（跨项目复用）
```

每轮决策依据 `docs/experimental_SOP.md` 的 Go/No-Go 表
（KD ≤ 2× 亲本、SEC > 95%、Tm1 不劣化、表位不变）。

---

## 8. 验证与基准

```bash
python3 tests/test_pipeline.py       # 单元+集成（含 4D5→曲妥珠基准）
python3 tests/backtest.py            # 深度回测（曲妥珠/贝伐珠）
python3 tests/backtest_scale.py      # HumAb25：25 个真实药物亲本规模化验证
```

金标准结论（docs/backtest_report.md）：结构关键位点召回率 1.0，无关键漏报。

---

## 9. 文档索引

| 文档 | 内容 |
|------|------|
| `docs/theory.md` | 回复突变理论（位置表+文献） |
| `docs/scoring.md` | 评分公式与权重 |
| `docs/improvement_roadmap.md` | 差距分析 + 行业对比 + 路线图 |
| `docs/backtest_report.md` | 回测方法与结论 |
| `docs/learning_loop.md` | 实验数据闭环 |
| `docs/ecosystem_analysis.md` | GitHub 同类项目全景 |
| `docs/tool_setup.md` | 服务器部署 |
| `docs/validation.md` | 验证方案 |
| `docs/experimental_SOP.md` | 实验 SOP（6-10 周） |
