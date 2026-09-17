# 抗体人源化流程指导文件

## 概述

本文件总结抗体人源化（CDR 移植 + 回突变设计）的完整流程，包括四步走标准流程、独立分析工具、结构验证、免疫原性评估等环节，以及各步骤的关键注意点和已知陷阱。

---

## 一、标准四步走流程

### Step 1：Germline 策略对比

**目标**：评估 9 种策略，选择最佳 germline 组合

**命令**：
```bash
python3 scripts/humanize/cli.py compare --input <seq.fasta>
```

**9 种策略**：
| 策略 | 说明 |
|------|------|
| `fr_best` | FR 同源性最高 |
| `cdr_best` | CDR 同源性最高 |
| `composite` | 0.7×FR + 0.3×CDR |
| `cvi_best` | canonical+vernier+interface 综合最优 |
| `min_backmutations` | 回突变最少 |
| `current` | top-30% FR with max CDR |
| `adimab_frequency` | Adimab 使用频率 |
| `pioneer_frequency` | Pioneer 使用频率 |
| `composite_3axis` | 0.5×CVI + 0.3×freq + 0.2×FR |

**注意点**：
- 每个策略必须返回 5 个候选
- `auto` 是默认路由策略（VH=`adimab_frequency`，VL=`current`），不需要单独跑
- 分数相同时按使用频率排序

---

### Step 2：按选定 germline 正式运行

**目标**：生成完整报告集（markdown/docx/json/csv/fasta）

**命令**：
```bash
python3 scripts/humanize/cli.py run --input <seq.fasta> \
  --germline-strategy <chosen> \
  --outdir outputs
```

**注意点**：
- 可用 `--force-germline VH=IGHV1-3*01 VL=IGKV1-39*01` 锁定具体基因
- 输出包括：变体序列、回突变评分、CDR 移植结果等

---

### Step 3：结构优化（可选）

**目标**：用实际结构增强回突变分级（接触保持/提升；暴露非接触降级为 T3；FR4 结构性回复突变）

**命令**：
```bash
python3 scripts/humanize/cli.py run --input <seq.fasta> \
  --donor-structure <pdb> \
  --outdir outputs
```

**进阶选项**：
```bash
# AF3 本地模式
--af3-mode local --af3-binary <run_alphafold.py> --antigen <seq>

# 结构验证
--af3-rmsd
```

**关键注意点**：
1. **必须使用 variant 特异性结构**：Donor 构象不充分，backmutation 改变局部 packing
2. **CDR-RMSD 评估**：使用 Kabat 编号范围划分 FR1/FR2/FR3/FR4
3. **FR4 为 J 区域**：不纳入框架叠合（fit），但叠合后单独报告 FR4 偏差（post-fit），量化 J 区替换对 CDR3 邻近骨架的影响
4. **pLDDT 过滤**：只使用高置信度区域进行叠合

---

### Step 4：可开发性优化（可选）

**目标**：检测高风险开发性基序（DD/NG/NS/DG/MW），运行 MPNN 优化

**命令**：
```bash
python3 scripts/humanize/cli.py run --input <seq.fasta> \
  --mpnn-mode local --mpnn-script <path> \
  --design-panel \
  --outdir outputs
```

**注意点**：
- `--design-panel` 生成 V_opt panel（分级回突变子集）
- 用于实验筛选

---

## 二、独立分析工具

### 1. 人源度评估 (humanness_score)

**功能**：Kabat + IMGT 双编号评估 FR/CDR 人源性

**用法**：
```bash
python3 tools/humanness_score/evaluate_humanness.py \
  --input <fasta> --donor <donor_fasta> --output <outdir>
```

**注意点**：
- 输出 14 个 Word 报告 + 合并报告 + Markdown 摘要
- 人源度分档：Human-like (>90%) / Humanized-like (80-90%) / Chimeric-like (70-80%) / Murine-like (<70%)

---

### 2. CDR 同源性分析 (cdr_homology)

**功能**：对比 TheraSAbDab 临床抗体数据库

**用法**：
```bash
python3 tools/cdr_homology/query.py \
  --input <fasta> --top 5 --output <outdir>
```

**⚠️ 关键陷阱**：
1. **FASTA 命名必须匹配配对逻辑**：
   - 正确：`>HV0_LV0_Heavy` + `>HV0_LV0_Light`
   - 错误：`>H_V0_L_V0_VH` + `>H_V0_L_V0_VL`
2. **配对逻辑**：移除 `_Heavy/_Light/_VH/_VL` 后缀后匹配 base name
3. **依赖 ANARCI**：需要安装 ANARCI 进行 CDR 提取
4. **依赖 BioPython**：需要安装 BioPython 进行序列比对

---

### 3. PTM 暴露度分析 (ptm_exposure)

**功能**：预测翻译后修饰位点（N-糖基化、脱酰胺、氧化等）

**用法**：
```bash
python3 tools/ptm_exposure/analyze_ptm.py \
  --input <fasta> --output <outdir>
```

**⚠️ 已知 Bug 及修复**：
1. **Cys 双重计数**：同一位置的 Cys 相关 PTM 会重复计数
   - 修复：添加 `seen_positions` 字典去重
2. **MAX_SASA 值不一致**：不同模块使用不同的参考值
   - 修复：统一使用 Tien et al. 2013 标准值

---

### 4. 免疫原性分析 (immunogenicity)

**功能**：MHC-II T 细胞表位预测 + 结构风险调整

**两步流程**：
1. **序列分析**：HLAIIPred 预测 MHC-II 结合
2. **结构确认**：FreeSASA 计算 relSASA，调整风险

**用法**：
```bash
# Step 1: 序列免疫原性
python3 tools/immunogenicity/immunogenicity_analyzer.py \
  --input <fasta> --output <outdir> --all-formats

# Step 2: 结构风险调整
python3 tools/immunogenicity/structural_risk.py \
  --immunogenicity-json <json> \
  --structure-dir <pdb_dir> \
  --backmutation-dir <csv_dir> \
  --output-dir <outdir> --method relSASA
# 注：工具会自动读取 <csv_dir>/variants_numbering.json 以精确对齐 Kabat 位置
# （含 FR4 与插入位点）；也可显式传 --numbering-json <path>
```

**⚠️ 关键注意点**：
1. **HLAIIPred 安装**：
   ```bash
   conda create -n hlapred python=3.11
   conda activate hlapred
   pip install torch --index-url https://download.pytorch.org/whl/cpu
   git clone https://github.com/zaferozk/HLAIIPred.git
   cd HLAIIPred && pip install -r requirements.txt
   ```
2. **Allele Padding**：runner 用 `0` 填充到 14 个 allele
3. **必须结合结构验证**：仅序列分析会高估免疫原性
4. **回突变 CSV 必须提供**：否则无法正确评估回突变位点埋藏/暴露状态

**⚠️ 已知 Bug 及修复**：
1. **per_residue dict key 类型不匹配**：JSON 字符串 key vs int lookup
   - 修复：key 归一化为 int
2. **无 donor 参数时 V0 未自动用作参考**
   - 修复：添加自动 fallback
3. **structural_risk.py 只考虑回突变位置**
   - 修复：无回突变数据时 fallback 到所有位置
4. **structural_risk.py 位置映射错误（Kabat vs 序列序号）**
   - 问题：回突变 CSV 是 Kabat 编号，工具却按序列序号匹配；插入位点（FR3 82A-C、CDR 27A-/100A-K）与 FR4 会被错配，且 `int("82A")` 直接崩溃
   - 修复：pipeline 输出 `variants_numbering.json`（每个变体的逐残基 Kabat 标签），工具自动发现并精确匹配；解析兼容插入字母
5. **FR4 是否被各工具覆盖**
   - `ptm_exposure`、`immunogenicity_analyzer` 序列扫描：**整链扫描，包含 FR4**
   - `structural_risk`：**包含 FR4**（T_FR4 权重 0.9），需 `variants_numbering.json` 才能精确对齐
   - `igfold`：只做结构预测（模型含 FR4），不做区域分析

---

### 5. TAP 可开发性分析 (tap_profiler)

**功能**：预测 TAP 转运效率

**用法**：
```bash
python3 tools/tap_profiler/tap_analyzer.py \
  --input <fasta> --output <outdir>
```

---

### 6. 稳定性预测 (stability_predictor)

**功能**：预测突变对稳定性的影响

**用法**：
```bash
python3 tools/stability_predictor/stability_analyzer.py \
  --input <fasta> --output <outdir>
```

---

## 三、排列组合结构验证

### IgFold 结构预测

**功能**：预测 49 个 H×L 组合的 3D 结构

**用法**：
```bash
python3 scripts/humanize/combinatorial_igfold.py \
  --input <fasta> --output <pdb_dir>
```

**关键优势**：
- 快速（~6秒/结构）
- 本地可复现
- 支持批量处理

### CDR-RMSD 评估

**功能**：评估变体结构与 donor 的 CDR 偏差

**用法**：
```bash
python3 scripts/humanize/cli.py rmsd \
  --input <fasta> --donor <donor.pdb> \
  --variants <variant1.pdb> <variant2.pdb> ... \
  --chain H
```

**注意点**：
- 使用 Kabat 编号划分 FR1/FR2/FR3/FR4
- FR4（J 区域）不纳入框架叠合，但会以 post-fit 方式单独报告 RMSD
- 输出 per-FR RMSD 值

---

## 四、关键经验教训

### 1. 结构验证是必须的

**问题**：仅序列分析会高估免疫原性

**解决方案**：
- 使用 IgFold/AF3 预测变体结构
- 计算 relSASA 评估表位暴露程度
- 埋藏表位（relSASA < 0.20）免疫原性风险降低

---

### 2. FASTA 命名规范

**问题**：工具间 FASTA 命名不一致导致配对失败

**解决方案**：
- 统一命名格式：`>HV0_LV0_Heavy` + `>HV0_LV0_Light`
- 避免使用 `>H_V0_L_V0_VH` 等复杂格式

---

### 3. 回突变数据必须提供

**问题**：无回突变 CSV 时，结构风险评估无法正确判断位点重要性

**解决方案**：
- Step 3 输出必须包含 `backmutations_*.csv`
- 结构风险评估必须使用 `--backmutation-dir` 参数

---

### 4. 依赖管理

**核心依赖**：
- ANARCI（CDR 编号）
- BioPython（序列比对）
- FreeSASA（SASA 计算）
- HLAIIPred（MHC-II 预测）
- python-docx（Word 报告）

**安装检查**：
```bash
python3 scripts/humanize/cli.py setup-check
```

---

### 5. 报告格式统一

**问题**：不同工具输出格式不一致

**解决方案**：
- 统一使用 Markdown + Word + JSON 三格式
- 回突变表述统一为 `位置编号 + 捐体氨基酸 > 人源氨基酸`（如 H5 L>V）
- 避免使用 IgK53I 等混淆写法

---

## 五、已知限制

1. **PTM 模块**：仅预测常见 PTM（N-糖基化、脱酰胺、氧化），不覆盖所有修饰
2. **免疫原性模块**：仅使用 9 个常用 DRB1 allele，不覆盖所有 HLA 类型
3. **结构预测**：IgFold 精度低于 AF3，但速度快 100 倍
4. **稳定性预测**：基于经验打分，不替代实验验证

---

## 六、推荐流程

```
1. Step 1: compare → 选择 germline
2. Step 2: run → 生成变体
3. Step 3: run --donor-structure → 结构增强
4. 独立工具: humanness + cdr_homology + ptm + immunogenicity
5. IgFold: 预测 49 个组合结构
6. RMSD: 评估结构偏差
7. 免疫原性 + 结构风险调整
8. 综合评估 → 选择候选
```

---

## 七、输出文件组织

```
outputs/
├── step1/
│   └── compare_results.md
├── step2/
│   ├── variants.fasta
│   └── backmutations_*.csv
├── step3/
│   ├── variants.fasta
│   ├── backmutations_*.csv
│   └── HUMANNESS_REPORT.md
├── cdr_homology/
│   ├── cdr_homology_report.md
│   └── cdr_homology_report.docx
├── immunogenicity/
│   ├── *_immunogenicity.json
│   ├── structural_risk_assessment.md
│   └── immunogenicity_with_structural_risk.json
└── combinatorial/
    ├── *.pdb
    ├── RMSD_REPORT.md
    └── COMBINED_PTM_RMSD_REPORT.md
```

---

*最后更新：2026-09-17*
*基于 AMG110 人源化项目经验总结*
