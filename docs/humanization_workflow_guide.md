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

**功能**：MHC-II T 细胞表位预测 + 结构风险调整（两步合一）

**核心指标（新定义）**：
- **Strong / Weak**：15-mer 滑窗计数，≥ strong_thr / [weak_thr, strong_thr)
- **Density**：(strong + weak) / 总窗口数
- **Score**：`100 * variant_mean / reference_mean`，**无 100 封顶**，>100 表示更高风险
- **参考序列**：提供 `--donor` 时永远用 donor 序列单独打分；否则用精确 `V0 / H_V0 / L_V0 / combo_V0`；可用 `--reference-variant` 显式指定
- **结构校正**：`--structural-risk` 开启后，**所有变体和参考均用各自结构算调整后均值**，再相对计算 Score；调整后指标覆盖 `mean_epitope_load / n_strong / n_weak / density / immunogenicity_score`

**推荐单条命令**：
```bash
python3 tools/immunogenicity/immunogenicity_analyzer.py   --input <fasta> --donor <donor.fasta>   --structural-risk --structure-dir <pdb_dir>   --backmutation-dir <step3_dir>   --allele-panel extended --pep-len 15   --strong-threshold 0.5 --weak-threshold 0.3   --output <outdir> --all-formats
```

**关键参数**：
| 参数 | 默认 | 说明 |
|------|------|------|
| `--allele-panel` | `9common` | `9common`(官方示例) / `extended`(24 alleles，含 DQ/DP) |
| `--pep-len` | `15` | 9-30 |
| `--strong-threshold` | `0.5` | Strong 窗口阈值 |
| `--weak-threshold` | `0.3` | Weak 窗口阈值 |
| `--reference-variant` | auto | 无 donor 时的精确参考变体名 |
| `--structural-risk` | off | 开启结构校分并直接改写最终 Score |

**⚠️ 关键注意点**：
1. **HLAIIPred 运行环境**：见 `tools/immunogenicity/README.md`；若不可用自动回退 heuristic（仅做近似）
2. **Allele 批次**：HLAIIPred 单次 ≤14 alleles；`extended` 24 个自动分批并用 `0` 填充
3. **结构匹配**：每个变体需独立 PDB/CIF；文件名支持 `V0.pdb` / `combo_V0_V1.pdb` / `igfold_V0.pdb` 等多种模式
4. **Position labels**：必须提供 `variants_numbering.json`（Kabat/IMGT 逐残基标签），否则回退到序列索引（插入位点会错位）
5. **Score 语义改变**：结构校正后 `immunogenicity_score` 为“相对结构校正后参考均值的百分比”，原始序列 Score 保留在 `raw_mean_epitope_load` 等字段

**已修复的历史 Bug**：
1. `per_residue` key 类型不匹配 → 归一化为 int
2. 无 donor 时未自动用 V0 → 已实现精确 `V0 / combo_V0` fallback
3. 结构校正只看回突变位点 → 现为 **所有肽窗口位置 core-weighted**
4. Kabat 标签映射错误 → `variants_numbering.json` + 精确 label 匹配（兼容 82A/100B 等插入）
5. FR4 权重硬编码 0.9 → **移除 tier 权重**，改为 core/flank 加权
6. 缺结构静默给 0.5 → **neutral 1.0**，并标记错误
7. Strong/Weak 为残基计数 → **肽窗口计数**
### 5. TAP 可开发性分析 (tap_profiler)

**功能**：Therapeutic Antibody Profiler（Raybould 2019/2024）五项可开发性指标：
Ltot（总 CDR 长度）、PSH（表面疏水斑块）、PPC/PNC（正/负电荷斑块）、SFvCSP（VH/VL 表面电荷乘积）。
结构模式需要 ANARCI + FreeSASA（或 JOY `psa`）；序列模式只算 Ltot。

**用法**：
```bash
# 结构（默认 ANARCI IMGT 编号 + FreeSASA + Fv 作用域 + 2024 阈值）
python3 tools/tap_profiler/tap_analyzer.py \
  --pdb structure.pdb --output outputs/tap/

# 序列模式（仅 Ltot）
python3 tools/tap_profiler/tap_analyzer.py --vh <VH_SEQ> --vl <VL_SEQ>
```
详见 `tools/tap_profiler/README.md`。

---

### 6. 稳定性预测 (stability_predictor)

**功能**：使用 ThermoMPNN 预测回复突变对折叠稳定性的影响 (ΔΔG)，并映射到真实 Kabat/IMGT 位点。推荐在 **V0 结构**（人源框架 + 供体 CDR + 供体位点回复突变）上评估，此时直接测量 human→donor 效应，无需符号翻转。

**用法**：
```bash
# 推荐：直接评估 pipeline 的 T1/T2 回复突变（V0 结构，双链）
python3 tools/stability_predictor/stability_analyzer.py \
  --pdb <V0.pdb> --chain all \
  --pipeline-json <step3>/humanization_result.json \
  --tiers T1,T2 \
  --thermompnn-dir ~/ThermoMPNN \
  --python ~/miniconda3/envs/thermoMPNN/bin/python \
  --output <outdir> --all-formats

# 或从 CSV 读取回复突变
python3 tools/stability_predictor/stability_analyzer.py \
  --pdb <V0.pdb> --chain H \
  --backmutations-csv <step3>/backmutations_<name>.H.csv \
  --thermompnn-dir ~/ThermoMPNN \
  --python ~/miniconda3/envs/thermoMPNN/bin/python \
  --output <outdir> --json

# 通用 SSM：--positions 是真实 Kabat 位点（不是 ThermoMPNN 的 0-based index）
python3 tools/stability_predictor/stability_analyzer.py \
  --pdb <structure.pdb> --chain H \
  --positions H37,H44,H45,H67,H69 \
  --thermompnn-dir ~/ThermoMPNN \
  --python ~/miniconda3/envs/thermoMPNN/bin/python \
  --output <outdir> --json
```

**关键点**：
- `ddG_backmutation <= -0.2` 支持回复突变；`>= 1.0` 反对；中间为中性/谨慎。
- V2/V3 结构（WT=donor）测量的是 donor→human，输出用 `-ddG` 近似；**决策优先用 V0**。
- 每条链按单体计算，丢失 VH/VL 界面；Fab 请同时看 H/L（`--chain all`），界面位点以结构接触/埋藏证据为准。
- checkpoint 必须位于 `ThermoMPNN/models/`；`vanilla_model_weights/` 是 ProteinMPNN 权重，不会被使用。
- 详见 `tools/stability_predictor/README.md`。

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
4. **稳定性预测**：使用 ThermoMPNN 预测 ΔΔG，仅反映相对突变效应，不预测绝对折叠稳定性 (ΔG)；单链计算丢失 VH/VL 界面，且 V2/V3 上的符号翻转是近似；不替代实验验证

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
