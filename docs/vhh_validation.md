# VHH (camelid nanobody) humanization validation（驼源抗体人源化验证）

本文件回答：**这套流程对驼源 VHH 的人源化效果是否验证过、验证了什么、还缺什么。**

## 1. 结论（TL;DR）

- ✅ **正确性验证已覆盖**：编号、hallmark 识别与保护、含二硫键 CDR3 的移植守恒、
  端到端运行（1BZQ / 1MEL / cAb-Lys3）。
- ✅ **效果回测（金标准）已新增**：cAb-Lys3 → hCAb-Lys3（Vincke 2009 通用人源化
  纳米抗体支架设计）——**precision 0.778 / recall 1.0**，且管线自动选择 VH3-23
  家族 germline（正是支架的框架基础）。
- ⚠️ **局限**：hCAb-Lys3 序列为按论文文档化设计**重建**（序列级细节需对照原文
  复核后再用于生产）；caplacizumab 等已上市人源化 VHH 因亲本序列未公开，
  无法做完整亲本→人源化逐位点回测。

## 2. 已验证的 VHH 能力

| 能力 | 验证方式 | 结果 |
|------|----------|------|
| VHH 检测（hallmark 37/44/45/47） | 1BZQ / 1MEL / cAb-Lys3 | 4/4 命中 |
| hallmark 全程保护（KEEP_DONOR） | 全部 VHH 案例（V0-V3/Vmin 不变） | ✅ |
| 含二硫键 CDR3 移植守恒 | 1MEL / cAb-Lys3（CDR1-Cys + CDR3-Cys 桥） | ✅ 逐位点一致 |
| 完整 CDR3 环（H93-102）守恒 | 4D5/1BZQ/1MEL 回归测试 | ✅ |
| 人源化 VHH 的编号（FR2=14+S49 边界） | hCAb-Lys3 重建序列 | ✅ 引擎正确处理 |
| 端到端管线（选架/变体/报告） | 1BZQ 示例 + HumAb25 泛化 | ✅ |

## 3. VHH 金标准回测：cAb-Lys3 → hCAb-Lys3

### 3.1 案例数据

| 项 | 内容 | 来源 |
|----|------|------|
| 亲本 | cAb-Lys3（抗溶菌酶 VHH） | PDB 1MEL（结构验证） |
| 人源化 | hCAb-Lys3 | Vincke et al., *Nat Methods* 6:343 (2009) 通用支架设计（**序列为按论文文档化设计重建**） |
| 人源化设计 | IGHV3-23 共识框架 + **camelid hallmark 保留** + 供体 CDR（含 CDR1↔CDR3 二硫键桥） | Vincke 2009 |
| Mode B germline | IGHV3-23*01（支架的框架基础） | — |

### 3.2 结果（`python3 tests/backtest.py`）

```
VHH hallmark detected: True (4/4: H37F H44E H45R H47G)
hallmark positions kept KEEP_DONOR: H37 H44 H45 H47
germline (auto): IGHV3-23*01        ← 管线自动选中支架的框架基础
precision 0.778 | recall 1.0
over-reverted: H71, H78             ← 保守方向（实际支架保留 germline）
```

- **recall 1.0**：hCAb-Lys3 中全部 7 个"供体保留位点"被覆盖——
  hallmark 4 个（KEEP_DONOR 机制）+ CDR1 茎部 H27/H29/H30（vernier T1/T2）。
- **precision 0.778**：多回 H71/H78（vernier 位点，实际支架保留 germline）——
  与 Fab 案例一致的保守方向（安全侧）。
- **CDR1-Cys + CDR3-Cys 二硫键对完整移植**（驼源 VHH 结合所必需）。

### 3.3 与 Fab 金标准的一致性

| 指标 | cAb-Lys3 (VHH) | 曲妥珠 VL (Fab) | 贝伐珠 VH (Fab) |
|------|----------------|------------------|-----------------|
| 关键位点召回率 | 1.0 | 1.0 | 0.75-1.0 |
| 过回复 | H71/H78（保守） | 无（L87 已降级） | H28/L4（保守） |
| germline 匹配历史 | IGHV3-23 自动命中 | 历史选法不同但 CVI 支持 | 基本一致 |

## 4. 已知局限与后续

1. **hCAb-Lys3 为重建序列**：建议从 Vincke 2009 原文/附件核对 FR2（S49）、
   FR3、CDR3 逐个残基后再作为正式金标准；`tests/backtest.py` 已标注该假设。
2. **缺上市 VHH 亲本对**：caplacizumab（已上市人源化 VHH）亲本序列未公开；
   若获得授权数据可补充。
3. **VHH 免疫原性**：人源化后的 T 细胞表位评估依赖外部工具
   （netMHCIIpan/BioPhi，路线图 P2），序列层面无法自证。
4. **VHH FR2 长度规则**：引擎对"13 残基 FR2 + hallmark 内容"的工程序列
   （部分人源化 VHH 设计）仍依赖长度校验回退；当前 hCAb-Lys3（FR2=14+S49）
   已正确处理。

## 5. 复现

```bash
python3 tests/test_pipeline.py     # 含 test_vhh_humanization_gold_standard
python3 tests/backtest.py          # 含 cAb-Lys3 -> hCAb-Lys3 回测
python3 tests/backtest_scale.py    # HumAb25 泛化（Fab）
```
