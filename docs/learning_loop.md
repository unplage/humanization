# Closed-loop learning from experimental data（实验数据闭环优化）

本流程支持用**实测实验数据持续校准评分模型**：每次人源化一轮实验后，
把（亲本序列、变体序列、亲和力 KD）记录成 JSON，运行 `humanize learn`
得到 `calibration.json`；后续所有 `humanize run` 自动加载它，用实测效应
修正回复突变分级与评分。

## 1. 数据格式（experiments.json）

```json
[
  {
    "name": "Ab1",
    "parent_vh": "EVQLQQSGPELVKPGASVKMSCKAS...",   // 亲本（非人源）VH
    "parent_vl": "DIQMTQTTSSLSASLGDRVTISC...",     // 亲本 VL（VHH 可省略）
    "parent_kd": 0.15,                              // 亲本 KD (nM)
    "variants": [
      { "name": "Ab1_hV0", "vh": "...", "vl": "...", "kd": 0.9 },
      { "name": "Ab1_hV2", "vh": "...", "vl": "...", "kd": 0.17 },
      { "name": "Ab1_hV2_H71A", "vh": "...", "vl": "...", "kd": 0.30 }
    ]
  },
  { "name": "Ab2", "parent_vh": "...", "parent_kd": 0.5,
    "variants": [ {"name": "Ab2_V2", "seq": "...", "kd": 0.6} ] }
]
```

- `seq` 用于单链（VHH）；Fab 用 `vh`+`vl` 两个键
- 可选字段：`tm`（DSF Tm1）、`expression`（mg/L）——后续版本可用于
  稳定性/表达量校准
- 一次实验内：**所有变体必须与同一亲本在同一 assay/格式下比较**
  （KD 比值才有意义）

## 2. 训练与使用

```bash
# 训练（输出 calibration.json）
python3 scripts/humanize/cli.py learn --experiments experiments.json \
  --out calibration.json

# 使用（后续所有人源化自动加载实测先验）
python3 scripts/humanize/cli.py run --input new_antibody.fasta \
  --outdir outputs --calibration calibration.json
```

报告自动新增 `empirical ddG` 列（kcal/mol，正值 = 保留 donor 有益）与
分级调整说明。

## 3. 统计模型（透明、可解释）

对每个框架位点 p（Kabat 坐标，FR1–FR3，CDR/J 除外）：

```
ΔΔG(变体) = RT·ln(KD变体 / KD亲本),   RT = 0.593 kcal/mol (25 °C)

effect(p) = mean(ΔΔG | 变体在该位点携带 donor 残基)
          − mean(ΔΔG | 变体在该位点携带 human 残基)
```

- **组间对比**（同一亲本内）：分离"该位点回复与否"对亲和力的净影响
- **收缩**：按变体数 n 向 0 收缩（`effect·n/(n+1)`），小样本不喧宾夺主
- 跨实验合并：多实验同一位点按样本量加权平均

分级调整阈值（`learning.py::effect_thresholds`）：

| ddG (kcal/mol) | 调整 |
|----------------|------|
| ≥ +0.41（KD 改善约 2 倍） | 强制 T1（保留 donor，score ≥ 70） |
| +0.20 ~ +0.41 | T3 → T2 提升 |
| −0.20 ~ +0.20 | 不变（n≥2 时记为"无效应"） |
| ≤ −0.20 | T1/T2 → T3 降级（score ≤ 40） |

## 4. 推荐实验设计（最大化学习效率）

| 轮次 | 变体面板 | 学习价值 |
|------|----------|----------|
| 第 1 轮 | V0（纯移植）、V1、V2、亲本 | 基线：V0 的活性损失量 = 整体框架效应；V2 确认 T1/T2 假设 |
| 第 2 轮 | V2 ± 单点回复（对 V2 中"可疑"位点各做一个单点变体，如 H71A、L87F） | **单点变体是逐位点效应的黄金标准**（其他位点固定） |
| 第 3 轮 | 组合面板（T3 位点子集，20–50 个，96 孔小试） | 覆盖 T3 空间，支持批量校准 |

要点：
- 至少包含 V0 与 V2；只有 V2 时学习到的效应是"批次联合效应"，区分度有限
- 单点变体每个成本低（1 个克隆），优先给 T1/T2 边缘位点（如 L87）
- KD 测量精度要求：同板、同批次、至少双复孔；KD 比值误差 < 1.5× 时
  ddG 分辨力约 ±0.24 kcal/mol —— 因此阈值定为 ±0.20
- 数据积累 ≥ 3 个抗体、≥ 30 个位点效应后，可把全库 merge 成
  "组织级校准文件"（跨项目通用先验）

## 5. 持续优化循环（组织级 SOP）

```
序列输入 ──► humanize run ──► 变体面板(V0-V3+Vmin+V_SDR)
                                    │ 基因合成/表达/BLI-SPR
                                    ▼
                        experiments.json（亲本+变体+KD）
                                    │
                                    ▼
                    humanize learn ──► calibration.json
                                    │  （每轮迭代覆盖）
                                    ▼
              下次 run 自动加载，评分随数据进化
```

- 每轮迭代后校验指标：V2 与亲本 KD 比值 ≤ 2×（若 >5× → 结构回炉）
- calibration.json 是**数据资产**：随项目归档，跨项目复用
- 用实测效应反向修正 `docs/theory.md` 中的位置表假设（例如：如果多个
  抗体都显示 L87 类位点无效应，从 interface 核心表降权）

## 6. 局限与注意事项

- 组间对比依赖变体间差异**主要在该位点**；高度耦合的位点（协同效应）
  需要单点/成对设计才能分离（组合面板轮次解决）
- 不同 assay（BLI vs SPR、不同捕获方式）的 KD 不可直接混合
- 效应是"该 germline 背景下的"效应；更换 germline 后部分效应会漂移
  （报告中注明 germline 版本）
- 免疫原性维度无法用 KD 学习——仍依赖 germline 匹配 + （可选）
  netMHCIIpan/BioPhi 类工具
