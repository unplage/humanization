"""Domain-knowledge tables and scoring configuration for back-mutation design.

Position sets are in **Kabat numbering** (the scheme of the portable engine
and of the original structural literature):

  * Vernier zone        - Foote & Winter, JMB 224:487 (1992)
  * VH/VL interface     - Chothia et al., JMB 186:651 (1985);
                          Vargas-Madrazo & Paz-Garcia, JMB 330:783 (2003)
  * Canonical residues  - Chothia & Lesk, JMB 196:901 (1987);
                          Al-Lazikani et al., JMB 273:927 (1997)
  * VHH hallmark        - Muyldermans, J Biotechnol 74:277 (2001)

Only framework positions (FR1-FR3) are ever considered for back-mutation;
CDR residues are grafted from the donor and FR4 comes from the human J gene.
"""

# ---------------------------------------------------------------------------
# Position sets (Kabat)
# ---------------------------------------------------------------------------

VERNIER_ZONE = {
    "H": {2, 27, 28, 29, 30, 47, 48, 49, 67, 69, 71, 73, 78, 93, 94, 103},
    "L": {2, 4, 35, 36, 37, 38, 43, 44, 45, 46, 48, 49, 58, 62, 63, 66, 67, 68, 69,
          71, 87, 88, 98},
}

# Core VH/VL packing interface (Chothia 1985). Extended set adds positions
# with moderate interface contribution (Vargas-Madrazo 2003).
INTERFACE_CORE = {
    "H": {37, 39, 45, 47, 91, 93, 95, 103},
    "L": {34, 36, 37, 38, 44, 45, 46, 87, 89, 91, 96, 98},
}
INTERFACE_EXTENDED = {
    "H": {32, 34, 50, 58, 104},
    "L": {32, 50, 53},
}

# Canonical-structure framework support residues (Chothia & Lesk 1987).
CANONICAL = {
    "H": {24, 26, 27, 29, 33, 34, 47, 48, 49, 57, 58, 71, 78, 94},
    "L": {2, 25, 33, 34, 36, 37, 45, 46, 48, 49, 58, 64, 71, 90, 94, 97, 98},
}

# Camelid VHH hallmark (Kabat): positions that must stay donor to keep the
# single-domain fold (hydrophilic FR2 patch, no VL partner).
VHH_HALLMARK = {37, 44, 45, 47}

# J-region anchors (never back-mutated: human J provides FR4).
J_ANCHOR = {"H": 103, "L": 98}

# ---------------------------------------------------------------------------
# FR4 structural analysis (Step3 only - requires structure data)
# ---------------------------------------------------------------------------
# FR4 positions may need reversion to donor if they have structural importance
# (CDR3 contact, antigen contact). This is an exception to the "never back-mutate
# FR4" rule, applied only when structure data reveals critical contacts.

FR4_STRUCTURAL_WEIGHTS = {
    "cdr3_contact": 0.45,      # Position contacts CDR3 loop (<4.5 Å)
    "antigen_contact": 0.35,   # Position contacts antigen (<4.5 Å)
    "buried": 0.15,            # Position is buried (relSASA < 0.20)
    "interface_core": 0.05,    # Position in VH/VL interface core
}

# Threshold for recommending FR4 reversion (sum of weights must exceed this)
# Lowered to 0.4 so positions with just CDR3 contact (0.45) can be recommended
FR4_REVERSION_THRESHOLD = 0.4

# ---------------------------------------------------------------------------
# Scoring weights (documented in docs/scoring.md)
# ---------------------------------------------------------------------------

WEIGHTS = {
    "structural": {
        "interface_core": 1.00,
        "vernier": 0.85,
        "canonical": 0.80,
        "interface_extended": 0.65,
        "buried": 0.70,          # needs structure (AF3): relSASA < 0.20
        "cdr_contact": 0.85,     # needs structure: heavy-atom < 4.5 A to CDR
        "antigen_contact": 0.60, # needs structure: heavy-atom < 4.5 A to antigen
        "vhh_hallmark": 1.00,
        "disulfide_cys": 1.00,   # donor Cys paired with CDR Cys (VHH CDR3)
    },
    "immunogenicity": {
        "exposed_mismatch": 1.00,   # exposed non-human residue
        "buried_mismatch": 0.35,    # buried non-human residue (weak benefit)
    },
    "chemical": {
        "removes_nglycan": 0.80,              # NxS/T motif broken by reversion (高风险)
        "removes_deamidation_ng": 0.55,       # NG 热点脱酰胺 (Lu 2018: CDR H2/L1)
        "removes_deamidation_ns": 0.50,       # NS 中等脱酰胺
        "removes_deamidation_nh": 0.40,       # NH 中低风险
        "removes_deamidation_nd": 0.35,       # ND 低风险
        "removes_isomerization_dg": 0.50,     # DG 最常见异构化位点 (FRIDA 2024)
        "removes_isomerization_ds": 0.40,     # DS 中等异构化
        "removes_isomerization_dt": 0.40,     # DT 中等异构化
        "removes_isomerization_dh": 0.35,     # DH 中低风险
        "removes_acid_hydrolysis": 0.45,      # D-X (X=small residue) 肽键断裂
        "removes_acid_hydrolysis_dd": 0.55,   # DD 高风险：异构化+酸性水解 (PubMed 2013)
        "removes_oxidation": 0.30,            # M / W / C 氧化
        "removes_base_hydrolysis": 0.25,      # K-X (X=D/E) 碱性水解 (罕见)
        "removes_met_lyscleavage": 0.25,      # MK 金属蛋白酶裂解
        "introduces_nglycan": -0.80,          # reversion would CREATE an N-glycan motif (惩罚)
    },
    "blend": (0.55, 0.30, 0.15),    # structural, immunogenicity, chemical
}

# ---------------------------------------------------------------------------
# Empirically no-effect positions (from gold-standard backtests)
# ---------------------------------------------------------------------------
# Positions where the approved humanized antibody kept the HUMAN germline
# residue while retaining affinity, although the literature position sets
# would classify them as structural (interface/vernier). Without structure
# data supporting burial/CDR-contact, these are demoted to T3 (see
# docs/backtest_report.md - 4D5 -> trastuzumab: VL L87 F->Y kept germline,
# KD ~0.1 nM retained).
EMPIRICAL_NO_EFFECT = {
    "H": set(),
    "L": {87},   # L87: interface+vernier by literature; germline Y87 fine
}
EMPIRICAL_NO_EFFECT_NOTE = {
    "L87": "backtest gold standard (trastuzumab): kept germline Y87, "
           "affinity retained; demoted to T3 unless structure shows "
           "burial/CDR-contact (AF3 mode overrides this demotion).",
}

# ---------------------------------------------------------------------------
# Tier rules (final recommendation)
#
# NOTE: tier assignment is implemented authoritatively in
# scripts/humanize/backmut.py::_assign_tier (hard-coded from the literature
# position sets above). The TIER_RULES table was a declarative duplicate that
# drifted out of sync and is intentionally removed -- keep _assign_tier as the
# single source of truth.
# ---------------------------------------------------------------------------

TIER_LABELS = {
    "T1": "must revert (structural pillar)",
    "T2": "strongly recommended (structural/CDR support)",
    "T3": "optional (immunogenicity/consensus)",
    "T_FR4": "FR4 structural reversion (CDR3/antigen contact)",
    "KEEP_HUMAN": "keep human (no structural role, exposed)",
    "KEEP_DONOR": "keep donor (VHH hallmark / disulfide / CDR3 anchor)",
}

# ---------------------------------------------------------------------------
# Developability risk motifs (sequence level)
# ---------------------------------------------------------------------------
# WARNING: 以下风险基团仅基于序列模式检测，未考虑三维结构暴露状态。
# 埋藏在蛋白核心的位点实际风险较低，表面暴露位点风险更高。
# 如需精确评估，请结合 AF3/结构数据判断 relSASA 或溶剂可及性。
# 参考：Rajan et al., MAbs 2021; Boosted et al., J Pharm Sci 2022.

RISK_MOTIFS = {
    # === 糖基化风险 ===
    "N-glycan (NxS/T)": r"N[^P][ST]",           # N-X-S/T (X≠P) 被糖基化

    # === 脱酰胺风险 (Asn → Asp/isoAsp) ===
    "deamidation (NG)": r"NG",                   # 高风险：Asn-Gly
    "deamidation (NS)": r"NS",                   # 中风险：Asn-Ser
    "deamidation (NH)": r"NH",                   # 中风险：Asn-His（新增）
    "deamidation (ND)": r"ND",                   # 低风险：Asn-Asp（新增）

    # === 异构化风险 (Asp → isoAsp) ===
    "isomerization (DG)": r"DG",                 # 高风险：Asp-Gly
    "isomerization (DS)": r"DS",                 # 中风险：Asp-Ser（新增）
    "isomerization (DT)": r"DT",                 # 中风险：Asp-Thr（新增）
    "isomerization (DH)": r"DH",                 # 低风险：Asp-His（新增）

    # === 酸性水解风险 (Asp-X 肽键断裂) ===
    # D-X 仅统计未被异构化类覆盖的残基 (G/S/T/H 已在上方单独计分，
    # D 在 DD 单独计分)，避免同一位点双计。
    "acid hydrolysis (D-X)": r"D[AVLIP]",        # Asp 后接小侧链残基（不含 G/S/T/H/D）
    "acid hydrolysis (DD)": r"DD",               # Asp-Asp 高风险：酸性水解（新增）

    # === 氧化风险 ===
    "oxidation (M)": r"M",                       # 甲硫氨酸 → 亚砜/砜
    "oxidation (W)": r"W",                       # 色氨酸 → 氧化产物
    "oxidation (C)": r"C",                       # 半胱氨酸 → 磺基丙氨酸（新增）

    # === 二硫键/聚集风险 ===
    "unpaired Cys": r"C",                        # 未配对半胱氨酸 → 错误配对/聚集

    # === 碱性水解风险（新增）===
    "base hydrolysis (K-X)": r"K[DE]",           # Lys-Asp/Glu 肽键在碱性条件断裂

    # === Met-Lys 裂解（新增）===
    "met-lyscleavage (MK)": r"MK",               # 金属蛋白酶裂解位点
}


# ---------------------------------------------------------------------------
# High-risk developability positions (for MPNN optimization)
# ---------------------------------------------------------------------------

# Risk levels: high-risk positions should be optimized if possible
HIGH_RISK_MOTIFS = {
    "N-glycan": {"pattern": r"N[^P][ST]", "risk": "high", "issue": "glycosylation"},
    "deamidation_NG": {"pattern": r"NG", "risk": "high", "issue": "deamidation"},
    "deamidation_NS": {"pattern": r"NS", "risk": "medium", "issue": "deamidation"},
    "isomerization_DG": {"pattern": r"DG", "risk": "high", "issue": "isomerization"},
    "acid_hydrolysis_DD": {"pattern": r"DD", "risk": "high", "issue": "acid hydrolysis"},
    "oxidation_MW": {"pattern": r"[MW]", "risk": "medium", "issue": "oxidation"},
}


def detect_high_risk_positions(sequence: str, chain_type: str = "H") -> list:
    """Detect high-risk developability positions in a sequence.
    
    Returns list of dicts: [{"pos": int, "motif": str, "risk": str, "issue": str, "context": str}]
    """
    import re
    results = []
    for name, info in HIGH_RISK_MOTIFS.items():
        pattern = info["pattern"]
        for match in re.finditer(pattern, sequence):
            start = match.start()
            # Map sequence position to Kabat position (approximate)
            # This is simplified - in practice use numbering module
            context = sequence[max(0, start-2):min(len(sequence), start+len(match.group())+2)]
            results.append({
                "seq_pos": start,
                "motif": match.group(),
                "motif_name": name,
                "risk": info["risk"],
                "issue": info["issue"],
                "context": context,
            })
    return results
