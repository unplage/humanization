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
    "L": {2, 4, 35, 36, 38, 43, 44, 46, 48, 49, 58, 62, 63, 66, 67, 68, 69,
          71, 87, 88, 98},
}

# Core VH/VL packing interface (Chothia 1985). Extended set adds positions
# with moderate interface contribution (Vargas-Madrazo 2003).
INTERFACE_CORE = {
    "H": {37, 39, 45, 47, 91, 93, 95, 103},
    "L": {34, 36, 38, 44, 46, 87, 89, 91, 96, 98},
}
INTERFACE_EXTENDED = {
    "H": {32, 34, 50, 58, 104},
    "L": {32, 50, 53},
}

# Canonical-structure framework support residues (Chothia & Lesk 1987).
CANONICAL = {
    "H": {24, 26, 27, 29, 33, 34, 47, 48, 49, 57, 58, 71, 78, 94},
    "L": {2, 25, 33, 34, 36, 46, 48, 49, 58, 64, 71, 90, 94, 97, 98},
}

# Camelid VHH hallmark (Kabat): positions that must stay donor to keep the
# single-domain fold (hydrophilic FR2 patch, no VL partner).
VHH_HALLMARK = {37, 44, 45, 47}

# J-region anchors (never back-mutated: human J provides FR4).
J_ANCHOR = {"H": 103, "L": 98}

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
        "removes_nglycan": 0.8,     # NxS/T motif broken by reversion
        "removes_deamidation": 0.5, # NG / NS
        "removes_isomerization": 0.4,  # DG
        "removes_oxidation": 0.3,   # M / W in exposed loops
        "introduces_nglycan": -0.8, # reversion would CREATE an N-glycan motif
    },
    "blend": (0.55, 0.30, 0.15),    # structural, immunogenicity, chemical
}

# ---------------------------------------------------------------------------
# Tier rules (final recommendation)
# ---------------------------------------------------------------------------

TIER_RULES = {
    # (feature, requires_structure, requires_buried) -> tier
    "T1_must": {
        "interface_core": (True, False),
        "disulfide_cys": (False, False),   # keep donor! handled separately
    },
    "T1_if_buried": {
        "vernier": (True, True),
        "canonical": (True, True),
    },
    "T2_recommended": {
        "vernier": (False, False),
        "canonical": (False, False),
        "cdr_contact": (False, False),
        "antigen_contact": (False, False),
        "buried": (False, False),
        "interface_extended": (False, False),
    },
}

TIER_LABELS = {
    "T1": "must revert (structural pillar)",
    "T2": "strongly recommended (structural/CDR support)",
    "T3": "optional (immunogenicity/consensus)",
    "KEEP_HUMAN": "keep human (no structural role, exposed)",
    "KEEP_DONOR": "keep donor (VHH hallmark / disulfide / CDR3 anchor)",
}

# ---------------------------------------------------------------------------
# Developability risk motifs (sequence level)
# ---------------------------------------------------------------------------

RISK_MOTIFS = {
    "N-glycan (NxS/T)": r"N[^P][ST]",
    "deamidation (NG)": r"NG",
    "deamidation (NS)": r"NS",
    "isomerization (DG)": r"DG",
    "oxidation (M)": r"M",
    "oxidation (W)": r"W",
    "unpaired Cys": r"C",
}
