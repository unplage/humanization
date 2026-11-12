#!/usr/bin/env python3
"""Structure-aware isoelectric point (pI) and Fv charge-patch analysis.

This module is the reusable core behind ``pi_analyzer.py``.  It has two
independent jobs:

1. **pI calculation**
   * sequence-only Henderson-Hasselbalch titration with selectable pKa sets
     (Bjellqvist/ExPASy, IPC_protein, EMBOSS, Sillero, Nozaki-Tanford);
   * structure-aware pI: per-residue pKa values from PROPKA when available,
     otherwise a transparent microenvironment correction based on solvent
     exposure, salt bridges and hydrogen bonds.

2. **Fv charge-patch / developability context**
   * reuses the validated TAP engine in ``tools/tap_profiler/tap_core.py``
     (PDB/mmCIF parsing, ANARCI IMGT numbering, FreeSASA, PSH/PPC/PNC/SFvCSP);
   * clusters exposed charged residues into 3D patches and proposes
     framework-local, developability-oriented substitutions.

Scientific notes
----------------
pI is the pH where the modelled net charge is zero.  The sequence model
assumes every ionizable group behaves like an isolated model compound.  In a
folded protein, buried groups and salt bridges shift individual pKa values;
the structure-aware path approximates that effect.  It is a screening
heuristic, not a replacement for experimental cIEF or a rigorous continuum
electrostatics calculation.  The PROPKA path is preferred when PROPKA is
installed.

pKa sets are taken from the comparison collected by Kozlowski (IPC,
Nucleic Acids Res 2016/2021) and the original references.  Bjellqvist also
uses residue-specific terminal corrections, as in ExPASy.

Key references
--------------
* Bjellqvist B et al. Electrophoresis 1993;14:1023-1031.
* Kozlowski LP. Nucleic Acids Res 2016;44:W226-W231 (IPC).
* Kozlowski LP. Nucleic Acids Res 2021;49:W285-W291 (IPC 2.0).
* EMBOSS iep documentation; Sillero & Ribeiro 1989; Nozaki & Tanford 1967.
* Olsson MHM et al. J Chem Theory Comput 2011;7:525-537 (PROPKA3).
* Thorsteinson N et al. mAbs 2021;13:1981805 (structure-based pI / Fv charge).
* Raybould MIJ et al. PNAS 2019;116:4025-4030 (TAP).
* Gupta P et al. mAbs 2022;14:2073057 (Fv pI 7.5-9 window).
"""

from __future__ import annotations

import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

# ---------------------------------------------------------------------------
# Repo / TAP engine import
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(os.path.dirname(_HERE))
_TAP_DIR = os.path.join(_REPO_ROOT, "tools", "tap_profiler")
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
if _TAP_DIR not in sys.path:
    sys.path.insert(0, _TAP_DIR)

try:
    import tap_core  # type: ignore  # noqa: E402
except Exception as exc:  # pragma: no cover - import guard
    raise ImportError(
        "pi_core requires tools/tap_profiler/tap_core.py. "
        "Run from the repository root and keep the tools/ tree intact."
    ) from exc


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ACIDIC = ("D", "E", "C", "Y")
BASIC = ("H", "K", "R")
IONIZABLE = set(ACIDIC) | set(BASIC)

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# pKa tables.  Keys: Nterm, Cterm, D, E, C, Y, H, K, R.
# Bjellqvist additionally carries residue-specific terminal corrections used
# by ExPASy's Compute pI/Mw.
PKA_SETS: Dict[str, Dict[str, object]] = {
    "Bjellqvist": {
        "Nterm": 7.5, "Cterm": 3.55,
        "D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0,
        "H": 5.98, "K": 10.0, "R": 12.0,
        "Nterm_overrides": {
            "A": 7.59, "M": 7.00, "S": 6.93, "P": 8.36,
            "T": 6.82, "V": 7.44, "E": 7.70,
        },
        "Cterm_overrides": {"D": 4.55, "E": 4.75},
    },
    "IPC_protein": {
        "Nterm": 9.094, "Cterm": 2.869,
        "D": 3.872, "E": 4.412, "C": 7.555, "Y": 10.85,
        "H": 5.637, "K": 9.052, "R": 11.84,
    },
    "EMBOSS": {
        "Nterm": 8.6, "Cterm": 3.6,
        "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1,
        "H": 6.5, "K": 10.8, "R": 12.5,
    },
    "Sillero": {
        "Nterm": 8.2, "Cterm": 3.2,
        "D": 4.0, "E": 4.5, "C": 9.0, "Y": 10.0,
        "H": 6.4, "K": 10.4, "R": 12.0,
    },
    "Nozaki_Tanford": {
        "Nterm": 7.5, "Cterm": 3.8,
        "D": 4.0, "E": 4.4, "C": 9.5, "Y": 9.6,
        "H": 6.3, "K": 10.4, "R": 12.0,
    },
}

PKA_SET_ALIASES = {
    "bjellqvist": "Bjellqvist",
    "expasy": "Bjellqvist",
    "exPASy": "Bjellqvist",
    "ipc": "IPC_protein",
    "ipc_protein": "IPC_protein",
    "emboss": "EMBOSS",
    "sillero": "Sillero",
    "nozaki": "Nozaki_Tanford",
    "nozaki_tanford": "Nozaki_Tanford",
    "nozaki-tanford": "Nozaki_Tanford",
}

# IMGT region boundaries used by the pipeline / TAP.
IMGT_REGIONS = (
    ("FR1", 1, 26),
    ("CDR1", 27, 38),
    ("FR2", 39, 55),
    ("CDR2", 56, 65),
    ("FR3", 66, 104),
    ("CDR3", 105, 117),
    ("FR4", 118, 128),
)
IMGT_CDR_RANGES = ((27, 38), (56, 65), (105, 117))
IMGT_ANCHOR_POS = (25, 26, 39, 40, 54, 55, 66, 67, 103, 104, 118, 119)

# Camelid VHH hallmark positions in IMGT numbering (Kabat 37/44/45/47).
VHH_HALLMARK_IMGT = {42, 49, 50, 52}

# Conserved intradomain disulfide Cys in Kabat numbering.
#
# The intradomain disulfide of an Ig V domain connects the FR1 1st-CYS with
# the FR3 2nd-CYS.  In Kabat numbering these are H22/H92 (heavy) and L23/L88
# (light).  Lambda light chains have their FR1 Cys at Kabat L22; the repo's
# numbering engine (scripts/humanize/numbering.py) places the lambda FR1 Cys
# on L23 by applying the canonical Kabat L10 deletion, so the *engine's*
# output for a lambda chain is L23 as well.  Both L22 and L23 are therefore
# accepted here, and the pairing requirement below decides whether a given
# Cys is really the disulfide partner.
CONSERVED_CYS_KABAT = {"H": {22, 92}, "L": {22, 23, 88}}

# pI interpretation windows (Fv).  The 7.5-9.0 window follows Gupta et al.
# 2022 (mAbs) and Thorsteinson et al. 2021 (mAbs).
PI_LOW = 6.5
PI_NEUTRAL = 7.5
PI_FAVORABLE_HIGH = 9.0

# Exposure / microenvironment thresholds.
SURFACE_REL_SASA = 7.5      # %, TAP definition
INTERFACE_CUTOFF = 4.5      # A, VH/VL interface heavy-atom cutoff
BURIED_REL_SASA = 10.0      # %, pKa microenvironment
SOLVATED_REL_SASA = 25.0    # %, above this no burial pKa penalty
SALT_BRIDGE_CUTOFF = 4.0    # A, donor N ... acceptor O (pKa microenvironment)
HBOND_CUTOFF = 3.5          # A, polar N/O ... N/O
HBOND_ANGLE_MIN = 120.0     # deg, antecedent-donor...acceptor angle
DISULFIDE_CUTOFF = 2.5      # A, SG ... SG
PATCH_CUTOFF = 8.0          # A, exposed charged-residue patch graph

# Terminus modelling.  A Fab/VHH V-domain C-terminus is normally fused to a
# constant domain or tag, so the "auto" default caps it for Fab and leaves it
# free for a single-domain VHH.
TERMINI_MODES = ("auto", "free", "capped", "blocked")
TRUE_TERMINUS_N_IMGT = 1
# The last FR4-IMGT position of the V domain depends on the chain type:
# a VH/VHH runs to IMGT 128 (FR4-IMGT 118-128) while a VL ends at IMGT 127
# (FR4-IMGT 118-127).  ANARCI reproduces exactly that: 4D5 VH -> H128,
# 4D5 VL -> L127, a VHH -> H128.  Using one threshold for both chain types
# either truncates a complete VL or accepts a VH that lost its final
# residue, so the threshold is resolved per chain type below.
TRUE_TERMINUS_C_IMGT = {"H": 128, "L": 127}

# Ionic strength / Debye-Hückel screening (approximate).
IONIC_STRENGTH_DEFAULT = 0.15   # M, physiological
DEBYE_KAPPA_COEFF = 0.3291      # 1/A at 298 K for I in mol/L

# Local electrostatic patch score bands.
#
# The local score is an *unordered* fractional-charge sum
# (sum of |q_i*q_j|/r^2 over like-charged pairs inside one patch), whereas
# the published TAP PPC/PNC values are *ordered* integer-charge sums over the
# whole CDR vicinity.  The two scales are therefore not interchangeable and
# the bands below are NOT the TAP thresholds.
#
# The default bands are published in the tool README together with the
# reasoning; they can be replaced at runtime from a JSON calibration file
# (--patch-calibration) so that a site can be re-fitted against its own
# experimental data instead of relying on a generic prior.
PATCH_RISK_POS = {"elevated": 0.50, "high": 2.00}
PATCH_RISK_NEG = {"elevated": 0.60, "high": 2.50}
PATCH_CDR_FACTOR = 1.25

# Structure quality.
PLDDT_WARN = 70.0           # AF3-style pLDDT stored in the B-factor column
CA_BREAK_CUTOFF = 4.5       # A, consecutive CA distance indicating a chain break
MISSING_SIDECHAIN_FRACTION = 0.5
# Crystallographic B-factors are NOT pLDDT.  A refined crystal structure has
# low B-factors for well-ordered residues and large ones for disordered
# loops; interpreting a low mean B-factor as "low confidence" inverts the
# meaning.  The threshold below is the conventional "poorly ordered"
# boundary for experimental structures.
CRYSTAL_BFACTOR_WARN = 60.0  # A, mean CA B-factor above which a refined
                              # structure is considered partly disordered

# Expected heavy side-chain atoms (used only for completeness warnings).
SIDECHAIN_HEAVY = {
    "R": ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "N": ("CB", "CG", "OD1", "ND2"),
    "D": ("CB", "CG", "OD1", "OD2"),
    "Q": ("CB", "CG", "CD", "OE1", "NE2"),
    "E": ("CB", "CG", "CD", "OE1", "OE2"),
    "H": ("CB", "CG", "ND1", "CD2", "CE1", "NE2"),
    "I": ("CB", "CG1", "CG2", "CD1"),
    "L": ("CB", "CG", "CD1", "CD2"),
    "K": ("CB", "CG", "CD", "CE", "NZ"),
    "M": ("CB", "CG", "SD", "CE"),
    "F": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "P": ("CB", "CG", "CD"),
    "S": ("CB", "OG"),
    "T": ("CB", "OG1", "CG2"),
    "W": ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "Y": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "V": ("CB", "CG1", "CG2"),
}

# Atoms used to define the angle at a polar contact (antecedent -> donor or
# acceptor).  Without explicit hydrogens this is the best available proxy:
# the antecedent is the heavy atom bonded to the polar atom.
HBOND_ANTECEDENT = {
    "D": {"OD1": "CG", "OD2": "CG"},
    "E": {"OE1": "CD", "OE2": "CD"},
    "C": {"SG": "CB"},
    "Y": {"OH": "CZ"},
    "H": {"ND1": "CG", "NE2": "CD2"},
    "K": {"NZ": "CE"},
    "R": {"NE": "CD", "NH1": "CZ", "NH2": "CZ"},
    "N": {"OD1": "CG", "ND2": "CG"},
    "Q": {"OE1": "CD", "NE2": "CD"},
    "S": {"OG": "CB"},
    "T": {"OG1": "CB"},
    "W": {"NE1": "CD1"},
}

# Hydrogen-bond donors and acceptors among the ionizable side chains.
#
# A hydrogen bond needs a donor (an electronegative atom bearing a hydrogen)
# and a distinct acceptor (an electronegative atom with a lone pair).  The
# previous implementation accepted any N/O pair, so it counted D-D and K-K
# contacts as hydrogen bonds.  These tables make the polarity explicit:
#
#   * carboxylate O (D/E) and thiolate/thiophenolate S/Y-O are acceptors;
#   * the guanidinium / ammonium / imidazole / amide N atoms are donors;
#   * His ND1/NE2 and Tyr OH are both donor and acceptor.
#
# Backbone N is a donor and backbone O is an acceptor, but backbone contacts
# are handled separately (they are excluded for sequential residues).
HBOND_DONOR_ATOMS = {
    "K": ("NZ",), "R": ("NE", "NH1", "NH2"), "H": ("ND1", "NE2"),
    "N": ("ND2",), "Q": ("NE2",), "S": ("OG",), "T": ("OG1",),
    "W": ("NE1",), "Y": ("OH",),
}
HBOND_ACCEPTOR_ATOMS = {
    "D": ("OD1", "OD2"), "E": ("OE1", "OE2"), "C": ("SG",),
    "H": ("ND1", "NE2"), "N": ("OD1",), "Q": ("OE1",),
    "S": ("OG",), "T": ("OG1",), "Y": ("OH",),
}

# Side-chain atoms used for microenvironment detection.
ACID_ATOMS = {
    "D": ("OD1", "OD2"),
    "E": ("OE1", "OE2"),
    "C": ("SG",),
    "Y": ("OH",),
}
BASE_ATOMS = {
    "H": ("ND1", "NE2"),
    "K": ("NZ",),
    "R": ("NE", "NH1", "NH2"),
}
SALT_DONORS = {"K": ("NZ",), "R": ("NE", "NH1", "NH2"), "H": ("ND1", "NE2")}
SALT_ACCEPTORS = {"D": ("OD1", "OD2"), "E": ("OE1", "OE2")}

# Conservative neutral-polar substitutions for charge-patch engineering.
POSITIVE_TO_NEUTRAL = ("Q", "N", "S", "T", "A")
NEGATIVE_TO_NEUTRAL = ("Q", "N", "S", "T", "A")
POSITIVE_TO_NEGATIVE = ("E", "D")
NEGATIVE_TO_POSITIVE = ("K", "R")

# Outputs that cannot be produced without a three-dimensional structure.
#
# In sequence-only mode these are absent from the result.  They are listed
# explicitly (in the JSON and in the reports) because an absent flag is
# easily mistaken for a clean one: an antibody with a real hydrophobic-patch
# or charge-asymmetry flag is indistinguishable from a clean antibody in a
# sequence-only run, so a silent omission would read as an all-clear.
#
# The reason given for each entry is the physical reason the metric is
# structure-dependent, not merely a statement that it is unavailable.
STRUCTURE_REQUIRED_METRICS: Tuple[Tuple[str, str, str], ...] = (
    ("structure_pI",
     "Structure-aware pI",
     "needs the solvent exposure and the salt-bridge / hydrogen-bond "
     "microenvironment of every ionizable group"),
    ("tap",
     "TAP flags (PSH, PPC, PNC, SFvCSP)",
     "all four metrics are defined on solvent-exposed charges and 3D "
     "distances; only Ltot (CDR length) is sequence-determinable"),
    ("patches",
     "3D charge patches and the mutation candidates derived from them",
     "patches are 3D proximity clusters of exposed charged residues"),
    ("ionizable_residues",
     "Per-residue microenvironment table",
     "per-residue exposure, salt-bridge and hydrogen-bond assignments are "
     "properties of the folded structure"),
    ("structure_quality",
     "Structure quality checks",
     "B-factor semantics, missing side chains and chain breaks describe "
     "the deposited model"),
)


def unavailable_metrics() -> List[Dict[str, str]]:
    """The structure-dependent outputs a sequence-only run cannot produce."""
    return [
        {"key": key, "label": label, "reason": reason}
        for key, label, reason in STRUCTURE_REQUIRED_METRICS
    ]


# ---------------------------------------------------------------------------
# Small data model
# ---------------------------------------------------------------------------

@dataclass
class ChainInfo:
    """One antibody V domain used for pI calculation."""
    name: str                    # display name, e.g. VH / VL / VHH
    chain_type: str              # H or L
    sequence: str
    structure_chain: Optional[str] = None
    excluded_cys: Set[int] = field(default_factory=set)  # 0-based seq indices

    @property
    def nterm_aa(self) -> str:
        return self.sequence[0].upper() if self.sequence else ""

    @property
    def cterm_aa(self) -> str:
        return self.sequence[-1].upper() if self.sequence else ""


@dataclass
class IonizableGroup:
    """One ionizable side chain or terminus."""
    chain: str
    resseq: int
    icode: str
    aa: str
    kind: str                    # "acid" or "base"
    pka_model: float
    pka_eff: float
    active: bool = True
    rel_sasa: Optional[float] = None   # fraction 0-1
    exposure: str = "n/a"
    salt_bridge: bool = False
    hbond: bool = False
    disulfide: bool = False
    region: str = ""
    label: str = ""
    structure_chain: Optional[str] = None
    shift_reason: str = ""

    def charge(self, ph: float) -> float:
        if not self.active:
            return 0.0
        if self.kind == "acid":
            return -1.0 / (1.0 + 10.0 ** (self.pka_eff - ph))
        return 1.0 / (1.0 + 10.0 ** (ph - self.pka_eff))


@dataclass
class Terminus:
    chain: str
    kind: str                    # "N" or "C"
    aa: str
    pka: float
    blocked: bool = False
    truncated: bool = False      # True when the modelled chain is not the true end
    source: str = "model"        # "model", "propka", or "truncated"

    def charge(self, ph: float) -> float:
        if self.blocked:
            return 0.0
        if self.kind == "N":
            return 1.0 / (1.0 + 10.0 ** (ph - self.pka))
        return -1.0 / (1.0 + 10.0 ** (self.pka - ph))


@dataclass
class Patch:
    """A 3D cluster of exposed charged residues."""
    patch_id: int
    sign: str                    # "positive" / "negative" / "mixed"
    ph: float
    residues: List[Dict[str, object]] = field(default_factory=list)
    total_charge: float = 0.0
    score: float = 0.0            # local unordered fractional-charge sum
    score_tap: float = 0.0        # TAP-style ordered integer-charge sum
    risk_basis: str = ""
    terminal_residues: List[str] = field(default_factory=list)
    n_cdr: int = 0
    n_framework: int = 0
    n_anchor: int = 0
    in_cdr_vicinity: bool = False
    chains: List[str] = field(default_factory=list)
    risk: str = "low"

    @property
    def labels(self) -> List[str]:
        return [str(r["label"]) for r in self.residues]


# ---------------------------------------------------------------------------
# Sequence / pKa helpers
# ---------------------------------------------------------------------------

def clean_sequence(sequence: str) -> Tuple[str, List[str]]:
    """Uppercase a sequence and map rare amino acids to titration analogues.

    U (selenocysteine) -> C, O (pyrrolysine) -> K, B/Z/X are ignored.
    Returns (clean_sequence, warnings).
    """
    warnings: List[str] = []
    raw = re.sub(r"[^A-Za-z]", "", sequence or "").upper()
    out: List[str] = []
    for aa in raw:
        if aa in "ACDEFGHIKLMNPQRSTVWY":
            out.append(aa)
        elif aa == "U":
            out.append("C")
            warnings.append("U (selenocysteine) treated as C for pI")
        elif aa == "O":
            out.append("K")
            warnings.append("O (pyrrolysine) treated as K for pI")
        else:
            warnings.append("ignored non-standard residue {!r}".format(aa))
    return "".join(out), warnings


def normalize_pka_set(name: str) -> str:
    if name in PKA_SETS:
        return name
    key = str(name or "").strip()
    if key in PKA_SET_ALIASES:
        return PKA_SET_ALIASES[key]
    low = key.lower().replace("-", "_")
    if low in PKA_SET_ALIASES:
        return PKA_SET_ALIASES[low]
    raise ValueError(
        "Unknown pKa set {!r}; choose from {}".format(name, ", ".join(PKA_SETS))
    )


def _pka_value(pka_set: str, aa: str, terminal: Optional[str] = None) -> float:
    """Return the model pKa for a side chain or a terminus."""
    table = PKA_SETS[pka_set]
    if terminal == "N":
        overrides = table.get("Nterm_overrides", {})
        return float(overrides.get(aa, table["Nterm"]))  # type: ignore[union-attr]
    if terminal == "C":
        overrides = table.get("Cterm_overrides", {})
        return float(overrides.get(aa, table["Cterm"]))  # type: ignore[union-attr]
    return float(table[aa])  # type: ignore[index]


def resolve_termini_mode(mode: str, kind: str) -> str:
    """Resolve the context-aware default terminus model.

    A Fab VH/VL C-terminus is fused to CH1/CL, so the neutral default is a
    capped C-terminus.  A single-domain VHH has a free C-terminus unless it is
    fused to a tag.
    """
    if mode not in TERMINI_MODES:
        raise ValueError(
            "--termini must be one of {}".format(", ".join(TERMINI_MODES)))
    if mode != "auto":
        return mode
    return "capped" if kind == "fab" else "free"


def debye_kappa(ionic_strength: float) -> float:
    """Approximate inverse Debye length in 1/A (298 K, I in mol/L)."""
    if ionic_strength is None or ionic_strength <= 0.0:
        return 0.0
    return DEBYE_KAPPA_COEFF * math.sqrt(float(ionic_strength))


def _screening_factor(distance: float, ionic_strength: float) -> float:
    """Distance-dependent electrostatic screening factor in (0, 1]."""
    kappa = debye_kappa(ionic_strength)
    if kappa <= 0.0:
        return 1.0
    return math.exp(-kappa * max(0.0, float(distance)))


def _angle_deg(a: Sequence[float], b: Sequence[float],
               c: Sequence[float]) -> float:
    """Angle ABC in degrees for three 3-D points."""
    v1 = (a[0] - b[0], a[1] - b[1], a[2] - b[2])
    v2 = (c[0] - b[0], c[1] - b[1], c[2] - b[2])
    n1 = math.sqrt(sum(v * v for v in v1))
    n2 = math.sqrt(sum(v * v for v in v2))
    if n1 <= 1e-9 or n2 <= 1e-9:
        return 0.0
    dot = sum(x * y for x, y in zip(v1, v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (n1 * n2)))))


def _region_imgt(pos: Optional[int]) -> str:
    """IMGT region for a position, or "unknown" if it is out of range.

    A position outside FR1..FR4 (1-128) is not a valid IMGT position;
    labelling it "FR4" would silently mis-classify it and, in the mutation
    filter, would wrongly exclude it from consideration.
    """
    if pos is None:
        return ""
    for name, lo, hi in IMGT_REGIONS:
        if lo <= pos <= hi:
            return name
    return "unknown"


def _is_cdr_imgt(pos: Optional[int]) -> bool:
    return pos is not None and any(lo <= pos <= hi for lo, hi in IMGT_CDR_RANGES)


def _is_anchor_imgt(pos: Optional[int]) -> bool:
    return pos is not None and pos in IMGT_ANCHOR_POS


def _make_groups_from_sequence(
    chains: Sequence[ChainInfo],
    pka_set: str,
) -> Tuple[List[IonizableGroup], List[Terminus]]:
    """Build per-residue groups and termini for sequence-only pI."""
    groups: List[IonizableGroup] = []
    termini: List[Terminus] = []
    for chain in chains:
        seq = chain.sequence.upper()
        for idx, aa in enumerate(seq):
            if aa not in IONIZABLE:
                continue
            if idx in chain.excluded_cys:
                continue
            kind = "acid" if aa in ACIDIC else "base"
            pka = _pka_value(pka_set, aa)
            groups.append(IonizableGroup(
                chain=chain.name, resseq=idx + 1, icode="", aa=aa,
                kind=kind, pka_model=pka, pka_eff=pka,
                region="", label="{}:{}".format(chain.name, idx + 1),
                structure_chain=chain.structure_chain,
            ))
        if seq:
            termini.append(Terminus(
                chain=chain.name, kind="N", aa=seq[0],
                pka=_pka_value(pka_set, seq[0], "N"),
            ))
            termini.append(Terminus(
                chain=chain.name, kind="C", aa=seq[-1],
                pka=_pka_value(pka_set, seq[-1], "C"),
            ))
    return groups, termini


def net_charge(groups: Iterable[IonizableGroup],
               termini: Iterable[Terminus], ph: float,
               capped_cterm: bool = False, blocked_nterm: bool = False) -> float:
    q = sum(g.charge(ph) for g in groups)
    for t in termini:
        blocked = t.blocked or (t.kind == "C" and capped_cterm) or (
            t.kind == "N" and blocked_nterm)
        if not blocked:
            q += t.charge(ph)
    return q


def isoelectric_point(groups: Sequence[IonizableGroup],
                      termini: Sequence[Terminus],
                      capped_cterm: bool = False,
                      blocked_nterm: bool = False,
                      lo: float = 0.0, hi: float = 14.0) -> Optional[float]:
    """Bisection on the monotone net-charge curve."""
    def q(ph: float) -> float:
        return net_charge(groups, termini, ph, capped_cterm, blocked_nterm)

    q_lo = q(lo)
    q_hi = q(hi)
    if q_lo <= 0.0:
        return lo
    if q_hi >= 0.0:
        return hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        q_mid = q(mid)
        if abs(q_mid) < 1e-7:
            return mid
        if q_mid > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sequence_pi_table(
    chains: Sequence[ChainInfo],
    pka_set: str,
    ph: float = 7.4,
    formulation_ph: float = 6.0,
    termini_mode: str = "auto",
    kind: Optional[str] = None,
) -> Dict[str, object]:
    """pI for the combined chains under free/capped/blocked/selected termini."""
    groups, termini = _make_groups_from_sequence(chains, pka_set)
    if kind is None:
        kind = "fab" if len(chains) >= 2 else "vhh"
    resolved = resolve_termini_mode(termini_mode, kind)
    capped = resolved in ("capped", "blocked")
    n_blocked = resolved == "blocked"
    return {
        "pI_free": isoelectric_point(groups, termini),
        "pI_capped": isoelectric_point(groups, termini, capped_cterm=True),
        "pI_blocked": isoelectric_point(groups, termini, capped_cterm=True,
                                        blocked_nterm=True),
        "pI_selected": isoelectric_point(
            groups, termini, capped_cterm=capped, blocked_nterm=n_blocked),
        "termini_mode": resolved,
        "net_charge_ph74_free": net_charge(groups, termini, ph),
        "net_charge_ph74_capped": net_charge(groups, termini, ph,
                                             capped_cterm=True),
        "net_charge_selected_ph": net_charge(
            groups, termini, ph, capped_cterm=capped, blocked_nterm=n_blocked),
        "net_charge_selected_formulation_ph": net_charge(
            groups, termini, formulation_ph,
            capped_cterm=capped, blocked_nterm=n_blocked),
    }


def chain_pi_table(
    chains: Sequence[ChainInfo],
    pka_set: str,
    ph: float = 7.4,
    formulation_ph: float = 6.0,
    termini_mode: str = "free",
) -> Dict[str, Dict[str, object]]:
    """Per-chain pI and charge at the resolved terminus model."""
    capped = termini_mode in ("capped", "blocked")
    n_blocked = termini_mode == "blocked"
    out: Dict[str, Dict[str, object]] = {}
    for chain in chains:
        groups, termini = _make_groups_from_sequence([chain], pka_set)
        out[chain.name] = {
            "chain_type": chain.chain_type,
            "length": len(chain.sequence),
            "sequence": chain.sequence,
            "pI_free": isoelectric_point(groups, termini),
            "pI_capped": isoelectric_point(groups, termini, capped_cterm=True),
            "pI_selected": isoelectric_point(
                groups, termini, capped_cterm=capped, blocked_nterm=n_blocked),
            "net_charge_ph74": net_charge(
                groups, termini, ph, capped_cterm=capped,
                blocked_nterm=n_blocked),
            "net_charge_ph74_free": net_charge(groups, termini, ph),
            "net_charge_formulation_ph": net_charge(
                groups, termini, formulation_ph, capped_cterm=capped,
                blocked_nterm=n_blocked),
            "n_ionizable": len(groups),
        }
    return out


def _true_c_terminus(chain_type: Optional[str]) -> int:
    """Last FR4-IMGT position that a complete V domain of this type reaches.

    ANARCI numbers a VH/VHH to H128 and a VL to L127 (FR4-IMGT 118-128 vs
    118-127).  A chain whose last modelled residue is below that boundary
    is truncated and its terminal group must not be treated as a genuine
    chain end.
    """
    if chain_type == "L":
        return TRUE_TERMINUS_C_IMGT["L"]
    # Default to the heavy-chain boundary: it is the stricter of the two, so
    # an unknown chain type is not silently accepted as complete.
    return TRUE_TERMINUS_C_IMGT["H"]


def default_ph_grid(start: float = 4.0, stop: float = 10.0,
                    step: float = 0.5) -> List[float]:
    """Inclusive pH grid used for the titration/charge curve."""
    if step <= 0:
        raise ValueError("pH grid step must be positive")
    values: List[float] = []
    ph = float(start)
    while ph <= stop + 1e-9:
        values.append(round(ph, 3))
        ph += step
    return values


def charge_curve(
    groups: Sequence[IonizableGroup],
    termini: Sequence[Terminus],
    ph_values: Sequence[float],
    capped_cterm: bool = False,
    blocked_nterm: bool = False,
) -> List[Dict[str, float]]:
    """Net charge sampled across pH (for formulation/IEX planning)."""
    return [
        {"ph": round(float(ph), 3),
         "net_charge": round(net_charge(
             groups, termini, float(ph),
             capped_cterm=capped_cterm, blocked_nterm=blocked_nterm), 4)}
        for ph in ph_values
    ]


def detect_structure_provenance(path: str) -> str:
    """Classify the B-factor semantics of a structure file.

    Returns ``"experimental"`` (B-factors are crystallographic displacement
    factors in A^2), ``"predicted"`` (B-factor column carries pLDDT, as in
    AlphaFold/AF3 models) or ``"unknown"``.

    The distinction matters: for a refined crystal structure a *low* mean
    B-factor means the model is well ordered, while for a predicted model a
    *low* mean B-factor means low confidence.  Treating one as the other
    inverts the message, so the provenance is read from the file header
    instead of guessed from the numeric range.
    """
    try:
        with open(path, "r", errors="replace") as fh:
            head = []
            for _ in range(4000):
                line = fh.readline()
                if not line:
                    break
                head.append(line)
                if len(head) >= 400:
                    break
    except OSError:
        return "unknown"
    text = "".join(head)
    upper = text.upper()
    if ("EXPDTA" in upper or "REMARK   2 RESOLUTION" in upper
            or "_EXPTL.METHOD" in upper or "REMARK   3" in upper):
        return "experimental"
    if ("ALPHAFOLD" in upper or "ALPHAFOLD3" in upper
            or "_MA_MAJOR_VERSION" in upper or "PLDDT" in upper):
        return "predicted"
    return "unknown"


def assess_structure_quality(
    residues: Sequence["tap_core.Residue"],
    kind: str = "fab",
    chain_types: Optional[Dict[str, str]] = None,
    provenance: str = "unknown",
) -> Dict[str, object]:
    """Report model-quality features that affect pI/SASA reliability.

    The B-factor column is only interpreted through ``provenance``:

    * ``experimental``: the column holds crystallographic displacement
      factors (A^2).  Low values mean *well ordered*, so only a high mean
      (disordered loops) is reported as a caveat.
    * ``predicted``: the column holds pLDDT (0-100), so a low mean means
      low confidence.
    * ``unknown``: no confidence claim is made; the numbers are reported
      without interpretation.

    ``chain_types`` maps a chain id to "H"/"L" so that the truncated-term
    check can use the correct FR4-IMGT end position per chain type.
    """
    chain_types = chain_types or {}
    bfactors = [float(r.bfactor) for r in residues
                if getattr(r, "bfactor", None) is not None]
    mean_bf = sum(bfactors) / len(bfactors) if bfactors else None
    plddt_like = bool(bfactors) and all(0.0 <= b <= 100.0 for b in bfactors)
    if provenance == "experimental":
        bfactor_meaning = "crystallographic B-factor (A^2)"
    elif provenance == "predicted":
        bfactor_meaning = "predicted-model pLDDT"
    else:
        bfactor_meaning = "unknown"

    missing_sidechain: List[str] = []
    missing_ionizable: List[str] = []
    for res in residues:
        expected = SIDECHAIN_HEAVY.get(res.aa)
        if not expected:
            continue
        names = {raw.strip() for raw, _e, _x, _y, _z in res.atoms}
        present = sum(1 for name in expected if name in names)
        frac = present / float(len(expected))
        if frac < MISSING_SIDECHAIN_FRACTION:
            missing_sidechain.append(_res_label(res))
        if res.aa in IONIZABLE:
            ion_atoms = ACID_ATOMS.get(res.aa, ()) + BASE_ATOMS.get(res.aa, ())
            if ion_atoms and not any(name in names for name in ion_atoms):
                missing_ionizable.append(_res_label(res))

    chain_breaks: Dict[str, int] = {}
    by_chain: Dict[str, List["tap_core.Residue"]] = defaultdict(list)
    for res in residues:
        by_chain[res.chain].append(res)
    for chain, chain_res in by_chain.items():
        breaks = 0
        prev_ca: Optional[Tuple[float, float, float]] = None
        for res in chain_res:
            ca = _atom_map(res).get("CA")
            if ca is None:
                prev_ca = None
                continue
            if prev_ca is not None:
                d2 = sum((a - b) ** 2 for a, b in zip(ca, prev_ca))
                if math.sqrt(d2) > CA_BREAK_CUTOFF:
                    breaks += 1
            prev_ca = ca
        if breaks:
            chain_breaks[chain] = breaks

    truncated_termini: List[str] = []
    for chain, chain_res in by_chain.items():
        if not chain_res:
            continue
        first, last = chain_res[0], chain_res[-1]
        n_pos = getattr(first, "imgt_pos", None)
        c_pos = getattr(last, "imgt_pos", None)
        if n_pos is not None and int(n_pos) != TRUE_TERMINUS_N_IMGT:
            truncated_termini.append("{}:N(IMGT {})".format(chain, n_pos))
        c_limit = _true_c_terminus(chain_types.get(chain))
        if c_pos is not None and int(c_pos) < c_limit:
            truncated_termini.append(
                "{}:C(IMGT {}, expected >= {})".format(chain, c_pos, c_limit))

    return {
        "kind": kind,
        "n_residues": len(residues),
        "n_chains": len(by_chain),
        "provenance": provenance,
        "bfactor_meaning": bfactor_meaning,
        "bfactor_mean": None if mean_bf is None else round(mean_bf, 1),
        "bfactor_min": None if not bfactors else round(min(bfactors), 1),
        "bfactor_max": None if not bfactors else round(max(bfactors), 1),
        "bfactor_plddt_like": plddt_like,
        "n_missing_sidechain": len(missing_sidechain),
        "missing_sidechain_fraction": (
            round(len(missing_sidechain) / float(len(residues)), 3)
            if residues else 0.0),
        "missing_sidechain": missing_sidechain[:50],
        "n_missing_ionizable_atoms": len(missing_ionizable),
        "missing_ionizable_atoms": missing_ionizable[:50],
        "chain_breaks": chain_breaks,
        "n_chain_breaks": sum(chain_breaks.values()),
        "ca_break_cutoff": CA_BREAK_CUTOFF,
        "truncated_termini": truncated_termini,
    }


def structure_quality_warnings(quality: Dict[str, object]) -> List[str]:
    """Human-readable warnings derived from ``assess_structure_quality``."""
    warnings: List[str] = []
    mean_raw = quality.get("bfactor_mean")
    mean_bf = float(mean_raw) if isinstance(mean_raw, (int, float)) else None
    provenance = str(quality.get("provenance", "unknown"))
    if provenance == "predicted":
        # The B-factor column holds pLDDT: low values mean low confidence.
        if mean_bf is not None and mean_bf < PLDDT_WARN:
            warnings.append(
                "Low mean CA pLDDT ({:.1f}); structure-aware pI and charge "
                "patches may be unreliable. Prefer a higher-confidence "
                "model.".format(mean_bf))
    elif provenance == "experimental":
        # The B-factor column holds crystallographic displacement factors:
        # low values mean *well ordered*, so only large values are a problem.
        if mean_bf is not None and mean_bf > CRYSTAL_BFACTOR_WARN:
            warnings.append(
                "High mean CA B-factor ({:.1f} A^2); part of the model is "
                "poorly ordered, so burial, salt-bridge and H-bond "
                "assignments in those regions are uncertain.".format(mean_bf))
    else:
        if mean_bf is not None:
            warnings.append(
                "B-factor semantics could not be determined from the file "
                "header; the B-factor column is reported but not "
                "interpreted (mean {:.1f}).".format(mean_bf))
    missing_raw = quality.get("n_missing_sidechain")
    n_missing = int(missing_raw) if isinstance(missing_raw, (int, float)) else 0
    missing_list = quality.get("missing_sidechain")
    if not isinstance(missing_list, list):
        missing_list = []
    if n_missing:
        warnings.append(
            "{} residue(s) have incomplete side chains; SASA/salt-bridge/"
            "H-bond microenvironment assignments may be unreliable "
            "(first: {}).".format(
                n_missing,
                ", ".join(str(x) for x in missing_list[:8])))
    breaks_raw = quality.get("n_chain_breaks")
    n_breaks = int(breaks_raw) if isinstance(breaks_raw, (int, float)) else 0
    if n_breaks:
        warnings.append(
            "{} chain break(s) detected (CA-CA > {:.1f} A); missing loops can "
            "change burial and patch connectivity.".format(
                n_breaks, CA_BREAK_CUTOFF))
    truncated = quality.get("truncated_termini")
    if not isinstance(truncated, list):
        truncated = []
    if truncated:
        warnings.append(
            "Truncated V-domain terminus/termini detected ({}); no artificial "
            "free terminus charge was added at the cut.".format(
                ", ".join(str(x) for x in truncated)))
    return warnings


# ---------------------------------------------------------------------------
# Disulfide handling
# ---------------------------------------------------------------------------

def find_disulfide_cys_structure(residues: Sequence["tap_core.Residue"]) -> Set[Tuple[str, int, str]]:
    """Geometric SG-SG disulfide detection (<= 2.5 A)."""
    sulfurs: List[Tuple[Tuple[str, int, str], float, float, float]] = []
    for res in residues:
        if res.aa != "C":
            continue
        for raw_name, _elem, x, y, z in res.atoms:
            if raw_name.strip() == "SG":
                sulfurs.append((res.key, x, y, z))
                break
    found: Set[Tuple[str, int, str]] = set()
    cutoff2 = DISULFIDE_CUTOFF * DISULFIDE_CUTOFF
    for i in range(len(sulfurs)):
        k1, x1, y1, z1 = sulfurs[i]
        for j in range(i + 1, len(sulfurs)):
            k2, x2, y2, z2 = sulfurs[j]
            if k1[0] == k2[0] and k1[1] == k2[1] and k1[2] == k2[2]:
                continue
            d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            if d2 <= cutoff2:
                found.add(k1)
                found.add(k2)
    return found


def find_disulfide_cys_sequence(chains: Sequence[ChainInfo]) -> Tuple[Set[Tuple[str, int]], List[str]]:
    """Use the portable Kabat numbering engine to locate conserved Cys.

    The intradomain disulfide links the FR1 Cys with the FR3 Cys, so a
    conserved Cys is only assigned to the disulfide when *both* partners are
    present.  A single unpaired Cys (e.g. a light chain whose FR3 partner
    has been mutated away) is a free thiol and stays ionizable, which is
    what the titration model requires.

    Returns (set of (chain_name, 0-based index), warnings).
    """
    warnings: List[str] = []
    out: Set[Tuple[str, int]] = set()
    try:
        from scripts.humanize.numbering import number_heavy, number_light
    except Exception as exc:  # pragma: no cover
        warnings.append("cannot import numbering engine for disulfide detection: {}".format(exc))
        return out, warnings

    for chain in chains:
        try:
            numbered = (number_heavy(chain.sequence)
                        if chain.chain_type == "H"
                        else number_light(chain.sequence))
        except Exception as exc:
            warnings.append(
                "numbering failed for {} ({}); conserved Cys not auto-excluded".format(
                    chain.name, exc))
            continue
        conserved = CONSERVED_CYS_KABAT.get(chain.chain_type, set())
        hits: Dict[int, int] = {}
        for res in numbered.residues:
            if res.aa.upper() != "C":
                continue
            m = re.match(r"^[A-Z](\d+)([A-Za-z]*)$", res.pos)
            if not m:
                continue
            pos = int(m.group(1))
            if pos in conserved:
                hits[pos] = res.index
        if not hits:
            continue
        # The conserved positions form the FR1/FR3 pair; both must be
        # present for the disulfide to exist.
        if len(hits) >= 2:
            out.update((chain.name, idx) for idx in hits.values())
        else:
            only_pos = sorted(hits)[0]
            warnings.append(
                "{}: only one conserved Cys (Kabat {}) was found; it is "
                "treated as a free (ionizable) thiol rather than a "
                "disulfide partner.".format(chain.name, only_pos))
    return out, warnings


# ---------------------------------------------------------------------------
# Structure microenvironment: salt bridges, H-bonds, pKa shifts
# ---------------------------------------------------------------------------

def _atom_map(res: "tap_core.Residue") -> Dict[str, Tuple[float, float, float]]:
    return {raw.strip(): (x, y, z) for raw, _elem, x, y, z in res.atoms}


def detect_salt_bridges(
    residues: Sequence["tap_core.Residue"],
) -> Tuple[Set[Tuple[str, int, str]], List[Dict[str, object]],
           Dict[Tuple[str, int, str], float]]:
    """Detect cation ... anion salt bridges (donor N ... acceptor O <= 4 A).

    This is the pKa-microenvironment definition.  It intentionally differs
    from TAP's stricter 3.2 A / exposed / one-partner rule; both are reported
    so the two views cannot be confused.

    Histidine is included as a donor: a neutral imidazole ring is
    protonated at low pH and does form charge-assisted hydrogen bonds to
    carboxylates, and ignoring it would systematically mis-assign the pKa
    shift of buried acidic groups that are paired with His.  Multiple
    partners per residue are allowed (the broader definition), and the
    recorded distance is the *minimum* over all contacting atom pairs.
    """
    donors: List[Tuple[int, str, float, float, float]] = []
    acceptors: List[Tuple[int, str, float, float, float]] = []
    for idx, res in enumerate(residues):
        atoms = _atom_map(res)
        for name in SALT_DONORS.get(res.aa, ()):
            if name in atoms:
                x, y, z = atoms[name]
                donors.append((idx, name, x, y, z))
        for name in SALT_ACCEPTORS.get(res.aa, ()):
            if name in atoms:
                x, y, z = atoms[name]
                acceptors.append((idx, name, x, y, z))

    cutoff2 = SALT_BRIDGE_CUTOFF * SALT_BRIDGE_CUTOFF
    partners: Set[Tuple[str, int, str]] = set()
    pairs: List[Dict[str, object]] = []
    min_dist: Dict[Tuple[str, int, str], float] = {}
    best: Dict[Tuple[int, int], float] = {}
    for d_idx, _d_name, x1, y1, z1 in donors:
        for a_idx, _a_name, x2, y2, z2 in acceptors:
            if d_idx == a_idx:
                continue
            d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            if d2 > cutoff2:
                continue
            dist = math.sqrt(d2)
            pair = (d_idx, a_idx) if d_idx < a_idx else (a_idx, d_idx)
            if pair not in best or dist < best[pair]:
                best[pair] = dist
    for (d_idx, a_idx), dist in sorted(best.items()):
        r1, r2 = residues[d_idx], residues[a_idx]
        partners.add(r1.key)
        partners.add(r2.key)
        for key in (r1.key, r2.key):
            if key not in min_dist or dist < min_dist[key]:
                min_dist[key] = dist
        pairs.append({
            "donor": _res_label(r1),
            "acceptor": _res_label(r2),
            "distance": round(dist, 2),
        })
    return partners, pairs, min_dist


def detect_hbonds(
    residues: Sequence["tap_core.Residue"],
    salt_bridge_keys: Set[Tuple[str, int, str]],
) -> Tuple[Set[Tuple[str, int, str]],
           Dict[Tuple[str, int, str], float]]:
    """Flag ionizable side chains engaged in a hydrogen bond.

    A hydrogen bond requires a donor (electronegative atom bearing H) and a
    *different* acceptor (electronegative atom with a lone pair).  The
    implementation therefore:

      * uses the ionizable side-chain atom itself, not any residue N/O;
      * requires donor/acceptor polarity complementarity, so a D-D or K-K
        contact is no longer counted as a hydrogen bond;
      * requires a plausible donor-H...acceptor geometry, i.e. the
        antecedent-donor...acceptor angle (a heavy-atom proxy for the
        D-H...A angle) must exceed ``HBOND_ANGLE_MIN``;
      * excludes covalent contacts and same-residue contacts;
      * excludes salt-bridge partners (handled separately);
      * ignores sequential backbone N/O contacts that are mostly peptide
        geometry rather than a side-chain hydrogen bond;
      * reports the *minimum* qualifying distance per residue instead of
        the first contact found.
    """
    atoms_by_res: Dict[int, Dict[str, Tuple[float, float, float]]] = {}
    backbone_by_res: Dict[int, List[Tuple[str, float, float, float]]] = defaultdict(list)
    for idx, res in enumerate(residues):
        atoms = _atom_map(res)
        atoms_by_res[idx] = atoms
        for raw_name, _elem, x, y, z in res.atoms:
            name = raw_name.strip()
            if name in ("N", "O", "OXT"):
                backbone_by_res[idx].append((name, x, y, z))

    # Enumerate polar contacts once, then classify polarity per residue.
    cutoff2 = HBOND_CUTOFF * HBOND_CUTOFF
    flagged: Set[Tuple[str, int, str]] = set()
    min_dist: Dict[Tuple[str, int, str], float] = {}

    donor_atoms: Dict[int, List[Tuple[str, float, float, float]]] = defaultdict(list)
    acceptor_atoms: Dict[int, List[Tuple[str, float, float, float]]] = defaultdict(list)
    for idx, res in enumerate(residues):
        atoms = atoms_by_res[idx]
        for name in HBOND_DONOR_ATOMS.get(res.aa, ()):
            if name in atoms:
                donor_atoms[idx].append((name, *atoms[name]))
        for name in HBOND_ACCEPTOR_ATOMS.get(res.aa, ()):
            if name in atoms:
                acceptor_atoms[idx].append((name, *atoms[name]))

    for idx, res in enumerate(residues):
        if res.aa not in IONIZABLE or res.key in salt_bridge_keys:
            continue
        atoms = atoms_by_res[idx]
        target_names = ACID_ATOMS.get(res.aa, ()) + BASE_ATOMS.get(res.aa, ())
        targets = [(n, atoms[n]) for n in target_names if n in atoms]
        if not targets:
            continue
        antecedents = HBOND_ANTECEDENT.get(res.aa, {})
        best: Optional[float] = None
        for jdx, other in enumerate(residues):
            if jdx == idx or other.key in salt_bridge_keys:
                continue
            sequential = abs(jdx - idx) == 1 and other.chain == res.chain
            other_atoms = atoms_by_res.get(jdx, {})
            other_donor = donor_atoms.get(jdx, [])
            other_acceptor = acceptor_atoms.get(jdx, [])
            # A backbone N is a donor, a backbone O an acceptor; these are
            # legitimate partners but sequential peptide contacts are not.
            for name2, x2, y2, z2 in backbone_by_res.get(jdx, ()):
                if sequential:
                    continue
                if name2 == "N":
                    other_donor = list(other_donor) + [(name2, x2, y2, z2)]
                elif name2 in ("O", "OXT"):
                    other_acceptor = list(other_acceptor) + [(name2, x2, y2, z2)]
            for name1, (x1, y1, z1) in targets:
                # This residue acts as donor (its atom can donate) and/or as
                # acceptor; the partner must play the complementary role.
                self_is_donor = name1 in HBOND_DONOR_ATOMS.get(res.aa, ())
                self_is_acceptor = name1 in HBOND_ACCEPTOR_ATOMS.get(res.aa, ())
                partners: List[Tuple[str, float, float, float]] = []
                if self_is_donor:
                    partners.extend(other_acceptor)
                if self_is_acceptor:
                    partners.extend(other_donor)
                for name2, x2, y2, z2 in partners:
                    if sequential and name2 in ("N", "O", "OXT"):
                        continue
                    d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    if d2 > cutoff2:
                        continue
                    dist = math.sqrt(d2)
                    if dist < 1.8:       # covalent / overlapping
                        continue
                    # Donor-side geometry: the hydrogen points from the
                    # donor heavy atom towards the acceptor, so the
                    # antecedent-donor...acceptor angle must be wide.
                    if self_is_donor:
                        ant1_name = antecedents.get(name1)
                        if ant1_name and ant1_name in atoms:
                            ang1 = _angle_deg(atoms[ant1_name],
                                              (x1, y1, z1), (x2, y2, z2))
                            if ang1 < HBOND_ANGLE_MIN:
                                continue
                    else:
                        # This residue is the acceptor; the *partner* donates.
                        ant2_name = HBOND_ANTECEDENT.get(other.aa, {}).get(name2)
                        if ant2_name and ant2_name in other_atoms:
                            ang2 = _angle_deg(other_atoms[ant2_name],
                                              (x2, y2, z2), (x1, y1, z1))
                            if ang2 < HBOND_ANGLE_MIN:
                                continue
                    if best is None or dist < best:
                        best = dist
        if best is not None:
            flagged.add(res.key)
            min_dist[res.key] = best
    return flagged, min_dist


def _res_label(res: "tap_core.Residue") -> str:
    """IMGT label like H42, H100A, L27."""
    pos = getattr(res, "imgt_pos", None)
    if pos is None:
        return "{}:{}{}".format(res.chain, res.resseq, res.icode.strip())
    ins = getattr(res, "imgt_ins", "") or ""
    return "{}{}{}".format(res.chain, pos, ins)


def _exposure_class(rel_sasa_frac: Optional[float]) -> str:
    if rel_sasa_frac is None:
        return "n/a"
    pct = rel_sasa_frac * 100.0
    if pct < BURIED_REL_SASA:
        return "buried"
    if pct < SOLVATED_REL_SASA:
        return "intermediate"
    return "exposed"


def _burial_factor(rel_sasa_frac: Optional[float]) -> float:
    """0 (fully solvated) .. 1 (deeply buried) using the 25 % solvation cut."""
    if rel_sasa_frac is None:
        return 0.0
    pct = rel_sasa_frac * 100.0
    if pct >= SOLVATED_REL_SASA:
        return 0.0
    return max(0.0, min(1.0, (SOLVATED_REL_SASA - pct) / SOLVATED_REL_SASA))


def _pka_shift(aa: str, rel_sasa_frac: Optional[float],
               salt_distance: Optional[float],
               hbond_distance: Optional[float],
               ionic_strength: float = IONIC_STRENGTH_DEFAULT) -> Tuple[float, str]:
    """Transparent empirical pKa shift for the structure-aware fallback.

    Acidic groups: burial raises pKa (charged state penalised); salt bridge /
    H-bond lower it (charged state stabilised).  Basic groups mirror this.
    Salt-bridge and H-bond stabilisation is screened by an approximate
    Debye-Hückel factor exp(-kappa*r), so a high-ionic-strength formulation
    does not over-reward long-range contacts.
    """
    burial = _burial_factor(rel_sasa_frac)
    if aa in ("D", "E"):
        burial_scale = 1.5
    elif aa in ("C", "Y"):
        burial_scale = 1.0
    elif aa == "H":
        burial_scale = 1.0
    else:  # K, R
        burial_scale = 1.5

    sign = 1.0 if aa in ACIDIC else -1.0
    delta = sign * burial_scale * burial
    reasons: List[str] = []
    if burial > 0.0:
        reasons.append("buried {:.0f}%".format(burial * 100))
    if salt_distance is not None:
        screen = _screening_factor(salt_distance, ionic_strength)
        delta += (-0.8 if aa in ACIDIC else 0.8) * screen
        reasons.append("salt bridge {:.1f} A (screen {:.2f})".format(
            salt_distance, screen))
    if hbond_distance is not None:
        screen = _screening_factor(hbond_distance, ionic_strength)
        delta += (-0.3 if aa in ACIDIC else 0.3) * screen
        reasons.append("H-bond {:.1f} A (screen {:.2f})".format(
            hbond_distance, screen))
    delta = max(-2.0, min(2.0, delta))
    return delta, ", ".join(reasons)


def structure_groups(
    residues: Sequence["tap_core.Residue"],
    pka_set: str,
    disulfide_keys: Set[Tuple[str, int, str]],
    salt_bridge_keys: Set[Tuple[str, int, str]],
    hbond_keys: Set[Tuple[str, int, str]],
    propka_pka: Optional[Dict[Tuple[str, int, str], float]] = None,
    propka_model_pka: Optional[Dict[Tuple[str, int, str], float]] = None,
    model: str = "heuristic",
    salt_distances: Optional[Dict[Tuple[str, int, str], float]] = None,
    hbond_distances: Optional[Dict[Tuple[str, int, str], float]] = None,
    ionic_strength: float = IONIC_STRENGTH_DEFAULT,
    propka_termini: Optional[Dict[Tuple[str, str], float]] = None,
    propka_termini_model: Optional[Dict[Tuple[str, str], float]] = None,
    chain_types: Optional[Dict[str, str]] = None,
) -> Tuple[List[IonizableGroup], List[Terminus]]:
    """Build per-residue pKa groups for a structure-aware pI.

    When PROPKA values are available they are used *as an environmental
    correction* on top of the selected pKa set:

        pKa_effective = pKa_set(model) + [pKa_PROPKA - pKa_model_PROPKA]

    The bracket is the desolvation / hydrogen-bond / salt-bridge / charge
    contribution that PROPKA computed for this residue; it is
    set-independent, so transferring it onto each bundled pKa set keeps the
    pKa-set sensitivity analysis meaningful instead of collapsing all five
    sets onto a single PROPKA value.  If PROPKA did not report a model pKa
    for a residue, the raw PROPKA value is used as a fallback.
    """
    groups: List[IonizableGroup] = []
    for res in residues:
        if res.aa not in IONIZABLE:
            continue
        rel = res.rel_sasa / 100.0 if getattr(res, "rel_sasa", None) is not None else None
        key = res.key
        disulfide = key in disulfide_keys
        salt = key in salt_bridge_keys
        hbond = key in hbond_keys
        salt_dist = salt_distances.get(key) if salt_distances else (
            SALT_BRIDGE_CUTOFF if salt else None)
        hbond_dist = hbond_distances.get(key) if hbond_distances else (
            HBOND_CUTOFF if hbond else None)
        kind = "acid" if res.aa in ACIDIC else "base"
        pka_model = _pka_value(pka_set, res.aa)
        reason = ""
        if disulfide:
            pka_eff = pka_model
            active = False
            reason = "disulfide (non-ionizable)"
        elif propka_pka is not None and key in propka_pka:
            propka_value = float(propka_pka[key])
            reference = None
            if propka_model_pka:
                reference = propka_model_pka.get(key)
            if reference is not None:
                env_delta = propka_value - float(reference)
                pka_eff = pka_model + env_delta
                reason = "PROPKA environment ({:+.2f} pH)".format(env_delta)
            else:
                pka_eff = propka_value
                reason = "PROPKA"
            active = pka_eff < 20.0
            if not active:
                reason = "PROPKA non-ionizable (pKa >= 20)"
        elif model == "none":
            pka_eff = pka_model
            active = True
            reason = "no microenvironment correction"
        else:
            delta, reason = _pka_shift(
                res.aa, rel, salt_dist, hbond_dist, ionic_strength)
            pka_eff = pka_model + delta
            active = True
        groups.append(IonizableGroup(
            chain=res.chain,
            resseq=res.resseq,
            icode=res.icode,
            aa=res.aa,
            kind=kind,
            pka_model=pka_model,
            pka_eff=pka_eff,
            active=active,
            rel_sasa=rel,
            exposure=_exposure_class(rel),
            salt_bridge=salt,
            hbond=hbond,
            disulfide=disulfide,
            region=_region_imgt(getattr(res, "imgt_pos", None)),
            label=_res_label(res),
            structure_chain=res.chain,
            shift_reason=reason,
        ))
    termini = _structure_termini(
        residues, pka_set, propka_termini,
        propka_termini_model=propka_termini_model,
        chain_types=chain_types)
    return groups, termini


def _structure_termini(
    residues: Sequence["tap_core.Residue"],
    pka_set: str,
    propka_termini: Optional[Dict[Tuple[str, str], float]] = None,
    propka_termini_model: Optional[Dict[Tuple[str, str], float]] = None,
    chain_types: Optional[Dict[str, str]] = None,
) -> List[Terminus]:
    """N/C termini for each Fv chain, using the first/last modelled residue.

    If a chain is truncated (the first/last modelled residue is not the true
    V-domain terminus), no artificial free terminus is created at the cut:
    the terminus is marked blocked so a missing loop does not become an
    unphysical charged group.  PROPKA N+/C- values are used when available.

    As for the side chains, a PROPKA terminal pKa is transferred as an
    environmental correction onto the pKa set's own terminal value, so that
    the terminal pKa also varies with the selected pKa set.

    The C-terminal completeness test uses the chain-type specific FR4-IMGT
    end (VH/VHH -> 128, VL -> 127); see ``_true_c_terminus``.
    """
    propka_termini = propka_termini or {}
    propka_termini_model = propka_termini_model or {}
    chain_types = chain_types or {}
    by_chain: Dict[str, List["tap_core.Residue"]] = defaultdict(list)
    for res in residues:
        by_chain[res.chain].append(res)
    termini: List[Terminus] = []
    for chain, chain_res in by_chain.items():
        if not chain_res:
            continue
        first, last = chain_res[0], chain_res[-1]
        n_pos = getattr(first, "imgt_pos", None)
        c_pos = getattr(last, "imgt_pos", None)
        n_true = n_pos is None or int(n_pos) == TRUE_TERMINUS_N_IMGT
        c_limit = _true_c_terminus(chain_types.get(chain))
        c_true = c_pos is None or int(c_pos) >= c_limit

        def _terminal_pka(kind: str, aa: str) -> Tuple[float, str]:
            """Effective terminal pKa and its provenance."""
            fallback = _pka_value(pka_set, aa, kind)
            value = propka_termini.get((chain, kind))
            if value is None:
                return fallback, "model"
            reference = propka_termini_model.get((chain, kind))
            if reference is None:
                return float(value), "propka"
            env_delta = float(value) - float(reference)
            return fallback + env_delta, "propka environment ({:+.2f} pH)".format(
                env_delta)

        n_pka, n_source = _terminal_pka("N", first.aa)
        c_pka, c_source = _terminal_pka("C", last.aa)
        termini.append(Terminus(
            chain=chain, kind="N", aa=first.aa,
            pka=n_pka,
            blocked=not n_true,
            truncated=not n_true,
            source=("truncated" if not n_true else n_source),
        ))
        termini.append(Terminus(
            chain=chain, kind="C", aa=last.aa,
            pka=c_pka,
            blocked=not c_true,
            truncated=not c_true,
            source=("truncated" if not c_true else c_source),
        ))
    return termini


# ---------------------------------------------------------------------------
# PROPKA integration (optional)
# ---------------------------------------------------------------------------

def find_propka() -> Optional[str]:
    for exe in ("propka3", "propka"):
        path = shutil.which(exe)
        if path:
            return path
    return None


def _parse_propka_rows(
    text: str,
) -> Tuple[Dict[Tuple[str, int, str], float],
           Dict[Tuple[str, int, str], float],
           Dict[Tuple[str, str], float],
           Dict[Tuple[str, str], float]]:
    """Parse the PROPKA summary section.

    Returns ``(pka, model_pka, termini_pka, termini_model_pka)``.

    Two details of the real PROPKA output matter:

    * the chain-ID column may be blank.  The residue number and the chain ID
      are then separated only by whitespace, so a naive ``split()`` reads the
      pKa as the chain ID and the model-pKa as the pKa.  The column is
      therefore disambiguated by type: the chain ID is never a number, the
      pKa always is.
    * PROPKA reports the *model* pKa it started from next to the final pKa.
      The difference between the two is the purely environmental
      contribution (desolvation, hydrogen bonds, salt bridges, charge
      interactions).  Keeping the model-pKa column allows that contribution
      to be transferred onto a different free-solution pKa set, which is what
      makes the pKa-set sensitivity analysis meaningful in PROPKA mode.
    """
    pka: Dict[Tuple[str, int, str], float] = {}
    model_pka: Dict[Tuple[str, int, str], float] = {}
    termini: Dict[Tuple[str, str], float] = {}
    termini_model: Dict[Tuple[str, str], float] = {}
    in_summary = False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        stripped = line.strip()
        if not stripped:
            if pka or termini:
                break
            continue
        if stripped.startswith(("Group", "SUMMARY", "---", "pKa")):
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        group = parts[0]

        def _as_float(token: str) -> Optional[float]:
            try:
                return float(token)
            except (TypeError, ValueError):
                return None

        def _looks_like_pka(token: str) -> bool:
            """PROPKA prints pKa with two decimals, a chain ID never has one."""
            return "." in token and _as_float(token) is not None

        # The chain-ID column is blank for some inputs.  It can be detected
        # unambiguously because a chain ID is never a number while the pKa
        # always is:  "ASP  25 A  5.07  3.80" (chain A) vs
        # "ASP  25    5.07  3.80" (no chain).
        chain_token = parts[2]
        blank_chain = _looks_like_pka(chain_token)
        chain = "" if blank_chain else chain_token
        offset = 2 if blank_chain else 3
        if len(parts) < offset + 2:
            continue
        pka_value = _as_float(parts[offset])
        model_value = _as_float(parts[offset + 1])
        if pka_value is None:
            continue

        if group in ("N+", "N", "C-", "C"):
            kind = "N" if group in ("N+", "N") else "C"
            termini[(chain, kind)] = pka_value
            if model_value is not None:
                termini_model[(chain, kind)] = model_value
            continue
        if group not in ("ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG"):
            continue
        m = re.match(r"^(\d+)([A-Za-z]*)$", parts[1])
        if not m:
            continue
        try:
            resseq = int(m.group(1))
            icode = m.group(2) or " "
        except ValueError:
            continue
        pka[(chain, resseq, icode)] = pka_value
        if model_value is not None:
            model_pka[(chain, resseq, icode)] = model_value
    return pka, model_pka, termini, termini_model


def parse_propka_summary(
    text: str,
) -> Tuple[Dict[Tuple[str, int, str], float], Dict[Tuple[str, str], float]]:
    """Parse PROPKA side-chain and N+/C- terminal pKa values."""
    pka, _model, termini, _termini_model = _parse_propka_rows(text)
    return pka, termini


def parse_propka_model_pka(
    text: str,
) -> Tuple[Dict[Tuple[str, int, str], float], Dict[Tuple[str, str], float]]:
    """Parse the PROPKA *model* pKa values (the free-solution reference)."""
    _pka, model, _termini, termini_model = _parse_propka_rows(text)
    return model, termini_model


def parse_propka_pka(text: str) -> Dict[Tuple[str, int, str], float]:
    """Backward-compatible wrapper returning only side-chain pKa values."""
    return parse_propka_summary(text)[0]


def run_propka_full(
    pdb_path: str,
    timeout: int = 600,
) -> Tuple[Optional[Dict[Tuple[str, int, str], float]],
           Dict[Tuple[str, str], float], Optional[str]]:
    """Run PROPKA and return (side-chain pKa, terminal pKa, error)."""
    pka, _model, termini, _termini_model, err = run_propka_full_detailed(
        pdb_path, timeout=timeout)
    return pka, termini, err


def run_propka_full_detailed(
    pdb_path: str,
    timeout: int = 600,
) -> Tuple[Optional[Dict[Tuple[str, int, str], float]],
           Dict[Tuple[str, int, str], float],
           Dict[Tuple[str, str], float],
           Dict[Tuple[str, str], float],
           Optional[str]]:
    """Run PROPKA and return pKa, model pKa, terminal pKa and terminal model pKa.

    The model pKa is the free-solution reference PROPKA started from; the
    difference between the two columns is the purely environmental
    contribution.  Keeping it lets the pKa-set sensitivity analysis be
    meaningful even when PROPKA supplies every pKa value.
    """
    exe = find_propka()
    if not exe:
        return None, {}, {}, {}, "PROPKA not found on PATH (install propka3 or use --structure-model heuristic)"
    with tempfile.TemporaryDirectory() as tmpdir:
        local = os.path.join(tmpdir, os.path.basename(pdb_path))
        try:
            shutil.copy(pdb_path, local)
        except OSError as exc:
            return None, {}, {}, {}, "cannot copy structure for PROPKA: {}".format(exc)
        try:
            proc = subprocess.run(
                [exe, os.path.basename(local)],
                cwd=tmpdir, capture_output=True, text=True, timeout=timeout,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return None, {}, {}, {}, "PROPKA failed: {}".format(exc)
        if proc.returncode != 0:
            return None, {}, {}, {}, "PROPKA exited {}: {}".format(
                proc.returncode, (proc.stderr or proc.stdout).strip()[:400])
        stem = os.path.splitext(os.path.basename(local))[0]
        pka_path = os.path.join(tmpdir, stem + ".pka")
        text = ""
        if os.path.exists(pka_path):
            try:
                with open(pka_path) as fh:
                    text = fh.read()
            except OSError as exc:
                return None, {}, {}, {}, "cannot read PROPKA output: {}".format(exc)
        else:
            text = proc.stdout
        parsed, model, termini, termini_model = _parse_propka_rows(text)
        if not parsed:
            return None, {}, {}, {}, "PROPKA produced no parsable pKa values"
        return parsed, model, termini, termini_model, None


def run_propka(
    pdb_path: str,
    timeout: int = 600,
) -> Tuple[Optional[Dict[Tuple[str, int, str], float]], Optional[str]]:
    """Backward-compatible PROPKA runner returning only side-chain pKa."""
    pka, _termini, err = run_propka_full(pdb_path, timeout=timeout)
    return pka, err


# ---------------------------------------------------------------------------
# Patch clustering
# ---------------------------------------------------------------------------

@dataclass
class PatchCalibration:
    """Risk bands for the local electrostatic patch score.

    The default bands are a *prior*, not a validated threshold set: the
    local score is an unordered fractional-charge sum and the published TAP
    PPC/PNC bands are ordered integer-charge sums over the whole CDR
    vicinity, so the two scales are not interchangeable.  A site with its own
    experimental data can therefore supply a calibration file that replaces
    the bands (see ``load_patch_calibration``).
    """

    pos: Dict[str, float] = field(
        default_factory=lambda: dict(PATCH_RISK_POS))
    neg: Dict[str, float] = field(
        default_factory=lambda: dict(PATCH_RISK_NEG))
    cdr_factor: float = PATCH_CDR_FACTOR
    source: str = "default"

    def bands_for(self, sign: str) -> Dict[str, float]:
        return self.pos if sign == "positive" else self.neg


def load_patch_calibration(path: Optional[str]) -> PatchCalibration:
    """Load patch risk bands from a JSON calibration file.

    Expected schema (all keys optional)::

        {
          "positive": {"elevated": 0.5, "high": 2.0},
          "negative": {"elevated": 0.6, "high": 2.5},
          "cdr_factor": 1.25
        }

    A missing file is not an error: the documented default bands are used
    and the source is reported as "default".
    """
    if not path:
        return PatchCalibration()
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(
            "Patch calibration file must contain a JSON object, got {}".format(
                type(data).__name__))
    cal = PatchCalibration(source=os.path.abspath(path))
    for key, target in (("positive", cal.pos), ("negative", cal.neg)):
        block = data.get(key)
        if isinstance(block, dict):
            for band in ("elevated", "high"):
                if band in block:
                    target[band] = float(block[band])
    if "cdr_factor" in data:
        cal.cdr_factor = float(data["cdr_factor"])
    for band in ("elevated", "high"):
        for target in (cal.pos, cal.neg):
            if target["high"] < target["elevated"]:
                raise ValueError(
                    "Patch calibration for {}: 'high' ({}) must not be below "
                    "'elevated' ({})".format(
                        "positive" if target is cal.pos else "negative",
                        target["high"], target["elevated"]))
    if cal.cdr_factor <= 0:
        raise ValueError("Patch calibration 'cdr_factor' must be positive")
    return cal


def _patch_risk(score: float, sign: str, in_cdr: bool,
                total_charge: float = 0.0,
                calibration: Optional[PatchCalibration] = None) -> Tuple[str, str]:
    """Classify a local electrostatic patch.

    The local score is an unordered fractional-charge sum, so the published
    TAP PPC/PNC bands do not apply directly.  The bands used here are
    expressed in the same physical units and are intentionally transparent.
    A patch in the CDR vicinity is boosted by the CDR factor because
    CDR-proximal charge clusters are more likely to affect antigen binding
    and polyspecificity.
    """
    cal = calibration or PatchCalibration()
    if sign == "mixed" and abs(total_charge) < 0.5:
        return "low", "balanced mixed patch"
    effective = score * (cal.cdr_factor if in_cdr else 1.0)
    dominant = sign
    if sign == "mixed":
        dominant = "positive" if total_charge > 0 else "negative"
    bands = cal.bands_for(dominant)
    if effective >= bands["high"]:
        risk = "high"
    elif effective >= bands["elevated"]:
        risk = "elevated"
    else:
        risk = "low"
    basis = ("local unordered fractional-charge score {:.3f}; effective "
             "{:.3f}{}".format(
                 score, effective,
                 " (x{:.2f} CDR vicinity)".format(cal.cdr_factor)
                 if in_cdr else ""))
    return risk, basis


def cluster_patches(
    residues: Sequence["tap_core.Residue"],
    pair_dist: Dict[Tuple[int, int], float],
    ph: float,
    groups_by_key: Optional[Dict[Tuple[Optional[str], int, str], IonizableGroup]] = None,
    min_rel_sasa_pct: float = SURFACE_REL_SASA,
    cutoff: float = PATCH_CUTOFF,
    interface_keys: Optional[Set[Tuple[str, int, str]]] = None,
    terminal_keys: Optional[Set[Tuple[str, int, str]]] = None,
    calibration: Optional[PatchCalibration] = None,
) -> List[Patch]:
    """Cluster exposed charged residues into 3D patches at a given pH."""
    interface_keys = interface_keys or set()
    terminal_keys = terminal_keys or set()
    nodes: List[int] = []
    node_charge: Dict[int, float] = {}
    for idx, res in enumerate(residues):
        if res.aa not in IONIZABLE:
            continue
        rel = getattr(res, "rel_sasa", 0.0) or 0.0
        if rel < min_rel_sasa_pct:
            continue
        key = res.key
        if groups_by_key and key in groups_by_key:
            group = groups_by_key[key]
            if not group.active:
                continue
            q = group.charge(ph)
            if group.salt_bridge:
                q = 0.0
        else:
            q = res.charge  # TAP charge at pH 7.4
        if abs(q) < 0.05:
            continue
        nodes.append(idx)
        node_charge[idx] = q

    if not nodes:
        return []

    # Union-find over the exposed charged residues.
    parent = {i: i for i in nodes}

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    node_set = set(nodes)
    for (i, j), dist in pair_dist.items():
        if i in node_set and j in node_set and dist <= cutoff:
            union(i, j)

    clusters: Dict[int, List[int]] = defaultdict(list)
    for idx in nodes:
        clusters[find(idx)].append(idx)

    patches: List[Patch] = []
    for patch_id, members in enumerate(sorted(clusters.values(), key=lambda m: (-len(m), m)), 1):
        member_set = set(members)
        total_q = sum(node_charge[i] for i in members)
        if total_q > 0.2:
            sign = "positive"
        elif total_q < -0.2:
            sign = "negative"
        else:
            sign = "mixed"

        score = 0.0
        for i in members:
            qi = node_charge[i]
            for j in members:
                if j <= i:
                    continue
                qj = node_charge[j]
                if (qi > 0 and qj > 0) or (qi < 0 and qj < 0):
                    d = pair_dist.get((i, j) if i < j else (j, i))
                    if d and d > 1e-6:
                        score += abs(qi * qj) / (d * d)

        # TAP-style ordered-pair sum with formal integer charges, reported
        # for cross-checking against the published PPC/PNC scale.  It is not
        # used for the local risk decision because TAP's thresholds apply to
        # the whole CDR vicinity, not to an isolated patch.
        score_tap = 0.0
        for i in members:
            qi = getattr(residues[i], "charge", 0.0) or 0.0
            for j in members:
                if j == i:
                    continue
                qj = getattr(residues[j], "charge", 0.0) or 0.0
                if not ((qi > 0 and qj > 0) or (qi < 0 and qj < 0)):
                    continue
                d = pair_dist.get((i, j) if i < j else (j, i))
                if d and d > 1e-6:
                    score_tap += abs(qi * qj) / (d * d)

        res_rows: List[Dict[str, object]] = []
        n_cdr = n_fw = n_anchor = 0
        in_vicinity = False
        chains: List[str] = []
        terminal_labels: List[str] = []
        for i in members:
            res = residues[i]
            label = _res_label(res)
            region = _region_imgt(getattr(res, "imgt_pos", None))
            is_terminal = res.key in terminal_keys
            if is_terminal:
                terminal_labels.append(label)
            if res.is_cdr:
                n_cdr += 1
            else:
                n_fw += 1
            if res.is_anchor:
                n_anchor += 1
            if getattr(res, "in_cdr_vicinity", False):
                in_vicinity = True
            if res.chain not in chains:
                chains.append(res.chain)
            res_rows.append({
                "label": label,
                "chain": res.chain,
                "aa": res.aa,
                "imgt_pos": getattr(res, "imgt_pos", None),
                "region": region,
                "rel_sasa": round((getattr(res, "rel_sasa", 0.0) or 0.0), 1),
                "charge": round(node_charge[i], 3),
                "is_cdr": bool(res.is_cdr),
                "is_anchor": bool(res.is_anchor),
                "is_terminus": bool(is_terminal),
                "in_cdr_vicinity": bool(getattr(res, "in_cdr_vicinity", False)),
                "is_interface": res.key in interface_keys,
            })

        risk, risk_basis = _patch_risk(
            score, sign, in_vicinity, total_q, calibration=calibration)
        patch = Patch(
            patch_id=patch_id,
            sign=sign,
            ph=ph,
            residues=sorted(res_rows, key=lambda r: (-abs(float(str(r["charge"]))), str(r["label"]))),
            total_charge=round(total_q, 3),
            score=round(score, 3),
            score_tap=round(score_tap, 3),
            risk_basis=risk_basis,
            terminal_residues=terminal_labels,
            n_cdr=n_cdr,
            n_framework=n_fw,
            n_anchor=n_anchor,
            in_cdr_vicinity=in_vicinity,
            chains=chains,
            risk=risk,
        )
        patches.append(patch)
    return patches


def _sequence_effect(
    chains: Sequence[ChainInfo],
    chain_index: int,
    position: int,
    mutant: str,
    pka_set: str,
    ph: float,
    formulation_ph: float,
    termini_mode: str,
    kind: str,
) -> Dict[str, Optional[float]]:
    """Estimate the sequence-level effect of one substitution."""
    baseline = sequence_pi_table(
        chains, pka_set, ph=ph, formulation_ph=formulation_ph,
        termini_mode=termini_mode, kind=kind)
    mutated: List[ChainInfo] = []
    for idx, chain in enumerate(chains):
        seq = chain.sequence
        if idx == chain_index and 0 <= position < len(seq):
            seq = seq[:position] + mutant + seq[position + 1:]
        mutated.append(ChainInfo(
            name=chain.name, chain_type=chain.chain_type, sequence=seq,
            structure_chain=chain.structure_chain))
    after = sequence_pi_table(
        mutated, pka_set, ph=ph, formulation_ph=formulation_ph,
        termini_mode=termini_mode, kind=kind)
    raw_base_pi = baseline.get("pI_selected")
    raw_new_pi = after.get("pI_selected")
    raw_base_q = baseline.get("net_charge_selected_formulation_ph")
    raw_new_q = after.get("net_charge_selected_formulation_ph")
    base_pi = float(raw_base_pi) if isinstance(raw_base_pi, (int, float)) else None
    new_pi = float(raw_new_pi) if isinstance(raw_new_pi, (int, float)) else None
    base_q = float(raw_base_q) if isinstance(raw_base_q, (int, float)) else None
    new_q = float(raw_new_q) if isinstance(raw_new_q, (int, float)) else None
    return {
        "delta_pI_estimate": (None if base_pi is None or new_pi is None
                              else round(new_pi - base_pi, 3)),
        "delta_charge_formulation_estimate": (
            None if base_q is None or new_q is None
            else round(new_q - base_q, 3)),
    }


def mutation_candidates(
    patch: Patch,
    residues: Sequence["tap_core.Residue"],
    interface_keys: Optional[Set[Tuple[str, int, str]]] = None,
    vhh_hallmark_keys: Optional[Set[Tuple[str, int, str]]] = None,
    terminal_keys: Optional[Set[Tuple[str, int, str]]] = None,
    pka_set: str = "IPC_protein",
    ph: float = 7.4,
    formulation_ph: float = 6.0,
    termini_mode: str = "auto",
    kind: str = "fab",
) -> List[Dict[str, object]]:
    """Suggest conservative framework-local substitutions for a charged patch.

    Only residues that move the patch toward neutrality are proposed.  CDR,
    FR4, VHH-hallmark, chain termini and (optionally) VH/VL-interface
    positions are excluded or flagged, because those changes carry
    affinity/pairing/processing risk.

    Each candidate also carries a sequence-level what-if estimate
    (``delta_pI_estimate`` / ``delta_charge_formulation_estimate``) so the
    user can rank hypotheses before building variants.  The estimate is
    sequence-based and does not model the structural relaxation of the
    mutation.
    """
    interface_keys = interface_keys or set()
    vhh_hallmark_keys = vhh_hallmark_keys or set()
    terminal_keys = terminal_keys or set()
    by_label = {_res_label(r): r for r in residues}
    target_sign = patch.sign
    if target_sign == "mixed":
        if abs(patch.total_charge) < 0.5:
            return []
        target_sign = "positive" if patch.total_charge > 0 else "negative"

    # Chain sequences and residue -> (chain index, sequence position) map.
    chain_order: List[str] = []
    chain_residues: Dict[str, List["tap_core.Residue"]] = defaultdict(list)
    for res in residues:
        if res.chain not in chain_residues:
            chain_order.append(res.chain)
        chain_residues[res.chain].append(res)
    chain_infos: List[ChainInfo] = []
    pos_by_key: Dict[Tuple[str, int, str], Tuple[int, int]] = {}
    for ci, cid in enumerate(chain_order):
        chain_res = chain_residues[cid]
        seq = "".join(r.aa for r in chain_res)
        chain_infos.append(ChainInfo(
            name="VH" if ci == 0 and kind == "fab" else (
                "VL" if kind == "fab" else "VHH"),
            chain_type="H" if ci == 0 else "L",
            sequence=seq, structure_chain=cid))
        for pi, res in enumerate(chain_res):
            pos_by_key[res.key] = (ci, pi)

    out: List[Dict[str, object]] = []
    for row in patch.residues:
        res = by_label.get(str(row["label"]))
        if res is None:
            continue
        aa = res.aa
        if aa not in IONIZABLE:
            continue
        if res.is_cdr:
            continue
        region = str(row.get("region", ""))
        if region == "FR4":
            continue
        if res.is_anchor:
            # IMGT anchors (26/39, 55/66, 104/118) delimit the CDR loops
            # and are part of the conserved paratope framework.  The README
            # forbids mutating them, so they are excluded rather than merely
            # flagged.
            continue
        if res.key in vhh_hallmark_keys:
            continue
        if res.key in terminal_keys:
            # Changing the first/last residue also changes the terminal pKa
            # and can create pyroglutamate/processing liabilities.
            continue
        if aa in ("C", "Y"):
            # Cys is handled as a disulfide; Tyr is an aromatic residue whose
            # mutation is rarely a charge-engineering solution.
            continue
        if target_sign == "negative" and aa in BASIC:
            continue
        if target_sign == "positive" and aa in ACIDIC:
            continue
        if aa in BASIC:
            neutral = list(POSITIVE_TO_NEUTRAL)
            reversal = list(POSITIVE_TO_NEGATIVE)
        else:
            neutral = list(NEGATIVE_TO_NEUTRAL)
            reversal = list(NEGATIVE_TO_POSITIVE)
        is_interface = res.key in interface_keys
        risk = "low"
        if res.is_anchor or row.get("in_cdr_vicinity"):
            risk = "medium"
        if is_interface:
            risk = "high"
        if risk == "high":
            rationale = ("VH/VL interface residue; mutation may perturb chain "
                         "pairing - validate pairing and affinity")
        elif risk == "medium":
            rationale = ("near CDR/anchors; confirm paratope independence "
                         "before mutating")
        else:
            rationale = ("framework-exposed charge; neutral-polar swap reduces "
                         "the patch without introducing an opposite charge")

        effects: Dict[str, Dict[str, Optional[float]]] = {}
        loc = pos_by_key.get(res.key)
        if loc is not None:
            ci, pos = loc
            for mutant in neutral[:3]:
                effects[mutant] = _sequence_effect(
                    chain_infos, ci, pos, mutant, pka_set, ph,
                    formulation_ph, termini_mode, kind)
        first = neutral[0] if neutral else None
        first_effect = effects.get(first, {}) if first else {}
        out.append({
            "label": str(row["label"]),
            "aa": aa,
            "region": region,
            "rel_sasa": row.get("rel_sasa"),
            "charge": row.get("charge"),
            "is_interface": is_interface,
            "is_terminus": False,
            "preferred_substitutions": neutral,
            "charge_reversal_options": reversal,
            "substitution_effects": effects,
            "delta_pI_estimate": first_effect.get("delta_pI_estimate"),
            "delta_charge_formulation_estimate": first_effect.get(
                "delta_charge_formulation_estimate"),
            "risk": risk,
            "rationale": rationale,
        })
    return out


# ---------------------------------------------------------------------------
# Recommendations
# ---------------------------------------------------------------------------

def _fmt(value: Optional[float], digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return ("{:.%df}" % digits).format(value)


def build_recommendations(result: Dict[str, object]) -> List[Dict[str, str]]:
    """Turn pI / TAP / patch results into actionable developability advice."""
    recs: List[Dict[str, str]] = []

    def add(category: str, severity: str, message: str, details: str = "") -> None:
        recs.append({
            "category": category,
            "severity": severity,
            "message": message,
            "details": details,
        })

    fv = result.get("fv", {})  # type: ignore[assignment]
    seq = result.get("sequence_pI", {})  # type: ignore[assignment]
    struct = result.get("structure_pI", {})  # type: ignore[assignment]
    pka_set = result.get("pka_set", "IPC_protein")

    pi_seq = None
    if isinstance(seq, dict):
        entry = seq.get(pka_set, {}) if pka_set in seq else {}
        if isinstance(entry, dict):
            pi_seq = entry.get("pI_selected", entry.get("pI_free"))

    pi_struct = None
    if isinstance(struct, dict):
        entry = struct.get(pka_set, {}) if pka_set in struct else {}
        if isinstance(entry, dict):
            pi_struct = entry.get("pI_selected", entry.get("pI_free"))

    # Prefer the structure-aware value when a structure was analysed.
    pi_selected = pi_struct if pi_struct is not None else pi_seq
    if pi_selected is None and isinstance(fv, dict):
        pi_selected = fv.get("pI_selected", fv.get("pI_free"))

    # --- pI window ---------------------------------------------------------
    if pi_selected is not None:
        pi = float(str(pi_selected))
        if pi < PI_LOW:
            add(
                "pI", "caution",
                "Fv pI is low ({}) - acidic Fv".format(_fmt(pi)),
                "Low-pI Fvs can be difficult to separate from acidic HCPs by "
                "anion exchange and may sit close to standard pH 5-6 "
                "formulation buffers. Consider raising pI toward 7.5-9.0 by "
                "neutralising exposed framework acids (D/E -> Q/N/S/T/A) or "
                "adding basic residues; check the Fc isotype charge match "
                "(IgG1 Fc is positive, IgG2/IgG4 Fc is negative).",
            )
        elif pi < PI_NEUTRAL:
            add(
                "pI", "caution",
                "Fv pI is near-neutral ({})".format(_fmt(pi)),
                "Avoid formulating within ~0.5 pH units of the pI. Confirm "
                "colloidal stability (AC-SINS/BVP/DLS) across pH 5.0-7.5 and "
                "check VH/VL charge asymmetry.",
            )
        elif pi <= PI_FAVORABLE_HIGH:
            add(
                "pI", "info",
                "Fv pI is in the favourable weakly-basic window ({})".format(_fmt(pi)),
                "The 7.5-9.0 Fv pI window balances repulsive self-interactions "
                "in acidic formulation with low non-specific binding at pH 7.4 "
                "(Gupta 2022; Thorsteinson 2021). Maintain this window during "
                "affinity/developability optimisation.",
            )
        else:
            add(
                "pI", "caution",
                "Fv pI is high ({}) - strongly basic Fv".format(_fmt(pi)),
                "High-pI Fvs can show heparin/glycocalyx binding, faster "
                "non-FcRn clearance and polyspecificity. Consider lowering pI "
                "toward 8.0-8.5 by neutralising exposed framework K/R "
                "(K/R -> Q/N/S/T/A); avoid charge reversal unless necessary. "
                "Check PPC and Fv_chml, and include a heparin-binding or "
                "polyreactivity assay.",
            )

    # --- structure vs sequence --------------------------------------------
    # The sequence table reports the pI of the *modelled* chains and the
    # structure table the pI with the microenvironment correction applied.
    # Only the selected-model delta is meaningful; the free/capped variants
    # differ only by the terminus model and are reported separately in the
    # JSON, so they must not produce a second, near-identical message.
    pi_seq_sel = None
    if isinstance(seq, dict):
        entry = seq.get(pka_set, {}) if pka_set in seq else {}
        if isinstance(entry, dict):
            pi_seq_sel = entry.get("pI_selected")
    if pi_struct is not None and pi_seq_sel is not None:
        delta_sel = float(str(pi_struct)) - float(str(pi_seq_sel))
        if abs(delta_sel) >= 0.5:
            add(
                "pI", "info",
                "Structure-aware pI differs from sequence pI by {:+.2f} pH units (selected model)".format(delta_sel),
                "Buried/salt-bridged ionizable groups shift the effective "
                "pI. Use the structure-aware value for formulation and "
                "chromatography planning, but confirm by cIEF. The "
                "fallback microenvironment model is a screening "
                "heuristic; install PROPKA for a stronger prediction.",
            )

    # --- charge asymmetry --------------------------------------------------
    if isinstance(fv, dict):
        # The published thresholds are defined on the *surface* descriptors:
        # SFvCSP is a product of solvent-exposed chain charges, and the
        # Fv_chml window refers to the same pair of surface chain charges.
        # The total-charge counterparts are reported for completeness but
        # must not be compared against surface thresholds, so the surface
        # values are used explicitly here (and are None without a
        # structure).
        chml = fv.get("Fv_chml_surface")
        if chml is None:
            chml = fv.get("Fv_chml_total")
        sfvcsp = fv.get("SFvCSP_surface")
        if chml is not None and abs(float(chml)) > 4.0:
            add(
                "charge_asymmetry", "caution",
                "VH-VL charge asymmetry is high (Fv_chml = {:+.2f})".format(float(chml)),
                "Large VH/VL charge differences are associated with attractive "
                "self-interactions and higher viscosity (Thorsteinson 2021). "
                "Balance the more charged chain using exposed framework "
                "residues; keep |Fv_chml| <= 4 where possible.",
            )
        if sfvcsp is not None and float(sfvcsp) < -6.0:
            add(
                "charge_asymmetry", "caution",
                "SFvCSP is outside the TAP green region ({:.2f})".format(float(sfvcsp)),
                "Oppositely charged VH and VL surfaces (negative product) are a "
                "TAP viscosity/self-association flag. Pair this with the "
                "per-chain charge and patch tables.",
            )

    # --- TAP flags ---------------------------------------------------------
    tap = result.get("tap", {})
    if isinstance(tap, dict):
        flags = tap.get("flags", {})
        if isinstance(flags, dict):
            if flags.get("PPC") in ("amber", "red"):
                add(
                    "charge_patch", "caution",
                    "Positive-charge patch flag (PPC={})".format(flags.get("PPC")),
                    "Reduce exposed positive charge in the CDR vicinity. "
                    "Prefer framework-local K/R -> Q/N/S/T/A substitutions; "
                    "avoid mutating CDR residues that contact antigen.",
                )
            if flags.get("PNC") in ("amber", "red"):
                add(
                    "charge_patch", "caution",
                    "Negative-charge patch flag (PNC={})".format(flags.get("PNC")),
                    "Reduce exposed negative charge in the CDR vicinity. "
                    "Prefer framework-local D/E -> Q/N/S/T/A substitutions; "
                    "check whether the negative patch is needed for antigen "
                    "binding before engineering.",
                )
            if flags.get("PSH") in ("amber", "red"):
                add(
                    "hydrophobic_patch", "caution",
                    "Hydrophobic patch flag (PSH={})".format(flags.get("PSH")),
                    "Reduce exposed CDR-vicinity hydrophobicity. Consider "
                    "S/T/N/Q/Y substitutions at framework positions, and use "
                    "ProteinMPNN/ThermoMPNN to filter designs for stability.",
                )

    # --- specific patches --------------------------------------------------
    patches = result.get("patches", [])
    if isinstance(patches, list):
        for patch in patches:
            if not isinstance(patch, dict):
                continue
            if patch.get("risk") not in ("elevated", "high"):
                continue
            cands = patch.get("mutation_candidates", [])
            labels = []
            if isinstance(cands, list):
                for c in cands[:6]:
                    if isinstance(c, dict):
                        labels.append("{} ({})".format(
                            c.get("label"), "/".join(c.get("preferred_substitutions", [])[:3])))
            if labels:
                add(
                    "charge_patch", "caution",
                    "{} {} patch (score {}) at pH {}".format(
                        str(patch.get("risk")).upper(),
                        patch.get("sign"), patch.get("score"), patch.get("ph")),
                    "Candidate framework positions: {}. Treat these as design "
                    "hypotheses: check humanness (BioPhi), stability "
                    "(ThermoMPNN), PTM motifs, and paratope contacts before "
                    "building variants.".format("; ".join(labels)),
                )

    # --- formulation -------------------------------------------------------
    if pi_selected is not None:
        pi = float(str(pi_selected))
        form_ph = result.get("formulation_ph", 6.0)
        try:
            dist = abs(float(str(form_ph)) - pi)
        except (TypeError, ValueError):
            dist = None
        if dist is not None and dist < 0.5:
            add(
                "formulation", "high",
                "Formulation pH {} is within 0.5 units of the Fv pI".format(form_ph),
                "Net charge is near zero at the pI, reducing electrostatic "
                "repulsion and increasing colloidal aggregation risk. Move the "
                "buffer pH at least ~1 unit away from the pI, or re-engineer "
                "the Fv charge.",
            )
        elif dist is not None and dist < 1.0:
            add(
                "formulation", "caution",
                "Formulation pH {} is within 1 unit of the Fv pI".format(form_ph),
                "Check concentration-dependent aggregation/viscosity at the "
                "intended formulation pH; a 0.5-1.0 unit pH shift may improve "
                "colloidal stability.",
            )
        else:
            add(
                "formulation", "info",
                "Formulation pH {} is {:.1f} units from the Fv pI".format(form_ph, dist),
                "The Fv carries a net charge at the intended formulation pH. "
                "Confirm the sign and magnitude against the Fc charge and the "
                "planned IEX step.",
            )

    # --- applicability scope ------------------------------------------------
    # The thresholds and reference distributions used here come from work on
    # mammalian IgG / Fab.  They are applied to a VHH or to a non-mammalian
    # construct only by analogy, and the analytical scope of the pI number
    # itself (Fv vs whole IgG) differs from the experimental readout.
    inp = result.get("input", {})
    kind = inp.get("kind") if isinstance(inp, dict) else None
    if result.get("analysis_type") == "sequence_only":
        add(
            "scope", "caution",
            "Structure-dependent developability metrics were NOT computed",
            "TAP flags (PSH/PPC/PNC/SFvCSP), the 3D charge patches and the "
            "structure-aware pI all require a three-dimensional structure; "
            "only Ltot (CDR length) and the pI / charge budget are "
            "sequence-determinable. An absent flag in this run is NOT a clean "
            "flag: an antibody with a real hydrophobic-patch or "
            "charge-asymmetry flag is indistinguishable from a clean one "
            "here. Re-run the same input through a predicted model "
            "(ESMFold/AF2) or an experimental structure to obtain them.",
        )
    if kind == "vhh":
        add(
            "scope", "caution",
            "VHH: the reference thresholds are Fab/IgG-derived and are not "
            "validated for camelid VHHs",
            "The TAP thresholds (Ltot, PSH, PPC, PNC, SFvCSP), the 7.5-9.0 "
            "pI window and the |Fv_chml| <= 4 limit were all derived on "
            "mammalian IgG/Fab. A VHH has a longer CDR3, a different solvent "
            "exposure profile and no VL partner, so the absolute flag values "
            "must be read as relative guidance only. Use the VHH's own "
            "historical data for calibration (--patch-calibration) before "
            "acting on a flag.",
        )
    if pi_selected is not None:
        add(
            "scope", "info",
            "This is the Fv (variable domain) pI, not the whole IgG pI",
            "cIEF and IEX measure the intact molecule, whose charge also "
            "contains the constant regions: an IgG1 Fc contributes roughly "
            "+7 at pH 7.4 (IgG2/IgG4 are close to neutral to slightly "
            "negative). The Fv pI reported here is therefore a variable "
            "domain quantity and is not directly comparable with an "
            "experimental IgG pI; use it for the relative charge engineering "
            "of the variable domains and confirm the assembled molecule "
            "experimentally.",
        )

    # --- validation / tooling ---------------------------------------------
    add(
        "validation", "info",
        "Experimental confirmation is required",
        "Use cIEF for pI, IEX for charge variants, HIC/AC-SINS/BVP for "
        "hydrophobicity and self-interaction, SEC for aggregation, DSF for Tm, "
        "and viscosity at the intended concentration. Computational pI is a "
        "screening estimate (typically within ~0.5-1 pH unit).",
    )
    add(
        "tooling", "info",
        "Recommended follow-up tools",
        "TAP: python3 tools/tap_profiler/tap_analyzer.py --pdb <pdb> "
        "--all-formats; PTM exposure: python3 tools/ptm_exposure/analyze_ptm.py "
        "--pdb <pdb> --all-formats; ThermoMPNN: python3 "
        "tools/stability_predictor/stability_analyzer.py --pdb <pdb> --chain all "
        "--positions <labels> --thermompnn-dir <repo>; BioPhi humanness: "
        "python3 tools/biophi/biophi_analyze.py sapiens --input <fasta> "
        "--outdir <dir>; ProteinMPNN: python3 "
        "tools/mpnn/proteinmpnn_adapter.py --pdb_path <pdb> --out_folder <dir>.",
    )
    return recs


# ---------------------------------------------------------------------------
# Top-level analysis
# ---------------------------------------------------------------------------

def analyze_sequence(
    chains: Sequence[ChainInfo],
    pka_set: str = "IPC_protein",
    all_pka_sets: bool = False,
    ph: float = 7.4,
    formulation_ph: float = 6.0,
    disulfide_mode: str = "auto",
    termini_mode: str = "auto",
    kind: Optional[str] = None,
    ph_grid: Optional[Sequence[float]] = None,
    name: str = "",
    input_path: str = "",
    patch_calibration: Optional[PatchCalibration] = None,
) -> Dict[str, object]:
    """Sequence-only pI analysis (no structure)."""
    pka_set = normalize_pka_set(pka_set)
    if kind is None:
        kind = "fab" if len(chains) >= 2 else "vhh"
    resolved_termini = resolve_termini_mode(termini_mode, kind)
    capped = resolved_termini in ("capped", "blocked")
    n_blocked = resolved_termini == "blocked"
    warnings: List[str] = []
    # Clean the sequences once and keep the warnings.  Cleaning twice (once to
    # build the chains, once to collect warnings) loses the warnings, because
    # the second pass sees the already-cleaned sequence.
    #
    # Non-standard residues are dropped or substituted, which shifts the
    # residue indices of the excluded (disulfide) Cys; the index mapping is
    # therefore tracked explicitly.
    cleaned_chains: List[ChainInfo] = []
    for c in chains:
        clean, warns = clean_sequence(c.sequence)
        warnings.extend("{}: {}".format(c.name, w) for w in warns)
        if not clean:
            raise ValueError(
                "Chain {} has no usable sequence after cleaning".format(c.name))
        # Map indices of the original sequence onto the cleaned sequence.
        # clean_sequence first strips every non-letter, then keeps exactly
        # one output character per retained residue, in order.  The retained
        # original indices can therefore be listed directly, and the
        # excluded-Cys indices (which refer to the original sequence) are
        # remapped onto cleaned-sequence coordinates.
        letters = [i for i, ch in enumerate(c.sequence) if ch.isalpha()]
        kept = [i for i in letters
                if c.sequence[i].upper() in ("A", "C", "D", "E", "F", "G",
                                            "H", "I", "K", "L", "M", "N",
                                            "P", "Q", "R", "S", "T", "V",
                                            "W", "Y", "U", "O")]
        original_to_cleaned = {orig: new for new, orig in enumerate(kept)}
        excluded = {original_to_cleaned[i] for i in c.excluded_cys
                    if i in original_to_cleaned}
        cleaned_chains.append(ChainInfo(
            name=c.name, chain_type=c.chain_type, sequence=clean,
            structure_chain=c.structure_chain, excluded_cys=excluded))
    chains = cleaned_chains

    if disulfide_mode == "exclude":
        for c in chains:
            c.excluded_cys = {i for i, aa in enumerate(c.sequence) if aa == "C"}
    elif disulfide_mode == "auto":
        found, warns = find_disulfide_cys_sequence(chains)
        warnings.extend(warns)
        for c in chains:
            c.excluded_cys = {i for (n, i) in found if n == c.name}
        if not found:
            warnings.append(
                "No conserved disulfide Cys could be assigned; free Cys "
                "treated as ionizable (use --disulfide exclude to force "
                "disulfide-only Cys).")
    elif disulfide_mode == "free":
        pass
    else:
        raise ValueError("--disulfide must be auto, exclude, or free")

    # Compute every pKa set once: the selected set is reported and the full
    # set provides the pKa-model uncertainty range.
    all_tables: Dict[str, object] = {}
    for s in PKA_SETS:
        all_tables[s] = sequence_pi_table(
            chains, s, ph=ph, formulation_ph=formulation_ph,
            termini_mode=termini_mode, kind=kind)
    seq_table: Dict[str, object] = (
        all_tables if all_pka_sets else {pka_set: all_tables[pka_set]})

    selected_pis: List[float] = []
    for entry in all_tables.values():
        if isinstance(entry, dict):
            value = entry.get("pI_selected")
            if isinstance(value, (int, float)):
                selected_pis.append(float(value))
    pI_range = ([round(min(selected_pis), 3), round(max(selected_pis), 3)]
                if selected_pis else None)

    chain_table = chain_pi_table(
        chains, pka_set, ph=ph, formulation_ph=formulation_ph,
        termini_mode=resolved_termini)
    groups, termini = _make_groups_from_sequence(chains, pka_set)
    pi_free = isoelectric_point(groups, termini)
    pi_capped = isoelectric_point(groups, termini, capped_cterm=True)
    pi_selected = isoelectric_point(
        groups, termini, capped_cterm=capped, blocked_nterm=n_blocked)
    q74 = net_charge(groups, termini, ph, capped_cterm=capped,
                     blocked_nterm=n_blocked)
    qform = net_charge(groups, termini, formulation_ph, capped_cterm=capped,
                       blocked_nterm=n_blocked)
    q74_free = net_charge(groups, termini, ph)
    q74_capped = net_charge(groups, termini, ph, capped_cterm=True)

    fv: Dict[str, object] = {
        "pI_free": pi_free,
        "pI_capped": pi_capped,
        "pI_selected": pi_selected,
        "pI_range": pI_range,
        "termini_mode": resolved_termini,
        "net_charge_ph": q74,
        "net_charge_ph74": q74,
        "net_charge_ph74_free": q74_free,
        "net_charge_formulation_ph": qform,
        "net_charge_ph74_capped": q74_capped,
        "Fv_chml": None,
        "Fv_chml_total": None,
        "Fv_chml_surface": None,
        "SFvCSP": None,
        "SFvCSP_total": None,
        "SFvCSP_surface": None,
        "SFvCSP_basis": None,
        "surface_charge_vh": None,
        "surface_charge_vl": None,
        "total_charge_vh": None,
        "total_charge_vl": None,
    }
    if len(chains) >= 2:
        q_by_chain = {}
        for c in chains:
            g, t = _make_groups_from_sequence([c], pka_set)
            q_by_chain[c.name] = net_charge(
                g, t, ph, capped_cterm=capped, blocked_nterm=n_blocked)
        if "VH" in q_by_chain and "VL" in q_by_chain:
            chml = q_by_chain["VH"] - q_by_chain["VL"]
            # In sequence mode there is no solvent-exposure information, so
            # the *surface* descriptors (the quantities the TAP thresholds
            # are defined on) cannot be computed.  The documented alias
            # therefore falls back to the total-charge product, and
            # "SFvCSP_basis" records which quantity the alias holds so that
            # the TAP thresholds are never applied to the wrong one.
            fv["Fv_chml"] = chml
            fv["Fv_chml_total"] = chml
            fv["Fv_chml_surface"] = None
            fv["SFvCSP"] = q_by_chain["VH"] * q_by_chain["VL"]
            fv["SFvCSP_total"] = fv["SFvCSP"]
            fv["SFvCSP_surface"] = None
            fv["SFvCSP_basis"] = "total"
            fv["total_charge_vh"] = q_by_chain["VH"]
            fv["total_charge_vl"] = q_by_chain["VL"]

    curve = charge_curve(
        groups, termini, ph_grid or default_ph_grid(),
        capped_cterm=capped, blocked_nterm=n_blocked)
    sensitivity = {
        "pka_sets": {
            s: (entry.get("pI_selected") if isinstance(entry, dict) else None)
            for s, entry in all_tables.items()},
        "pI_selected_range": pI_range,
        "termini_mode": resolved_termini,
    }

    result: Dict[str, object] = {
        "name": name or "sequence",
        "input": {
            "type": "sequence",
            "path": input_path,
            "kind": kind,
            "chains": [{"name": c.name, "chain_type": c.chain_type,
                        "length": len(c.sequence)} for c in chains],
        },
        "pka_set": pka_set,
        "pka_sets_available": list(PKA_SETS),
        "ph": ph,
        "formulation_ph": formulation_ph,
        "disulfide_mode": disulfide_mode,
        "termini_mode": resolved_termini,
        "analysis_type": "sequence_only",
        "unavailable_metrics": unavailable_metrics(),
        "sequence_pI": seq_table,
        "structure_pI": {},
        "sensitivity": sensitivity,
        "charge_curve": curve,
        "chains": chain_table,
        "fv": fv,
        "tap": {},
        "patches": [],
        "ionizable_residues": [],
        "structure_quality": None,
        "recommendations": [],
        "warnings": warnings,
    }
    result["recommendations"] = build_recommendations(result)
    return result


def _write_fv_pdb(
    residues: Sequence["tap_core.Residue"],
    path: str,
) -> Dict[str, str]:
    """Write the analysed Fv chains to a PDB file, preserving numbering.

    The original chain ID, residue number and insertion code are kept so
    that residue keys stay identical to the ones used everywhere else in
    this module, and so that PROPKA (which reports pKa values keyed by
    chain/residue number) can be matched back without a renumbering step.

    PDB format only allows a single-character chain ID, so longer IDs
    (possible in mmCIF) are remapped to A, B, C, ...  The returned mapping
    translates the written label back to the original chain ID.

    Returns the mapping ``written label -> original chain ID``.
    """
    chain_ids: List[str] = []
    for res in residues:
        if res.chain not in chain_ids:
            chain_ids.append(res.chain)
    labels: Dict[str, str] = {}
    used: Set[str] = set()
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
    for cid in chain_ids:
        label = cid if len(cid) == 1 and cid not in used else ""
        if not label:
            for cand in alphabet:
                if cand not in used:
                    label = cand
                    break
            if not label:
                label = "X"
        used.add(label)
        labels[label] = cid

    serial = 1
    with open(path, "w") as fh:
        for res in residues:
            label = None
            for cand, orig in labels.items():
                if orig == res.chain:
                    label = cand
                    break
            icode = (res.icode or " ")[:1]
            for raw_name, element, x, y, z in res.atoms:
                name = raw_name if len(raw_name) == 4 else " %-3s" % raw_name.strip()
                fh.write(
                    "ATOM  {:5d} {:4s} {:>3s} {:1s}{:4d}{:1s}   "
                    "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format(
                        serial, name, res.resname, label, res.resseq, icode,
                        x, y, z, 1.00, 0.00, element))
                serial += 1
        fh.write("END\n")
    return labels


def _run_propka_on_fv(
    residues: Sequence["tap_core.Residue"],
    timeout: int = 600,
) -> Tuple[Optional[Dict[Tuple[str, int, str], float]],
           Dict[Tuple[str, int, str], float],
           Dict[Tuple[str, str], float],
           Dict[Tuple[str, str], float],
           Optional[str]]:
    """Run PROPKA on the trimmed Fv and translate its chain labels back.

    PROPKA is executed on a PDB file that contains exactly the analysed Fv
    residues, with their original numbering preserved.  This keeps the
    PROPKA environment identical to the SASA-based environment used by the
    heuristic model (no constant regions, no antigen) and keeps the residue
    keys directly comparable.
    """
    exe = find_propka()
    if not exe:
        return None, {}, {}, {}, (
            "PROPKA not found on PATH (install propka3 or use "
            "--structure-model heuristic)")
    with tempfile.TemporaryDirectory() as tmpdir:
        local = os.path.join(tmpdir, "fv.pdb")
        try:
            labels = _write_fv_pdb(residues, local)
        except OSError as exc:
            return None, {}, {}, {}, "cannot write Fv for PROPKA: {}".format(exc)
        try:
            proc = subprocess.run(
                [exe, os.path.basename(local)],
                cwd=tmpdir, capture_output=True, text=True, timeout=timeout,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return None, {}, {}, {}, "PROPKA failed: {}".format(exc)
        if proc.returncode != 0:
            return None, {}, {}, {}, "PROPKA exited {}: {}".format(
                proc.returncode, (proc.stderr or proc.stdout).strip()[:400])
        pka_path = os.path.join(tmpdir, "fv.pka")
        text = ""
        if os.path.exists(pka_path):
            try:
                with open(pka_path) as fh:
                    text = fh.read()
            except OSError as exc:
                return None, {}, {}, {}, "cannot read PROPKA output: {}".format(exc)
        else:
            text = proc.stdout
        pka, model, termini, termini_model = _parse_propka_rows(text)
        if not pka:
            return None, {}, {}, {}, "PROPKA produced no parsable pKa values"

    # Translate the written chain labels back to the original chain IDs.
    def _remap(key: Tuple[str, int, str]) -> Tuple[str, int, str]:
        chain, resseq, icode = key
        return (labels.get(chain, chain), resseq, icode)

    def _remap_term(key: Tuple[str, str]) -> Tuple[str, str]:
        chain, kind = key
        return (labels.get(chain, chain), kind)

    return (
        {_remap(k): v for k, v in pka.items()},
        {_remap(k): v for k, v in model.items()},
        {_remap_term(k): v for k, v in termini.items()},
        {_remap_term(k): v for k, v in termini_model.items()},
        None,
    )


def _select_chains_structure(
    residues: Sequence["tap_core.Residue"],
    vh_chain: Optional[str] = None,
    vl_chain: Optional[str] = None,
) -> Tuple[Optional[str], Optional[str], List[str], str]:
    """Return (vh_chain, vl_chain, warnings, kind) for a structure.

    kind is "fab" or "vhh".  Falls back to the TAP selector for two-chain Fvs.
    """
    warnings: List[str] = []
    chain_ids = tap_core._chain_order(list(residues))
    if vh_chain and vh_chain not in chain_ids:
        raise ValueError("VH chain {!r} not found (chains: {})".format(
            vh_chain, ", ".join(chain_ids)))
    if vl_chain and vl_chain not in chain_ids:
        raise ValueError("VL chain {!r} not found (chains: {})".format(
            vl_chain, ", ".join(chain_ids)))

    if vh_chain and vl_chain:
        if vh_chain == vl_chain:
            raise ValueError(
                "--vh-chain and --vl-chain must be different chains "
                "(both were {!r}); an Fv needs one heavy and one light "
                "chain.".format(vh_chain))
        return vh_chain, vl_chain, warnings, "fab"
    if vh_chain and not vl_chain:
        return vh_chain, None, warnings, "vhh"

    # Single chain -> VHH (or a truncated heavy chain).
    if len(chain_ids) == 1:
        cid = chain_ids[0]
        seq = "".join(r.aa for r in residues if r.chain == cid)
        ctype = tap_core._anarci_chain_type(seq)
        if ctype in ("K", "L"):
            raise ValueError(
                "Only one light chain found; an Fv requires VH+VL. "
                "Pass --vh-chain/--vl-chain or provide a paired structure.")
        return cid, None, warnings, "vhh"

    # Two or more chains: use the validated TAP selector.
    try:
        vh, vl, w = tap_core.select_fv_chains(list(residues), vh_chain, vl_chain)
        warnings.extend(w)
        return vh, vl, warnings, "fab"
    except Exception as exc:
        # Last-resort: pick the longest H-like and L-like chains.
        warnings.append("TAP chain selection failed: {}".format(exc))
        types: Dict[str, Optional[str]] = {}
        for cid in chain_ids:
            seq = "".join(r.aa for r in residues if r.chain == cid)
            types[cid] = tap_core._anarci_chain_type(seq)
        h = [c for c in chain_ids if types.get(c) == "H"]
        l = [c for c in chain_ids if types.get(c) in ("K", "L")]
        if h and l:
            vh = max(h, key=lambda c: len([r for r in residues if r.chain == c]))
            vl = max(l, key=lambda c: len([r for r in residues if r.chain == c]))
            return vh, vl, warnings, "fab"
        # Single H-like chain among several (e.g. antigen) -> VHH.
        if h and not l:
            vh = max(h, key=lambda c: len([r for r in residues if r.chain == c]))
            return vh, None, warnings, "vhh"
        raise ValueError(
            "Cannot identify Fv chains automatically (chains: {}). "
            "Pass --vh-chain/--vl-chain.".format(", ".join(chain_ids)))


def analyze_structure(
    pdb_path: str,
    pka_set: str = "IPC_protein",
    all_pka_sets: bool = False,
    ph: float = 7.4,
    formulation_ph: float = 6.0,
    vh_chain: Optional[str] = None,
    vl_chain: Optional[str] = None,
    sasa_engine: str = "auto",
    psa_binary: Optional[str] = None,
    disulfide_mode: str = "auto",
    structure_model: str = "propka",
    tap_thresholds: str = "2024",
    termini_mode: str = "auto",
    ionic_strength: float = IONIC_STRENGTH_DEFAULT,
    ph_grid: Optional[Sequence[float]] = None,
    name: str = "",
    patch_calibration: Optional[PatchCalibration] = None,
) -> Dict[str, object]:
    """Structure-aware pI + TAP + charge-patch analysis for an Fv/VHH."""
    pka_set = normalize_pka_set(pka_set)
    if structure_model not in ("heuristic", "propka", "none"):
        raise ValueError("--structure-model must be heuristic, propka, or none")
    if disulfide_mode not in ("auto", "exclude", "free"):
        raise ValueError("--disulfide must be auto, exclude, or free")
    if ionic_strength < 0.0:
        raise ValueError("--ionic-strength must be >= 0")

    warnings: List[str] = []
    all_residues = tap_core.load_residues(pdb_path)
    vh, vl, chain_warnings, kind = _select_chains_structure(
        all_residues, vh_chain, vl_chain)
    warnings.extend(chain_warnings)
    resolved_termini = resolve_termini_mode(termini_mode, kind)
    capped = resolved_termini in ("capped", "blocked")
    n_blocked = resolved_termini == "blocked"

    # Number and trim the Fv chains to the V domain.
    fvi: List["tap_core.Residue"] = []
    chain_order: List[str] = []
    chain_types: Dict[str, str] = {}
    for cid, ctype in ((vh, "H"), (vl, "L")):
        if cid is None:
            continue
        chain_res = [r for r in all_residues if r.chain == cid]
        try:
            numbered = tap_core.number_chain_imgt(chain_res, ctype)
        except Exception as exc:
            raise RuntimeError(
                "IMGT numbering failed for chain {}: {}. ANARCI is required "
                "for structure mode.".format(cid, exc))
        fvi.extend(numbered)
        chain_order.append(cid)
        chain_types[cid] = ctype

    if not fvi:
        raise RuntimeError("No Fv residues could be numbered")
    if len(set(chain_order)) != len(chain_order):
        raise RuntimeError(
            "The same chain was selected twice ({}); an Fv needs one heavy "
            "and one light chain.".format(chain_order))

    # Surface exposure and TAP annotations.
    rel_sasa, engine, sasa_warnings = tap_core.compute_rel_sasa(
        fvi, engine=sasa_engine, psa_binary=psa_binary)
    warnings.extend(sasa_warnings)
    for idx, res in enumerate(fvi):
        res.rel_sasa = rel_sasa.get(idx, 0.0)
    counts = tap_core.annotate_residues(fvi)

    # Two distance maps are needed:
    #  * the TAP map (7.5 A) reproduces the published neighbour definition
    #    used for PSH/PPC/PNC and must not change;
    #  * the patch map uses the documented 8 A patch-graph cutoff.
    pair_dist = tap_core.build_pair_distances(fvi)
    patch_pair_dist = (pair_dist if tap_core.NEIGHBOUR_CUTOFF >= PATCH_CUTOFF
                       else tap_core.build_pair_distances(
                           fvi, cutoff=PATCH_CUTOFF))
    vhh_hallmark_keys: Set[Tuple[str, int, str]] = set()
    if kind == "vhh" and vh is not None:
        for r in fvi:
            if r.chain == vh and getattr(r, "imgt_pos", None) in VHH_HALLMARK_IMGT:
                vhh_hallmark_keys.add(r.key)
    interface_keys: Set[Tuple[str, int, str]] = set()
    for (i, j), dist in pair_dist.items():
        if dist <= INTERFACE_CUTOFF and fvi[i].chain != fvi[j].chain:
            interface_keys.add(fvi[i].key)
            interface_keys.add(fvi[j].key)
    terminal_keys: Set[Tuple[str, int, str]] = set()
    for cid in chain_order:
        chain_res = [r for r in fvi if r.chain == cid]
        if chain_res:
            terminal_keys.add(chain_res[0].key)
            terminal_keys.add(chain_res[-1].key)
    psh, ppc, pnc = tap_core.compute_surface_scores(fvi, pair_dist)
    ltot = sum(1 for r in fvi if r.is_cdr)

    sfvcsp = None
    qh = ql = None
    if vh is not None and vl is not None:
        sfvcsp, qh, ql = tap_core.compute_sfvcsp(fvi, vh, vl)

    values: Dict[str, Optional[float]] = {
        "Ltot": float(ltot),
        "PSH": psh,
        "PPC": ppc,
        "PNC": pnc,
        "SFvCSP": sfvcsp,
    }
    flags: Dict[str, str] = {}
    spec = tap_core.THRESHOLD_SETS[tap_thresholds]
    for metric in tap_core.METRIC_ORDER:
        value = values.get(metric)
        flags[metric] = ("unknown" if value is None
                         else tap_core.flag_value(value, spec[metric]))
    available_flags = {k: v for k, v in flags.items() if v != "unknown"}
    tap_summary = (tap_core.risk_summary(available_flags)
                   if available_flags else "UNKNOWN")

    # Disulfides.
    if disulfide_mode == "exclude":
        disulfide_keys = {r.key for r in fvi if r.aa == "C"}
    elif disulfide_mode == "free":
        disulfide_keys = set()
    else:
        disulfide_keys = find_disulfide_cys_structure(fvi)
        if not disulfide_keys:
            # Structures without SG atoms (e.g. backbone-only models) fall
            # back to the conserved-Cys numbering rule.
            seq_chains = []
            for cid, ctype in ((vh, "H"), (vl, "L")):
                if cid is None:
                    continue
                seq = "".join(r.aa for r in fvi if r.chain == cid)
                seq_chains.append(ChainInfo(
                    name=cid, chain_type=ctype, sequence=seq,
                    structure_chain=cid))
            found, warns = find_disulfide_cys_sequence(seq_chains)
            warnings.extend(warns)
            by_chain = {c.name: c for c in seq_chains}
            for cname, idx in found:
                c = by_chain.get(cname)
                if c is None:
                    continue
                chain_res = [r for r in fvi if r.chain == cname]
                if idx < len(chain_res):
                    disulfide_keys.add(chain_res[idx].key)
            if not disulfide_keys:
                warnings.append(
                    "No disulfide bond detected; Cys treated as ionizable "
                    "(use --disulfide exclude to force exclusion).")

    # Microenvironment interactions.  The pKa salt-bridge definition is
    # intentionally broader than TAP's exposed 3.2 A rule; both are reported.
    salt_keys, salt_pairs, salt_dists = detect_salt_bridges(fvi)
    hbond_keys, hbond_dists = detect_hbonds(fvi, salt_keys)

    # Optional PROPKA (side-chain and terminal pKa).
    #
    # PROPKA is run on the *trimmed Fv* rather than on the input file, so
    # that the PROPKA environment and the heuristic environment describe
    # the same molecular system.  Running it on the full input would let
    # constant regions or a bound antigen shift the Fv pKa values while the
    # SASA-based (heuristic) path still sees only the Fv, making the two
    # structure models non-comparable and the structure-vs-sequence delta
    # meaningless.
    propka_pka = None
    propka_model_pka: Dict[Tuple[str, int, str], float] = {}
    propka_termini: Dict[Tuple[str, str], float] = {}
    propka_termini_model: Dict[Tuple[str, str], float] = {}
    if structure_model == "propka":
        propka_pka, propka_model_pka, propka_termini, propka_termini_model, err = \
            _run_propka_on_fv(fvi)
        if propka_pka is None:
            warnings.append(
                "PROPKA is unavailable ({}), falling back to the heuristic "
                "pKa model. Install propka3 for PROPKA pKa values.".format(err))
            structure_model = "heuristic"
            propka_pka = None
            propka_model_pka = {}
            propka_termini = {}
            propka_termini_model = {}
    elif structure_model == "heuristic":
        if find_propka():
            warnings.append(
                "PROPKA is installed but --structure-model heuristic was "
                "selected; use --structure-model propka for PROPKA pKa values.")

    # Structure quality (B-factor semantics, missing side chains, breaks).
    provenance = detect_structure_provenance(pdb_path)
    quality = assess_structure_quality(
        fvi, kind, chain_types=chain_types, provenance=provenance)
    warnings.extend(structure_quality_warnings(quality))

    # Sequence pI table (using the modelled V-domain sequences).
    chains_info: List[ChainInfo] = []
    for cid, ctype, label in ((vh, "H", "VH"), (vl, "L", "VL")):
        if cid is None:
            continue
        seq = "".join(r.aa for r in fvi if r.chain == cid)
        chains_info.append(ChainInfo(
            name=label, chain_type=ctype, sequence=seq,
            structure_chain=cid))
    if len(chains_info) == 1:
        chains_info[0].name = "VHH"

    propka_arg = propka_pka if structure_model == "propka" else None
    propka_model_arg = (propka_model_pka if structure_model == "propka"
                        else None)
    propka_term_arg = propka_termini if structure_model == "propka" else None
    propka_term_model_arg = (propka_termini_model if structure_model == "propka"
                             else None)

    # Compute every pKa set: the selected set is reported and the full set
    # provides the pKa-model uncertainty range.
    all_struct: Dict[str, Dict[str, object]] = {}
    all_seq: Dict[str, Dict[str, object]] = {}
    selected_pis: List[float] = []
    for s in PKA_SETS:
        groups_s, termini_s = structure_groups(
            fvi, s, disulfide_keys, salt_keys, hbond_keys,
            propka_pka=propka_arg,
            propka_model_pka=propka_model_arg,
            model=structure_model,
            salt_distances=salt_dists,
            hbond_distances=hbond_dists,
            ionic_strength=ionic_strength,
            propka_termini=propka_term_arg,
            propka_termini_model=propka_term_model_arg,
            chain_types=chain_types,
        )
        pi_free = isoelectric_point(groups_s, termini_s)
        pi_capped = isoelectric_point(groups_s, termini_s, capped_cterm=True)
        pi_blocked = isoelectric_point(
            groups_s, termini_s, capped_cterm=True, blocked_nterm=True)
        pi_selected = isoelectric_point(
            groups_s, termini_s, capped_cterm=capped, blocked_nterm=n_blocked)
        q74 = net_charge(groups_s, termini_s, ph, capped_cterm=capped,
                         blocked_nterm=n_blocked)
        qform = net_charge(groups_s, termini_s, formulation_ph,
                           capped_cterm=capped, blocked_nterm=n_blocked)
        if pi_selected is not None:
            selected_pis.append(float(pi_selected))
        seq_entry = sequence_pi_table(
            chains_info, s, ph=ph, formulation_ph=formulation_ph,
            termini_mode=termini_mode, kind=kind)
        all_seq[s] = seq_entry
        delta = None
        delta_free = None
        if isinstance(seq_entry, dict):
            seq_sel = seq_entry.get("pI_selected")
            seq_free = seq_entry.get("pI_free")
            if isinstance(seq_sel, (int, float)) and pi_selected is not None:
                delta = pi_selected - float(seq_sel)
            if isinstance(seq_free, (int, float)) and pi_free is not None:
                delta_free = pi_free - float(seq_free)
        all_struct[s] = {
            "pI_free": pi_free,
            "pI_capped": pi_capped,
            "pI_blocked": pi_blocked,
            "pI_selected": pi_selected,
            "net_charge_ph": q74,
            "net_charge_ph74": q74,
            "net_charge_formulation_ph": qform,
            "delta_vs_sequence": delta,
            "delta_vs_sequence_free": delta_free,
            "model": structure_model,
            "n_groups": len(groups_s),
        }

    seq_table: Dict[str, Dict[str, object]] = (
        all_seq if all_pka_sets else {pka_set: all_seq[pka_set]})
    struct_table: Dict[str, Dict[str, object]] = (
        all_struct if all_pka_sets else {pka_set: all_struct[pka_set]})
    pI_range = ([round(min(selected_pis), 3), round(max(selected_pis), 3)]
                if selected_pis else None)

    # Per-chain table and charge descriptors at the selected pKa set.
    groups, termini = structure_groups(
        fvi, pka_set, disulfide_keys, salt_keys, hbond_keys,
        propka_pka=propka_arg,
        propka_model_pka=propka_model_arg,
        model=structure_model,
        salt_distances=salt_dists,
        hbond_distances=hbond_dists,
        ionic_strength=ionic_strength,
        propka_termini=propka_term_arg,
        propka_termini_model=propka_term_model_arg,
        chain_types=chain_types,
    )
    groups_by_key = {(g.structure_chain, g.resseq, g.icode): g for g in groups}
    chain_table: Dict[str, Dict[str, object]] = {}
    total_q_by_chain: Dict[str, float] = {}
    for cid, ctype, label in ((vh, "H", "VH"), (vl, "L", "VL")):
        if cid is None:
            continue
        chain_groups = [g for g in groups if g.structure_chain == cid]
        chain_termini = [t for t in termini if t.chain == cid]
        seq = "".join(r.aa for r in fvi if r.chain == cid)
        q_sel = net_charge(chain_groups, chain_termini, ph,
                           capped_cterm=capped, blocked_nterm=n_blocked)
        total_q_by_chain[label] = q_sel
        chain_table[label] = {
            "chain_type": ctype,
            "structure_chain": cid,
            "length": len(seq),
            "sequence": seq,
            "pI_free": isoelectric_point(chain_groups, chain_termini),
            "pI_capped": isoelectric_point(chain_groups, chain_termini,
                                           capped_cterm=True),
            "pI_selected": isoelectric_point(
                chain_groups, chain_termini, capped_cterm=capped,
                blocked_nterm=n_blocked),
            "net_charge_ph74": q_sel,
            "net_charge_ph74_free": net_charge(chain_groups, chain_termini, ph),
            "net_charge_formulation_ph": net_charge(
                chain_groups, chain_termini, formulation_ph,
                capped_cterm=capped, blocked_nterm=n_blocked),
            "n_ionizable": len(chain_groups),
        }

    pi_free = float(str(struct_table[pka_set]["pI_free"]))
    pi_capped = float(str(struct_table[pka_set]["pI_capped"]))
    pi_selected = float(str(struct_table[pka_set]["pI_selected"]))
    q74 = float(str(struct_table[pka_set]["net_charge_ph74"]))
    qform = float(str(struct_table[pka_set]["net_charge_formulation_ph"]))

    total_vh = total_q_by_chain.get("VH", total_q_by_chain.get("VHH"))
    total_vl = total_q_by_chain.get("VL")
    fv_chml_total = (total_vh - total_vl
                     if total_vh is not None and total_vl is not None else None)
    fv_chml_surface = (qh - ql
                       if qh is not None and ql is not None else None)
    sfvcsp_total = (total_vh * total_vl
                    if total_vh is not None and total_vl is not None else None)

    # Surface charges use TAP's exposed-residue definition at pH 7.4.
    fv: Dict[str, object] = {
        "pI_free": pi_free,
        "pI_capped": pi_capped,
        "pI_selected": pi_selected,
        "pI_range": pI_range,
        "termini_mode": resolved_termini,
        "net_charge_ph": q74,
        "net_charge_ph74": q74,
        "net_charge_formulation_ph": qform,
        "Fv_chml": (fv_chml_surface if fv_chml_surface is not None
                    else fv_chml_total),
        "Fv_chml_total": fv_chml_total,
        "Fv_chml_surface": fv_chml_surface,
        "SFvCSP": sfvcsp,
        "SFvCSP_surface": sfvcsp,
        "SFvCSP_total": sfvcsp_total,
        "SFvCSP_basis": "surface",
        "surface_charge_vh": qh,
        "surface_charge_vl": ql,
        "total_charge_vh": total_vh,
        "total_charge_vl": total_vl,
    }

    patch_calibration = patch_calibration or PatchCalibration()

    # Patch analysis at formulation and physiological pH.
    #
    # The patch graph uses the documented 8 A cutoff (PATCH_CUTOFF), which
    # is wider than the 7.5 A neighbour map used for the TAP metrics.
    patches: List[Dict[str, object]] = []
    patch_counter = 0
    for patch_ph in sorted({float(formulation_ph), float(ph)}):
        for patch in cluster_patches(
                fvi, patch_pair_dist, patch_ph, groups_by_key=groups_by_key,
                interface_keys=interface_keys, terminal_keys=terminal_keys,
                calibration=patch_calibration):
            patch_counter += 1
            row = {
                "patch_id": patch_counter,
                "sign": patch.sign,
                "ph": patch.ph,
                "total_charge": patch.total_charge,
                "score": patch.score,
                "score_tap": patch.score_tap,
                "risk": patch.risk,
                "risk_basis": patch.risk_basis,
                "n_cdr": patch.n_cdr,
                "n_framework": patch.n_framework,
                "n_anchor": patch.n_anchor,
                "in_cdr_vicinity": patch.in_cdr_vicinity,
                "terminal_residues": patch.terminal_residues,
                "chains": patch.chains,
                "residues": patch.residues,
                "mutation_candidates": mutation_candidates(
                    patch, fvi, interface_keys=interface_keys,
                    vhh_hallmark_keys=vhh_hallmark_keys,
                    terminal_keys=terminal_keys,
                    pka_set=pka_set, ph=ph,
                    formulation_ph=formulation_ph,
                    termini_mode=resolved_termini, kind=kind),
            }
            patches.append(row)

    ionizable_rows = []
    for g in sorted(groups, key=lambda x: (x.structure_chain or "", x.resseq, x.icode)):
        key = (None if g.structure_chain is None
               else (g.structure_chain, g.resseq, g.icode))
        ionizable_rows.append({
            "label": g.label,
            "chain": g.structure_chain,
            "aa": g.aa,
            "region": g.region,
            "kind": g.kind,
            "pka_model": round(g.pka_model, 2),
            "pka_eff": round(g.pka_eff, 2),
            "delta_pka": round(g.pka_eff - g.pka_model, 2),
            "active": g.active,
            "rel_sasa": None if g.rel_sasa is None else round(g.rel_sasa * 100.0, 1),
            "exposure": g.exposure,
            "salt_bridge": g.salt_bridge,
            "hbond": g.hbond,
            "disulfide": g.disulfide,
            "salt_bridge_distance": (None if key is None or key not in salt_dists
                                     else round(salt_dists[key], 2)),
            "hbond_distance": (None if key is None or key not in hbond_dists
                               else round(hbond_dists[key], 2)),
            "charge_ph74": round(g.charge(ph), 3),
            "charge_formulation_ph": round(g.charge(formulation_ph), 3),
            "reason": g.shift_reason,
        })

    curve = charge_curve(
        groups, termini, ph_grid or default_ph_grid(),
        capped_cterm=capped, blocked_nterm=n_blocked)
    sensitivity = {
        "pka_sets": {
            s: (entry.get("pI_selected") if isinstance(entry, dict) else None)
            for s, entry in all_struct.items()},
        "pI_selected_range": pI_range,
        "termini_mode": resolved_termini,
        "structure_model": structure_model,
        "sasa_engine": engine,
        "ionic_strength": ionic_strength,
    }

    if pI_range and (pI_range[1] - pI_range[0]) >= 0.5:
        warnings.append(
            "pI is sensitive to the pKa model (selected pI range {:.2f}-{:.2f} "
            "across bundled sets); report a range rather than a single value."
            .format(pI_range[0], pI_range[1]))
    if pi_capped is not None and pi_free is not None and abs(pi_capped - pi_free) >= 0.5:
        warnings.append(
            "Terminus model changes Fv pI by {:+.2f} units (free vs capped). "
            "The selected '{}' value is used for recommendations; confirm the "
            "construct's C-terminal fusion state.".format(
                pi_capped - pi_free, resolved_termini))

    result: Dict[str, object] = {
        "name": name or os.path.splitext(os.path.basename(pdb_path))[0],
        "input": {
            "type": "structure",
            "path": pdb_path,
            "kind": kind,
            "vh_chain": vh,
            "vl_chain": vl,
            "chains": chain_order,
        },
        "pka_set": pka_set,
        "pka_sets_available": list(PKA_SETS),
        "ph": ph,
        "formulation_ph": formulation_ph,
        "disulfide_mode": disulfide_mode,
        "structure_model": structure_model,
        "sasa_engine": engine,
        "termini_mode": resolved_termini,
        "ionic_strength": ionic_strength,
        "analysis_type": "structure",
        "unavailable_metrics": [],
        "sequence_pI": seq_table,
        "structure_pI": struct_table,
        "sensitivity": sensitivity,
        "charge_curve": curve,
        "chains": chain_table,
        "fv": fv,
        "tap": {
            "metrics": values,
            "flags": flags,
            "risk_summary": tap_summary,
            "thresholds": tap_thresholds,
            "n_surface": counts["n_surface"],
            "n_cdr_vicinity": counts["n_cdr_vicinity"],
            "n_cdr_residues": counts["n_cdr"],
            "n_salt_bridges": counts["n_salt_bridges"],
            "n_salt_bridges_tap": counts["n_salt_bridges"],
            "n_salt_bridges_pka": len(salt_pairs),
            "n_salt_bridge_residues_pka": len(salt_keys),
            "salt_bridges": salt_pairs,
            "salt_bridges_pka": salt_pairs,
            "salt_bridge_definitions": {
                "tap": "3.2 A, both partners exposed, one partner per residue",
                "pka_microenvironment": "4.0 A, any burial, multiple partners",
            },
            "n_residues": counts["n_residues"],
        },
        "patches": patches,
        "patch_calibration": {
            "source": patch_calibration.source,
            "positive": dict(patch_calibration.pos),
            "negative": dict(patch_calibration.neg),
            "cdr_factor": patch_calibration.cdr_factor,
        },
        "ionizable_residues": ionizable_rows,
        "structure_quality": quality,
        "recommendations": [],
        "warnings": warnings,
    }
    result["recommendations"] = build_recommendations(result)
    return result
