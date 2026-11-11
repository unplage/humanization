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
SALT_BRIDGE_CUTOFF = 4.0    # A, donor N ... acceptor O
HBOND_CUTOFF = 3.5          # A, polar N/O ... N/O
DISULFIDE_CUTOFF = 2.5      # A, SG ... SG
PATCH_CUTOFF = 8.0          # A, exposed charged-residue patch graph

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
SALT_DONORS = {"K": ("NZ",), "R": ("NE", "NH1", "NH2")}
SALT_ACCEPTORS = {"D": ("OD1", "OD2"), "E": ("OE1", "OE2")}

# Conservative neutral-polar substitutions for charge-patch engineering.
POSITIVE_TO_NEUTRAL = ("Q", "N", "S", "T", "A")
NEGATIVE_TO_NEUTRAL = ("Q", "N", "S", "T", "A")
POSITIVE_TO_NEGATIVE = ("E", "D")
NEGATIVE_TO_POSITIVE = ("K", "R")


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
    score: float = 0.0
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


def _region_imgt(pos: Optional[int]) -> str:
    if pos is None:
        return ""
    for name, lo, hi in IMGT_REGIONS:
        if lo <= pos <= hi:
            return name
    return "FR4"


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


def sequence_pi_table(chains: Sequence[ChainInfo], pka_set: str) -> Dict[str, object]:
    """pI for the combined chains under free/capped/blocked termini."""
    groups, termini = _make_groups_from_sequence(chains, pka_set)
    return {
        "pI_free": isoelectric_point(groups, termini),
        "pI_capped": isoelectric_point(groups, termini, capped_cterm=True),
        "pI_blocked": isoelectric_point(groups, termini, capped_cterm=True,
                                        blocked_nterm=True),
        "net_charge_ph74_free": net_charge(groups, termini, 7.4),
        "net_charge_ph74_capped": net_charge(groups, termini, 7.4,
                                             capped_cterm=True),
    }


def chain_pi_table(
    chains: Sequence[ChainInfo],
    pka_set: str,
    ph: float = 7.4,
    formulation_ph: float = 6.0,
) -> Dict[str, Dict[str, object]]:
    """Per-chain pI and charge."""
    out: Dict[str, Dict[str, object]] = {}
    for chain in chains:
        groups, termini = _make_groups_from_sequence([chain], pka_set)
        out[chain.name] = {
            "chain_type": chain.chain_type,
            "length": len(chain.sequence),
            "sequence": chain.sequence,
            "pI_free": isoelectric_point(groups, termini),
            "pI_capped": isoelectric_point(groups, termini, capped_cterm=True),
            "net_charge_ph74": net_charge(groups, termini, ph),
            "net_charge_formulation_ph": net_charge(groups, termini, formulation_ph),
            "n_ionizable": len(groups),
        }
    return out


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
        for res in numbered.residues:
            if res.aa.upper() != "C":
                continue
            m = re.match(r"^[A-Z](\d+)([A-Za-z]*)$", res.pos)
            if not m:
                continue
            pos = int(m.group(1))
            if pos in CONSERVED_CYS_KABAT.get(chain.chain_type, set()):
                out.add((chain.name, res.index))
    return out, warnings


# ---------------------------------------------------------------------------
# Structure microenvironment: salt bridges, H-bonds, pKa shifts
# ---------------------------------------------------------------------------

def _atom_map(res: "tap_core.Residue") -> Dict[str, Tuple[float, float, float]]:
    return {raw.strip(): (x, y, z) for raw, _elem, x, y, z in res.atoms}


def detect_salt_bridges(residues: Sequence["tap_core.Residue"]) -> Tuple[Set[Tuple[str, int, str]], List[Dict[str, object]]]:
    """Detect K/R ... D/E salt bridges (donor N ... acceptor O <= 4 A)."""
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
    seen: Set[Tuple[int, int]] = set()
    for d_idx, d_name, x1, y1, z1 in donors:
        for a_idx, a_name, x2, y2, z2 in acceptors:
            if d_idx == a_idx:
                continue
            pair = (d_idx, a_idx) if d_idx < a_idx else (a_idx, d_idx)
            if pair in seen:
                continue
            d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            if d2 <= cutoff2:
                seen.add(pair)
                r1, r2 = residues[d_idx], residues[a_idx]
                partners.add(r1.key)
                partners.add(r2.key)
                pairs.append({
                    "donor": _res_label(r1),
                    "acceptor": _res_label(r2),
                    "distance": round(math.sqrt(d2), 2),
                })
    return partners, pairs


def detect_hbonds(residues: Sequence["tap_core.Residue"],
                  salt_bridge_keys: Set[Tuple[str, int, str]]) -> Set[Tuple[str, int, str]]:
    """Flag ionizable side chains with a polar contact (N/O ... N/O <= 3.5 A).

    Backbone and side-chain N/O atoms are considered.  Salt-bridge partners
    are excluded here because they are already handled separately.
    """
    polar: List[Tuple[int, str, float, float, float]] = []
    for idx, res in enumerate(residues):
        for raw_name, elem, x, y, z in res.atoms:
            name = raw_name.strip()
            if elem in ("N", "O") or name.startswith(("N", "O")):
                polar.append((idx, name, x, y, z))
    if not polar:
        return set()

    cutoff2 = HBOND_CUTOFF * HBOND_CUTOFF
    flagged: Set[Tuple[str, int, str]] = set()
    # Pre-group by residue for a cheap neighbour scan.
    by_res: Dict[int, List[Tuple[str, float, float, float]]] = defaultdict(list)
    for idx, name, x, y, z in polar:
        by_res[idx].append((name, x, y, z))

    for idx, res in enumerate(residues):
        if res.aa not in IONIZABLE or res.key in salt_bridge_keys:
            continue
        atoms = _atom_map(res)
        target_names = ACID_ATOMS.get(res.aa, ()) + BASE_ATOMS.get(res.aa, ())
        targets = [(n, atoms[n]) for n in target_names if n in atoms]
        if not targets:
            continue
        found = False
        for jdx, other in enumerate(residues):
            if jdx == idx:
                continue
            if other.key in salt_bridge_keys:
                continue
            for name2, x2, y2, z2 in by_res.get(jdx, ()):
                for _n1, (x1, y1, z1) in targets:
                    d2 = (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    if d2 <= cutoff2:
                        found = True
                        break
                if found:
                    break
            if found:
                break
        if found:
            flagged.add(res.key)
    return flagged


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
               salt_bridge: bool, hbond: bool) -> Tuple[float, str]:
    """Transparent empirical pKa shift for the structure-aware fallback.

    Acidic groups: burial raises pKa (charged state penalised); salt bridge /
    H-bond lower it (charged state stabilised).  Basic groups mirror this.
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
    if salt_bridge:
        delta += -0.8 if aa in ACIDIC else 0.8
        reasons.append("salt bridge")
    if hbond:
        delta += -0.3 if aa in ACIDIC else 0.3
        reasons.append("H-bond")
    delta = max(-2.0, min(2.0, delta))
    return delta, ", ".join(reasons)


def structure_groups(
    residues: Sequence["tap_core.Residue"],
    pka_set: str,
    disulfide_keys: Set[Tuple[str, int, str]],
    salt_bridge_keys: Set[Tuple[str, int, str]],
    hbond_keys: Set[Tuple[str, int, str]],
    propka_pka: Optional[Dict[Tuple[str, int, str], float]] = None,
    model: str = "heuristic",
) -> Tuple[List[IonizableGroup], List[Terminus]]:
    """Build per-residue pKa groups for a structure-aware pI."""
    groups: List[IonizableGroup] = []
    for res in residues:
        if res.aa not in IONIZABLE:
            continue
        rel = res.rel_sasa / 100.0 if getattr(res, "rel_sasa", None) is not None else None
        key = res.key
        disulfide = key in disulfide_keys
        salt = key in salt_bridge_keys
        hbond = key in hbond_keys
        kind = "acid" if res.aa in ACIDIC else "base"
        pka_model = _pka_value(pka_set, res.aa)
        reason = ""
        if disulfide:
            pka_eff = pka_model
            active = False
            reason = "disulfide (non-ionizable)"
        elif propka_pka is not None and key in propka_pka:
            pka_eff = float(propka_pka[key])
            active = True
            reason = "PROPKA"
        elif model == "none":
            pka_eff = pka_model
            active = True
            reason = "no microenvironment correction"
        else:
            delta, reason = _pka_shift(res.aa, rel, salt, hbond)
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
    termini = _structure_termini(residues, pka_set)
    return groups, termini


def _structure_termini(residues: Sequence["tap_core.Residue"],
                       pka_set: str) -> List[Terminus]:
    """N/C termini for each Fv chain, using the first/last modelled residue."""
    by_chain: Dict[str, List["tap_core.Residue"]] = defaultdict(list)
    for res in residues:
        by_chain[res.chain].append(res)
    termini: List[Terminus] = []
    for chain, chain_res in by_chain.items():
        if not chain_res:
            continue
        first, last = chain_res[0], chain_res[-1]
        termini.append(Terminus(
            chain=chain, kind="N", aa=first.aa,
            pka=_pka_value(pka_set, first.aa, "N"),
        ))
        termini.append(Terminus(
            chain=chain, kind="C", aa=last.aa,
            pka=_pka_value(pka_set, last.aa, "C"),
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


def parse_propka_pka(text: str) -> Dict[Tuple[str, int, str], float]:
    """Parse the SUMMARY OF THIS PREDICTION section of a PROPKA ``.pka`` file."""
    out: Dict[Tuple[str, int, str], float] = {}
    in_summary = False
    for line in text.splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        stripped = line.strip()
        if not stripped:
            if out:
                break
            continue
        if stripped.startswith(("Group", "SUMMARY", "---", "pKa")):
            continue
        parts = stripped.split()
        if len(parts) < 4 or parts[0] not in (
                "ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG"):
            continue
        m = re.match(r"^(\d+)([A-Za-z]*)$", parts[1])
        if not m:
            continue
        try:
            resseq = int(m.group(1))
            icode = m.group(2) or " "
            chain = parts[2]
            pka = float(parts[3])
        except (ValueError, IndexError):
            continue
        out[(chain, resseq, icode)] = pka
    return out


def run_propka(pdb_path: str, timeout: int = 600) -> Tuple[Optional[Dict[Tuple[str, int, str], float]], Optional[str]]:
    """Run PROPKA on a structure and return {(chain, resseq, icode): pKa}."""
    exe = find_propka()
    if not exe:
        return None, "PROPKA not found on PATH (install propka3 or use --structure-model heuristic)"
    with tempfile.TemporaryDirectory() as tmpdir:
        local = os.path.join(tmpdir, os.path.basename(pdb_path))
        try:
            shutil.copy(pdb_path, local)
        except OSError as exc:
            return None, "cannot copy structure for PROPKA: {}".format(exc)
        try:
            proc = subprocess.run(
                [exe, os.path.basename(local)],
                cwd=tmpdir, capture_output=True, text=True, timeout=timeout,
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            return None, "PROPKA failed: {}".format(exc)
        if proc.returncode != 0:
            return None, "PROPKA exited {}: {}".format(
                proc.returncode, (proc.stderr or proc.stdout).strip()[:400])
        stem = os.path.splitext(os.path.basename(local))[0]
        pka_path = os.path.join(tmpdir, stem + ".pka")
        text = ""
        if os.path.exists(pka_path):
            try:
                with open(pka_path) as fh:
                    text = fh.read()
            except OSError as exc:
                return None, "cannot read PROPKA output: {}".format(exc)
        else:
            text = proc.stdout
        parsed = parse_propka_pka(text)
        if not parsed:
            return None, "PROPKA produced no parsable pKa values"
        return parsed, None


# ---------------------------------------------------------------------------
# Patch clustering
# ---------------------------------------------------------------------------

def _patch_risk(score: float, sign: str, in_cdr: bool) -> str:
    # TAP's published PPC/PNC green/red boundaries are a useful scale even
    # though a local patch score is not identical to the CDR-vicinity TAP sum.
    if sign == "positive":
        if score >= 4.22:
            return "high"
        if score >= 1.32:
            return "elevated"
    else:
        if score >= 4.42:
            return "high"
        if score >= 2.00:
            return "elevated"
    return "low"


def cluster_patches(
    residues: Sequence["tap_core.Residue"],
    pair_dist: Dict[Tuple[int, int], float],
    ph: float,
    groups_by_key: Optional[Dict[Tuple[Optional[str], int, str], IonizableGroup]] = None,
    min_rel_sasa_pct: float = SURFACE_REL_SASA,
    cutoff: float = PATCH_CUTOFF,
    interface_keys: Optional[Set[Tuple[str, int, str]]] = None,
) -> List[Patch]:
    """Cluster exposed charged residues into 3D patches at a given pH."""
    interface_keys = interface_keys or set()
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

        res_rows: List[Dict[str, object]] = []
        n_cdr = n_fw = n_anchor = 0
        in_vicinity = False
        chains: List[str] = []
        for i in members:
            res = residues[i]
            label = _res_label(res)
            region = _region_imgt(getattr(res, "imgt_pos", None))
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
                "in_cdr_vicinity": bool(getattr(res, "in_cdr_vicinity", False)),
                "is_interface": res.key in interface_keys,
            })

        patch = Patch(
            patch_id=patch_id,
            sign=sign,
            ph=ph,
            residues=sorted(res_rows, key=lambda r: (-abs(float(str(r["charge"]))), str(r["label"]))),
            total_charge=round(total_q, 3),
            score=round(score, 3),
            n_cdr=n_cdr,
            n_framework=n_fw,
            n_anchor=n_anchor,
            in_cdr_vicinity=in_vicinity,
            chains=chains,
            risk=_patch_risk(score, sign, n_cdr > 0),
        )
        patches.append(patch)
    return patches


def mutation_candidates(
    patch: Patch,
    residues: Sequence["tap_core.Residue"],
    interface_keys: Optional[Set[Tuple[str, int, str]]] = None,
    vhh_hallmark_keys: Optional[Set[Tuple[str, int, str]]] = None,
) -> List[Dict[str, object]]:
    """Suggest conservative framework-local substitutions for a charged patch.

    Only residues that move the patch toward neutrality are proposed.  CDR,
    FR4, VHH-hallmark and (optionally) VH/VL-interface positions are excluded
    or flagged, because those changes carry affinity/pairing risk.
    """
    interface_keys = interface_keys or set()
    vhh_hallmark_keys = vhh_hallmark_keys or set()
    by_label = {_res_label(r): r for r in residues}
    target_sign = patch.sign
    if target_sign == "mixed":
        if abs(patch.total_charge) < 0.5:
            return []
        target_sign = "positive" if patch.total_charge > 0 else "negative"

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
        if res.key in vhh_hallmark_keys:
            continue
        if aa in ("C",):
            # Cys is not a charge-patch residue in practice; keep safe.
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
        out.append({
            "label": str(row["label"]),
            "aa": aa,
            "region": region,
            "rel_sasa": row.get("rel_sasa"),
            "charge": row.get("charge"),
            "is_interface": is_interface,
            "preferred_substitutions": neutral,
            "charge_reversal_options": reversal,
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
            pi_seq = entry.get("pI_free")

    pi_struct = None
    if isinstance(struct, dict):
        entry = struct.get(pka_set, {}) if pka_set in struct else {}
        if isinstance(entry, dict):
            pi_struct = entry.get("pI_free")

    # Prefer the structure-aware value when a structure was analysed.
    pi_free = pi_struct if pi_struct is not None else pi_seq
    if pi_free is None and isinstance(fv, dict):
        pi_free = fv.get("pI_free")

    # --- pI window ---------------------------------------------------------
    if pi_free is not None:
        pi = float(str(pi_free))
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
    if pi_struct is not None and pi_seq is not None:
        delta = float(str(pi_struct)) - float(str(pi_seq))
        if abs(delta) >= 0.5:
            add(
                "pI", "info",
                "Structure-aware pI differs from sequence pI by {:+.2f} pH units".format(delta),
                "Buried/salt-bridged ionizable groups shift the effective "
                "pI. Use the structure-aware value for formulation and "
                "chromatography planning, but confirm by cIEF. The "
                "fallback microenvironment model is a screening "
                "heuristic; install PROPKA for a stronger prediction.",
            )

    # --- charge asymmetry --------------------------------------------------
    if isinstance(fv, dict):
        chml = fv.get("Fv_chml")
        sfvcsp = fv.get("SFvCSP")
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
    if pi_free is not None:
        pi = float(str(pi_free))
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
    name: str = "",
    input_path: str = "",
) -> Dict[str, object]:
    """Sequence-only pI analysis (no structure)."""
    pka_set = normalize_pka_set(pka_set)
    warnings: List[str] = []
    chains = [
        ChainInfo(name=c.name, chain_type=c.chain_type,
                  sequence=clean_sequence(c.sequence)[0],
                  structure_chain=c.structure_chain,
                  excluded_cys=set(c.excluded_cys))
        for c in chains
    ]
    for c in chains:
        _seq, warns = clean_sequence(c.sequence)
        warnings.extend("{}: {}".format(c.name, w) for w in warns)

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

    sets = list(PKA_SETS) if all_pka_sets else [pka_set]
    seq_table: Dict[str, object] = {}
    for s in sets:
        seq_table[s] = sequence_pi_table(chains, s)

    chain_table = chain_pi_table(chains, pka_set, ph=ph, formulation_ph=formulation_ph)
    groups, termini = _make_groups_from_sequence(chains, pka_set)
    pi_free = isoelectric_point(groups, termini)
    pi_capped = isoelectric_point(groups, termini, capped_cterm=True)
    q74 = net_charge(groups, termini, ph)
    qform = net_charge(groups, termini, formulation_ph)
    q74_capped = net_charge(groups, termini, ph, capped_cterm=True)

    fv: Dict[str, object] = {
        "pI_free": pi_free,
        "pI_capped": pi_capped,
        "net_charge_ph": q74,
        "net_charge_ph74": q74,
        "net_charge_formulation_ph": qform,
        "net_charge_ph74_capped": q74_capped,
        "Fv_chml": None,
        "SFvCSP": None,
        "surface_charge_vh": None,
        "surface_charge_vl": None,
    }
    if len(chains) >= 2:
        q_by_chain = {}
        for c in chains:
            g, t = _make_groups_from_sequence([c], pka_set)
            q_by_chain[c.name] = net_charge(g, t, ph)
        if "VH" in q_by_chain and "VL" in q_by_chain:
            fv["Fv_chml"] = q_by_chain["VH"] - q_by_chain["VL"]
            fv["SFvCSP"] = q_by_chain["VH"] * q_by_chain["VL"]
            fv["surface_charge_vh"] = q_by_chain["VH"]
            fv["surface_charge_vl"] = q_by_chain["VL"]

    result: Dict[str, object] = {
        "name": name or "sequence",
        "input": {
            "type": "sequence",
            "path": input_path,
            "chains": [{"name": c.name, "chain_type": c.chain_type,
                        "length": len(c.sequence)} for c in chains],
        },
        "pka_set": pka_set,
        "pka_sets_available": list(PKA_SETS),
        "ph": ph,
        "formulation_ph": formulation_ph,
        "disulfide_mode": disulfide_mode,
        "sequence_pI": seq_table,
        "structure_pI": {},
        "chains": chain_table,
        "fv": fv,
        "tap": {},
        "patches": [],
        "ionizable_residues": [],
        "recommendations": [],
        "warnings": warnings,
    }
    result["recommendations"] = build_recommendations(result)
    return result


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
    structure_model: str = "heuristic",
    tap_thresholds: str = "2024",
    name: str = "",
) -> Dict[str, object]:
    """Structure-aware pI + TAP + charge-patch analysis for an Fv/VHH."""
    pka_set = normalize_pka_set(pka_set)
    if structure_model not in ("heuristic", "propka", "none"):
        raise ValueError("--structure-model must be heuristic, propka, or none")
    if disulfide_mode not in ("auto", "exclude", "free"):
        raise ValueError("--disulfide must be auto, exclude, or free")

    warnings: List[str] = []
    all_residues = tap_core.load_residues(pdb_path)
    vh, vl, chain_warnings, kind = _select_chains_structure(
        all_residues, vh_chain, vl_chain)
    warnings.extend(chain_warnings)

    # Number and trim the Fv chains to the V domain.
    fvi: List["tap_core.Residue"] = []
    chain_order: List[str] = []
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

    if not fvi:
        raise RuntimeError("No Fv residues could be numbered")

    # Surface exposure and TAP annotations.
    rel_sasa, engine, sasa_warnings = tap_core.compute_rel_sasa(
        fvi, engine=sasa_engine, psa_binary=psa_binary)
    warnings.extend(sasa_warnings)
    for idx, res in enumerate(fvi):
        res.rel_sasa = rel_sasa.get(idx, 0.0)
    counts = tap_core.annotate_residues(fvi)

    pair_dist = tap_core.build_pair_distances(fvi)
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

    # Microenvironment interactions.
    salt_keys, salt_pairs = detect_salt_bridges(fvi)
    hbond_keys = detect_hbonds(fvi, salt_keys)

    # Optional PROPKA.
    propka_pka = None
    if structure_model == "propka":
        propka_pka, err = run_propka(pdb_path)
        if propka_pka is None:
            raise RuntimeError(
                "PROPKA requested but unavailable: {}. Install propka3 or use "
                "--structure-model heuristic.".format(err))
    elif structure_model == "heuristic":
        if find_propka():
            warnings.append(
                "PROPKA is installed but --structure-model heuristic was "
                "selected; use --structure-model propka for PROPKA pKa values.")

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

    sets = list(PKA_SETS) if all_pka_sets else [pka_set]
    seq_table: Dict[str, object] = {}
    struct_table: Dict[str, object] = {}
    for s in sets:
        seq_table[s] = sequence_pi_table(chains_info, s)
        groups, termini = structure_groups(
            fvi, s, disulfide_keys, salt_keys, hbond_keys,
            propka_pka=propka_pka if structure_model == "propka" else None,
            model=structure_model,
        )
        pi_free = isoelectric_point(groups, termini)
        pi_capped = isoelectric_point(groups, termini, capped_cterm=True)
        q74 = net_charge(groups, termini, ph)
        qform = net_charge(groups, termini, formulation_ph)
        seq_entry = seq_table[s]  # type: ignore[assignment]
        delta = None
        if isinstance(seq_entry, dict) and seq_entry.get("pI_free") is not None and pi_free is not None:
            delta = pi_free - float(seq_entry["pI_free"])
        struct_table[s] = {
            "pI_free": pi_free,
            "pI_capped": pi_capped,
            "net_charge_ph": q74,
            "net_charge_ph74": q74,
            "net_charge_formulation_ph": qform,
            "delta_vs_sequence": delta,
            "model": structure_model,
            "n_groups": len(groups),
        }

    # Per-chain table and charge descriptors at the selected pKa set.
    groups, termini = structure_groups(
        fvi, pka_set, disulfide_keys, salt_keys, hbond_keys,
        propka_pka=propka_pka if structure_model == "propka" else None,
        model=structure_model,
    )
    groups_by_key = {(g.structure_chain, g.resseq, g.icode): g for g in groups}
    chain_table: Dict[str, Dict[str, object]] = {}
    for cid, ctype, label in ((vh, "H", "VH"), (vl, "L", "VL")):
        if cid is None:
            continue
        chain_groups = [g for g in groups if g.structure_chain == cid]
        chain_termini = [t for t in termini if t.chain == cid]
        seq = "".join(r.aa for r in fvi if r.chain == cid)
        chain_table[label] = {
            "chain_type": ctype,
            "structure_chain": cid,
            "length": len(seq),
            "sequence": seq,
            "pI_free": isoelectric_point(chain_groups, chain_termini),
            "pI_capped": isoelectric_point(chain_groups, chain_termini,
                                           capped_cterm=True),
            "net_charge_ph74": net_charge(chain_groups, chain_termini, ph),
            "net_charge_formulation_ph": net_charge(
                chain_groups, chain_termini, formulation_ph),
            "n_ionizable": len(chain_groups),
        }

    pi_free = struct_table[pka_set]["pI_free"]  # type: ignore[index]
    pi_capped = struct_table[pka_set]["pI_capped"]  # type: ignore[index]
    q74 = struct_table[pka_set]["net_charge_ph74"]  # type: ignore[index]
    qform = struct_table[pka_set]["net_charge_formulation_ph"]  # type: ignore[index]

    # Surface charges use TAP's exposed-residue definition at pH 7.4.
    fv: Dict[str, object] = {
        "pI_free": pi_free,
        "pI_capped": pi_capped,
        "net_charge_ph": q74,
        "net_charge_ph74": q74,
        "net_charge_formulation_ph": qform,
        "Fv_chml": (qh - ql) if (qh is not None and ql is not None) else None,
        "SFvCSP": sfvcsp,
        "surface_charge_vh": qh,
        "surface_charge_vl": ql,
    }

    # Patch analysis at formulation and physiological pH.
    patches: List[Dict[str, object]] = []
    for patch_ph in sorted({float(formulation_ph), float(ph)}):
        for patch in cluster_patches(
                fvi, pair_dist, patch_ph, groups_by_key=groups_by_key,
                interface_keys=interface_keys):
            row = {
                "patch_id": patch.patch_id,
                "sign": patch.sign,
                "ph": patch.ph,
                "total_charge": patch.total_charge,
                "score": patch.score,
                "risk": patch.risk,
                "n_cdr": patch.n_cdr,
                "n_framework": patch.n_framework,
                "n_anchor": patch.n_anchor,
                "in_cdr_vicinity": patch.in_cdr_vicinity,
                "chains": patch.chains,
                "residues": patch.residues,
                "mutation_candidates": mutation_candidates(
                    patch, fvi, interface_keys=interface_keys,
                    vhh_hallmark_keys=vhh_hallmark_keys),
            }
            patches.append(row)

    ionizable_rows = []
    for g in sorted(groups, key=lambda x: (x.structure_chain or "", x.resseq, x.icode)):
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
            "charge_ph74": round(g.charge(ph), 3),
            "charge_formulation_ph": round(g.charge(formulation_ph), 3),
            "reason": g.shift_reason,
        })

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
        "sequence_pI": seq_table,
        "structure_pI": struct_table,
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
            "salt_bridges": salt_pairs,
            "n_residues": counts["n_residues"],
        },
        "patches": patches,
        "ionizable_residues": ionizable_rows,
        "recommendations": [],
        "warnings": warnings,
    }
    result["recommendations"] = build_recommendations(result)
    return result
