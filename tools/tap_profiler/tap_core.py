#!/usr/bin/env python3
"""TAP (Therapeutic Antibody Profiler) core engine.

Faithful, dependency-light re-implementation of the five TAP developability
guidelines published in:

  * Raybould MIJ et al. "Five computational developability guidelines for
    therapeutic antibody profiling." PNAS 2019;116(10):4025-4030.
  * Raybould MIJ et al. "Contextualising the developability risk of antibodies
    with lambda light chains using enhanced therapeutic antibody profiling."
    Commun Biol 2024;7:62.

Metric definitions implemented here (identical to the published Methods):

  Ltot     total IMGT CDR length of the Fv (all six CDR loops).
  PSH      patches of surface hydrophobicity over the CDR vicinity:
           sum over ordered surface residue pairs (R1,R2) with a closest
           heavy-atom distance r12 < 7.5 A of  H(R1)*H(R2)/r12**2,
           H = Kyte-Doolittle hydrophobicity normalised to [1, 2].
  PPC/PNC  patches of positive/negative charge: the same sum with
           |Q(R1)|*|Q(R2)|/r12**2 and both residues of the requested sign.
  SFvCSP   structural Fv charge symmetry parameter = Q(VH) * Q(VL), the
           *product* of the net surface charges of the two chains
           (oppositely charged chains -> negative product -> higher risk).

Supporting definitions (also from the published Methods):

  surface-exposed : relative accessible surface area >= 7.5 %
                    (Ala-X-Ala reference when JOY `psa` is used)
  CDR vicinity    : surface-exposed IMGT CDR + anchor residues, plus every
                    surface-exposed residue whose closest heavy atom lies
                    within 4 A of one of those seeds
  salt bridge     : K/R N  vs  D/E O within 3.2 A, both residues exposed;
                    partners get charge 0 and glycine hydrophobicity
  IMGT CDRs       : 27-38, 56-65, 105-117 (insertions included)
  IMGT anchors    : 25, 26, 39, 40, 54, 55, 66, 67, 103, 104, 118, 119
  charges @pH 7.4 : D/E -1, K/R +1, H +0.1, everything else 0

Implementation notes (fixes over the legacy tap_analyzer.py):

  * metrics are computed with the published pairwise 1/r**2 formulas, not
    with sequence-order patch counting; PSH/PPC/PNC are therefore directly
    comparable to the published thresholds (which are areas/scores, not
    residue counts)
  * SFvCSP is a product, not a difference of chain charges
  * only the Fv chains participate (antigen/constant chains are excluded
    unless --scope complex is requested), and chains are identified by
    ANARCI rather than by chain letter or chain length
  * every structure-dependent step runs on a self-written work PDB with
    deterministic residue numbering, so insertion codes, missing residues
    and multi-model files cannot desynchronise the SASA <-> residue mapping
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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

# ---------------------------------------------------------------------------
# Constants (published TAP methodology)
# ---------------------------------------------------------------------------

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Kyte & Doolittle hydropathy, normalised to [1, 2] exactly as used by the
# reference TAP implementations ((value + 4.5) / 4.5 with KD in [-4.5, 4.5]).
NORMALISED_HYDROPHOBICITY = {
    "ILE": 2.0,
    "VAL": 1.9666666666666666,
    "LEU": 1.9222222222222223,
    "PHE": 1.8111111111111111,
    "CYS": 1.7777777777777777,
    "MET": 1.7111111111111112,
    "ALA": 1.7,
    "GLY": 1.4555555555555555,
    "THR": 1.4222222222222223,
    "SER": 1.4111111111111112,
    "TRP": 1.4,
    "TYR": 1.3555555555555556,
    "PRO": 1.3222222222222222,
    "HIS": 1.1444444444444444,
    "GLU": 1.1111111111111112,
    "GLN": 1.1111111111111112,
    "ASP": 1.1111111111111112,
    "ASN": 1.1111111111111112,
    "LYS": 1.0666666666666667,
    "ARG": 1.0,
}

# Formal charges at pH 7.4 (Henderson-Hasselbalch, pKa 6.0 -> +0.1).
RESIDUE_CHARGES = {
    "ALA": 0.0, "ARG": 1.0, "ASN": 0.0, "ASP": -1.0, "CYS": 0.0,
    "GLN": 0.0, "GLU": -1.0, "GLY": 0.0, "HIS": 0.1, "ILE": 0.0,
    "LEU": 0.0, "LYS": 1.0, "MET": 0.0, "PHE": 0.0, "PRO": 0.0,
    "SER": 0.0, "THR": 0.0, "TRP": 0.0, "TYR": 0.0, "VAL": 0.0,
}

# Theoretical maximum solvent-accessible area per residue (Tien et al. 2013),
# used for the FreeSASA fallback surface definition.
REF_SASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLU": 223.0, "GLN": 225.0, "GLY": 104.0,
    "HIS": 224.0, "ILE": 197.0, "LEU": 201.0, "LYS": 236.0,
    "MET": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

IMGT_CDR_RANGES = ((27, 38), (56, 65), (105, 117))
IMGT_ANCHOR_POS = (25, 26, 39, 40, 54, 55, 66, 67, 103, 104, 118, 119)

NEIGHBOUR_CUTOFF = 7.5       # A, closest heavy-atom distance for PSH/PPC/PNC pairs
SALT_BRIDGE_CUTOFF = 3.2     # A, donor N ... acceptor O
VICINITY_CUTOFF = 4.0        # A, CDR vicinity expansion
SURFACE_MIN_REL_SASA = 7.5   # %, side-chain exposure threshold

# Salt-bridge atom names (JOY/TAP convention).
SALT_DONORS = {"LYS": ("NZ",), "ARG": ("NH1", "NH2")}
SALT_ACCEPTORS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}

# ---------------------------------------------------------------------------
# Thresholds: green/amber regions; anything outside them is red.
# Flag precedence is amber -> green -> red (mirrors the reference code).
# ---------------------------------------------------------------------------

THRESHOLD_SETS: Dict[str, Dict[str, dict]] = {
    # Raybould 2024 (Commun Biol 7:62), ABodyBuilder2 models of 664 CSTs.
    "2024": {
        "Ltot": {"green": ((43, 54),), "amber": ((37, 42), (55, 63))},
        "PSH": {"green": ((110.11, 168.06),), "amber": ((95.58, 110.11), (168.06, 201.59))},
        "PPC": {"green": ((0.0, 1.32),), "amber": ((1.32, 4.22),)},
        "PNC": {"green": ((0.0, 2.00),), "amber": ((2.00, 4.42),)},
        "SFvCSP": {"green": ((-6.00, 1.0e9),), "amber": ((-30.60, -6.00),)},
    },
    # Raybould 2019 (PNAS) thresholds as tabulated in the 2024 update
    # (ABodyBuilder1 column) - kept for comparison with older results.
    "2019": {
        "Ltot": {"green": ((43, 54),), "amber": ((37, 42), (55, 63))},
        "PSH": {"green": ((114.61, 174.26),), "amber": ((90.88, 114.61), (174.26, 219.20))},
        "PPC": {"green": ((0.0, 1.27),), "amber": ((1.27, 4.21),)},
        "PNC": {"green": ((0.0, 1.93),), "amber": ((1.93, 4.36),)},
        "SFvCSP": {"green": ((-6.00, 1.0e9),), "amber": ((-30.00, -6.00),)},
    },
    # Thresholds hard-coded in Exscientia/ab-characterisation (MIT), which
    # recalibrated the bands on a private therapeutic set.
    "exscientia": {
        "Ltot": {"green": ((43, 55),), "amber": ((37, 43), (55, 63))},
        "PSH": {"green": ((137.61, 200.71),), "amber": ((106.44, 137.61), (200.71, 225.85))},
        "PPC": {"green": ((0.0, 1.19),), "amber": ((1.19, 3.58),)},
        "PNC": {"green": ((0.0, 1.67),), "amber": ((1.67, 3.50),)},
        "SFvCSP": {"green": ((-4.20, 1.0e9),), "amber": ((-20.50, -4.20),)},
    },
    # Bands used by the legacy tool (kept so old reports stay interpretable).
    # NOT the published values - prefer "2024".
    "legacy": {
        "Ltot": {"green": ((43, 54),), "amber": ((37, 42), (55, 65))},
        "PSH": {"green": ((111.41, 167.63),), "amber": ((95.77, 111.40), (167.64, 211.65))},
        "PPC": {"green": ((0.0, 1.33),), "amber": ((1.34, 4.20),)},
        "PNC": {"green": ((0.0, 1.98),), "amber": ((1.99, 4.43),)},
        "SFvCSP": {"green": ((-6.00, 1.0e9),), "amber": ((-30.60, -6.00),)},
    },
}

METRIC_ORDER = ("Ltot", "PSH", "PPC", "PNC", "SFvCSP")


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class Residue:
    """A residue of the analysed structure (heavy atoms only)."""
    chain: str
    resseq: int
    icode: str
    resname: str
    aa: str
    atoms: List[Tuple[str, str, float, float, float]] = field(default_factory=list)
    # (raw 4-char name, element, x, y, z)

    imgt_pos: Optional[int] = None
    imgt_ins: str = ""
    rel_sasa: float = 0.0
    is_surface: bool = False
    is_cdr: bool = False
    is_anchor: bool = False
    in_cdr_vicinity: bool = False
    hydrophobicity: float = 0.0
    charge: float = 0.0
    salt_bridge_partner: Optional[Tuple[str, int, str]] = None

    @property
    def key(self) -> Tuple[str, int, str]:
        return (self.chain, self.resseq, self.icode)


@dataclass
class TAPResult:
    """Everything the report layer needs."""
    name: str
    metrics: Dict[str, float]           # Ltot / PSH / PPC / PNC / SFvCSP (None if n/a)
    flags: Dict[str, str]               # green / amber / red / unknown
    risk_summary: str
    details: Dict


# ---------------------------------------------------------------------------
# Structure parsing (pure stdlib)
# ---------------------------------------------------------------------------

def _element_of(raw_name: str) -> str:
    """Best-effort element symbol from a 4-char PDB atom-name field."""
    name = raw_name.strip()
    if not name:
        return "?"
    if name[0].isdigit():
        return name[1] if len(name) > 1 else name[0]
    return name[0]


def parse_pdb(path: str) -> List[Residue]:
    """Parse ATOM records of a PDB file (first model only).

    Alternate locations are de-duplicated per atom name (first wins), waters
    and other HETATM groups are ignored, and residues are returned in file
    order.
    """
    residues: List[Residue] = []
    index: Dict[Tuple[str, int, str], int] = {}
    seen_atoms = set()
    in_model = False
    model_seen = False
    with open(path) as fh:
        for line in fh:
            rec = line[:6]
            if rec == "MODEL ":
                if model_seen:
                    break
                model_seen = True
                in_model = True
                continue
            if rec == "ENDMDL" and in_model:
                break
            if rec != "ATOM  ":          # waters/ligands are not part of the Fv
                continue
            if len(line) < 54:
                continue
            if line[16] not in (" ", "A"):      # altLoc
                continue
            resname = line[17:20].strip()
            if resname not in AA3TO1:
                continue
            try:
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            chain = line[21]
            icode = line[26] if len(line) > 26 and line[26] != " " else " "
            raw_name = line[12:16] if len(line) >= 16 else line[12:].ljust(4)
            name = raw_name.strip()
            if not name or name.startswith("H") or _element_of(raw_name) == "H":
                continue
            element = _element_of(raw_name)
            if element not in ("C", "N", "O", "S"):
                element = name[0]
            key = (chain, resseq, icode)
            atom_key = (key, name)
            if atom_key in seen_atoms:
                continue
            seen_atoms.add(atom_key)
            if key not in index:
                index[key] = len(residues)
                residues.append(Residue(
                    chain=chain, resseq=resseq, icode=icode,
                    resname=resname, aa=AA3TO1[resname],
                ))
            residues[index[key]].atoms.append((raw_name, element, x, y, z))
    return residues


def parse_cif(path: str) -> List[Residue]:
    """Parse the _atom_site loop of an mmCIF file (PDBx/mmCIF)."""
    residues: List[Residue] = []
    index: Dict[Tuple[str, int, str], int] = {}
    seen_atoms = set()
    cols: List[str] = []
    in_loop = False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                in_loop = True
                cols.append(s[len("_atom_site."):].strip())
                continue
            if in_loop and (s.startswith("#") or s == "" or s.startswith("loop_")):
                in_loop = False
                cols = []
                continue
            if not in_loop or s.startswith(("_", "#")):
                continue
            parts = s.split()
            if len(parts) < len(cols):
                continue
            d = dict(zip(cols, parts))
            if d.get("group_PDB") not in ("ATOM", "HETATM"):
                continue
            if d.get("label_atom_id") in ("H", "D"):
                continue
            if d.get("pdbx_PDB_model_num", "1") != "1":
                continue
            resname = d.get("auth_comp_id", d.get("label_comp_id", ""))
            if resname not in AA3TO1:
                continue
            altloc = d.get("label_alt_id", ".")
            if altloc not in (".", "?", "A"):
                continue
            try:
                resseq = int(d.get("auth_seq_id", d.get("label_seq_id", "0")))
                x = float(d["Cartn_x"])
                y = float(d["Cartn_y"])
                z = float(d["Cartn_z"])
            except (ValueError, KeyError):
                continue
            chain = d.get("auth_asym_id", d.get("label_asym_id", " "))
            ic = d.get("pdbx_PDB_ins_code", "?")
            icode = " " if ic in (".", "?") else ic
            atom_field = d.get("auth_atom_id", d.get("label_atom_id", "?"))
            raw_name = " %-3s" % atom_field if len(atom_field) < 4 else atom_field[:4]
            name = raw_name.strip()
            if not name or name[0] == "H" or (name[0].isdigit() and name[1] == "H"):
                continue
            element = d.get("type_symbol", _element_of(raw_name)).capitalize()
            if element not in ("C", "N", "O", "S"):
                element = name[0]
            key = (chain, resseq, icode)
            atom_key = (key, name)
            if atom_key in seen_atoms:
                continue
            seen_atoms.add(atom_key)
            if key not in index:
                index[key] = len(residues)
                residues.append(Residue(
                    chain=chain, resseq=resseq, icode=icode,
                    resname=resname, aa=AA3TO1[resname],
                ))
            residues[index[key]].atoms.append((raw_name, element, x, y, z))
    return residues


def load_residues(path: str) -> List[Residue]:
    lower = path.lower()
    if lower.endswith((".pdb", ".ent")):
        residues = parse_pdb(path)
    elif lower.endswith((".cif", ".mmcif")):
        residues = parse_cif(path)
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    if not residues:
        raise ValueError(f"No protein residues found in {path}")
    # Drop residues whose atoms all vanished (e.g. hydrogen-only records).
    return [r for r in residues if r.atoms]


# ---------------------------------------------------------------------------
# Chain selection and IMGT numbering
# ---------------------------------------------------------------------------

def _chain_order(residues: List[Residue]) -> List[str]:
    seen: Dict[str, None] = {}
    for r in residues:
        seen.setdefault(r.chain, None)
    return list(seen)


def _chain_residues(residues: List[Residue], chain: str) -> List[Residue]:
    return [r for r in residues if r.chain == chain]


def _anarci_chain_type(sequence: str) -> Optional[str]:
    """Return 'H'/'K'/'L' when ANARCI recognises the sequence as an antibody
    domain, otherwise None (non-antibody chain)."""
    try:
        from anarci import anarci as _anarci
    except Exception:
        return None
    try:
        result = _anarci([("query", sequence)], scheme="imgt", output=False)
    except Exception:
        return None
    try:
        alignments = result[1]
        if alignments and alignments[0]:
            for hit in alignments[0]:
                if isinstance(hit, dict):
                    ctype = hit.get("chain_type", "")
                    if ctype in ("H", "K", "L"):
                        return ctype
    except (IndexError, TypeError, KeyError):
        return None
    return None


def select_fv_chains(residues: List[Residue],
                     vh_chain: Optional[str] = None,
                     vl_chain: Optional[str] = None) -> Tuple[str, str, List[str]]:
    """Identify the VH and VL chains of the structure.

    Preference order: explicit CLI values, ANARCI chain-type detection,
    chain-name heuristics (H/L then A/B). Raises when nothing sensible is
    found.
    """
    warnings: List[str] = []
    chains = _chain_order(residues)
    if vh_chain and vl_chain:
        for c in (vh_chain, vl_chain):
            if c not in chains:
                raise ValueError(f"Requested chain {c!r} not found in structure "
                                 f"(chains present: {', '.join(chains)})")
        return vh_chain, vl_chain, warnings

    def seq_of(chain: str) -> str:
        return "".join(r.aa for r in _chain_residues(residues, chain))

    types: Dict[str, Optional[str]] = {}
    anarci_ok = False
    try:
        from scripts.humanize.anarci_adapter import is_anarci_available
        anarci_ok = is_anarci_available()
    except Exception:
        anarci_ok = False

    if anarci_ok:
        for chain in chains:
            types[chain] = _anarci_chain_type(seq_of(chain))
        heaviest = max(chains, key=lambda c: len(_chain_residues(residues, c)))
        vh_candidates = [c for c in chains if types[c] == "H"] or []
        vl_candidates = [c for c in chains if types[c] in ("K", "L")] or []
        if vh_candidates and vl_candidates:
            vh = max(vh_candidates, key=lambda c: len(_chain_residues(residues, c)))
            vl = max(vl_candidates, key=lambda c: len(_chain_residues(residues, c)))
            return vh, vl, warnings
        warnings.append(
            "ANARCI could not classify both chains as antibody domains; "
            "falling back to chain-name heuristics")
    else:
        warnings.append("ANARCI not installed - chain selection uses name heuristics "
                        "(install anarci for reliable detection)")

    # Heuristic fallback.
    by_name = {c.upper(): c for c in chains}
    vh = by_name.get("H") or by_name.get("A")
    vl = by_name.get("L") or by_name.get("B")
    if vh and vl and vh != vl:
        return vh, vl, warnings
    if len(chains) == 2:
        a, b = chains
        longer, shorter = (a, b) if len(_chain_residues(residues, a)) >= len(
            _chain_residues(residues, b)) else (b, a)
        warnings.append(f"Chain names are not H/L - assumed VH={longer}, VL={shorter}")
        return longer, shorter, warnings
    raise ValueError(
        f"Cannot identify the Fv chains automatically (chains: {', '.join(chains)}). "
        "Pass --vh-chain/--vl-chain explicitly.")


def _parse_pos_label(label: str) -> Tuple[Optional[int], str]:
    m = re.match(r"^[A-Z](\d+)([A-Za-z]*)$", label)
    if not m:
        return None, ""
    return int(m.group(1)), m.group(2)


def number_chain_imgt(residues: List[Residue], chain_type: str) -> List[Residue]:
    """Number one chain with ANARCI/IMGT and annotate every residue.

    Residues beyond the V domain (constant region) are stripped, because TAP
    is defined on the Fv. Raises when ANARCI is unavailable or fails.
    """
    from scripts.humanize.anarci_adapter import is_anarci_available, number_with_anarci_imgt

    if not is_anarci_available():
        raise RuntimeError(
            "ANARCI is required for IMGT numbering but is not installed "
            "(pip install anarci). Alternatively use --numbering native on a "
            "structure that is already IMGT-numbered.")
    sequence = "".join(r.aa for r in residues)
    chain = number_with_anarci_imgt(sequence, chain_type)
    if chain is None or not chain.residues:
        raise RuntimeError(f"ANARCI failed to number a {chain_type} chain "
                           f"({len(sequence)} residues)")
    n = len(chain.residues)
    if n > len(residues):
        raise RuntimeError("ANARCI returned more residues than the chain has")
    if n < len(residues):
        # Constant region appended after the V domain: keep only the V part.
        residues = residues[:n]
    for res, nr in zip(residues, chain.residues):
        pos, ins = _parse_pos_label(nr.pos)
        res.imgt_pos = pos
        res.imgt_ins = ins
        res.is_cdr = pos is not None and any(lo <= pos <= hi for lo, hi in IMGT_CDR_RANGES)
        res.is_anchor = pos is not None and pos in IMGT_ANCHOR_POS
    return residues


def number_chain_native(residues: List[Residue]) -> List[Residue]:
    """Interpret existing residue numbers as IMGT positions (no ANARCI)."""
    for res in residues:
        res.imgt_pos = res.resseq
        res.imgt_ins = res.icode.strip()
        res.is_cdr = any(lo <= res.resseq <= hi for lo, hi in IMGT_CDR_RANGES)
        res.is_anchor = res.resseq in IMGT_ANCHOR_POS
    return residues


# ---------------------------------------------------------------------------
# Solvent accessibility
# ---------------------------------------------------------------------------

def write_work_pdb(residues: List[Residue], path: str
                   ) -> List[Tuple[str, int]]:
    """Write a minimal, canonically numbered PDB used for SASA.

    Residues are re-numbered 1..N per chain so that the SASA engine returns
    exactly the residue keys we expect (insertion codes and original gaps can
    never desynchronise the mapping). Returns the (chain, resseq) keys in
    file order.
    """
    keys: List[Tuple[str, int]] = []
    serial = 1
    counters: Dict[str, int] = defaultdict(int)
    with open(path, "w") as fh:
        for res in residues:
            counters[res.chain] += 1
            resseq = counters[res.chain]
            keys.append((res.chain, resseq))
            for raw_name, element, x, y, z in res.atoms:
                name = raw_name if len(raw_name) == 4 else " %-3s" % raw_name.strip()
                fh.write("ATOM  {:5d} {:4s} {:>3s} {:1s}{:4d}    "
                         "{:8.3f}{:8.3f}{:8.3f}  1.00 20.00          {:>2s}\n"
                         .format(serial, name, res.resname, res.chain, resseq,
                                 x, y, z, element))
                serial += 1
        fh.write("END\n")
    return keys


def _find_psa(psa_binary: Optional[str] = None) -> Optional[str]:
    if psa_binary:
        path = os.path.expanduser(psa_binary)
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
        raise FileNotFoundError(f"--psa-binary not executable: {psa_binary}")
    return shutil.which("psa")


def _rel_sasa_psa(work_pdb: str, keys: List[Tuple[str, int]],
                  residues: List[Residue], psa_path: str) -> Dict[int, float]:
    """Side-chain relative ASA (%) exactly as JOY psa reports it.

    Each ACCESS line carries five (area, percent) pairs; the reference
    implementations read the fourth pair's percentage (fixed column 61:67),
    which is the side-chain value relative to the Ala-X-Ala reference.
    """
    proc = subprocess.run([psa_path, "-t", work_pdb],
                          capture_output=True, text=True)
    lines = [ln for ln in proc.stdout.splitlines() if ln.startswith("ACCESS")]
    if len(lines) != len(keys):
        raise RuntimeError(
            f"JOY psa returned {len(lines)} residue lines for {len(keys)} residues: "
            f"{(proc.stderr or proc.stdout).strip()[:300]}")
    rel: Dict[int, float] = {}
    for idx, ((chain, resseq), line) in enumerate(zip(keys, lines)):
        res = residues[idx]
        reported_num = line[6:12].strip()
        reported_res = line[14:17]
        if reported_num and (reported_num != str(resseq) or reported_res != res.resname):
            raise RuntimeError(
                f"JOY psa residue mismatch at line {idx}: expected "
                f"{resseq} {res.resname}, got {reported_num} {reported_res}")
        try:
            value = float(line[61:67])
        except ValueError:
            nums = re.findall(r"-?\d+\.\d+", line[17:])
            if len(nums) < 8:
                raise RuntimeError(f"Cannot parse psa ACCESS line: {line!r}")
            value = float(nums[7])      # fourth percent value
        rel[idx] = value
    return rel


def _rel_sasa_freesasa(work_pdb: str, keys: List[Tuple[str, int]],
                       residues: List[Residue]) -> Dict[int, float]:
    """Whole-residue relative ASA (%) via FreeSASA + Tien et al. maxima.

    Measured agreement with the JOY psa surface/buried decision: 92.9 % of
    residues on a Fab+antigen test structure (see README).
    """
    try:
        import freesasa
    except ImportError as exc:
        raise RuntimeError(
            "No SASA engine available: install FreeSASA (pip install freesasa) "
            "or point --psa-binary at the JOY psa executable") from exc
    structure = freesasa.Structure(work_pdb)
    calc = freesasa.calc(structure)
    per_res: Dict[Tuple[str, int], float] = defaultdict(float)
    for i in range(structure.nAtoms()):
        key = (structure.chainLabel(i), int(str(structure.residueNumber(i)).strip()))
        per_res[key] += calc.atomArea(i)
    rel: Dict[int, float] = {}
    for idx, key in enumerate(keys):
        ref = REF_SASA.get(residues[idx].resname, 150.0)
        area = per_res.get(key, 0.0)
        rel[idx] = 100.0 * area / ref if ref > 0 else 0.0
    return rel


def compute_rel_sasa(residues: List[Residue], engine: str = "auto",
                     psa_binary: Optional[str] = None
                     ) -> Tuple[Dict[int, float], str, List[str]]:
    """Compute relative SASA for every residue.

    Returns ({residue_index: rel_sasa_percent}, engine_used, warnings).
    """
    warnings: List[str] = []
    psa_path = None
    if engine in ("auto", "psa"):
        try:
            psa_path = _find_psa(psa_binary)
        except FileNotFoundError:
            raise
    if psa_path:
        with tempfile.TemporaryDirectory() as tmpdir:
            work = os.path.join(tmpdir, "fv.pdb")
            keys = write_work_pdb(residues, work)
            return _rel_sasa_psa(work, keys, residues, psa_path), "psa (JOY/Ala-X-Ala)", warnings
    if engine == "psa":
        raise RuntimeError("JOY psa requested (--sasa psa) but not found on PATH")
    warnings.append(
        "FreeSASA fallback: surface exposure uses whole-residue rSASA vs Tien "
        "2013 maxima instead of the JOY psa Ala-X-Ala side-chain reference "
        "(~7 % of residues flip surface/buried decision); install JOY psa for "
        "exact TAP parity")
    with tempfile.TemporaryDirectory() as tmpdir:
        work = os.path.join(tmpdir, "fv.pdb")
        keys = write_work_pdb(residues, work)
        return _rel_sasa_freesasa(work, keys, residues), "FreeSASA (Tien rSASA)", warnings


# ---------------------------------------------------------------------------
# Geometry: neighbour distances and salt bridges
# ---------------------------------------------------------------------------

def build_pair_distances(residues: List[Residue],
                         cutoff: float = NEIGHBOUR_CUTOFF) -> Dict[Tuple[int, int], float]:
    """Closest heavy-atom distance for every residue pair within `cutoff`.

    A uniform spatial hash (cell size == cutoff) keeps this O(n) in practice
    while staying pure Python.
    """
    cell = cutoff
    grid: Dict[Tuple[int, int, int], List[Tuple[int, float, float, float]]] = defaultdict(list)
    for idx, res in enumerate(residues):
        for _raw, _elem, x, y, z in res.atoms:
            grid[(int(math.floor(x / cell)), int(math.floor(y / cell)),
                  int(math.floor(z / cell)))].append((idx, x, y, z))
    # Only "forward" neighbour offsets so each cell pair is visited once.
    offsets = [(0, 0, 0)]
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) == (0, 0, 0):
                    continue
                if dx > 0 or (dx == 0 and dy > 0) or (dx == 0 and dy == 0 and dz > 0):
                    offsets.append((dx, dy, dz))
    limit = cutoff * cutoff
    best: Dict[Tuple[int, int], float] = {}
    for (gx, gy, gz), items in grid.items():
        for ox, oy, oz in offsets:
            others = grid.get((gx + ox, gy + oy, gz + oz))
            if not others:
                continue
            same_cell = (ox, oy, oz) == (0, 0, 0)
            for i, x1, y1, z1 in items:
                for j, x2, y2, z2 in others:
                    if same_cell and j <= i:
                        continue          # each unordered pair once per cell
                    dx = x1 - x2
                    dy = y1 - y2
                    dz = z1 - z2
                    d2 = dx * dx + dy * dy + dz * dz
                    if d2 > limit:
                        continue
                    pair = (i, j) if i < j else (j, i)
                    cur = best.get(pair)
                    if cur is None or d2 < cur:
                        best[pair] = d2
    return {pair: math.sqrt(d2) for pair, d2 in best.items()}


def assign_salt_bridges(residues: List[Residue]) -> int:
    """Flag K/R ... D/E salt bridges (<= 3.2 A, both partners exposed).

    Each residue participates at most once; assignment order is residue file
    order, which makes the result deterministic. Returns the number of
    bridges assigned.
    """
    donors: List[Tuple[int, str, float, float, float]] = []
    acceptors: List[Tuple[int, str, float, float, float]] = []
    for idx, res in enumerate(residues):
        for raw_name, _elem, x, y, z in res.atoms:
            name = raw_name.strip()
            if name in SALT_DONORS.get(res.resname, ()):
                donors.append((idx, name, x, y, z))
            elif name in SALT_ACCEPTORS.get(res.resname, ()):
                acceptors.append((idx, name, x, y, z))
    if not donors or not acceptors:
        return 0
    cell = SALT_BRIDGE_CUTOFF
    grid: Dict[Tuple[int, int, int], List[Tuple[int, str, float, float, float]]] = defaultdict(list)
    for rec in acceptors:
        _idx, _name, x, y, z = rec
        grid[(int(math.floor(x / cell)), int(math.floor(y / cell)),
              int(math.floor(z / cell)))].append(rec)
    offsets = [(0, 0, 0)]
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) != (0, 0, 0):
                    offsets.append((dx, dy, dz))
    limit = cell * cell
    assigned = 0
    for d_idx, d_name, x1, y1, z1 in donors:
        d_res = residues[d_idx]
        if d_res.salt_bridge_partner or not d_res.is_surface:
            continue
        gx, gy, gz = int(math.floor(x1 / cell)), int(math.floor(y1 / cell)), int(math.floor(z1 / cell))
        best_acceptor: Optional[Tuple[int, str, float]] = None
        for ox, oy, oz in offsets:
            for a_idx, a_name, x2, y2, z2 in grid.get((gx + ox, gy + oy, gz + oz), ()):
                a_res = residues[a_idx]
                if a_idx == d_idx or a_res.salt_bridge_partner or not a_res.is_surface:
                    continue
                dx = x1 - x2
                dy = y1 - y2
                dz = z1 - z2
                d2 = dx * dx + dy * dy + dz * dz
                if d2 > limit:
                    continue
                cand = (a_idx, a_name, d2)
                if best_acceptor is None or d2 < best_acceptor[2]:
                    best_acceptor = cand
        if best_acceptor is None:
            continue
        a_idx = best_acceptor[0]
        d_res.salt_bridge_partner = residues[a_idx].key
        residues[a_idx].salt_bridge_partner = d_res.key
        assigned += 1
    return assigned


def build_adjacency(pair_dist: Dict[Tuple[int, int], float],
                    n: int, cutoff: float) -> Dict[int, List[Tuple[int, float]]]:
    adj: Dict[int, List[Tuple[int, float]]] = defaultdict(list)
    for (i, j), d in pair_dist.items():
        if d < cutoff:
            adj[i].append((j, d))
            adj[j].append((i, d))
    for key in adj:
        adj[key].sort()
    return adj


# ---------------------------------------------------------------------------
# Annotation and metrics
# ---------------------------------------------------------------------------

def annotate_residues(residues: List[Residue]) -> Dict[str, int]:
    """Apply surface, salt-bridge, property and CDR-vicinity annotations.

    Must run after IMGT numbering. Returns counts for reporting.
    """
    n = len(residues)
    for res in residues:                     # surface flag first: salt bridges
        res.is_surface = res.rel_sasa >= SURFACE_MIN_REL_SASA  # need it
    pair_dist = build_pair_distances(residues)
    n_bridges = assign_salt_bridges(residues)

    for res in residues:
        if res.salt_bridge_partner:
            res.hydrophobicity = NORMALISED_HYDROPHOBICITY["GLY"]
            res.charge = 0.0
        else:
            res.hydrophobicity = NORMALISED_HYDROPHOBICITY.get(
                res.resname, NORMALISED_HYDROPHOBICITY["GLY"])
            res.charge = RESIDUE_CHARGES.get(res.resname, 0.0)

    # CDR vicinity: seeds are exposed CDR/anchor residues; expand to every
    # exposed residue within 4 A (closest heavy-atom distance).
    adjacency = build_adjacency(pair_dist, n, VICINITY_CUTOFF)
    seeds = [i for i, r in enumerate(residues)
             if r.is_surface and (r.is_cdr or r.is_anchor)]
    for i in seeds:
        residues[i].in_cdr_vicinity = True
    for i in seeds:
        for j, dist in adjacency.get(i, ()):
            if dist < VICINITY_CUTOFF and residues[j].is_surface:
                residues[j].in_cdr_vicinity = True

    return {
        "n_residues": n,
        "n_surface": sum(1 for r in residues if r.is_surface),
        "n_cdr": sum(1 for r in residues if r.is_cdr),
        "n_cdr_vicinity": sum(1 for r in residues if r.in_cdr_vicinity),
        "n_salt_bridges": n_bridges,
        "n_neighbour_pairs": len(pair_dist),
    }


def compute_surface_scores(residues: List[Residue],
                           pair_dist: Dict[Tuple[int, int], float]
                           ) -> Tuple[float, float, float]:
    """PSH, PPC, PNC as published: ordered-pair sums of x1*x2/r**2."""
    vicinity = [i for i, r in enumerate(residues) if r.in_cdr_vicinity]
    psh = ppc = pnc = 0.0
    for i in vicinity:
        ri = residues[i]
        for j in vicinity:
            if j == i:
                continue
            pair = (i, j) if i < j else (j, i)
            dist = pair_dist.get(pair)
            if dist is None or dist <= 1e-6:
                continue          # missing or degenerate (overlapping) geometry
            rj = residues[j]
            if ri.charge > 0 and rj.charge > 0:
                ppc += (ri.charge * rj.charge) / (dist * dist)
            if ri.charge < 0 and rj.charge < 0:
                pnc += (ri.charge * rj.charge) / (dist * dist)
            psh += (ri.hydrophobicity * rj.hydrophobicity) / (dist * dist)
    return psh, ppc, pnc


def compute_sfvcsp(residues: List[Residue], vh_chain: str, vl_chain: str) -> Tuple[float, float, float]:
    """SFvCSP = (net surface charge of VH) * (net surface charge of VL)."""
    qh = sum(r.charge for r in residues
             if r.chain == vh_chain and r.is_surface)
    ql = sum(r.charge for r in residues
             if r.chain == vl_chain and r.is_surface)
    return qh * ql, qh, ql


def flag_value(value: float, spec: dict) -> str:
    """amber -> green -> red precedence, as in the reference implementations."""
    for lo, hi in spec.get("amber", ()):
        if lo <= value <= hi:
            return "amber"
    for lo, hi in spec.get("green", ()):
        if lo <= value <= hi:
            return "green"
    return "red"


def risk_summary(flags: Dict[str, str]) -> str:
    red = sum(1 for f in flags.values() if f == "red")
    amber = sum(1 for f in flags.values() if f == "amber")
    if red >= 2:
        return "HIGH RISK"
    if red == 1:
        return "MEDIUM-HIGH RISK"
    if amber >= 2:
        return "MEDIUM RISK"
    if amber == 1:
        return "LOW-MEDIUM RISK"
    return "LOW RISK"


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def analyze_structure(pdb_path: str, name: str = "",
                      vh_chain: Optional[str] = None,
                      vl_chain: Optional[str] = None,
                      sasa_engine: str = "auto",
                      psa_binary: Optional[str] = None,
                      numbering: str = "anarci",
                      scope: str = "fv",
                      thresholds: str = "2024") -> TAPResult:
    """Run the full TAP analysis on a PDB/mmCIF structure.

    scope='fv' (default) analyses the Fv chains only, matching the official
    TAP protocol; scope='complex' keeps every chain so surface burial by
    antigen is taken into account.
    """
    if thresholds not in THRESHOLD_SETS:
        raise ValueError(f"Unknown threshold set {thresholds!r}; "
                         f"choose from {sorted(THRESHOLD_SETS)}")
    if numbering not in ("anarci", "native"):
        raise ValueError("--numbering must be 'anarci' or 'native'")
    if scope not in ("fv", "complex"):
        raise ValueError("--scope must be 'fv' or 'complex'")
    if not name:
        name = Path(pdb_path).stem

    all_residues = load_residues(pdb_path)
    vh, vl, chain_warnings = select_fv_chains(all_residues, vh_chain, vl_chain)
    warnings = list(chain_warnings)

    # IMGT numbering + V-domain trimming for the Fv chains.
    fvi: List[Residue] = []
    for chain in (vh, vl):
        chain_res = _chain_residues(all_residues, chain)
        if numbering == "native":
            fvi.extend(number_chain_native(chain_res))
        else:
            ctype = "H" if chain == vh else "L"
            fvi.extend(number_chain_imgt(chain_res, ctype))
    fvi_ids = set(id(r) for r in fvi)

    if scope == "fv":
        analysed = fvi
    else:
        analysed = [r for r in all_residues if id(r) in fvi_ids or r.chain not in (vh, vl)]

    rel_sasa, engine, sasa_warnings = compute_rel_sasa(
        analysed, engine=sasa_engine, psa_binary=psa_binary)
    warnings.extend(sasa_warnings)
    for idx, res in enumerate(analysed):
        res.rel_sasa = rel_sasa.get(idx, 0.0)

    counts = annotate_residues(analysed)

    pair_dist = build_pair_distances(analysed)
    psh, ppc, pnc = compute_surface_scores(analysed, pair_dist)
    sfdcsp, qh, ql = compute_sfvcsp(analysed, vh, vl)
    ltot = sum(1 for r in fvi if r.is_cdr)

    cdr_lengths = {}
    for chain in (vh, vl):
        rows = {}
        for cdr_idx, (lo, hi) in enumerate(IMGT_CDR_RANGES, start=1):
            rows[f"CDR{cdr_idx}"] = sum(
                1 for r in fvi
                if r.chain == chain and r.imgt_pos is not None
                and lo <= r.imgt_pos <= hi)
        cdr_lengths[chain] = rows

    values = {"Ltot": ltot, "PSH": psh, "PPC": ppc, "PNC": pnc, "SFvCSP": sfdcsp}
    spec = THRESHOLD_SETS[thresholds]
    flags = {m: flag_value(values[m], spec[m]) for m in METRIC_ORDER}
    summary = risk_summary(flags)

    def seq_of(chain: str) -> str:
        return "".join(r.aa for r in _chain_residues(all_residues, chain))

    details = {
        "pdb_file": pdb_path,
        "numbering": "IMGT (ANARCI)" if numbering == "anarci" else "IMGT (structure-native)",
        "sasa_engine": engine,
        "scope": scope,
        "thresholds": thresholds,
        "vh_chain": vh,
        "vl_chain": vl,
        "vh_sequence": seq_of(vh),
        "vl_sequence": seq_of(vl),
        "cdr_lengths": cdr_lengths,
        "n_surface": counts["n_surface"],
        "n_cdr_vicinity": counts["n_cdr_vicinity"],
        "n_cdr_residues": counts["n_cdr"],
        "n_salt_bridges": counts["n_salt_bridges"],
        "net_charge_vh": qh,
        "net_charge_vl": ql,
        "n_residues_analysed": counts["n_residues"],
        "warnings": warnings,
    }
    return TAPResult(name=name, metrics=values, flags=flags,
                     risk_summary=summary, details=details)


def analyze_sequences(vh_sequence: str, vl_sequence: str, name: str = "",
                      thresholds: str = "2024") -> TAPResult:
    """Sequence-only TAP: Ltot is exact, the four structure metrics are n/a."""
    if thresholds not in THRESHOLD_SETS:
        raise ValueError(f"Unknown threshold set {thresholds!r}")
    if not name:
        name = "sequence"

    from scripts.humanize.anarci_adapter import is_anarci_available, number_with_anarci_imgt
    if not is_anarci_available():
        raise RuntimeError("ANARCI is required for sequence-only TAP "
                           "(pip install anarci)")

    cdr_lengths: Dict[str, Dict[str, int]] = {}
    ltot = 0
    warnings = ["Sequence-only mode: PSH/PPC/PNC/SFvCSP require a structure"]
    for label, seq in (("VH", vh_sequence), ("VL", vl_sequence)):
        chain = number_with_anarci_imgt(seq.upper().strip(),
                                        "H" if label == "VH" else "L")
        if chain is None:
            raise RuntimeError(f"ANARCI failed to number the {label} chain")
        rows = {}
        for cdr_idx, (lo, hi) in enumerate(IMGT_CDR_RANGES, start=1):
            n = 0
            for nr in chain.residues:
                pos, _ins = _parse_pos_label(nr.pos)
                if pos is not None and lo <= pos <= hi:
                    n += 1
            rows[f"CDR{cdr_idx}"] = n
            ltot += n
        cdr_lengths[label] = rows

    values: Dict[str, float] = {"Ltot": float(ltot), "PSH": None, "PPC": None,
                                "PNC": None, "SFvCSP": None}  # type: ignore
    flags: Dict[str, str] = {"Ltot": flag_value(ltot, THRESHOLD_SETS[thresholds]["Ltot"])}
    for m in METRIC_ORDER[1:]:
        flags[m] = "unknown"

    return TAPResult(
        name=name, metrics=values, flags=flags,
        risk_summary="UNKNOWN (sequence-only analysis)",
        details={
            "analysis_type": "sequence_only",
            "numbering": "IMGT (ANARCI)",
            "thresholds": thresholds,
            "vh_sequence": vh_sequence.upper().strip(),
            "vl_sequence": vl_sequence.upper().strip(),
            "cdr_lengths": cdr_lengths,
            "warnings": warnings,
        })
