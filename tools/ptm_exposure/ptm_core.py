#!/usr/bin/env python3
"""Core engine for the PTM exposure tool (structure parsing, exposure, risk).

Key design decisions:

* Every residue is addressed by (chain, resSeq, iCode) - never by a running
  index - so structures with insertion codes or missing residues keep an
  exact one-to-one mapping between sequence position, residue identity and
  solvent accessibility.
* Missing (disordered) residues are reconstructed as 'X' placeholders before
  numbering, so Kabat/IMGT labels stay correct for every residue that IS
  present (the numbering engine otherwise shifts every label after the gap).
* Waters and ions are removed before the SASA calculation (they are solvent,
  not structure); glycans and other covalently attached hetero groups are
  kept because they genuinely occlude the surface.
* Each PTM site is reported once per residue (the most severe motif is
  primary; other motifs are listed as additional_motifs).
* Conservative Cys (intradomain disulfides) are excluded entirely; all
  other Cys are checked for real disulfide partners (SG-SG <= 2.5 A).
* Canonical VH N-glycosylation site (Kabat H28-H33 loop) is identified
  structurally (glycan proximity) or by position and flagged as
  "canonical-glyco-site" rather than a de novo risk.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
import tempfile
import textwrap
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Any

try:
    import freesasa
except ImportError:
    freesasa = None  # type: ignore

# Repo root for imports
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from scripts.humanize.numbering import number_heavy, number_light
from scripts.humanize.anarci_adapter import number_with_anarci_imgt
from scripts.humanize.config import RISK_MOTIFS

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# 3-letter -> 1-letter (standard residues)
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# 1-letter -> 3-letter (reverse mapping for maxASA lookup)
AA_1TO3 = {v: k for k, v in AA_3TO1.items()}

# Common modified residues -> standard AA
MODIFIED_AA = {
    'MSE': 'M',   # selenomethionine
    'SEP': 'S',   # phosphoserine
    'TPO': 'T',   # phosphothreonine
    'PTR': 'Y',   # phosphotyrosine
    'HYP': 'P',   # hydroxyproline
    'CSO': 'C',   # cysteine sulfinic acid
    'CME': 'C',   # S-methylcysteine
    'MLY': 'K',   # methyllysine
    'KCX': 'K',   # carboxylated lysine
    'LLP': 'K',   # pyrrolidone carboxylic acid
    'M3L': 'K',   # trimethyllysine
    'ASX': 'N',   # Asp/Asn ambiguous
    'GLX': 'Q',   # Glu/Gln ambiguous
    'PCA': 'E',   # pyroglutamic acid
    'FME': 'M',   # formylmethionine
    'SAH': None,  # S-adenosylhomocysteine (non-polymer)
}

# Residues to REMOVE from SASA calculation (solvent/buffer)
WATER_ION_RESNAMES = {
    'HOH', 'DOD', 'WAT', 'TIP', 'H2O',  # waters
    'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CO', 'NI', 'CU', 'CD',  # ions
    'SO4', 'PO4', 'NO3', 'ACT', 'EDO', 'GOL', 'PEG', 'MPD',  # common buffer/salts (small)
}

# Max SASA tables
MAX_SASA_TIEN2013 = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0,
    'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
    'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
    'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
}

MAX_SASA_MILLER1987 = {
    'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0,
    'CYS': 140.0, 'GLU': 183.0, 'GLN': 179.0, 'GLY': 85.0,
    'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0,
    'MET': 203.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0,
    'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0,
}

# Conserved intradomain disulfide Cys (Kabat numbering)
CONSERVED_CYS = {
    "H": {22, 92},      # VH
    "L": {22, 23, 88},  # VL (kappa 23/88, lambda 22/88)
}

# Risk levels order
_RISK_ORDER = {"high": 2, "medium": 1, "low": 0}

# Default per-motif exposure thresholds and base risk
# Motif names MUST match those in RISK_MOTIFS (config.py)
MOTIF_SPECS = {
    # deamidation
    "deamidation (NG)": ("high", 0.30),
    "deamidation (NS)": ("medium", 0.40),
    "deamidation (NH)": ("medium", 0.40),
    "deamidation (ND)": ("low", 0.50),
    # isomerization
    "isomerization (DG)": ("high", 0.30),
    "isomerization (DS)": ("medium", 0.40),
    "isomerization (DT)": ("medium", 0.40),
    "isomerization (DH)": ("low", 0.50),
    # acid hydrolysis
    "acid hydrolysis (D-X)": ("medium", 0.40),
    "acid hydrolysis (DD)": ("high", 0.30),
    # oxidation
    "oxidation (M)": ("medium", 0.20),
    "oxidation (W)": ("medium", 0.20),
    "oxidation (C)": ("medium", 0.20),
    # unpaired Cys (real pairing checked separately)
    "unpaired Cys": ("medium", 0.20),
    # N-glycosylation
    "N-glycan (NxS/T)": ("high", 0.30),
    # Met-Lys cleavage
    "met-lyscleavage (MK)": ("low", 0.50),
    # base hydrolysis (from config.py)
    "base hydrolysis (K-X)": ("low", 0.50),
}

# ---------------------------------------------------------------------------
# Parameter dataclass
# ---------------------------------------------------------------------------

@dataclass
class PTMParameters:
    """Runtime parameters; all values are also recorded in the JSON report."""
    motif_thresholds: Dict[str, float] = field(default_factory=dict)  # per-motif relSASA gate
    buried_cutoff: float = 0.10          # relSASA below this -> not a reportable site
    class_buried: float = 0.20           # display classification boundaries
    class_exposed: float = 0.50
    borderline_band: float = 0.03        # +/- around motif threshold
    maxasa_table: str = "tien2013"       # "tien2013" | "miller1987"
    plddt_min: float = 50.0              # below -> not evaluable (requires --plddt-from-bfactor)
    plddt_low: float = 70.0              # below -> low confidence
    use_plddt: bool = False
    sasa_algorithm: str = "sr"           # "sr" (ShrakeRupley) | "lr" (LeeRichards)
    sasa_probe_radius: float = 1.4
    sasa_n_points: int = 960
    include_buried: bool = False         # list buried sites in main table
    include_antigen: bool = False        # include non-antibody chains in summary
    min_reportable_rel_sasa: float = 0.10  # synonym for buried_cutoff

    def __post_init__(self):
        # If motif_thresholds not provided, use defaults from MOTIF_SPECS
        if not self.motif_thresholds:
            self.motif_thresholds = {k: v[1] for k, v in MOTIF_SPECS.items()}

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        # Add risk levels for each motif (from MOTIF_SPECS)
        d["motif_risk_levels"] = {k: v[0] for k, v in MOTIF_SPECS.items()}
        d["motif_thresholds"] = self.motif_thresholds
        return d

    @classmethod
    def from_args(cls, args) -> "PTMParameters":
        return cls(
            motif_thresholds={},  # defaults will fill
            buried_cutoff=getattr(args, "min_relsasa", 0.10),
            borderline_band=0.03,
            maxasa_table=getattr(args, "maxasa_table", "tien2013"),
            plddt_min=getattr(args, "plddt_min", 50.0),
            plddt_low=getattr(args, "plddt_low", 70.0),
            use_plddt=getattr(args, "plddt_from_bfactor", False),
            sasa_algorithm=getattr(args, "freesasa_algorithm", "sr"),
            sasa_probe_radius=getattr(args, "freesasa_probe", 1.4),
            sasa_n_points=getattr(args, "freesasa_npoints", 960),
            include_buried=getattr(args, "include_buried", False),
            include_antigen=getattr(args, "include_antigen", False),
            min_reportable_rel_sasa=getattr(args, "min_relsasa", 0.10),
        )


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ResidueKey:
    """Immutable key for a residue in a specific chain."""
    chain: str
    res_seq: int
    i_code: str

    def __str__(self) -> str:
        if self.i_code:
            return f"{self.chain}{self.res_seq}{self.i_code}"
        return f"{self.chain}{self.res_seq}"


@dataclass
class ResidueRef:
    """Reference to a present polymer residue in the structure."""
    key: ResidueKey
    res_name: str          # 3-letter code (original)
    aa: str                # 1-letter code ('X' if unknown)
    bfactor_mean: float = 0.0
    sg_pos: Optional[Tuple[float, float, float]] = None  # for Cys
    n_pos: Optional[Tuple[float, float, float]] = None   # backbone N
    ca_pos: Optional[Tuple[float, float, float]] = None  # C-alpha (fallback)


@dataclass
class ChainModel:
    """Reconstructed chain model with numbering and annotations."""
    chain_id: str
    residues: List[ResidueRef]               # present polymer residues, in order
    sequence: str                            # X-reconstructed sequence for numbering
    seq_to_residue: Dict[int, ResidueRef]    # reconstructed seq index -> residue (present only)
    missing_resnums: List[int]               # internal missing resSeqs
    role: str                                # "antibody" | "other"
    chain_type: str                          # "H" | "L" | ""
    kabat_regions: Dict[int, str]            # struct_index (1-based) -> region
    imgt_regions: Dict[int, str]
    kabat_labels: Dict[int, str]             # struct_index -> "H31"
    imgt_labels: Dict[int, str]
    numbering_reliable: bool
    warnings: List[str] = field(default_factory=list)
    conserved_cys_indices: Set[int] = field(default_factory=set)  # struct indices


@dataclass
class PTMSite:
    chain: str
    res_seq: int
    i_code: str
    res_name: str
    aa: str
    kabat_pos: str
    imgt_pos: str
    region: str
    region_imgt: str
    primary_motif: str
    additional_motifs: List[str]
    motif_risk: str
    adjusted_risk: str
    sequence_context: str
    abs_sasa: float
    rel_sasa: float
    sasa_classification: str
    exposure_threshold: float
    is_exposed: bool
    borderline: bool
    plddt: Optional[float]
    low_confidence: bool
    priority: str
    canonical_glyco: bool
    note: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _parse_residue_label(label: str) -> Tuple[int, str]:
    """Parse FreeSASA/PDB residue label like '  52A' -> (52, 'A')."""
    s = label.strip()
    if not s:
        
        return 0, ""
    # Find the last non-digit character (iCode starts after the last digit)
    i = len(s) - 1
    while i >= 0 and s[i].isdigit():
        i -= 1
    if i == len(s) - 1:
        # Last char is non-digit -> single-char iCode
        if len(s) > 1:
            return int(s[:-1]), s[-1]
        else:
            return 0, s  # single letter
    elif i == -1:
        # All digits
        return int(s), ""
    # i is index of last non-digit (somewhere in middle)
    res_seq = int(s[:i+1])
    i_code = s[i+1:]
    return res_seq, i_code


def _resname_to_aa(resname: str) -> str:
    """Map 3-letter residue name to 1-letter AA; 'X' if unknown."""
    if resname in AA_3TO1:
        return AA_3TO1[resname]
    if resname in MODIFIED_AA:
        aa = MODIFIED_AA[resname]
        return aa if aa is not None else 'X'
    return 'X'


def _maxasa_table(name: str) -> Dict[str, float]:
    if name == "tien2013":
        return MAX_SASA_TIEN2013
    if name == "miller1987":
        return MAX_SASA_MILLER1987
    raise ValueError(f"Unknown maxASA table: {name}")


def _freesasa_params(params: PTMParameters) -> Any:
    if freesasa is None:
        raise RuntimeError("freesasa not installed")
    alg = freesasa.LeeRichards if params.sasa_algorithm == "lr" else freesasa.ShrakeRupley
    return freesasa.Parameters({
        'algorithm': alg,
        'probe-radius': params.sasa_probe_radius,
        'n-points': params.sasa_n_points,
    })


def _rel_sasa(abs_sasa: float, aa: str, params: PTMParameters) -> float:
    table = _maxasa_table(params.maxasa_table)
    aa3 = AA_1TO3.get(aa, aa)
    max_sasa = table.get(aa3, 0.0)
    return abs_sasa / max_sasa if max_sasa > 0 else 0.0


def _exposure_class(rel_sasa: float, params: PTMParameters) -> str:
    if rel_sasa < params.class_buried:
        return "BURIED"
    if rel_sasa < params.class_exposed:
        return "INTERMEDIATE"
    return "EXPOSED"


def _risk_rank(risk: str) -> int:
    return _RISK_ORDER.get(risk, -1)


def _format_seq_label(res_seq: int, i_code: str) -> str:
    return f"{res_seq}{i_code}" if i_code else str(res_seq)


# ---------------------------------------------------------------------------
# Structure parsing
# ---------------------------------------------------------------------------

def _is_water_or_ion(resname: str) -> bool:
    return resname in WATER_ION_RESNAMES


def _is_polymer_residue(resname: str) -> bool:
    """Is this residue part of the polypeptide chain (including modified AA)?"""
    return resname in AA_3TO1 or resname in MODIFIED_AA


def _is_glycan_residue(resname: str) -> bool:
    return resname in {'NAG', 'NAB', 'NDG', 'MAN', 'BMA', 'FUC', 'SIA', 'GAL', 'GLC', 'SGC'}


def parse_structure(pdb_path: str) -> Tuple[Dict[str, ChainModel], Dict[str, List[dict]], List[Tuple[str, str, str, str]], int]:
    """
    Parse a PDB file into ChainModels per antibody/antigen chain.
    Returns (chain_models, hetero_atoms, ssbonds, n_models).
    """
    if freesasa is None:
        raise RuntimeError("freesasa not installed")

    # First pass: read all records
    chains_atoms: Dict[str, Dict[ResidueKey, List[dict]]] = defaultdict(lambda: defaultdict(list))
    ssbonds: List[Tuple[str, str, str, str]] = []
    models_seen = 0
    current_model = 0

    with open(pdb_path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "MODEL":
                models_seen += 1
                current_model += 1
                if current_model > 1:
                    continue  # we only parse model 1
            elif rec == "ENDMDL":
                current_model -= 1
                continue
            if current_model > 1:
                continue

            if rec == "SSBOND":
                # SSBOND columns: 12-14 res1, 15 chain1, 17-19 seq1, 20 iCode1, 21-23 res2, 24 chain2, 26-28 seq2, 29 iCode2
                r1 = line[11:14].strip()
                c1 = line[15].strip()
                s1 = int(line[16:20].strip())
                i1 = line[20].strip()
                r2 = line[21:24].strip()
                c2 = line[24].strip()
                s2 = int(line[25:29].strip())
                i2 = line[29].strip()
                ssbonds.append((c1, _format_seq_label(s1, i1), c2, _format_seq_label(s2, i2)))
            elif line.startswith(("ATOM", "HETATM")):
                chain = line[21].strip() or ' '
                res_seq = int(line[22:26].strip())
                i_code = line[26].strip()
                resname = line[17:20].strip()
                atom_name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                bfactor = float(line[60:66]) if len(line) >= 66 else 0.0
                key = ResidueKey(chain, res_seq, i_code)
                chains_atoms[chain][key].append({
                    'name': atom_name,
                    'resname': resname,
                    'x': x, 'y': y, 'z': z,
                    'bfactor': bfactor,
                })

    n_models = max(1, models_seen)

    # Build ChainModels per chain
    chain_models: Dict[str, ChainModel] = {}
    hetero_atoms: Dict[str, List[dict]] = defaultdict(list)  # non-water HETATM for glycan/SS checks

    for chain, residues_atoms in chains_atoms.items():
        # Separate polymer vs non-polymer residues
        polymer_keys = [k for k, atoms in residues_atoms.items() if _is_polymer_residue(atoms[0]['resname'])]
        if not polymer_keys:
            # No polymer residues - likely water/ion only or non-protein
            for k, atoms in residues_atoms.items():
                for a in atoms:
                    if not _is_water_or_ion(a['resname']):
                        hetero_atoms[chain].append(a)
            continue

        # Sort polymer residues by (res_seq, i_code)
        polymer_keys.sort(key=lambda k: (k.res_seq, k.i_code))

        # Build ResidueRefs for polymer residues
        residues: List[ResidueRef] = []
        bf_sum = defaultdict(float)
        bf_count = defaultdict(int)
        for k in polymer_keys:
            atoms_list = residues_atoms[k]
            # Mean B-factor
            bf_sum[k] = sum(a['bfactor'] for a in atoms_list)
            bf_count[k] = len(atoms_list)
            mean_bf = bf_sum[k] / bf_count[k] if bf_count[k] else 0.0

            resname = atoms_list[0]['resname']
            aa = _resname_to_aa(resname)

            # Special atoms for Cys/N
            sg_pos = None
            n_pos = None
            ca_pos = None
            for a in atoms_list:
                if a['name'] in ('SG',):
                    sg_pos = (a['x'], a['y'], a['z'])
                if a['name'] in ('N',):
                    n_pos = (a['x'], a['y'], a['z'])
                if a['name'] in ('CA',):
                    ca_pos = (a['x'], a['y'], a['z'])

            residues.append(ResidueRef(
                key=k,
                res_name=resname,
                aa=aa,
                bfactor_mean=mean_bf,
                sg_pos=sg_pos,
                n_pos=n_pos,
                ca_pos=ca_pos,
            ))

        # Collect hetero atoms (non-water, non-polymer) for this chain
        for k, atoms in residues_atoms.items():
            if k in polymer_keys:
                continue
            rn = atoms[0]['resname']
            if _is_water_or_ion(rn):
                continue
            for a in atoms:
                hetero_atoms[chain].append(a)

        # Detect internal gaps in resSeq (accounting for iCodes)
        missing = []
        prev_res_seq = None
        for r in residues:
            if prev_res_seq is not None:
                gap = r.key.res_seq - prev_res_seq - 1
                if gap > 0:
                    missing.extend(range(prev_res_seq + 1, r.key.res_seq))
            prev_res_seq = r.key.res_seq

        # Build X-reconstructed sequence and mapping
        seq_chars = []
        seq_to_residue = {}
        prev_res = None
        for idx, r in enumerate(residues):
            if prev_res is not None:
                gap = r.key.res_seq - prev_res.key.res_seq - 1
                if gap > 0:
                    for _ in range(gap):
                        seq_chars.append('X')
            seq_chars.append(r.aa)
            seq_to_residue[len(seq_chars) - 1] = r
            prev_res = r
        sequence = "".join(seq_chars)

        # Number the X-reconstructed sequence
        kabat_regions = {}
        imgt_regions = {}
        kabat_labels = {}
        imgt_labels = {}
        chain_type = ""
        numbering_reliable = True
        warnings = []

        # Kabat numbering
        kabat_numbered = None
        try:
            kabat_numbered = number_heavy(sequence) if sequence else None
        except Exception:
            pass
        if kabat_numbered is None:
            try:
                kabat_numbered = number_light(sequence)
            except Exception:
                kabat_numbered = None

        if kabat_numbered is not None:
            chain_type = kabat_numbered.chain_type
            for r in kabat_numbered.residues:
                si = r.index + 1  # 1-based structural index among PRESENT residues? No: r.index is index in X-reconstructed sequence.
                # Need to map X-reconstructed index -> present residue
                # If r.index maps to an X placeholder, it won't be in seq_to_residue
                if r.index in seq_to_residue:
                    ref = seq_to_residue[r.index]
                    struct_idx = residues.index(ref) + 1  # 1-based over present residues
                    kabat_labels[struct_idx] = r.pos
                    kabat_regions[struct_idx] = r.region
        else:
            warnings.append("Kabat numbering failed")

        # IMGT numbering (requires ANARCI)
        imgt_numbered = None
        try:
            imgt_numbered = number_with_anarci_imgt(sequence, chain_type if chain_type in ("H", "L") else "H")
        except Exception:
            pass
        if imgt_numbered is None:
            warnings.append("IMGT numbering failed")
        else:
            for r in imgt_numbered.residues:
                if r.index in seq_to_residue:
                    ref = seq_to_residue[r.index]
                    struct_idx = residues.index(ref) + 1
                    imgt_labels[struct_idx] = r.pos
                    imgt_regions[struct_idx] = r.region

        # Determine chain role
        role = "other"
        if chain_type in ("H", "L") and len(residues) >= 80:
            # Verify conserved Cys anchors
            conserved_expected = CONSERVED_CYS.get(chain_type, set())
            cys_at_anchor = 0
            for pos in conserved_expected:
                label = f"{'H' if chain_type == 'H' else 'L'}{pos}"
                # Check if any residue has this kabat label and is Cys
                for si, lbl in kabat_labels.items():
                    if lbl == label and si <= len(residues) and residues[si - 1].aa == 'C':
                        cys_at_anchor += 1
                        break
            if cys_at_anchor >= 1:  # at least one anchor Cys confirmed
                role = "antibody"
            else:
                warnings.append(f"Conserved Cys anchors not found at expected Kabat labels ({conserved_expected}); numbering may be shifted")
                numbering_reliable = False
        else:
            warnings.append("Chain too short or not recognized as antibody V domain")

        # Identify conserved Cys indices (structural indices)
        conserved_cys_indices = set()
        for si, lbl in kabat_labels.items():
            if si <= len(residues) and residues[si - 1].aa == 'C':
                m = re.match(r'^[HL](\d+)', lbl)
                if m and int(m.group(1)) in CONSERVED_CYS.get(chain_type, set()):
                    conserved_cys_indices.add(si)

        chain_models[chain] = ChainModel(
            chain_id=chain,
            residues=residues,
            sequence=sequence,
            seq_to_residue=seq_to_residue,
            missing_resnums=missing,
            role=role,
            chain_type=chain_type,
            kabat_regions=kabat_regions,
            imgt_regions=imgt_regions,
            kabat_labels=kabat_labels,
            imgt_labels=imgt_labels,
            numbering_reliable=numbering_reliable and (len(warnings) == 0 or all("numbering" not in w for w in warnings)),
            warnings=warnings,
            conserved_cys_indices=conserved_cys_indices,
        )

    return chain_models, dict(hetero_atoms), ssbonds, n_models


# ---------------------------------------------------------------------------
# SASA calculation (on filtered structure)
# ---------------------------------------------------------------------------

def _write_filtered_pdb(pdb_path: str, tmp_path: str) -> int:
    """Write a filtered PDB (model 1 only, altLoc A/blank, no water/ions) and return model count."""
    model_count = 0
    in_model = True  # True for model 1 (or if no MODEL records)
    model_seen = False
    with open(pdb_path) as fin, open(tmp_path, 'w') as fout:
        for line in fin:
            rec = line[:6].strip()
            if rec == "MODEL":
                model_seen = True
                model_count += 1
                in_model = (model_count == 1)
                if in_model:
                    fout.write(line)
                continue
            if rec == "ENDMDL":
                in_model = False
                continue
            if not in_model:
                continue
            if line.startswith(("ATOM", "HETATM")):
                alt = line[16]
                if alt not in (' ', 'A'):
                    continue
                resname = line[17:20].strip()
                if _is_water_or_ion(resname):
                    continue
                fout.write(line)
            else:
                fout.write(line)
    return max(1, model_count)


def calculate_sasa(pdb_path: str, params: PTMParameters) -> Tuple[Dict[str, Dict[ResidueKey, float]], int]:
    """Compute per-residue SASA on a filtered structure."""
    if freesasa is None:
        raise RuntimeError("freesasa not installed")

    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        n_models = _write_filtered_pdb(pdb_path, tmp_path)
        # Debug: read back the filtered file
        with open(tmp_path) as f:
            lines = f.readlines()
        
        structure = freesasa.Structure(tmp_path)
        
        for i in range(structure.nAtoms()):
            chain = structure.chainLabel(i) or ' '
            label = str(structure.residueNumber(i)).strip()
            if not label:
                print(f"DEBUG: EMPTY label for atom {i} chain={repr(chain)} resname={structure.residueName(i)}", file=sys.stderr)
                break
        fs_params = _freesasa_params(params)
        result = freesasa.calc(structure, fs_params)

        sasa_by_res: Dict[str, Dict[ResidueKey, float]] = defaultdict(lambda: defaultdict(float))
        n_atoms = structure.nAtoms()
        for i in range(n_atoms):
            chain = structure.chainLabel(i) or ' '
            label = str(structure.residueNumber(i)).strip()
            if not label:
                
                continue
            res_seq, i_code = _parse_residue_label(label)
            
            area = result.atomArea(i)
            key = ResidueKey(chain, res_seq, i_code)
            sasa_by_res[chain][key] += area
        return dict(sasa_by_res), n_models
    finally:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# Disulfide pairing and glycan occupancy
# ---------------------------------------------------------------------------

def _find_disulfide_partners(chain_models: Dict[str, ChainModel]) -> Dict[Tuple[str, int], Set[Tuple[str, int]]]:
    """
    Find all disulfide partners for each Cys (by structural index within its chain).
    Returns dict: (chain, struct_idx) -> set of (partner_chain, partner_struct_idx).
    """
    partners: Dict[Tuple[str, int], Set[Tuple[str, int]]] = defaultdict(set)
    # Collect all Cys with SG position
    cys_list: List[Tuple[str, int, Tuple[float, float, float]]] = []
    for chain, model in chain_models.items():
        for si, ref in enumerate(model.residues, start=1):
            if ref.aa == 'C' and ref.sg_pos is not None:
                cys_list.append((chain, si, ref.sg_pos))
    # Pairwise distance check
    for i, (c1, si1, pos1) in enumerate(cys_list):
        for c2, si2, pos2 in cys_list[i + 1:]:
            dx = pos1[0] - pos2[0]
            dy = pos1[1] - pos2[1]
            dz = pos1[2] - pos2[2]
            dist = (dx * dx + dy * dy + dz * dz) ** 0.5
            if dist <= 2.5:  # disulfide bond cutoff
                partners[(c1, si1)].add((c2, si2))
                partners[(c2, si2)].add((c1, si1))
    return partners


def _check_glycan_occupancy(chain: str, si: int, chain_models: Dict[str, ChainModel],
                             hetero_atoms: Dict[str, List[dict]]) -> bool:
    """Check if a sequon N has a glycan attached within 2.5 A."""
    if si <= 0 or si > len(chain_models[chain].residues):
        return False
    ref = chain_models[chain].residues[si - 1]
    if ref.n_pos is None:
        return False
    nx, ny, nz = ref.n_pos
    # Check hetero atoms on the same chain or any chain (glycans can be on different chain label)
    for h_chain, atoms in hetero_atoms.items():
        for a in atoms:
            if not _is_glycan_residue(a['resname']):
                continue
            dx = a['x'] - nx
            dy = a['y'] - ny
            dz = a['z'] - nz
            if (dx * dx + dy * dy + dz * dz) ** 0.5 <= 2.5:
                return True
    return False


# ---------------------------------------------------------------------------
# Site detection
# ---------------------------------------------------------------------------

def _detect_motif_hits(sequence: str, motif_patterns: Dict[str, Tuple[str, float]]) -> Dict[int, List[Tuple[str, str, float]]]:
    """
    Scan sequence for all motifs.
    Returns: struct_index (1-based in reconstructed sequence) -> list of (motif_name, risk, threshold).
    """
    hits: Dict[int, List[Tuple[str, str, float]]] = defaultdict(list)
    for motif_name, (risk, threshold) in motif_patterns.items():
        pattern = RISK_MOTIFS.get(motif_name)
        if pattern is None:
            continue
        for m in re.finditer(pattern, sequence):
            start = m.start()
            struct_idx = start + 1  # 1-based in reconstructed sequence
            hits[struct_idx].append((motif_name, risk, threshold))
    return hits


def detect_sites(
    chain_models: Dict[str, ChainModel],
    sasa: Dict[str, Dict[ResidueKey, float]],
    params: PTMParameters,
    hetero_atoms: Dict[str, List[dict]],
    partners: Dict[Tuple[str, int], Set[Tuple[str, int]]],
) -> Tuple[List[PTMSite], List[PTMSite]]:
    """
    Detect PTM sites across all chains.
    Returns (reportable_sites, buried_sites).
    """
    motif_specs = MOTIF_SPECS
    reportable = []
    buried = []

    for chain, model in chain_models.items():
        if model.role == "other" and not params.include_antigen:
            continue
        if chain not in sasa:
            continue

        # Scan motifs on X-reconstructed sequence
        hits = _detect_motif_hits(model.sequence, motif_specs)

        # Disulfide partners for this chain
        chain_partners = {si: prs for (c, si), prs in partners.items() if c == chain}

        for seq_idx, motif_list in hits.items():
            # Map reconstructed index -> present residue
            if seq_idx not in model.seq_to_residue:
                # Hit falls on an X placeholder (missing residue) -> not evaluable
                continue
            ref = model.seq_to_residue[seq_idx]
            si = model.residues.index(ref) + 1  # structural index (1-based)

            # Skip conserved Cys
            if si in model.conserved_cys_indices:
                continue

            # Get SASA
            key = ref.key
            abs_sasa = sasa.get(chain, {}).get(key, 0.0)
            rel_sasa = _rel_sasa(abs_sasa, ref.aa, params)

            # pLDDT
            plddt = ref.bfactor_mean if params.use_plddt else None
            low_conf = False
            not_evaluable = False
            if params.use_plddt and plddt is not None:
                if plddt < params.plddt_min:
                    not_evaluable = True
                elif plddt < params.plddt_low:
                    low_conf = True

            # Determine exposure per motif, pick primary
            best = None
            best_rank = -1
            all_motifs = []
            for motif_name, risk, threshold in motif_list:
                is_exposed = rel_sasa >= threshold
                adjusted_risk = risk if is_exposed else {"high": "medium", "medium": "low", "low": "low"}.get(risk, "low")
                all_motifs.append((motif_name, adjusted_risk, is_exposed, threshold))
                rank = _risk_rank(adjusted_risk)
                if rank > best_rank:
                    best_rank = rank
                    best = (motif_name, adjusted_risk, is_exposed, threshold)

            if best is None:
                continue

            primary_motif, adjusted_risk, is_exposed, primary_threshold = best
            additional = [m for m, ar, ie, th in all_motifs if m != primary_motif]

            # Cys pairing: if primary is oxidation (C) or unpaired Cys, check pairing
            if primary_motif in ("oxidation (C)", "unpaired Cys"):
                has_partner = (si in chain_partners) and len(chain_partners[si]) > 0
                if has_partner:
                    # Not unpaired; keep oxidation (C) if exposed, else skip unpaired
                    if primary_motif == "unpaired Cys":
                        continue  # skip this hit entirely; a paired Cys is not an unpaired risk
                    # oxidation (C) stays
                else:
                    # No partner -> if primary was oxidation (C), promote unpaired Cys
                    if primary_motif == "oxidation (C)":
                        primary_motif = "unpaired Cys"
                        adjusted_risk = "medium" if is_exposed else "low"
                        # Note: additional will still include oxidation (C) if it was in the list
                        # (we'll handle that below)

            # Glycan occupancy check for N-glycan motifs
            canonical_glyco = False
            if primary_motif == "N-glycan (NxS/T)":
                if _check_glycan_occupancy(chain, si, chain_models, {}):  # hetero_atoms passed in caller
                    canonical_glyco = True
                    adjusted_risk = "low"  # canonical site is informational, not a de novo risk
                # Also check if the N is in the canonical VH loop (Kabat H28-H33)
                kabat_lbl = model.kabat_labels.get(seq_idx, "")
                m = re.match(r'H(\d+)', kabat_lbl)
                if m and 28 <= int(m.group(1)) <= 33:
                    canonical_glyco = True
                    adjusted_risk = "low"

            # Buried filter
            sasa_class = _exposure_class(rel_sasa, params)
            is_buried = rel_sasa < params.buried_cutoff

            # Borderline
            borderline = abs(rel_sasa - primary_threshold) <= params.borderline_band

            # Sequence context
            ctx_start = max(0, seq_idx - 4)
            ctx_end = min(len(model.sequence), seq_idx + 3)
            context = model.sequence[ctx_start:ctx_end]

            # Regions
            region = model.kabat_regions.get(si, "-")
            region_imgt = model.imgt_regions.get(si, "-")
            kabat_pos = model.kabat_labels.get(si, "-")
            imgt_pos = model.imgt_labels.get(si, "-")

            # Priority
            if canonical_glyco:
                priority = "canonical-glyco-site"
            elif is_exposed and (region in ("CDR1", "CDR2", "CDR3") or region_imgt in ("CDR1", "CDR2", "CDR3")):
                priority = "CDR-hotspot"
            elif is_exposed:
                priority = "framework-hotspot"
            else:
                priority = "framework-buried"

            # Note
            note_parts = []
            if not_evaluable:
                note_parts.append(f"pLDDT {plddt:.1f} below minimum ({params.plddt_min}); not evaluable")
            if low_conf:
                note_parts.append(f"low confidence (pLDDT {plddt:.1f})")
            if borderline:
                note_parts.append(f"borderline (|relSASA-threshold|={abs(rel_sasa - primary_threshold):.3f} <= {params.borderline_band})")
            if not model.numbering_reliable:
                note_parts.append("numbering unreliable (missing residues or anchor check failed)")

            site = PTMSite(
                chain=chain,
                res_seq=ref.key.res_seq,
                i_code=ref.key.i_code,
                res_name=ref.res_name,
                aa=ref.aa,
                kabat_pos=kabat_pos,
                imgt_pos=imgt_pos,
                region=region,
                region_imgt=region_imgt,
                primary_motif=primary_motif,
                additional_motifs=additional,
                motif_risk=MOTIF_SPECS.get(primary_motif, ("medium", 0.0))[0],
                adjusted_risk=adjusted_risk,
                sequence_context=context,
                abs_sasa=abs_sasa,
                rel_sasa=rel_sasa,
                sasa_classification=sasa_class,
                exposure_threshold=primary_threshold,
                is_exposed=is_exposed,
                borderline=borderline,
                plddt=plddt,
                low_confidence=low_conf,
                priority=priority,
                canonical_glyco=canonical_glyco,
                note="; ".join(note_parts),
            )

            if not_evaluable:
                continue  # skip from both lists
            if is_buried:
                if params.include_buried:
                    reportable.append(site)
                else:
                    buried.append(site)
            else:
                reportable.append(site)

    return reportable, buried


# ---------------------------------------------------------------------------
# Main analysis entry
# ---------------------------------------------------------------------------

def analyze_structure(pdb_path: str, params: Optional[PTMParameters] = None) -> Dict[str, Any]:
    """Full analysis of a single structure."""
    if params is None:
        params = PTMParameters()

    # Parse structure
    chain_models, hetero_atoms, ssbonds, n_models = parse_structure(pdb_path)

    # SASA
    sasa, n_models = calculate_sasa(pdb_path, params)

    # Disulfide partners
    partners = _find_disulfide_partners(chain_models)

    # Detect sites
    all_sites, buried_sites = detect_sites(chain_models, sasa, params, hetero_atoms, partners)

    # File hash
    with open(pdb_path, 'rb') as f:
        file_sha256 = hashlib.sha256(f.read()).hexdigest()

    # Build result
    result = {
        "pdb_file": pdb_path,
        "file_sha256": file_sha256,
        "analysis_date": datetime.now().isoformat(),
        "tool_version": "2.0",
        "parameters": params.to_dict(),
        "structure": {
            "n_models": n_models,
            "model_used": 1,
            "chains": {},
        },
        "chains": {},
        "summary": {
            "total_ptm_sites": len(all_sites),
            "high_risk": 0,
            "medium_risk": 0,
            "low_risk": 0,
            "not_evaluable": 0,
            "by_priority": defaultdict(int),
        },
    }

    for chain, model in chain_models.items():
        chain_sites = [s for s in all_sites if s.chain == chain]
        chain_buried = [s for s in buried_sites if s.chain == chain]

        result["structure"]["chains"][chain] = {
            "role": model.role,
            "chain_type": model.chain_type,
            "n_residues": len(model.residues),
            "missing_residues": model.missing_resnums,
            "numbering_reliable": model.numbering_reliable,
            "warnings": model.warnings,
        }

        chain_sites_data = [s.to_dict() for s in chain_sites]
        result["chains"][chain] = {
            "sequence": model.sequence,
            "ptm_sites": chain_sites_data,
            "buried_sites": [s.to_dict() for s in chain_buried] if params.include_buried else [],
            "residue_count": len(model.residues),
        }

        for s in chain_sites:
            result["summary"][s.adjusted_risk + "_risk"] += 1
            result["summary"]["by_priority"][s.priority] += 1

    result["summary"]["total_ptm_sites"] = len(all_sites)
    return result


# For backward compatibility with existing tests
def get_chain_regions(sequence: str, chain: str) -> Tuple[Dict[int, str], Dict[int, str], str, Set[int]]:
    """Legacy API used by tests/test_tool_numbering.py."""
    # This mimics the old behavior: number the given sequence (no X reconstruction)
    # and return kabat/imgt regions by 1-based index, chain type, and conserved cys indices.
    kabat_regions = {}
    imgt_regions = {}
    chain_type = ""
    conserved = set()

    try:
        kb = number_heavy(sequence) if chain == "H" else number_light(sequence)
        chain_type = kb.chain_type
        for r in kb.residues:
            kabat_regions[r.index + 1] = r.region
            m = re.match(r'^[HL](\d+)', r.pos)
            if m and int(m.group(1)) in CONSERVED_CYS.get(chain_type, set()):
                conserved.add(r.index + 1)
    except Exception:
        pass

    try:
        imgt = number_with_anarci_imgt(sequence, chain_type if chain_type in ("H", "L") else "H")
        if imgt:
            for r in imgt.residues:
                imgt_regions[r.index + 1] = r.region
    except Exception:
        pass

    return kabat_regions, imgt_regions, chain_type, conserved


from datetime import datetime