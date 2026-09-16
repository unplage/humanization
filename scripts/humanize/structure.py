"""AlphaFold3 adapter and structure-based scoring.

* Builds AF3 request JSONs (Fv or Fv + antigen) and submits them to a local
  AF3 binary or an AF3-compatible API endpoint.
* Parses PDB output (AF3 `--output_format=pdb`) with a tiny dependency-free
  PDB reader and computes the structure hints consumed by back-mutation
  scoring: per-residue buriedness (SASA proxy), CDR contacts and
  antigen contacts, plus CDR-loop RMSD vs a reference model.
* `--mock` mode produces no structure data (all hints None); the pipeline
  then relies on sequence/structure-free rules.

Requires biopython only for the full SASA computation; without it, a
cheap heavy-atom-neighbourhood proxy is used.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .backmut import StructureHints

AA_ATOMS = {
    "A": 5, "R": 11, "N": 8, "D": 8, "C": 6, "Q": 9, "E": 9, "G": 4, "H": 10,
    "I": 8, "L": 8, "K": 9, "M": 8, "F": 11, "P": 7, "S": 6, "T": 7, "W": 14,
    "Y": 12, "V": 7,
}


@dataclass
class PDBAtom:
    name: str
    resname: str
    chain: str
    resseq: int
    x: float
    y: float
    z: float
    plddt: float = 0.0  # AF3 confidence score (B-factor column)


@dataclass
class PDBModel:
    atoms: List[PDBAtom] = field(default_factory=list)

    def residues(self) -> Dict[Tuple[str, int], Dict[str, Tuple[float, float, float]]]:
        out = {}
        for a in self.atoms:
            out.setdefault((a.chain, a.resseq), {})[a.name] = (a.x, a.y, a.z)
        return out

    def heavy_atoms(self) -> List[Tuple[Tuple[str, int], str, Tuple[float, float, float]]]:
        return [
            ((a.chain, a.resseq), a.name, (a.x, a.y, a.z))
            for a in self.atoms if a.name[0] != "H"
        ]


def parse_pdb(path: str) -> Optional[PDBModel]:
    model = PDBModel()
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Hydrogen filtering is deferred to heavy_atoms() which uses
                # the PDB element column (line[76:78]) for robustness; the
                # atom name prefix check there handles most cases already.
                try:
                    # AF3 stores pLDDT in the B-factor column (60:66)
                    plddt = float(line[60:66]) if len(line) > 66 else 0.0
                    model.atoms.append(PDBAtom(
                        name=line[12:16].strip(),
                        resname=line[17:20].strip(),
                        chain=line[21],
                        resseq=int(line[22:26].strip()),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        plddt=plddt,
                    ))
                except (ValueError, IndexError):
                    continue
    return model if model.atoms else None


def parse_cif_atom_site(path: str) -> Optional[PDBModel]:
    """Parse the _atom_site loop of an mmCIF file (AF3 fallback)."""
    model = PDBModel()
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
                # End of loop: a blank line, a comment, or a new loop_
                # keyword all terminate the current _atom_site loop.
                in_loop = False
                cols = []
                continue
            if in_loop and not s.startswith(("_", "#")):
                parts = s.split()
                if len(parts) >= len(cols):
                    d = dict(zip(cols, parts))
                    if d.get("group_PDB") in ("ATOM", "HETATM"):
                        try:
                            b_iso = d.get("B_iso_or_equiv", "0")
                            plddt = float(b_iso) if b_iso not in (".", "?") else 0.0
                            model.atoms.append(PDBAtom(
                                name=d.get("auth_atom_id", d.get("label_atom_id", "?")),
                                resname=d.get("auth_comp_id", d.get("label_comp_id", "?")),
                                chain=d.get("auth_asym_id", d.get("label_asym_id", "?")),
                                resseq=int(d.get("auth_seq_id", d.get("label_seq_id", "0"))),
                                x=float(d["Cartn_x"]),
                                y=float(d["Cartn_y"]),
                                z=float(d["Cartn_z"]),
                                plddt=plddt,
                            ))
                        except (ValueError, KeyError):
                            continue
    return model if model.atoms else None


def load_model(path: str) -> Optional[PDBModel]:
    if path.endswith(".pdb"):
        m = parse_pdb(path)
        if m:
            return m
    if path.endswith(".cif"):
        return parse_cif_atom_site(path)
    return None


# ---------------------------------------------------------------------------
# chain assignment (PDB chain -> input chain)
# ---------------------------------------------------------------------------

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def longest_common_run(a: str, b: str) -> int:
    """Length of the longest contiguous residue run shared by two sequences
    (classic DP; O(len_a * len_b), fine for ~120-aa V domains). Robust to
    terminal tags, truncations and numbering offsets in structures."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for ca in a:
        cur = [0] * (len(b) + 1)
        for j, cb in enumerate(b, 1):
            if ca == cb:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best = cur[j]
        prev = cur
    return best


def match_pdb_chain(pdb_chains: Dict[str, list], donor_seq: str,
                    min_run: int = 20) -> Optional[str]:
    """Pick the PDB chain whose residue sequence shares the longest exact
    run with the donor sequence. First-residue matching is unreliable:
    unrelated VH and VL chains frequently share the same N-terminal residue,
    and two heavy domains even share the QVQL.. prefix."""
    best_cid, best_run = None, 0
    for cid, atoms in pdb_chains.items():
        cas = sorted((a.resseq, a.resname) for a in atoms if a.name == "CA")
        s = "".join(_AA3TO1.get(rn, "X") for _, rn in cas)
        run = longest_common_run(donor_seq, s)
        if run > best_run:
            best_cid, best_run = cid, run
    return best_cid if best_run >= min_run else None


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _try_freesasa_buriedness(
    pdb_path: str,
    chain_label: str,
    resseq_to_pos: Dict[int, str],
) -> Optional[Dict[str, Tuple[bool, float]]]:
    """Try to compute buriedness using FreeSASA. Returns None if unavailable.
    
    Returns: {pos: (is_buried, rel_sasa)} - includes relSASA for confidence interval.
    """
    try:
        import freesasa
        from collections import defaultdict
        
        structure = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure)
        
        # Reference SASA values (Tien et al. 2013)
        ref_sasa = {
            'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
            'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0,
            'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
            'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
            'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
        }
        
        # Calculate per-residue SASA
        sasa_by_residue = defaultdict(float)
        residue_names = {}
        for i in range(structure.nAtoms()):
            chain = structure.chainLabel(i)
            res_num = structure.residueNumber(i)
            if isinstance(res_num, str):
                res_num = int(res_num.strip())
            area = result.atomArea(i)
            sasa_by_residue[(chain, res_num)] += area
            residue_names[(chain, res_num)] = structure.residueName(i)
        
        # Convert to buried/exposed with relSASA
        buried = {}
        for resseq, pos in resseq_to_pos.items():
            key = (chain_label, resseq)
            if key in sasa_by_residue:
                abs_sasa = sasa_by_residue[key]
                resname = residue_names.get(key)
                if resname and resname in ref_sasa:
                    rel_sasa = abs_sasa / ref_sasa[resname]
                    buried[pos] = (rel_sasa < 0.20, rel_sasa)
        
        return buried if buried else None
    except ImportError:
        return None


def compute_hints(
    model: PDBModel,
    chain_label: str,
    position_map: Dict[str, int],
    cdr_positions: Dict[str, int],
    antigen_chains: Optional[List[str]] = None,
    pdb_path: Optional[str] = None,
    literature_positions: Optional[set] = None,
    use_plddt: bool = True,
) -> StructureHints:
    """Compute per-position hints for a chain in the model.

    chain_label:  PDB chain id of the target chain (e.g. "H").
    position_map: {Kabat pos: residue number} for ALL residues of the chain.
    cdr_positions:{Kabat pos: residue number} for CDR residues only.
    pdb_path:     Optional path to PDB file for FreeSASA calculation.
    literature_positions: Optional set of positions with literature evidence
                         (vernier/canonical/interface). If provided, structural
                         evidence is only applied to these positions.
    use_plddt:    Interpret the B-factor column as an AF3 pLDDT score. Must be
                  False for experimental PDB/CIF files, whose B-factor is a
                  crystallographic temperature factor, not a confidence score.
    """
    heavy = model.heavy_atoms()
    by_res: Dict[Tuple[str, int], List] = {}
    for (ch, resseq), name, xyz in heavy:
        by_res.setdefault((ch, resseq), []).append((name, xyz))

    # which residues are CDR residues
    cdr_resseqs = {v for v in cdr_positions.values() if isinstance(v, int)}
    # chain-aware filtering is essential: AF3 numbers every chain from 1,
    # so resseq-only filtering would mix the target chain's own atoms into
    # the CDR/antigen pools and corrupt every contact hint.
    cdr_atoms = [a for a in heavy
                 if a[0][0] == chain_label and a[0][1] in cdr_resseqs]
    ag_atoms = ([a for a in heavy if a[0][0] in antigen_chains]
                if antigen_chains else [])

    # map resseq -> kabat pos
    resseq_to_pos = {v: k for k, v in position_map.items()}

    # Extract pLDDT from model (AF3 B-factor column). Skipped for experimental
    # structures where the B-factor is a temperature factor, not pLDDT.
    plddt_by_pos: Dict[str, float] = {}
    if use_plddt:
        for a in model.atoms:
            if a.chain == chain_label and a.plddt > 0:
                pos = resseq_to_pos.get(a.resseq)
                if pos and pos not in plddt_by_pos:
                    plddt_by_pos[pos] = a.plddt

    # Try FreeSASA for buriedness calculation
    # Returns: {pos: (is_buried, rel_sasa)} or None
    freesasa_data = None
    if pdb_path:
        freesasa_data = _try_freesasa_buriedness(pdb_path, chain_label, resseq_to_pos)
    
    # Build buried dict with confidence interval logic
    buried: Dict[str, Optional[bool]] = {}
    rel_sasa: Dict[str, float] = {}
    if freesasa_data:
        for pos, (is_buried, sasa_val) in freesasa_data.items():
            rel_sasa[pos] = sasa_val
            # 方案3: relSASA 置信区间 0.15-0.25 → 标记为 uncertain (None)
            if 0.15 <= sasa_val <= 0.25:
                buried[pos] = None  # uncertain, treated as exposed
            else:
                buried[pos] = is_buried
    else:
        # Fallback to contact-based heuristic
        for (ch, resseq), atoms in by_res.items():
            if ch != chain_label:
                continue
            pos = resseq_to_pos.get(resseq)
            if pos is None:
                continue
            contacts = 0
            for (oc, oseq), oa in by_res.items():
                if (oc, oseq) == (ch, resseq):
                    continue
                for _, xyz in oa:
                    for _, axyz in atoms[: min(3, len(atoms))]:
                        if _dist(axyz, xyz) < 6.0:
                            contacts += 1
                            break
            buried[pos] = contacts >= 45

    cdr_contact: Dict[str, bool] = {}
    cdr_contact_sc: Dict[str, bool] = {}
    ag_contact: Dict[str, bool] = {}
    cdr_partners: Dict[str, List[str]] = {}
    cdr_partners_sc: Dict[str, List[str]] = {}
    for (ch, resseq), atoms in by_res.items():
        if ch != chain_label:
            continue
        pos = resseq_to_pos.get(resseq)
        if pos is None:
            continue
        # CDR contact: any heavy atom within 4.5 A of any CDR heavy atom.
        # All atoms of the residue are tested (not just the first one), so a
        # side-chain-mediated contact is not missed. We additionally record
        # whether the contact involves the framework SIDE CHAIN: contacts
        # mediated only by backbone atoms (N/CA/C/O) are fixed beta-sheet
        # geometry and are unchanged by a side-chain substitution, so they
        # must not be treated as functional CDR contacts during Step 3.
        backbone = {"N", "CA", "C", "O"}
        sidechain = [(n1, a1) for (n1, a1) in atoms if n1 not in backbone]
        partners = set()
        partners_sc = set()
        for (r2, _n2, a2) in cdr_atoms:
            ppos = resseq_to_pos.get(r2[1])
            if not ppos:
                continue
            if any(_dist(a2, a1) < 4.5 for (_n1, a1) in atoms):
                partners.add(ppos)
            if sidechain and any(_dist(a2, a1) < 4.5 for (_n1, a1) in sidechain):
                partners_sc.add(ppos)
        cdr_contact[pos] = bool(partners)
        if partners:
            cdr_partners[pos] = sorted(partners)
        cdr_contact_sc[pos] = bool(partners_sc)
        if partners_sc:
            cdr_partners_sc[pos] = sorted(partners_sc)
        if ag_atoms:
            ag_contact[pos] = any(
                _dist(a2, a1) < 4.5
                for (r2, _n2, a2) in ag_atoms for (_n1, a1) in atoms
            )

    # 方案5: 文献交叉验证 - 仅限文献位置集合
    # 如果提供了文献位置集合，仅在这些位置应用结构证据
    if literature_positions is not None:
        for pos in list(buried.keys()):
            if pos not in literature_positions:
                del buried[pos]
        for pos in list(cdr_contact.keys()):
            if pos not in literature_positions:
                del cdr_contact[pos]
        for pos in list(cdr_contact_sc.keys()):
            if pos not in literature_positions:
                del cdr_contact_sc[pos]
        for pos in list(ag_contact.keys()):
            if pos not in literature_positions:
                del ag_contact[pos]

    return StructureHints({
        "buried": buried,
        "cdr_contact": cdr_contact,
        "cdr_contact_sc": cdr_contact_sc,
        "antigen_contact": ag_contact,
        "cdr_partners": cdr_partners,
        "cdr_partners_sc": cdr_partners_sc,
        "plddt": plddt_by_pos,
        "rel_sasa": rel_sasa,
    })


def compute_multi_model_consensus(
    pdb_paths: List[str],
    chain_label: str,
    position_map: Dict[str, int],
    cdr_positions: Dict[str, int],
    antigen_chains: Optional[List[str]] = None,
    min_consensus: int = 3,
    use_plddt: bool = True,
) -> StructureHints:
    """Compute structure hints from multiple models (e.g., AF3 rank_1-N).

    Buried / CDR-contact / antigen-contact calls require at least
    ``min_consensus`` models to agree (auto-relaxed when fewer models are
    available). CDR partner lists are unioned across models, and pLDDT /
    relSASA are averaged, so the structure-derived features consumed by
    Step 3 (minimal reversion, paratope grafting) are all preserved.

    Args:
        pdb_paths: PDB file paths (e.g. rank_1.pdb ... rank_N.pdb)
        chain_label: PDB chain id of the target chain
        position_map: {Kabat pos: residue number}
        cdr_positions: {Kabat pos: residue number} for CDR residues
        antigen_chains: Optional antigen chain ids
        min_consensus: Minimum number of agreeing models (default: 3)
        use_plddt: Interpret the B-factor column as an AF3 pLDDT score.

    Returns:
        Consensus StructureHints (all six fields populated).
    """
    if not pdb_paths:
        return StructureHints()

    buried_votes: Dict[str, List[bool]] = {}
    cdr_votes: Dict[str, List[bool]] = {}
    cdr_sc_votes: Dict[str, List[bool]] = {}
    ag_votes: Dict[str, List[bool]] = {}
    partners_all: Dict[str, set] = {}
    partners_sc_all: Dict[str, set] = {}
    plddt_vals: Dict[str, List[float]] = {}
    sasa_vals: Dict[str, List[float]] = {}
    n_models = 0

    for pdb_path in pdb_paths:
        if not os.path.exists(pdb_path):
            continue
        model = load_model(pdb_path)
        if not model:
            continue
        hints = compute_hints(
            model, chain_label, position_map, cdr_positions,
            antigen_chains, pdb_path=pdb_path, use_plddt=use_plddt,
        )
        n_models += 1
        for pos, v in (hints.data.get("buried") or {}).items():
            if v is not None:
                buried_votes.setdefault(pos, []).append(v)
        for pos, v in (hints.data.get("cdr_contact") or {}).items():
            cdr_votes.setdefault(pos, []).append(bool(v))
        for pos, v in (hints.data.get("cdr_contact_sc") or {}).items():
            cdr_sc_votes.setdefault(pos, []).append(bool(v))
        for pos, v in (hints.data.get("antigen_contact") or {}).items():
            ag_votes.setdefault(pos, []).append(bool(v))
        for pos, ps in (hints.data.get("cdr_partners") or {}).items():
            partners_all.setdefault(pos, set()).update(ps)
        for pos, ps in (hints.data.get("cdr_partners_sc") or {}).items():
            partners_sc_all.setdefault(pos, set()).update(ps)
        for pos, v in (hints.data.get("plddt") or {}).items():
            plddt_vals.setdefault(pos, []).append(v)
        for pos, v in (hints.data.get("rel_sasa") or {}).items():
            sasa_vals.setdefault(pos, []).append(v)

    if n_models == 0:
        return StructureHints()
    thr = min(min_consensus, n_models) if min_consensus else 1

    consensus_buried: Dict[str, Optional[bool]] = {}
    for pos, vs in buried_votes.items():
        if len(vs) >= thr:
            t = sum(1 for v in vs if v)
            f = len(vs) - t
            consensus_buried[pos] = True if t >= thr else (False if f >= thr else None)

    consensus_cdr = {
        pos: sum(1 for v in vs if v) >= thr
        for pos, vs in cdr_votes.items() if len(vs) >= thr
    }
    consensus_cdr_sc = {
        pos: sum(1 for v in vs if v) >= thr
        for pos, vs in cdr_sc_votes.items() if len(vs) >= thr
    }
    consensus_ag = {
        pos: sum(1 for v in vs if v) >= thr
        for pos, vs in ag_votes.items() if len(vs) >= thr
    }
    consensus_partners = {
        pos: sorted(ps) for pos, ps in partners_all.items() if ps
    }
    consensus_partners_sc = {
        pos: sorted(ps) for pos, ps in partners_sc_all.items() if ps
    }
    avg_plddt = {p: sum(v) / len(v) for p, v in plddt_vals.items() if v}
    avg_sasa = {p: sum(v) / len(v) for p, v in sasa_vals.items() if v}

    return StructureHints({
        "buried": consensus_buried,
        "cdr_contact": consensus_cdr,
        "cdr_contact_sc": consensus_cdr_sc,
        "antigen_contact": consensus_ag,
        "cdr_partners": consensus_partners,
        "cdr_partners_sc": consensus_partners_sc,
        "plddt": avg_plddt,
        "rel_sasa": avg_sasa,
    })


def cdr_rmsd(model_a: PDBModel, model_b: PDBModel, resseqs: List[int]) -> Optional[float]:
    """CA-based RMSD of a residue subset between two models in the SAME frame
    (no superposition). Kept for callers that already share a coordinate
    frame; use ``structure_rmsd`` for independently predicted models."""
    def ca(model):
        out = {}
        for a in model.atoms:
            if a.name == "CA":
                out[(a.chain, a.resseq)] = (a.x, a.y, a.z)
        return out

    ca_a, ca_b = ca(model_a), ca(model_b)
    pairs = []
    for ch in set(k[0] for k in ca_a) & set(k[0] for k in ca_b):
        for resseq in resseqs:
            if (ch, resseq) in ca_a and (ch, resseq) in ca_b:
                pairs.append((ca_a[(ch, resseq)], ca_b[(ch, resseq)]))
    if not pairs:
        return None
    s = sum(_dist(a, b) ** 2 for a, b in pairs)
    return math.sqrt(s / len(pairs))


# ---------------------------------------------------------------------------
# Superposition + CDR-RMSD validation (AF3 variant structures)
# ---------------------------------------------------------------------------

def ca_map(model: PDBModel, chain: str) -> Dict[int, Tuple[float, float, float]]:
    """{resseq: CA xyz} for one chain."""
    return {a.resseq: (a.x, a.y, a.z)
            for a in model.atoms if a.chain == chain and a.name == "CA"}


def kabsch_superpose(
    mobile: List[Tuple[float, float, float]],
    target: List[Tuple[float, float, float]],
):
    """Least-squares superposition of ``mobile`` onto ``target`` (Horn 1987).

    Returns ``(rmsd, R, t)`` where R is a 3x3 rotation matrix and t a
    translation, such that ``R @ mobile_i + t ~= target_i``. Pure stdlib
    (quaternion power iteration); no numpy required.
    """
    n = len(mobile)
    if n == 0 or n != len(target):
        return None
    cm = [sum(p[i] for p in mobile) / n for i in range(3)]
    ct = [sum(p[i] for p in target) / n for i in range(3)]
    P = [[p[i] - cm[i] for i in range(3)] for p in mobile]
    Q = [[p[i] - ct[i] for i in range(3)] for p in target]

    # covariance S = sum p_i q_i^T
    S = [[sum(P[k][i] * Q[k][j] for k in range(n)) for j in range(3)]
         for i in range(3)]
    Sxx, Sxy, Sxz = S[0]
    Syx, Syy, Syz = S[1]
    Szx, Szy, Szz = S[2]
    N = [
        [Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx],
        [Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz],
        [Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy],
        [Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz],
    ]
    # shift eigenvalues positive so power iteration finds the rotation
    shift = sum(N[i][i] for i in range(4))
    for i in range(4):
        N[i][i] += shift
    v = [1.0, 0.0, 0.0, 0.0]
    for _ in range(200):
        w = [sum(N[i][j] * v[j] for j in range(4)) for i in range(4)]
        norm = math.sqrt(sum(x * x for x in w))
        if norm < 1e-12:
            break
        v = [x / norm for x in w]
    qw, qx, qy, qz = v
    R = [
        [qw * qw + qx * qx - qy * qy - qz * qz, 2 * (qx * qy - qw * qz),
         2 * (qx * qz + qw * qy)],
        [2 * (qx * qy + qw * qz), qw * qw - qx * qx + qy * qy - qz * qz,
         2 * (qy * qz - qw * qx)],
        [2 * (qx * qz - qw * qy), 2 * (qy * qz + qw * qx),
         qw * qw - qx * qx - qy * qy + qz * qz],
    ]
    t = [ct[i] - sum(R[i][j] * cm[j] for j in range(3)) for i in range(3)]
    s = 0.0
    for p, q in zip(mobile, target):
        rp = [sum(R[i][j] * p[j] for j in range(3)) + t[i] for i in range(3)]
        s += _dist(rp, q) ** 2
    return math.sqrt(s / n), R, tuple(t)


def _apply(R, t, pt):
    return tuple(sum(R[i][j] * pt[j] for j in range(3)) + t[i] for i in range(3))


def count_clashes(model: PDBModel, min_dist: float = 2.0) -> Dict:
    """Count inter-residue heavy-atom clashes below ``min_dist`` (grid-accelerated).

    Returns ``{n_clashes, worst, examples}``; ``examples`` are
    ``(chain,resseq,atom,chain,resseq,atom,distance)`` tuples (up to 5).
    """
    heavy = model.heavy_atoms()  # [((ch, resseq), name, xyz)]
    if len(heavy) < 2:
        return {"n_clashes": 0, "worst": None, "examples": []}
    cell = max(min_dist, 0.1)
    grid: Dict[Tuple[int, int, int], List[int]] = {}
    key_of = [(int(math.floor(xyz[0] / cell)),
               int(math.floor(xyz[1] / cell)),
               int(math.floor(xyz[2] / cell))) for _, _, xyz in heavy]
    for i, k in enumerate(key_of):
        grid.setdefault(k, []).append(i)
    n_clashes = 0
    worst = None
    examples = []
    for i, ((res_i, name_i, xyz_i), k) in enumerate(zip(heavy, key_of)):
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for j in grid.get((k[0] + dx, k[1] + dy, k[2] + dz), ()):
                        if j <= i:
                            continue
                        res_j, name_j, xyz_j = heavy[j]
                        if res_i == res_j:
                            continue
                        d = _dist(xyz_i, xyz_j)
                        if d < min_dist:
                            n_clashes += 1
                            if worst is None or d < worst:
                                worst = d
                            if len(examples) < 5:
                                examples.append(
                                    (res_i[0], res_i[1], name_i,
                                     res_j[0], res_j[1], name_j, round(d, 2)))
    return {"n_clashes": n_clashes,
            "worst": round(worst, 2) if worst is not None else None,
            "examples": examples}


def _pos_key(pos: str):
    return (int("".join(c for c in pos if c.isdigit())), pos)


def _get_fr_region(pos: str, chain_type: str) -> str:
    """Get FR region (FR1/FR2/FR3/FR4) from Kabat position string.
    
    VH: FR1 1-30, CDR1 31-35, FR2 36-49, CDR2 50-65, FR3 66-94, CDR3 95-102, FR4 103-113
    VL: FR1 1-23, CDR1 24-34, FR2 35-49, CDR2 50-56, FR3 57-88, CDR3 89-97, FR4 98-107
    """
    # Extract numeric position from strings like "H31", "H100A", "L27E"
    num_str = ""
    for c in pos[1:]:  # Skip chain letter
        if c.isdigit():
            num_str += c
        elif num_str:  # Stop at first non-digit after digits
            break
    
    if not num_str:
        return "FR4"  # Default for unknown positions
    
    num = int(num_str)
    
    if chain_type == "H":
        if num <= 30:
            return "FR1"
        elif num <= 49:
            return "FR2"
        elif num <= 94:
            return "FR3"
        else:
            return "FR4"
    else:  # L chain
        if num <= 23:
            return "FR1"
        elif num <= 49:
            return "FR2"
        elif num <= 88:
            return "FR3"
        else:
            return "FR4"


def structure_rmsd(
    donor_model: PDBModel,
    donor_chain: str,
    donor_pos_to_resseq: Dict[str, int],
    variant_model: PDBModel,
    variant_chain: str,
    variant_pos_to_resseq: Dict[str, int],
    cdr_sets: Dict[str, set],
    framework_positions: set,
    chain_type: str = "H",
) -> Dict:
    """Superpose a variant on the donor framework, then measure CDR CA-RMSD.

    ``*_pos_to_resseq`` are {Kabat pos: PDB residue number} maps. Returns
    ``{cdr_rmsd, fr_rmsd, n_cdr, n_fr, per_cdr, per_fr}`` (RMSD in Angstrom, None when
    there is not enough aligned structure).
    """
    dca = ca_map(donor_model, donor_chain)
    vca = ca_map(variant_model, variant_chain)

    def xyz_map(pos_to_resseq, ca):
        out = {}
        for pos, resseq in pos_to_resseq.items():
            if resseq in ca:
                out[pos] = ca[resseq]
        return out

    dxyz = xyz_map(donor_pos_to_resseq, dca)
    vxyz = xyz_map(variant_pos_to_resseq, vca)

    fr_common = [p for p in sorted(framework_positions, key=_pos_key)
                 if p in dxyz and p in vxyz]
    if len(fr_common) < 3:
        return {"cdr_rmsd": None, "fr_rmsd": None, "n_cdr": 0,
                "n_fr": len(fr_common), "per_cdr": {}, "per_fr": {}}
    fit = kabsch_superpose([dxyz[p] for p in fr_common],
                           [vxyz[p] for p in fr_common])
    if fit is None:
        return {"cdr_rmsd": None, "fr_rmsd": None, "n_cdr": 0,
                "n_fr": len(fr_common), "per_cdr": {}, "per_fr": {}}
    fr_rmsd, R, t = fit

    per_cdr: Dict[str, Optional[float]] = {}
    all_pairs: List[float] = []
    for name, posset in cdr_sets.items():
        common = [p for p in sorted(posset, key=_pos_key)
                  if p in dxyz and p in vxyz]
        if not common:
            per_cdr[name] = None
            continue
        ss = sum(_dist(_apply(R, t, dxyz[p]), vxyz[p]) ** 2 for p in common)
        per_cdr[name] = round(math.sqrt(ss / len(common)), 3)
        all_pairs.extend(
            _dist(_apply(R, t, dxyz[p]), vxyz[p]) for p in common)

    cdr_rmsd = (round(math.sqrt(sum(d * d for d in all_pairs) / len(all_pairs)), 3)
                if all_pairs else None)
    
    # Calculate per-FR RMSD (FR1, FR2, FR3, FR4)
    per_fr: Dict[str, Optional[float]] = {}
    fr_regions = {"FR1": [], "FR2": [], "FR3": [], "FR4": []}
    for p in fr_common:
        region = _get_fr_region(p, chain_type)
        fr_regions[region].append(p)
    
    for region, positions in fr_regions.items():
        if not positions:
            per_fr[region] = None
            continue
        ss = sum(_dist(_apply(R, t, dxyz[p]), vxyz[p]) ** 2 for p in positions)
        per_fr[region] = round(math.sqrt(ss / len(positions)), 3)
    
    return {
        "cdr_rmsd": cdr_rmsd,
        "fr_rmsd": round(fr_rmsd, 3),
        "n_cdr": len(all_pairs),
        "n_fr": len(fr_common),
        "per_cdr": per_cdr,
        "per_fr": per_fr,
    }


# ---------------------------------------------------------------------------
# AF3 runner
# ---------------------------------------------------------------------------

@dataclass
class AF3Config:
    mode: str = "off"           # off | local | api
    binary: str = "run_alphafold.py"
    docker: bool = False
    workdir: str = "af3"
    model_dir: str = ""
    db_dir: str = ""
    api_url: str = ""
    api_token: str = ""
    n_structs: int = 2          # seeds per job


def make_af3_input(
    sequences: Dict[str, str],
    name: str,
    outdir: str,
) -> str:
    """Write an AF3 request JSON: sequences = {chain_id: aa_sequence}."""
    os.makedirs(outdir, exist_ok=True)
    payload = {
        "name": name,
        "modelSeeds": [1, 2, 3][:3],
        "sequences": [{"proteinChain": {"sequence": seq, "count": 1}} for seq in sequences.values()],
    }
    path = os.path.join(outdir, f"{name}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh)
    return path


def run_af3(cfg: AF3Config, input_json: str, outdir: str, tag: str = "") -> Optional[str]:
    """Submit to local AF3 (or API) and return the output PDB path."""
    if cfg.mode == "off":
        return None
    cmd: List[str] = []
    if cfg.mode == "local":
        if cfg.docker:
            cmd = ["docker", "run", "--gpus", "all", "-v", f"{outdir}:/app/data",
                   "-v", f"{cfg.model_dir}:/app/models", cfg.binary]
        else:
            cmd = [cfg.binary]
        cmd += [
            "--json_path", input_json,
            "--output_dir", outdir,
            "--output_format", "pdb",
        ]
        if cfg.model_dir:
            cmd += ["--model_dir", cfg.model_dir]
        if cfg.db_dir:
            cmd += ["--db_dir", cfg.db_dir]
        env = dict(os.environ)
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
        if proc.returncode != 0:
            raise RuntimeError(f"AF3 failed: {proc.stderr[-2000:]}")
        # locate the output pdb
        import glob
        pdbs = sorted(glob.glob(os.path.join(outdir, "**", "*_model_*.pdb"), recursive=True))
        return pdbs[0] if pdbs else None
    if cfg.mode == "api":
        import urllib.request
        from urllib.error import URLError, HTTPError
        with open(input_json, "rb") as fh:
            payload = fh.read()
        req = urllib.request.Request(
            cfg.api_url, data=payload,
            headers={"Content-Type": "application/json",
                     "Authorization": f"Bearer {cfg.api_token}"},
        )
        try:
            with urllib.request.urlopen(req, timeout=600) as r:
                out = json.loads(r.read())
        except (URLError, HTTPError, json.JSONDecodeError, OSError) as e:
            import warnings as _w
            _w.warn(f"[AF3 API] request failed: {e}")
            return None
        return out.get("pdb_path") or out.get("structure")
    return None


def predict_fv(
    cfg: AF3Config,
    vh_seq: str,
    vl_seq: Optional[str],
    antigen_seq: Optional[str],
    tag: str,
) -> Optional[str]:
    """Predict Fv (and optional complex) structure. Returns PDB path or None."""
    if cfg.mode == "off":
        return None
    chains: Dict[str, str] = {"H": vh_seq}
    if vl_seq:
        chains["L"] = vl_seq
    if antigen_seq:
        chains["A"] = antigen_seq
    os.makedirs(cfg.workdir, exist_ok=True)
    job = os.path.join(cfg.workdir, tag)
    inp = make_af3_input(chains, tag, job)
    return run_af3(cfg, inp, job, tag)
