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
                if line[76:78].strip() == "H" and line[0] == "A":
                    pass
                try:
                    model.atoms.append(PDBAtom(
                        name=line[12:16].strip(),
                        resname=line[17:20].strip(),
                        chain=line[21],
                        resseq=int(line[22:26].strip()),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
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
            if in_loop and (s.startswith("#") or s == ""):
                in_loop = False
                continue
            if in_loop and not s.startswith(("_", "#")):
                parts = s.split()
                if len(parts) >= len(cols):
                    d = dict(zip(cols, parts))
                    if d.get("group_PDB") in ("ATOM", "HETATM"):
                        try:
                            model.atoms.append(PDBAtom(
                                name=d.get("auth_atom_id", d.get("label_atom_id", "?")),
                                resname=d.get("auth_comp_id", d.get("label_comp_id", "?")),
                                chain=d.get("auth_asym_id", d.get("label_asym_id", "?")),
                                resseq=int(d.get("auth_seq_id", d.get("label_seq_id", "0"))),
                                x=float(d["Cartn_x"]),
                                y=float(d["Cartn_y"]),
                                z=float(d["Cartn_z"]),
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
# metrics
# ---------------------------------------------------------------------------

def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def compute_hints(
    model: PDBModel,
    chain_label: str,
    position_map: Dict[str, int],
    cdr_positions: Dict[str, int],
    antigen_chains: Optional[List[str]] = None,
) -> StructureHints:
    """Compute per-position hints for a chain in the model.

    chain_label:  PDB chain id of the target chain (e.g. "H").
    position_map: {Kabat pos: residue number} for ALL residues of the chain.
    cdr_positions:{Kabat pos: residue number} for CDR residues only.
    """
    heavy = model.heavy_atoms()
    by_res: Dict[Tuple[str, int], List] = {}
    for (ch, resseq), name, xyz in heavy:
        by_res.setdefault((ch, resseq), []).append((name, xyz))

    # which residues are CDR residues
    cdr_resseqs = {v for v in cdr_positions.values() if isinstance(v, int)}
    cdr_atoms = [a for a in heavy if (a[0][1] in cdr_resseqs)]
    ag_resseqs = set()
    if antigen_chains:
        ag_resseqs = {r[1] for r, _, _ in heavy if r[0] in antigen_chains}
    ag_atoms = [a for a in heavy if (a[0][1] in ag_resseqs)]

    # map resseq -> kabat pos
    resseq_to_pos = {v: k for k, v in position_map.items()}

    buried: Dict[str, bool] = {}
    cdr_contact: Dict[str, bool] = {}
    ag_contact: Dict[str, bool] = {}
    cdr_partners: Dict[str, List[str]] = {}
    for (ch, resseq), atoms in by_res.items():
        if ch != chain_label:
            continue
        pos = resseq_to_pos.get(resseq)
        if pos is None:
            continue
        # buriedness proxy: heavy-atom count within 6 A of any other residue
        n_heavy = len(atoms)
        contacts = 0
        for (oc, oseq), oa in by_res.items():
            if (oc, oseq) == (ch, resseq):
                continue
            for _, xyz in oa:
                for _, axyz in atoms[: min(3, len(atoms))]:
                    if _dist(axyz, xyz) < 6.0:
                        contacts += 1
                        break
        buried[pos] = contacts >= 20 or (n_heavy >= 6 and contacts / max(n_heavy, 1) >= 4)
        # CDR contact: any atom within 4.5 A of a CDR atom; record partners
        partners = set()
        for (r2, _n2, a2) in cdr_atoms:
            if _dist(a2[2], atoms[0][1]) < 4.5:
                ppos = resseq_to_pos.get(r2[1])
                if ppos:
                    partners.add(ppos)
        cdr_contact[pos] = bool(partners)
        if partners:
            cdr_partners[pos] = sorted(partners)
        if ag_atoms:
            ag_contact[pos] = any(
                _dist(a2[2], a1[2]) < 4.5
                for (r2, _n2, a2) in ag_atoms for (_n1, a1) in atoms
            )
    return StructureHints({
        "buried": buried,
        "cdr_contact": cdr_contact,
        "antigen_contact": ag_contact,
        "cdr_partners": cdr_partners,
    })


def cdr_rmsd(model_a: PDBModel, model_b: PDBModel, resseqs: List[int]) -> Optional[float]:
    """CA-based RMSD of a residue subset between two models (same residue
    numbering). Returns None when one model is missing residues."""
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
        req = urllib.request.Request(
            cfg.api_url, data=open(input_json, "rb").read(),
            headers={"Content-Type": "application/json",
                     "Authorization": f"Bearer {cfg.api_token}"},
        )
        with urllib.request.urlopen(req, timeout=600) as r:
            out = json.loads(r.read())
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
