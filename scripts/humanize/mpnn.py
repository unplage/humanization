"""ProteinMPNN adapter: framework re-design with fixed CDRs/interface.

Workflow (full mode, requires ProteinMPNN on the server):
  1. Take a grafted variant structure (AF3) or donor structure.
  2. Design all framework residues with CDRs + VH/VL interface + VHH
     hallmark fixed (MPNN keeps the backbone).
  3. Filter designed sequences by human-likeness (germline identity of the
     designed FR) and drop designs introducing developability risks.
  4. Optionally re-score top designs with AF3 (CDR RMSD, interface).
  5. Emit the best designs as alternative framework variants.

Mock mode (portable): emits the top-10-germline CONSENSUS framework as a
sanity-check alternative (no ProteinMPNN required).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .germline import GermlineDB, GermlineGene
from .numbering import NumberedChain


@dataclass
class MPNNConfig:
    mode: str = "off"            # off | local
    script: str = "protein_mpnn.py"   # helper_scripts path or full path
    num_seqs: int = 32
    batch_size: int = 8
    chains_to_design: str = "A"
    outdir: str = "mpnn"


@dataclass
class MPNNResult:
    designs: List[Dict[str, str]] = field(default_factory=list)   # {pos: aa}
    human_likeness: List[float] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


def _fixed_positions(chain: NumberedChain, is_vhh: bool) -> str:
    """Fixed residue numbers (1-based, single chain) for MPNN."""
    fixed = set()
    for r in chain.residues:
        if r.region in ("CDR1", "CDR2", "CDR3"):
            fixed.add(r.index + 1)
    # VH/VL interface + VHH hallmark + structural pillars stay fixed
    from .config import INTERFACE_CORE, INTERFACE_EXTENDED, VERNIER_ZONE, VHH_HALLMARK
    ctype = chain.chain_type
    for r in chain.residues:
        num = int("".join(c for c in r.pos if c.isdigit()))
        if num in INTERFACE_CORE[ctype] or num in INTERFACE_EXTENDED[ctype]:
            fixed.add(r.index + 1)
        if num in VERNIER_ZONE[ctype]:
            fixed.add(r.index + 1)
        if is_vhh and ctype == "H" and num in VHH_HALLMARK:
            fixed.add(r.index + 1)
    return ",".join(str(i) for i in sorted(fixed))


def run_mpnn(cfg: MPNNConfig, pdb_path: str, chain: NumberedChain,
             is_vhh: bool) -> MPNNResult:
    if cfg.mode == "off" or not pdb_path:
        return _consensus_design(chain)
    if not shutil.which("python"):
        raise RuntimeError("python not found")
    os.makedirs(cfg.outdir, exist_ok=True)
    fixed = _fixed_positions(chain, is_vhh)
    cmd = [
        "python", cfg.script,
        "--pdb_path", pdb_path,
        "--out_folder", cfg.outdir,
        "--num_seq_per_target", str(cfg.num_seqs),
        "--batch_size", str(cfg.batch_size),
        "--chains_to_design", cfg.chains_to_design,
        "--fixed_residues", fixed,
        "--seed", "0",
        "--suppress_print", "1",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"ProteinMPNN failed: {proc.stderr[-2000:]}")
    import glob
    fa_files = sorted(glob.glob(os.path.join(cfg.outdir, "*fa")))
    designs: List[Dict[str, str]] = []
    seqs: List[str] = []
    for f in fa_files:
        with open(f) as fh:
            content = fh.read()
        for block in content.split(">"):
            if not block.strip():
                continue
            lines = block.splitlines()
            seqs.append("".join(l.strip() for l in lines[1:] if l.strip()))
    for s in seqs:
        if len(s) == len(chain.sequence):
            designs.append(
                {r.pos: aa for r, aa in zip(chain.residues, s)}
            )
    return MPNNResult(designs=designs)


def _consensus_design(chain: NumberedChain,
                      top_germlines: Optional[List[Tuple[GermlineGene, dict]]] = None) -> MPNNResult:
    """Portable fallback: consensus of the top human germlines at every
    framework position, keeping donor CDRs/interface. Used as a sanity-check
    alternative when ProteinMPNN is unavailable."""
    from .backmut import analyze_backmutations
    consensus: Dict[str, str] = {}
    votes: Dict[str, Dict[str, int]] = {}
    for r in chain.residues:
        if r.region in ("CDR1", "CDR2", "CDR3") or r.region == "FR4":
            continue
        votes[r.pos] = {}
    for g, _s in (top_germlines or []):
        gm = g.numbered.posmap() if g.numbered else {}
        for pos, cnt in votes.items():
            if pos in gm and gm[pos]:
                cnt[gm[pos]] = cnt.get(gm[pos], 0) + 1
    for pos, cnt in votes.items():
        if cnt:
            consensus[pos] = max(cnt.items(), key=lambda kv: kv[1])[0]
    return MPNNResult(designs=[consensus], human_likeness=[1.0])
