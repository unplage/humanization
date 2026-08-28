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


@dataclass
class DevelopabilityRisk:
    """A single developability risk found in a sequence."""
    seq_pos: int
    kabat_pos: str
    motif: str
    motif_name: str
    risk: str  # high/medium
    issue: str  # deamidation/isomerization/etc
    context: str
    donor_aa: str = ""
    human_aa: str = ""
    skip_reason: Optional[str] = None  # Why this position was skipped for optimization
    rel_sasa: Optional[float] = None  # Relative SASA if available


@dataclass
class DevelopabilityOptimizationResult:
    """Results from developability optimization."""
    risks: List[DevelopabilityRisk] = field(default_factory=list)
    skipped_risks: List[DevelopabilityRisk] = field(default_factory=list)  # Risks skipped due to structural constraints
    optimized_designs: List[Dict[str, str]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    summary: str = ""


def detect_developability_risks(sequence: str, chain_type: str = "H") -> List[DevelopabilityRisk]:
    """Detect high-risk developability positions in a sequence.
    
    Returns list of DevelopabilityRisk objects.
    """
    import re
    from .config import HIGH_RISK_MOTIFS
    
    risks = []
    for name, info in HIGH_RISK_MOTIFS.items():
        pattern = info["pattern"]
        for match in re.finditer(pattern, sequence):
            start = match.start()
            context = sequence[max(0, start-2):min(len(sequence), start+len(match.group())+2)]
            risks.append(DevelopabilityRisk(
                seq_pos=start,
                kabat_pos="",  # Will be filled by numbering
                motif=match.group(),
                motif_name=name,
                risk=info["risk"],
                issue=info["issue"],
                context=context,
            ))
    return risks


def enrich_risks_with_structure(
    risks: List[DevelopabilityRisk],
    chain: NumberedChain,
    structure_hints: Optional["StructureHints"] = None,
) -> List[DevelopabilityRisk]:
    """Enrich developability risks with Kabat positions and relSASA from structure.
    
    Args:
        risks: List of detected risks
        chain: Numbered chain object for position mapping
        structure_hints: StructureHints from Step 3 (optional)
    
    Returns:
        List of enriched risks with kabat_pos and rel_sasa populated
    """
    if not structure_hints:
        return risks
    
    chain_type = chain.chain_type
    
    for risk in risks:
        # Map sequence position to Kabat position
        for r in chain.residues:
            if r.index == risk.seq_pos:  # 0-based index match
                risk.kabat_pos = r.pos
                # Get relSASA from structure hints
                rel_sasa = structure_hints.rel_sasa(r.pos)
                risk.rel_sasa = rel_sasa
                break
    
    return risks


def run_developability_optimization(
    cfg: MPNNConfig,
    pdb_path: str,
    chain: NumberedChain,
    risks: List[DevelopabilityRisk],
    structure_hints: Optional["StructureHints"] = None,  # From Step 3
    is_vhh: bool = False,
    top_germlines: Optional[List[Tuple[GermlineGene, dict]]] = None,
) -> DevelopabilityOptimizationResult:
    """Optimize high-risk developability positions using ProteinMPNN.
    
    This function:
    1. Identifies which positions have high-risk motifs (DD, NG, DG, etc.)
    2. Classifies positions based on structural data (buried/CDR-contact/exposed)
    3. Only optimizes surface-exposed positions (conservative strategy)
    4. Runs MPNN to generate alternative sequences
    5. Filters designs that eliminate the risk motifs
    6. Returns optimized designs with human-likeness scores
    
    Args:
        cfg: MPNN configuration
        pdb_path: Path to PDB structure
        chain: Numbered chain object
        risks: List of detected developability risks
        structure_hints: StructureHints from Step 3 (optional, for structural classification)
        is_vhh: Whether this is a VHH antibody
        top_germlines: Top germline candidates for human-likeness scoring
    
    Returns:
        DevelopabilityOptimizationResult with optimized designs
    """
    import re
    result = DevelopabilityOptimizationResult()
    
    if not risks:
        result.summary = "No high-risk developability positions detected."
        return result
    
    # Classify risks based on structural data
    optimizable = []  # Surface-exposed, safe to optimize
    skipped = []      # Structural constraints, skip optimization
    
    chain_type = chain.chain_type
    
    for risk in risks:
        # Map sequence position to Kabat position
        # Find the residue at this sequence position
        kabat_pos = ""
        for r in chain.residues:
            if r.index == risk.seq_pos:  # 0-based index match
                kabat_pos = r.pos
                risk.kabat_pos = r.pos
                break
        
        if not kabat_pos:
            # Cannot map to Kabat, skip with warning
            risk.skip_reason = "cannot map to Kabat position"
            skipped.append(risk)
            continue
        
        # Check structural constraints if available
        if structure_hints:
            # Check CDR/antigen contact - skip (affinity risk)
            if structure_hints.cdr_contact(chain_type, kabat_pos):
                risk.skip_reason = "CDR contact (affinity risk)"
                skipped.append(risk)
                continue
            
            if structure_hints.antigen_contact(chain_type, kabat_pos):
                risk.skip_reason = "antigen contact (binding risk)"
                skipped.append(risk)
                continue
            
            # Check buried - skip (stability risk)
            buried = structure_hints.buried(chain_type, kabat_pos)
            if buried is True:
                risk.skip_reason = "buried (stability risk)"
                skipped.append(risk)
                continue
            
            # Check relSASA - skip if too low
            rel_sasa = structure_hints.rel_sasa(kabat_pos)
            risk.rel_sasa = rel_sasa
            if rel_sasa is not None and rel_sasa < 0.20:
                risk.skip_reason = f"low SASA ({rel_sasa:.2f})"
                skipped.append(risk)
                continue
        
        # Surface-exposed or no structure data - safe to optimize
        risk.skip_reason = None
        optimizable.append(risk)
    
    result.risks = optimizable
    result.skipped_risks = skipped
    
    if not optimizable:
        result.summary = f"Found {len(risks)} high-risk positions, but all are structurally constrained (buried/CDR-contact). No optimization performed."
        return result
    
    if cfg.mode == "off" or not pdb_path:
        # Use consensus design as fallback
        mpnn_result = _consensus_design(chain, top_germlines)
        result.optimized_designs = mpnn_result.designs
        result.warnings.append("ProteinMPNN not available, using germline consensus")
    else:
        # Build fixed positions, excluding ONLY optimizable high-risk positions
        fixed = set()
        for r in chain.residues:
            if r.region in ("CDR1", "CDR2", "CDR3"):
                fixed.add(r.index + 1)
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
        
        # Remove ONLY optimizable positions from fixed set (not skipped ones)
        optimizable_seq_positions = {r.seq_pos + 1 for r in optimizable}  # Convert to 1-based
        fixed = fixed - optimizable_seq_positions
        
        fixed_str = ",".join(str(i) for i in sorted(fixed))
        
        # Run MPNN
        os.makedirs(cfg.outdir, exist_ok=True)
        cmd = [
            "python", cfg.script,
            "--pdb_path", pdb_path,
            "--out_folder", cfg.outdir,
            "--num_seq_per_target", str(cfg.num_seqs),
            "--batch_size", str(cfg.batch_size),
            "--chains_to_design", cfg.chains_to_design,
            "--fixed_residues", fixed_str,
            "--seed", "0",
            "--suppress_print", "1",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            result.warnings.append(f"ProteinMPNN failed: {proc.stderr[-500:]}")
            return result
        
        # Parse MPNN output
        import glob
        fa_files = sorted(glob.glob(os.path.join(cfg.outdir, "*fa")))
        designs = []
        for f in fa_files:
            with open(f) as fh:
                content = fh.read()
            for block in content.split(">"):
                if not block.strip():
                    continue
                lines = block.splitlines()
                seq = "".join(l.strip() for l in lines[1:] if l.strip())
                if len(seq) == len(chain.sequence):
                    designs.append(seq)
        
        # Filter designs that eliminate risk motifs
        from .config import HIGH_RISK_MOTIFS
        optimized = []
        for seq in designs:
            has_risk = False
            for name, info in HIGH_RISK_MOTIFS.items():
                if re.search(info["pattern"], seq):
                    has_risk = True
                    break
            if not has_risk:
                optimized.append({r.pos: aa for r, aa in zip(chain.residues, seq)})
        
        result.optimized_designs = optimized[:10]  # Keep top 10
    
    # Generate summary
    risk_summary = {}
    for r in risks:
        key = f"{r.risk}_{r.issue}"
        risk_summary[key] = risk_summary.get(key, 0) + 1
    
    result.summary = f"Found {len(risks)} high-risk positions: "
    result.summary += ", ".join(f"{v}x {k}" for k, v in risk_summary.items())
    result.summary += f". {len(optimizable)} optimizable, {len(skipped)} skipped (structural constraints)."
    if result.optimized_designs:
        result.summary += f" Generated {len(result.optimized_designs)} optimized designs."
    
    return result
