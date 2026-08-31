"""End-to-end humanization pipeline (portable mode: no external binaries).

Flow:
  input FASTA -> classify (VH/VL/VHH) -> number (Kabat engine)
  -> germline selection (human V/J per chain)
  -> CDR grafting (per scheme) -> back-mutation analysis & scoring
  -> variant ladder (V0-V3) -> optional AF3/MPNN (structure mode)
  -> reports (md/csv/json)

Structure mode (AF3/ProteinMPNN) is enabled via config on the server;
everything else runs with the Python standard library only.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional
from .minimal import MinimalReversion, MatrixEntry

from .backmut import BackMutationResult, StructureHints, analyze_backmutations
from .germline import (
    GermlineChoice,
    GermlineDB,
    choose_germlines,
    load_germline_db,
)
from .graft import GraftResult, graft_chain
from .mpnn import MPNNConfig
from .numbering import NumberedChain
from .sequences import InputChain, parse_input
from .structure import AF3Config, StructureHints as _SH, predict_fv
from .variants import Variant, assemble_variants


@dataclass
class PipelineConfig:
    germline_dir: str = ""
    germline_strategy: str = "auto"    # fr_best|cdr_best|composite|cvi_best|min_backmutations|current|auto
    forced_germlines: Dict[str, str] = field(default_factory=dict)  # {"H": "IGHV1-3*01", "L": "IGKV1-39*01"}
    cdr_scheme: str = "kabat"          # kabat|chothia|abm|imgt (graft def)
    report_schemes: List[str] = field(default_factory=lambda: ["kabat", "chothia", "imgt"])
    format: str = "auto"               # auto|fab|vhh
    af3: AF3Config = field(default_factory=AF3Config)
    mpnn: MPNNConfig = field(default_factory=MPNNConfig)
    antigen_seq: Optional[str] = None  # for AF3 complex prediction
    donor_structure: Optional[str] = None  # PDB/CIF of donor (for CDR RMSD)
    calibration_path: Optional[str] = None  # calibration.json from `humanize learn`
    biophi_env: Optional[str] = None        # conda env with biophi (server)
    oasis_db: Optional[str] = None          # OASis 9-mer DB path (server)
    mock_structures: bool = True       # run without AF3/MPNN

    def __post_init__(self):
        if self.cdr_scheme not in ("kabat", "chothia", "abm", "imgt"):
            raise ValueError(f"unknown CDR scheme: {self.cdr_scheme}")


@dataclass
class ChainReport:
    input_chain: InputChain
    germline: GermlineChoice
    backmut: BackMutationResult
    grafts: Dict[str, GraftResult] = field(default_factory=dict)
    variants: List[Variant] = field(default_factory=list)
    human_likeness: Dict[str, float] = field(default_factory=dict)
    structure_hints: StructureHints = field(default_factory=StructureHints)
    cvi_homology: float = 0.0
    minimal_reversion: Optional["MinimalReversion"] = None
    sdr_graft: Optional[GraftResult] = None
    matrix: List = field(default_factory=list)
    humanness: Dict[str, dict] = field(default_factory=dict)
    developability_optimization: Optional[object] = None


@dataclass
class RunResult:
    format: str
    chains: List[ChainReport]
    germline_db: GermlineDB
    warnings: List[str] = field(default_factory=list)
    config: Optional[PipelineConfig] = None


def human_likeness_percent(seq_graft: NumberedChain, v_gene) -> float:
    """% of framework residues matching the human germline."""
    if v_gene.numbered is None:
        return 0.0
    d = seq_graft.posmap()
    g = v_gene.numbered.posmap()
    fr = [p for p in d if p in g and seq_graft.region_of(p) in ("FR1", "FR2", "FR3")]
    if not fr:
        return 0.0
    return 100.0 * sum(1 for p in fr if d[p] == g[p]) / len(fr)


def run_pipeline(
    input_path: str,
    config: PipelineConfig,
    outdir: str = "outputs",
) -> RunResult:
    os.makedirs(outdir, exist_ok=True)
    fmt, chains = parse_input(input_path, config.format)
    warnings: List[str] = []
    for c in chains:
        warnings.extend(c.warnings)
    if fmt == "vhh_suspect":
        warnings.append(
            "single heavy chain without camelid hallmark detected - "
            "treating as VHH; verify species"
        )

    # germline DB (NCBI FASTA preferred, bundled JSON fallback)
    try:
        db = load_germline_db(config.germline_dir or _default_germline_dir())
    except FileNotFoundError as e:
        raise RuntimeError(
            f"germline database unavailable: {e}\n"
            "  run:  python3 scripts/humanize/cli.py setup-germline  (server)"
        ) from e

    reports: List[ChainReport] = []
    for chain in chains:
        rep = _process_chain(
            chain, fmt, db, config, outdir,
            antigen=config.antigen_seq,
            all_chains=chains,
        )
        reports.append(rep)

    return RunResult(
        format=fmt,
        chains=reports,
        germline_db=db,
        warnings=warnings,
        config=config,
    )


def _default_germline_dir() -> str:
    return os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "data", "germline",
    )


def _top_homologous_germlines(
    donor: NumberedChain,
    db: GermlineDB,
    n: int = 20,
    min_fr: float = 0.60,
) -> List[tuple]:
    """Top-N most homologous germlines (by FR identity), used as the
    reference panel for donor-residue conservation scoring.
    Distinct from the per-strategy winners: this panel reflects the full
    human germline repertoire, so conservation estimates are unbiased."""
    from .germline import compare_to_germline
    ctype = donor.chain_type
    scored = []
    for g in db.human(ctype):
        if g.numbered is None:
            continue
        s = compare_to_germline(donor, g)
        if s["fr_identity"] < min_fr:
            continue
        scored.append((g, s))
    scored.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))
    return scored[:n]


def _matrix_alternatives(
    donor: NumberedChain,
    main_v_gene,
    db: GermlineDB,
    n: int = 3,
    min_fr: float = 0.60,
) -> List[tuple]:
    """Deduplicated germline panel for the framework-matrix variants:
    keep the best allele per gene family (by FR identity), exclude the
    main choice's family, then take the top-N by FR identity."""
    from .germline import compare_to_germline
    ctype = donor.chain_type
    best_by_family: Dict[str, tuple] = {}
    for g in db.human(ctype):
        if g.numbered is None:
            continue
        s = compare_to_germline(donor, g)
        if s["fr_identity"] < min_fr:
            continue
        fam = g.gene_id.split("*")[0]
        cur = best_by_family.get(fam)
        if cur is None or s["fr_identity"] > cur[1]["fr_identity"]:
            best_by_family[fam] = (g, s)
    main_fam = main_v_gene.gene_id.split("*")[0]
    alts = [(g, s) for fam, (g, s) in best_by_family.items() if fam != main_fam]
    alts.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))
    return alts[:n]


def _process_chain(
    chain: InputChain,
    fmt: str,
    db: GermlineDB,
    config: PipelineConfig,
    outdir: str,
    antigen: Optional[str],
    all_chains: Optional[List[InputChain]] = None,
) -> ChainReport:
    donor = chain.numbered
    if donor is None:
        raise RuntimeError(f"[{chain.name}] no numbering available")
    ctype = donor.chain_type
    is_vhh = chain.is_vhh and ctype == "H"

    # ---- germline selection ----
    from .multi_strategy_germline import choose_germlines_multi_strategy
    
    # Check for forced germline
    forced_gene_id = config.forced_germlines.get(ctype)
    if forced_gene_id:
        # Use forced germline
        v_genes = db.v_for(ctype)
        v_gene = None
        for g in v_genes:
            if g.gene_id == forced_gene_id:
                v_gene = g
                break
        if v_gene is None:
            raise RuntimeError(f"[{chain.name}] forced germline {forced_gene_id} not found in database")
        
        # Select best J gene
        from .germline import score_j_match
        j_genes = db.j_for(ctype)
        j_scored = []
        for jg in j_genes:
            idn, n = score_j_match(donor, jg)
            j_scored.append((jg, idn, n))
        j_scored.sort(key=lambda t: (-t[1], -t[2]))
        j_gene = j_scored[0][0] if j_scored else None
        
        # Calculate identity for the forced gene
        from .germline import compare_to_germline
        scores = compare_to_germline(donor, v_gene)
        
        # Create choice object
        from .germline import GermlineChoice
        choice = GermlineChoice(
            v_gene=v_gene,
            j_gene=j_gene,
            scores=scores,
            alternatives=[],
        )
        
        # Also create multi_result for report compatibility
        multi_result = choose_germlines_multi_strategy(donor, db, is_vhh=is_vhh)
    else:
        # Use strategy-based selection
        strategy = config.germline_strategy
        if strategy == "auto":
            # 自动策略：VH 使用 adimab_frequency（回测最优），VL 使用 current（VL 最优）
            strategy = "adimab_frequency" if ctype == "H" else "current"
        
        multi_result = choose_germlines_multi_strategy(donor, db, is_vhh=is_vhh)
        candidate = multi_result.get_best(strategy)
        
        if candidate is None:
            # 回退到默认策略
            choice = choose_germlines(donor, db)
            v_gene, j_gene = choice.v_gene, choice.j_gene
        else:
            v_gene = candidate.gene
            # 选择 J 基因
            from .germline import score_j_match
            j_genes = db.j_for(ctype)
            j_scored = []
            for jg in j_genes:
                idn, n = score_j_match(donor, jg)
                j_scored.append((jg, idn, n))
            j_scored.sort(key=lambda t: (-t[1], -t[2]))
            j_gene = j_scored[0][0] if j_scored else None
            
            # 创建兼容的 choice 对象
            from .germline import GermlineChoice
            # Extract first candidate from each strategy's list
            all_alts = []
            for cand_list in multi_result.candidates.values():
                if isinstance(cand_list, list):
                    for c in cand_list[:3]:  # top 3 from each strategy
                        all_alts.append((c.gene, {"fr_identity": c.fr_identity, "cdr_identity": c.cdr_identity}))
                else:
                    all_alts.append((cand_list.gene, {"fr_identity": cand_list.fr_identity, "cdr_identity": cand_list.cdr_identity}))
            choice = GermlineChoice(
                v_gene=v_gene,
                j_gene=j_gene,
                scores={"fr_identity": candidate.fr_identity, "cdr_identity": candidate.cdr_identity},
                alternatives=all_alts,
            )
    
    if v_gene is None or j_gene is None:
        raise RuntimeError(f"[{chain.name}] no viable human germline found")

    # ---- structure hints (AF3 or donor structure) ----
    hints = StructureHints()
    af3_pdb = None
    if config.af3.mode != "off":
        os.makedirs(config.af3.workdir, exist_ok=True)
        af3_pdb = predict_fv(
            config.af3,
            donor.sequence,
            None if fmt == "vhh" else _partner_sequence(chain, all_chains),
            antigen,
            f"{chain.name}_{ctype}",
        )
    
    # Load structure from AF3 prediction or donor structure
    structure_path = af3_pdb if af3_pdb and os.path.exists(af3_pdb) else config.donor_structure
    if structure_path and os.path.exists(structure_path):
        from .structure import load_model, match_pdb_chain, compute_multi_model_consensus
        model = load_model(structure_path)
        if model:
            # Assign the PDB chain by sequence identity (longest common
            # residue run). First-residue matching is unreliable: unrelated
            # VH/VL chains often share the same N-terminal residue.
            pdb_chains = {}
            for atom in model.atoms:
                pdb_chains.setdefault(atom.chain, []).append(atom)

            label = match_pdb_chain(pdb_chains, donor.sequence)
            if label is None:
                # Fallback: conventional Fab chain ids
                label = "H" if ctype == "H" else "L"
            
            all_pos = {r.pos: r.index + 1 for r in donor.residues}
            from .graft import is_cdr_loop_position
            cdrs = {p: n for p, n in all_pos.items()
                    if is_cdr_loop_position(ctype, int("".join(c for c in p if c.isdigit())))}
            ag_chains = ["A"] if antigen else None
            
            # 方案2: 多模型共识 - 查找 rank_1-5 PDB 文件
            import glob
            pdb_dir = os.path.dirname(structure_path)
            pdb_base = os.path.basename(structure_path)
            # Match rank_*.pdb pattern
            all_pdbs = sorted(glob.glob(os.path.join(pdb_dir, "rank_*.pdb")))
            
            if len(all_pdbs) >= 3:
                # Use multi-model consensus (3+ models)
                hints = compute_multi_model_consensus(
                    all_pdbs, label, all_pos, cdrs, ag_chains,
                    min_consensus=3,
                )
            else:
                # Single model
                hints = _compute_hints_with_model(model, label, all_pos, cdrs, ag_chains, pdb_path=structure_path)

    # ---- back-mutation analysis ----
    # Conservation reference: top-N homologous germlines (unbiased panel),
    # not the per-strategy winners (which can repeat the same gene).
    top = _top_homologous_germlines(donor, db)
    calibration = None
    if config.calibration_path and os.path.exists(config.calibration_path):
        from .learning import load_calibration
        calibration = load_calibration(config.calibration_path)
    backmut = analyze_backmutations(
        donor, v_gene, is_vhh=is_vhh, structure=hints, top_germlines=top,
        calibration=calibration,
    )

    # ---- minimal-reversion & precision design ----
    from .minimal import (
        build_paratope_variant,
        cvi_homology,
        matrix_alternatives,
        minimal_reversion_set,
    )
    minrev = minimal_reversion_set(donor, backmut, structure=hints)
    sdr_graft = None
    if config.antigen_seq:
        sdr_graft = build_paratope_variant(
            donor, v_gene, j_gene, config.cdr_scheme, backmut, hints,
            is_vhh=is_vhh,
        )
    cvi = cvi_homology(donor, v_gene)
    matrix: List[MatrixEntry] = []
    matrix_alts = _matrix_alternatives(donor, v_gene, db)
    if matrix_alts:
        matrix = matrix_alternatives(
            donor, matrix_alts, j_gene, config.cdr_scheme,
            is_vhh=is_vhh, n=min(3, len(matrix_alts)),
        )

    # ---- grafts (all requested schemes, for reporting) ----
    grafts = {}
    for scheme in config.report_schemes:
        try:
            grafts[scheme] = graft_chain(donor, v_gene, j_gene, scheme, is_vhh)
        except ValueError as e:
            # InputChain always carries a warnings list; run_pipeline extends
            # RunResult.warnings from it, so the failure surfaces in the CLI
            # output and reports. (Do not use getattr-with-default here: a
            # fresh list would silently swallow the message.)
            chain.warnings.append(f"[{chain.name}] graft({scheme}) failed: {e}")

    # ---- variant ladder (main scheme) ----
    variants = assemble_variants(
        donor, v_gene, j_gene, config.cdr_scheme, backmut, is_vhh=is_vhh,
    )
    # append the affinity-preserving minimal variant (V2 + V_min comparison)
    if minrev.positions and set(minrev.positions) != set(backmut.revert_positions(("T1", "T2"))):
        from .variants import Variant
        from .graft import graft_variant
        g_min = graft_variant(
            donor, v_gene, j_gene, config.cdr_scheme, minrev.positions, is_vhh=is_vhh,
        )
        variants.append(Variant(
            name=f"{ctype}_Vmin",
            description=f"minimal reversion set ({minrev.method}, "
                        f"{minrev.covered_contacts}/{minrev.total_contacts} contacts kept)",
            graft=g_min,
            backmutations=minrev.positions,
        ))
    if sdr_graft is not None:
        from .variants import Variant
        variants.append(Variant(
            name=f"{ctype}_V_SDR",
            description="paratope-only grafting (antigen-contacting CDR "
                        "residues + structural pillars)",
            graft=sdr_graft,
            backmutations=[],
        ))

    # ---- human-likeness ----
    hl = {}
    for scheme, graft in grafts.items():
        hl[scheme] = round(human_likeness_percent(graft.numbered, v_gene), 1)

    # ---- BioPhi/Sapiens humanness cross-check (server, optional) ----
    humanness = {}
    if config.biophi_env:
        from .humanness import run_oasis, run_sapiens
        seqs = {f"{chain.name}|{v.name}": v.sequence for v in variants}
        sap = run_sapiens(seqs, env=config.biophi_env)
        for k, r in sap.items():
            humanness[k] = {"sapiens_mean": r.sapiens_mean, "note": r.note}
        if config.oasis_db:
            oas = run_oasis(seqs, config.oasis_db, env=config.biophi_env)
            for k, r in oas.items():
                humanness.setdefault(k, {})["oasis_identity"] = r.oasis_identity

    # ---- Step 4: Developability optimization (if high-risk motifs found AND mpnn enabled) ----
    from .mpnn import detect_developability_risks, enrich_risks_with_structure, run_developability_optimization, DevelopabilityOptimizationResult
    dev_opt_result = DevelopabilityOptimizationResult()
    
    # Check V2 variant (standard production candidate) for high-risk motifs
    v2_variant = None
    for v in variants:
        if v.name.endswith("_V2"):
            v2_variant = v
            break
    
    if v2_variant and v2_variant.graft and v2_variant.graft.numbered:
        v2_sequence = v2_variant.graft.numbered.sequence
        risks = detect_developability_risks(v2_sequence, ctype)
        
        if risks:
            # Enrich risks with Kabat positions and relSASA from structure
            risks = enrich_risks_with_structure(risks, donor, hints)
            
            if config.mpnn.mode != "off":
                # MPNN enabled - run optimization with structural classification
                dev_opt_result = run_developability_optimization(
                    config.mpnn, structure_path or "", donor, risks,
                    structure_hints=hints,  # Pass Step 3 structure data
                    is_vhh=is_vhh, top_germlines=top,
                )
            else:
                # MPNN mode off - still classify positions for reporting
                dev_opt_result = run_developability_optimization(
                    MPNNConfig(mode="off"), "", donor, risks,
                    structure_hints=hints,  # Pass Step 3 structure data
                    is_vhh=is_vhh, top_germlines=top,
                )

    return ChainReport(
        input_chain=chain,
        germline=choice,
        backmut=backmut,
        grafts=grafts,
        variants=variants,
        human_likeness=hl,
        structure_hints=hints,
        cvi_homology=cvi,
        minimal_reversion=minrev,
        sdr_graft=sdr_graft,
        matrix=matrix,
        humanness=humanness,
        developability_optimization=dev_opt_result,
    )


def _partner_sequence(chain: InputChain, _all_chains) -> Optional[str]:
    """Find the partner chain (VL for VH, VH for VL) for AF3 complex prediction.
    Previously returned None unconditionally, causing AF3 to predict a
    monomer instead of an Fv complex for Fab chains."""
    if not _all_chains:
        return None
    partner_type = "L" if chain.chain_type == "H" else "H"
    for c in _all_chains:
        if c.chain_type == partner_type and c.sequence:
            return c.sequence
    return None


def _compute_hints_with_model(model, label, all_pos, cdrs, ag_chains, pdb_path=None):
    from .structure import compute_hints
    return compute_hints(model, label, all_pos, cdrs, ag_chains, pdb_path=pdb_path)
