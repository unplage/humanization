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
    germline_strategy: str = "auto"    # fr_best|cdr_best|composite|cvi_best|min_backmutations|current|adimab_frequency|pioneer_frequency|composite_3axis|auto
    forced_germlines: Dict[str, str] = field(default_factory=dict)  # {"H": "IGHV1-3*01", "L": "IGKV1-39*01"}
    cdr_scheme: str = "kabat"          # kabat|chothia|abm|imgt (graft def)
    numbering: str = "kabat"           # kabat|imgt (coordinate space; imgt requires ANARCI)
    report_schemes: List[str] = field(default_factory=lambda: ["kabat", "chothia", "imgt"])
    format: str = "auto"               # auto|fab|vhh
    af3: AF3Config = field(default_factory=AF3Config)
    mpnn: MPNNConfig = field(default_factory=MPNNConfig)
    antigen_seq: Optional[str] = None  # for AF3 complex prediction
    donor_structure: Optional[str] = None  # PDB/CIF of donor (for CDR RMSD)
    # Treat the B-factor column of --donor-structure as a pLDDT score. True
    # for AF3/Boltz predictions (their B-factor IS the pLDDT), false for
    # experimental crystal/NMR structures (temperature factor). Without this
    # the pLDDT confidence gate and the structural confidence factor are
    # inactive for donor models.
    donor_structure_plddt: bool = False
    calibration_path: Optional[str] = None  # calibration.json from `humanize learn`
    biophi_env: Optional[str] = None        # conda env with biophi (server)
    oasis_db: Optional[str] = None          # OASis 9-mer DB path (server)
    mock_structures: bool = True       # unused legacy field; AF3/MPNN are controlled by af3.mode/mpnn.mode
    interactive_indel: bool = False    # enable interactive indel selection
    af3_validate_variants: bool = False  # predict each variant and compute CDR-RMSD
    design_panel: bool = False         # emit structure-guided V_opt panel
    design_panel_steps: int = 3
    # Optional per-position immunogenicity (MHC-II epitope) refinement
    immunogenicity_path: Optional[str] = None   # precomputed {kabat_pos: score}
    netmhciipan_binary: Optional[str] = None    # NetMHCIIpan executable
    netmhciipan_alleles: str = ""
    netmhciipan_env: Optional[str] = None       # conda env for NetMHCIIpan

    def __post_init__(self):
        if self.cdr_scheme not in ("kabat", "chothia", "abm", "imgt"):
            raise ValueError(f"unknown CDR scheme: {self.cdr_scheme}")
        if self.numbering not in ("kabat", "imgt"):
            raise ValueError(f"unknown numbering: {self.numbering}")
        if self.numbering == "imgt":
            # IMGT numbering uses the native IMGT CDR windows only; the other
            # three schemes are Kabat-backbone windows and would mix spaces.
            if self.cdr_scheme != "imgt":
                raise ValueError(
                    f"--numbering imgt requires CDR scheme 'imgt' "
                    f"(got {self.cdr_scheme!r})")
            self.report_schemes = ["imgt"]


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
    structure_validation: Dict[str, dict] = field(default_factory=dict)
    donor_structure: Optional[object] = None   # (model, chain_label, {pos: resseq})


@dataclass
class RunResult:
    format: str
    chains: List[ChainReport]
    germline_db: GermlineDB
    warnings: List[str] = field(default_factory=list)
    config: Optional[PipelineConfig] = None


def human_likeness_percent(seq_graft: NumberedChain, v_gene, corr=None) -> float:
    """% of framework residues matching the human germline.

    When a donor<->germline correspondence is supplied, aligned columns are
    compared (donor insertions are excluded), so an FR insertion does not count
    as a framework mismatch.
    """
    if v_gene.numbered is None:
        return 0.0
    d = seq_graft.posmap()
    g = v_gene.numbered.posmap()
    if corr is not None and corr.germline_to_donor:
        pairs = [
            (dpos, gpos)
            for gpos, dpos in corr.germline_to_donor.items()
            if dpos in d and gpos in g
            and seq_graft.region_of(dpos) in ("FR1", "FR2", "FR3")
        ]
        if not pairs:
            return 0.0
        return 100.0 * sum(1 for dpos, gpos in pairs if d[dpos] == g[gpos]) / len(pairs)
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
    fmt, chains = parse_input(input_path, config.format,
                              numbering=config.numbering)
    warnings: List[str] = []
    for c in chains:
        warnings.extend(c.warnings)
    if fmt == "vhh_suspect":
        warnings.append(
            "single heavy chain without camelid hallmark detected - "
            "treating as VHH; verify species"
        )

    # germline DB (bundled Kabat-space JSON preferred, NCBI FASTA fallback)
    try:
        db = load_germline_db(config.germline_dir or _default_germline_dir(),
                              numbering=config.numbering)
    except FileNotFoundError as e:
        raise RuntimeError(
            f"germline database unavailable: {e}\n"
            "  run:  python3 scripts/humanize/cli.py setup-germline  (server)"
        ) from e
    _validate_numbering_space(chains, db, config)

    reports: List[ChainReport] = []
    for chain in chains:
        rep = _process_chain(
            chain, fmt, db, config, outdir,
            antigen=config.antigen_seq,
            all_chains=chains,
        )
        reports.append(rep)

    result = RunResult(
        format=fmt,
        chains=reports,
        germline_db=db,
        warnings=warnings,
        config=config,
    )
    # Step 3b: variant structure validation needs both chains (Fab), so it runs
    # once here rather than per chain.
    _validate_all_variants(result, config)
    return result


def _default_germline_dir() -> str:
    return os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "data", "germline",
    )


def _validate_numbering_space(chains: List[InputChain], db: GermlineDB,
                              config: PipelineConfig) -> None:
    """Fail fast when donor chains or the germline DB are not in the
    configured coordinate space (IMGT mode only; Kabat keeps historical
    lenient behaviour)."""
    if config.numbering != "imgt":
        return
    from .schemespec import cdr_windows, get_spec
    if "imgt" not in (db.source_dir or ""):
        raise RuntimeError(
            f"IMGT numbering requested but the germline DB is not IMGT-space "
            f"({db.source_dir!r})")
    spec = get_spec(config.numbering)
    windows = cdr_windows(config.numbering, config.cdr_scheme)
    for c in chains:
        n = c.numbered
        if n is None:
            continue
        ct = n.chain_type
        probes = (f"{ct}{windows[ct]['CDR1'][0]}",
                  f"{ct}{spec.j_anchor_for(ct)}")
        missing = [lab for lab in probes if n.residue(lab) is None]
        if missing:
            raise RuntimeError(
                f"[{c.name}] numbering space mismatch: --numbering imgt but "
                f"residues {missing} absent (chain is not IMGT-numbered)")


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
        from .germline_frequency import get_j_frequency
        j_genes = db.j_for(ctype)
        j_scored = []
        for jg in j_genes:
            idn, n = score_j_match(donor, jg, config.numbering)
            j_scored.append((jg, idn, n, get_j_frequency(jg.gene_id)))
        # identity desc, then usage-frequency desc (industry convention tiebreak), then matched length
        j_scored.sort(key=lambda t: (-t[1], -t[3], -t[2]))
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
        multi_result = choose_germlines_multi_strategy(
            donor, db, is_vhh=is_vhh, numbering=config.numbering)
    else:
        # Use strategy-based selection
        strategy = config.germline_strategy
        if strategy == "auto":
            # 自动策略：VH 使用 adimab_frequency（回测最优），VL 使用 current（VL 最优）
            strategy = "adimab_frequency" if ctype == "H" else "current"
        
        multi_result = choose_germlines_multi_strategy(
            donor, db, is_vhh=is_vhh, numbering=config.numbering)
        candidate = multi_result.get_best(strategy)
        
        if candidate is None:
            # 回退到默认策略
            choice = choose_germlines(donor, db, numbering=config.numbering)
            v_gene, j_gene = choice.v_gene, choice.j_gene
        else:
            v_gene = candidate.gene
            # 选择 J 基因
            from .germline import score_j_match
            j_genes = db.j_for(ctype)
            j_scored = []
            for jg in j_genes:
                idn, n = score_j_match(donor, jg, config.numbering)
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
    donor_struct = None  # (model, chain_label, {kabat_pos: resseq})
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
    from_af3 = bool(af3_pdb and os.path.exists(af3_pdb))
    structure_path = af3_pdb if from_af3 else config.donor_structure
    if structure_path and os.path.exists(structure_path):
        from .structure import (
            load_model, match_pdb_chain, match_pdb_chain_loose,
            compute_multi_model_consensus,
        )
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
                # Exact-run match failed (renumbered/tagged structure, or a
                # sequence that drifted from the model). Fall back to global
                # sequence identity; refuse to guess if it is too low, so we
                # never annotate structural features on the wrong chain.
                label = match_pdb_chain_loose(pdb_chains, donor.sequence)
                if label is not None:
                    chain.warnings.append(
                        f"[{chain.name}] exact chain match failed; using chain "
                        f"'{label}' by sequence identity")
                else:
                    chain.warnings.append(
                        f"[{chain.name}] donor structure has no chain matching "
                        f"the input sequence; structure hints skipped (check "
                        f"that --donor-structure belongs to this antibody)")

            if label is not None:
                # Antigen chains: every chain that is neither the target chain
                # nor its partner. AF3 renumbers output chains (A/B/C), so a
                # fixed chain id must never be assumed.
                ag_chains = None
                if antigen:
                    used = {label}
                    partner_seq = None if fmt == "vhh" else _partner_sequence(chain, all_chains)
                    if partner_seq:
                        plabel = match_pdb_chain(pdb_chains, partner_seq)
                        if plabel:
                            used.add(plabel)
                    ag_chains = [c for c in pdb_chains if c not in used] or None

                # Map donor Kabat positions to PDB residue numbers. When the
                # PDB chain has exactly one CA per donor residue (the usual AF3
                # case, and tagged experimental PDBs), use the *actual* PDB
                # numbering instead of assuming resseq == index+1.
                ca_sorted = sorted(
                    (a.resseq, a.resname) for a in pdb_chains[label] if a.name == "CA")
                if len(ca_sorted) == len(donor.residues):
                    all_pos = {r.pos: ca_sorted[i][0]
                               for i, r in enumerate(donor.residues)}
                else:
                    all_pos = {r.pos: r.index + 1 for r in donor.residues}
                    chain.warnings.append(
                        f"[{chain.name}] structure chain has {len(ca_sorted)} CA "
                        f"vs {len(donor.residues)} donor residues; assuming "
                        f"sequential numbering")
                donor_struct = (model, label, all_pos)
                from .graft import is_cdr_loop_position, loop_of_position
                cdrs = {p: n for p, n in all_pos.items()
                        if is_cdr_loop_position(
                            ctype, int("".join(c for c in p if c.isdigit())),
                            config.cdr_scheme, config.numbering)}
                # Per-loop map (CDR1/CDR2/CDR3) for the loop-class correction:
                # the literature canonical/vernier sets are loop-specific, so
                # the correction needs to know which loop each position claims.
                cdr_loops = {}
                for p in cdrs:
                    loop = loop_of_position(
                        ctype, int("".join(c for c in p if c.isdigit())),
                        config.cdr_scheme, config.numbering)
                    if loop:
                        cdr_loops[p] = loop
                # B-factor is only a pLDDT score for AF3 predictions; an
                # experimental `--donor-structure` carries a temperature factor.
                use_plddt = from_af3 or config.donor_structure_plddt

                # Multi-model consensus over sibling models if available
                import glob
                pdb_dir = os.path.dirname(structure_path)
                found = set()
                for pat in ("rank_*.pdb", "*_model*.pdb"):
                    found.update(glob.glob(os.path.join(pdb_dir, pat)))
                all_pdbs = sorted(found)

                if len(all_pdbs) >= 3:
                    hints = compute_multi_model_consensus(
                        all_pdbs, label, all_pos, cdrs, ag_chains,
                        min_consensus=3, use_plddt=use_plddt,
                        cdr_loops=cdr_loops,
                    )
                else:
                    # Single model
                    hints = _compute_hints_with_model(
                        model, label, all_pos, cdrs, ag_chains,
                        pdb_path=structure_path, use_plddt=use_plddt,
                        cdr_loops=cdr_loops)

    # ---- FR indel detection and interactive selection ----
    from .fr_indel import detect_fr_indels, select_insertion_interactive, update_indel_selection
    
    fr_indels = detect_fr_indels(donor, v_gene, config.numbering)
    indel_overrides: Dict[str, str] = {}  # region -> selected position
    
    if config.interactive_indel and fr_indels:
        for indel in fr_indels:
            if indel.candidates and len(indel.candidates) > 1:
                selected_pos = select_insertion_interactive(indel)
                if selected_pos and selected_pos != indel.position:
                    indel_overrides[indel.fr_region] = selected_pos
                    # Update the indel object
                    indel = update_indel_selection(indel, selected_pos)
                    print(f"[{ctype}] 已更新 {indel.fr_region} 插入位置: {selected_pos}")

    # ---- back-mutation analysis ----
    # Conservation reference: the FULL human repertoire (frequency-weighted in
    # backmut._conservation), not just the closest homologs.
    top = _top_homologous_germlines(donor, db, n=100000, min_fr=0.0)
    calibration = None
    if config.calibration_path and os.path.exists(config.calibration_path):
        from .learning import load_calibration
        calibration = load_calibration(config.calibration_path)

    # Optional per-position immunogenicity (MHC-II epitope) refinement
    immunogenicity = None
    if config.immunogenicity_path and os.path.exists(config.immunogenicity_path):
        if config.numbering != "kabat":
            # Precomputed maps are keyed by Kabat positions; IMGT labels with
            # the same name are different residues. Live NetMHCIIpan below is
            # numbering-safe (mapped via sequence index).
            chain.warnings.append(
                f"[{chain.name}] precomputed immunogenicity map is Kabat-keyed "
                f"and ignored in IMGT numbering mode (use --netmhciipan-binary "
                f"for a live, numbering-safe run)")
        else:
            from .immunogenicity import load_scores_json
            immunogenicity = load_scores_json(config.immunogenicity_path)
    elif config.netmhciipan_binary:
        from .immunogenicity import netmhciipan_scores
        immunogenicity = netmhciipan_scores(
            donor, config.netmhciipan_binary,
            alleles=config.netmhciipan_alleles or
            ("HLA-DRB1*01:01,HLA-DRB1*04:01,HLA-DRB1*07:01,HLA-DRB1*15:01,"
             "HLA-DRB3*01:01,HLA-DRB3*02:02,HLA-DRB4*01:01,HLA-DRB5*01:01"),
            env=config.netmhciipan_env,
            workdir=os.path.join(outdir, "immunogenicity", ctype),
        )

    backmut = analyze_backmutations(
        donor, v_gene, is_vhh=is_vhh, structure=hints, top_germlines=top,
        calibration=calibration, indel_overrides=indel_overrides,
        j_gene=j_gene, immunogenicity=immunogenicity,
        numbering=config.numbering,
    )
    
    # P1: Annotate with functional conservation scores
    from .backmut import _annotate_with_functional_conservation
    backmut = _annotate_with_functional_conservation(
        backmut, db, donor, structure_data=hints
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
            is_vhh=is_vhh, numbering=config.numbering,
        )
    cvi = cvi_homology(donor, v_gene, config.numbering)
    matrix: List[MatrixEntry] = []
    matrix_alts = _matrix_alternatives(donor, v_gene, db)
    if matrix_alts:
        matrix = matrix_alternatives(
            donor, matrix_alts, j_gene, config.cdr_scheme,
            is_vhh=is_vhh, n=min(3, len(matrix_alts)),
            numbering=config.numbering,
        )

    # ---- grafts (all requested schemes, for reporting) ----
    grafts = {}
    for scheme in config.report_schemes:
        try:
            grafts[scheme] = graft_chain(
                donor, v_gene, j_gene, scheme, is_vhh,
                fr_indels=backmut.fr_indels, numbering=config.numbering)
        except ValueError as e:
            # InputChain always carries a warnings list; run_pipeline extends
            # RunResult.warnings from it, so the failure surfaces in the CLI
            # output and reports. (Do not use getattr-with-default here: a
            # fresh list would silently swallow the message.)
            chain.warnings.append(f"[{chain.name}] graft({scheme}) failed: {e}")

    # ---- variant ladder (main scheme) ----
    variants = assemble_variants(
        donor, v_gene, j_gene, config.cdr_scheme, backmut, is_vhh=is_vhh,
        numbering=config.numbering,
    )
    # append the affinity-preserving minimal variant (V2 + V_min comparison)
    if minrev.positions and set(minrev.positions) != set(backmut.revert_positions(("T1", "T2"))):
        from .variants import Variant
        from .graft import graft_variant
        g_min = graft_variant(
            donor, v_gene, j_gene, config.cdr_scheme, minrev.positions,
            is_vhh=is_vhh, fr_indels=backmut.fr_indels,
            numbering=config.numbering,
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
    # structure-guided multi-objective design panel (opt-in)
    if config.design_panel:
        from .variants import structure_guided_panel
        variants.extend(structure_guided_panel(
            donor, v_gene, j_gene, config.cdr_scheme, backmut,
            structure=hints, is_vhh=is_vhh, n_extra=config.design_panel_steps,
            numbering=config.numbering,
        ))

    # ---- human-likeness ----
    hl = {}
    for scheme, graft in grafts.items():
        hl[scheme] = round(human_likeness_percent(
            graft.numbered, v_gene, graft.fr_correspondence), 1)

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

    # ---- Step 4: Developability optimization (runs when high-risk motifs are
    # found; mpnn.mode controls real MPNN vs germline-consensus fallback) ----
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
                    numbering=config.numbering,
                )
            else:
                # MPNN mode off - still classify positions for reporting
                dev_opt_result = run_developability_optimization(
                    MPNNConfig(mode="off"), "", donor, risks,
                    structure_hints=hints,  # Pass Step 3 structure data
                    is_vhh=is_vhh, top_germlines=top,
                    numbering=config.numbering,
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
        structure_validation={},
        donor_structure=donor_struct,
    )


def build_cdr_framework_sets(
    numbered: NumberedChain, scheme: str = "kabat", numbering: str = "kabat"
):
    """Return (cdr_sets, framework_positions) as position-label sets.

    ``cdr_sets`` maps the scheme CDR name (CDR1/CDR2/CDR3) to its position
    labels; ``framework_positions`` is FR1-FR3 excluding any CDR-loop residue
    (e.g. H93/H94 are labelled FR3 in the strict numbering scheme but belong
    to the CDR3 loop set used here (graft windows), so they must not
    enter the fit).
    """
    from .schemespec import cdr_windows
    ctype = numbered.chain_type
    bounds = cdr_windows(numbering, scheme)[ctype]
    cdr_sets: Dict[str, set] = {name: set() for name in bounds}
    for r in numbered.residues:
        num = int("".join(c for c in r.pos if c.isdigit()) or -1)
        for name, (lo, hi) in bounds.items():
            if lo <= num <= hi:
                cdr_sets[name].add(r.pos)
                break
    cdr_all = set().union(*cdr_sets.values()) if cdr_sets else set()
    framework = {
        r.pos for r in numbered.residues
        if r.region in ("FR1", "FR2", "FR3") and r.pos not in cdr_all
    }
    return cdr_sets, framework


def pdb_chain_pos_map(model, chain: str, numbered: NumberedChain) -> Dict[str, int]:
    """Map numbered Kabat positions to PDB residue numbers for one chain.

    Prefers the actual PDB numbering when the chain has exactly one CA per
    numbered residue; otherwise falls back to sequential 1..N.
    """
    ca = sorted((a.resseq, a.resname)
                for a in model.atoms if a.chain == chain and a.name == "CA")
    if len(ca) == len(numbered.residues):
        return {r.pos: ca[i][0] for i, r in enumerate(numbered.residues)}
    return {r.pos: r.index + 1 for r in numbered.residues}


def align_pos_map(donor_numbered, vseq: str, model, chain: str,
                  min_identity: float = 0.5) -> Optional[Dict[str, int]]:
    """Map donor Kabat positions to variant PDB resseqs via sequence alignment.

    Two independent numbering runs (donor FASTA vs a variant PDB sequence) can
    place CDR insertions at different anchor residues (e.g. donor L30A-F vs
    variant L27A-F for the same 17-residue VL CDR1); pairing residues by label
    then mismatches them and inflates CDR-RMSD. Align donor sequence to the
    variant CA sequence and key the variant map by donor labels; gapped
    residues are omitted. Returns None when the alignment is too poor to trust
    (caller should fall back to label-based pairing).
    """
    from .fr_indel import _needleman_wunsch

    ca = sorted((a.resseq, a.resname)
                for a in model.atoms if a.chain == chain and a.name == "CA")
    if len(ca) != len(vseq):
        return None
    aln = _needleman_wunsch(list(donor_numbered.sequence), list(vseq))
    matches = sum(1 for a, b in aln if a != "-" and b != "-" and a == b)
    aligned = sum(1 for a, b in aln if a != "-" and b != "-")
    if not aligned or matches / aligned < min_identity:
        return None
    out: Dict[str, int] = {}
    di = vi = 0
    for a, b in aln:
        if a != "-" and b != "-":
            out[donor_numbered.residues[di].pos] = ca[vi][0]
        if a != "-":
            di += 1
        if b != "-":
            vi += 1
    return out


def _plddt_by_pos(model, label: str, numbered: NumberedChain) -> Dict[str, float]:
    """{Kabat pos: AF3 pLDDT} for one chain (empty for experimental PDBs)."""
    pos2res = pdb_chain_pos_map(model, label, numbered)
    by_res = {a.resseq: a.plddt for a in model.atoms
              if a.chain == label and a.name == "CA" and a.plddt > 0}
    return {p: by_res[res] for p, res in pos2res.items() if res in by_res}


def _measure_variant(donor_struct, config: PipelineConfig, donor_numbered,
                     variant, vmodel, vlabel: str, pdb: str) -> Dict:
    """CDR-RMSD + clashes + CDR pLDDT for one variant chain vs its donor."""
    from .structure import count_clashes, structure_rmsd, _AA3TO1
    donor_model, donor_label, donor_pos = donor_struct
    cdr_sets, framework = build_cdr_framework_sets(
        donor_numbered, config.cdr_scheme, config.numbering)
    # Pair donor positions to variant PDB resseqs via sequence alignment.
    # Two independent numbering runs (donor FASTA vs variant PDB sequence)
    # can place CDR insertions at different anchor residues; pairing by
    # label alone then mismatches residues and inflates CDR-RMSD.
    ca = sorted((a.resseq, a.resname)
                for a in vmodel.atoms if a.chain == vlabel and a.name == "CA")
    vseq = "".join(_AA3TO1.get(rn, "X") for _, rn in ca)
    vpos = align_pos_map(donor_numbered, vseq, vmodel, vlabel)
    if vpos is None:
        vpos = pdb_chain_pos_map(vmodel, vlabel, variant.graft.numbered)
    # Use only high-confidence framework residues for the superposition so a
    # flexible/low-pLDDT loop or terminus cannot drag the whole fit.
    plddt = _plddt_by_pos(vmodel, vlabel, variant.graft.numbered)
    fit_positions = framework
    if plddt:
        high_conf = {p for p in framework if plddt.get(p, 100.0) >= 70.0}
        if len(high_conf) >= 8:
            fit_positions = high_conf
    donor_regions = {r.pos: r.region for r in donor_numbered.residues}
    res = structure_rmsd(donor_model, donor_label, donor_pos,
                         vmodel, vlabel, vpos, cdr_sets, fit_positions,
                         pos_regions=donor_regions)
    clashes = count_clashes(vmodel)
    res["n_clashes"] = clashes["n_clashes"]
    res["worst_clash"] = clashes["worst"]
    cdr_all = set().union(*cdr_sets.values()) if cdr_sets else set()
    cdr_vals = [plddt[p] for p in cdr_all if p in plddt]
    res["cdr_plddt"] = round(sum(cdr_vals) / len(cdr_vals), 1) if cdr_vals else None
    res["pdb"] = pdb
    return res


def _validate_all_variants(result: "RunResult", config: PipelineConfig) -> None:
    """Predict variants and fill each chain's ``structure_validation``.

    For Fab inputs the VH and VL variants sharing a suffix (e.g. H_V2 / L_V2)
    are predicted TOGETHER as a real variant Fv (variant VH + variant VL); if
    only one side has that variant, the donor partner chain is used. VHH and
    single-chain inputs are predicted as monomers (plus antigen when given).
    """
    if not (config.af3_validate_variants and config.af3.mode != "off"):
        return
    from .structure import load_model, match_pdb_chain

    by_type: Dict[str, "ChainReport"] = {}
    for rep in result.chains:
        by_type.setdefault(rep.input_chain.chain_type, rep)
    antigen = config.antigen_seq

    def _predict(vh_seq, vl_seq, tag):
        try:
            pdb = predict_fv(config.af3, vh_seq, vl_seq, antigen, tag)
        except Exception as e:  # AF3 failure must not abort the run
            result.warnings.append(f"[AF3 RMSD] {tag} failed: {e}")
            return None, {}
        if not pdb or not os.path.exists(pdb):
            return None, {}
        model = load_model(pdb)
        if not model:
            return None, {}
        chains = {}
        for atom in model.atoms:
            chains.setdefault(atom.chain, []).append(atom)
        return model, chains

    # single-chain (VHH / unpaired): monomer + optional antigen
    if result.format == "vhh" or not ({"H", "L"} <= set(by_type)):
        for rep in result.chains:
            if rep.donor_structure is None:
                continue
            for v in rep.variants:
                model, chains = _predict(v.sequence, None, f"{rep.input_chain.name}_{v.name}")
                if model is None:
                    continue
                lab = match_pdb_chain(chains, v.sequence)
                if lab is None:
                    continue
                rep.structure_validation[v.name] = _measure_variant(
                    rep.donor_structure, config, rep.input_chain.numbered,
                    v, model, lab, "")
        return

    hrep, lrep = by_type["H"], by_type["L"]
    dh, dl = hrep.donor_structure, lrep.donor_structure
    if dh is None and dl is None:
        return

    def _suffix(name: str) -> str:
        return name.split("_", 1)[1] if "_" in name else name

    hmap = {_suffix(v.name): v for v in hrep.variants}
    lmap = {_suffix(v.name): v for v in lrep.variants}
    donor_h_seq = hrep.input_chain.sequence
    donor_l_seq = lrep.input_chain.sequence

    # Fab: predict each variant pair ONCE as a real variant Fv
    for suffix in sorted(set(hmap) | set(lmap)):
        hv, lv = hmap.get(suffix), lmap.get(suffix)
        vh_seq = hv.sequence if hv is not None else donor_h_seq
        vl_seq = lv.sequence if lv is not None else donor_l_seq
        tag = f"{hrep.input_chain.name}_{suffix}"
        model, chains = _predict(vh_seq, vl_seq, tag)
        if model is None:
            continue
        if hv is not None and dh is not None:
            lab = match_pdb_chain(chains, hv.sequence)
            if lab is not None:
                hrep.structure_validation[hv.name] = _measure_variant(
                    dh, config, hrep.input_chain.numbered, hv, model, lab, "")
        if lv is not None and dl is not None:
            lab = match_pdb_chain(chains, lv.sequence)
            if lab is not None:
                lrep.structure_validation[lv.name] = _measure_variant(
                    dl, config, lrep.input_chain.numbered, lv, model, lab, "")


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


def _compute_hints_with_model(model, label, all_pos, cdrs, ag_chains,
                              pdb_path=None, use_plddt=True, cdr_loops=None):
    from .structure import compute_hints
    return compute_hints(model, label, all_pos, cdrs, ag_chains,
                         pdb_path=pdb_path, use_plddt=use_plddt,
                         cdr_loops=cdr_loops)
