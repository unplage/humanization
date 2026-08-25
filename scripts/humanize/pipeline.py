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
) -> ChainReport:
    donor = chain.numbered
    if donor is None:
        raise RuntimeError(f"[{chain.name}] no numbering available")
    ctype = donor.chain_type
    is_vhh = chain.is_vhh and ctype == "H"

    # ---- germline selection ----
    from .multi_strategy_germline import choose_germlines_multi_strategy
    
    strategy = config.germline_strategy
    if strategy == "auto":
        # 自动策略：VH 使用 cvi_best，VL 使用 cdr_best
        strategy = "cvi_best" if ctype == "H" else "cdr_best"
    
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
        choice = GermlineChoice(
            v_gene=v_gene,
            j_gene=j_gene,
            scores={"fr_identity": candidate.fr_identity, "cdr_identity": candidate.cdr_identity},
            alternatives=[(c.gene, {"fr_identity": c.fr_identity, "cdr_identity": c.cdr_identity}) 
                         for c in multi_result.candidates.values()],
        )
    
    if v_gene is None or j_gene is None:
        raise RuntimeError(f"[{chain.name}] no viable human germline found")

    # ---- structure hints (AF3) ----
    hints = StructureHints()
    af3_pdb = None
    if config.af3.mode != "off":
        os.makedirs(config.af3.workdir, exist_ok=True)
        af3_pdb = predict_fv(
            config.af3,
            donor.sequence,
            None if fmt == "vhh" else _partner_sequence(chain, None),
            antigen,
            f"{chain.name}_{ctype}",
        )
    if af3_pdb and os.path.exists(af3_pdb):
        from .structure import load_model
        model = load_model(af3_pdb)
        if model:
            # map Kabat positions to model residue numbers (chain label H/L/A)
            label = "H" if ctype == "H" else "L"
            all_pos = {r.pos: r.index + 1 for r in donor.residues}
            cdrs = {p: n for p, n in all_pos.items() if (donor.region_of(p) or "").startswith("CDR")}
            ag_chains = ["A"] if antigen else None
            hints = _compute_hints_with_model(model, label, all_pos, cdrs, ag_chains)

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
            warnings = getattr(chain, "warnings", [])
            warnings.append(f"[{chain.name}] graft({scheme}) failed: {e}")

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
    )


def _partner_sequence(chain: InputChain, _all_chains) -> Optional[str]:
    return None  # Fv-only partner chain prediction handled at pipeline level


def _compute_hints_with_model(model, label, all_pos, cdrs, ag_chains):
    from .structure import compute_hints
    return compute_hints(model, label, all_pos, cdrs, ag_chains)
