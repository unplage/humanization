#!/usr/bin/env python3
"""Antibody humanization pipeline CLI.

Usage:
  humanize run    --input seq.fasta [--outdir outputs] [--format fab|vhh|auto]
                  [--scheme kabat|chothia|abm|imgt] [--germline-dir DIR]
                  [--germline-strategy fr_best|cdr_best|composite|cvi_best|min_backmutations|current|auto]
                  [--af3-mode off|local|api] [--mpnn-mode off|local]
                  [--antigen SEQ] [--donor-structure PDB]
  humanize setup-germline [--dir DIR]      # download NCBI IgBLAST germline
  humanize setup-check                     # report tool availability
  humanize stability --pdb FILE            # predict ΔG stability
  humanize stability --wt FILE --mutant FILE  # predict ΔΔG
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from humanize.pipeline import PipelineConfig, run_pipeline
from humanize.report import write_all
from humanize.structure import AF3Config
from humanize.mpnn import MPNNConfig


def cmd_run(args):
    # Parse forced germlines
    forced_germlines = {}
    if args.force_germline:
        for item in args.force_germline:
            if "=" in item:
                chain, gene = item.split("=", 1)
                # Normalize chain key: VH->H, VL->L, H->H, L->L
                chain_key = chain.upper().replace("VH", "H").replace("VL", "L")
                forced_germlines[chain_key] = gene
    
    cfg = PipelineConfig(
        germline_dir=args.germline_dir or "",
        germline_strategy=args.germline_strategy or "auto",
        forced_germlines=forced_germlines,
        cdr_scheme=args.scheme,
        report_schemes=["kabat", "chothia", "abm", "imgt"],
        format=args.format,
        antigen_seq=args.antigen,
        donor_structure=args.donor_structure,
        calibration_path=args.calibration,
        biophi_env=args.biophi_env,
        oasis_db=args.oasis_db,
        af3=AF3Config(
            mode=args.af3_mode,
            binary=args.af3_binary or "",
            workdir=os.path.join(args.outdir, "af3"),
            api_url=args.af3_api or "",
            api_token=os.environ.get("AF3_TOKEN", ""),
        ),
        mpnn=MPNNConfig(
            mode=args.mpnn_mode,
            script=args.mpnn_script or "",
            outdir=os.path.join(args.outdir, "mpnn"),
        ),
        interactive_indel=getattr(args, 'interactive_indel', False),
        af3_validate_variants=getattr(args, 'af3_rmsd', False),
        design_panel=getattr(args, 'design_panel', False),
        design_panel_steps=getattr(args, 'design_panel_steps', 3),
        immunogenicity_path=getattr(args, 'immunogenicity_json', None),
        netmhciipan_binary=getattr(args, 'netmhciipan', None),
        netmhciipan_alleles=getattr(args, 'netmhciipan_alleles', "") or "",
        netmhciipan_env=getattr(args, 'netmhciipan_env', None),
    )
    
    # Determine which step is being run and create appropriate subdirectory
    has_structure = args.donor_structure is not None and args.donor_structure != ""
    has_mpnn = args.mpnn_mode != "off"
    
    if has_structure and has_mpnn:
        step_dir = os.path.join(args.outdir, "step4")
    elif has_structure:
        step_dir = os.path.join(args.outdir, "step3")
    else:
        step_dir = os.path.join(args.outdir, "step2")
    
    os.makedirs(step_dir, exist_ok=True)
    
    # Update MPNN outdir to use step subdirectory
    cfg.mpnn.outdir = os.path.join(step_dir, "mpnn")
    
    result = run_pipeline(args.input, cfg, outdir=step_dir)
    paths = write_all(step_dir, result)
    print(f"\n[humanize] format: {result.format.upper()}")
    print(f"[humanize] germline strategy: {args.germline_strategy}")
    for rep in result.chains:
        c = rep.input_chain
        v = rep.germline.v_gene
        j = rep.germline.j_gene
        s = rep.germline.scores
        print(f"\n  {c.name} ({c.chain_type}{'|VHH' if c.is_vhh else ''})")
        print(f"    germline: {v.gene_id if v else '?'} + {j.gene_id if j else '?'} "
              f"(FR id {s.get('fr_identity')}, CDR id {s.get('cdr_identity')})")
        tiers = {}
        for b in rep.backmut.candidates:
            tiers[b.tier] = tiers.get(b.tier, 0) + 1
        print(f"    back-mutations: {sum(tiers.values())} "
              f"({', '.join(f'{k}:{v}' for k, v in sorted(tiers.items()))})")
        hl = rep.human_likeness
        if hl:
            print(f"    human-likeness (graft): {', '.join(f'{k} {v}%' for k, v in hl.items())}")
        for v in rep.variants:
            print(f"    {v.name}: {len(v.backmutations)} back-mutations")
    if result.warnings:
        print("\n  warnings:")
        for w in result.warnings:
            print(f"    - {w}")
    print("\n[humanize] outputs:")
    for k, p in paths.items():
        print(f"  {k}: {p}")
    return 0


def cmd_rmsd(args):
    """Superpose externally predicted variant structures on a donor structure
    and compute framework-superposed CDR CA-RMSD.

    The donor reference comes from ``--donor``; the variant structures may be
    any PDB/CIF that contains the chain being evaluated (they do NOT need to
    contain a donor chain). The donor sequence/numbering comes from ``--input``.
    """
    from humanize.sequences import parse_input
    from humanize.numbering import number_heavy, number_light
    from humanize.structure import (
        _AA3TO1, count_clashes, load_model, match_pdb_chain, structure_rmsd,
    )
    from humanize.pipeline import build_cdr_framework_sets, pdb_chain_pos_map

    _fmt, chains = parse_input(args.input, args.format)
    ctype = args.chain
    donor_chain = next((c for c in chains if c.chain_type == ctype), None)
    if donor_chain is None or donor_chain.numbered is None:
        print(f"[humanize] no {ctype} chain in {args.input}", file=sys.stderr)
        return 1
    donor = donor_chain.numbered

    dmodel = load_model(args.donor)
    if dmodel is None:
        print(f"[humanize] could not read donor structure {args.donor}",
              file=sys.stderr)
        return 1
    dchains = {}
    for a in dmodel.atoms:
        dchains.setdefault(a.chain, []).append(a)
    dlabel = args.donor_chain or match_pdb_chain(dchains, donor.sequence)
    if dlabel is None or dlabel not in dchains:
        print("[humanize] could not match the donor chain by sequence; "
              "pass --donor-chain", file=sys.stderr)
        return 1
    donor_pos = pdb_chain_pos_map(dmodel, dlabel, donor)
    cdr_sets, framework = build_cdr_framework_sets(donor, args.scheme)

    results = {}
    for vpath in args.variants:
        vmodel = load_model(vpath)
        if vmodel is None:
            print(f"[humanize] could not read {vpath}", file=sys.stderr)
            continue
        vchains = {}
        for a in vmodel.atoms:
            vchains.setdefault(a.chain, []).append(a)
        vlabel = args.variant_chain or match_pdb_chain(vchains, donor.sequence)
        if vlabel is None or vlabel not in vchains:
            print(f"[humanize] could not match a chain in {vpath}",
                  file=sys.stderr)
            continue
        ca = sorted((a.resseq, a.resname)
                    for a in vchains[vlabel] if a.name == "CA")
        vseq = "".join(_AA3TO1.get(rn, "X") for _, rn in ca)
        try:
            vnum = number_heavy(vseq) if ctype == "H" else number_light(vseq)
        except ValueError as e:
            print(f"[humanize] numbering failed for {vpath}: {e}",
                  file=sys.stderr)
            continue
        vpos = pdb_chain_pos_map(vmodel, vlabel, vnum)
        res = structure_rmsd(
            dmodel, dlabel, donor_pos, vmodel, vlabel, vpos, cdr_sets, framework)
        clashes = count_clashes(vmodel)
        res["n_clashes"] = clashes["n_clashes"]
        res["worst_clash"] = clashes["worst"]
        by_res = {a.resseq: a.plddt for a in vmodel.atoms
                  if a.chain == vlabel and a.name == "CA" and a.plddt > 0}
        cdr_all = set().union(*cdr_sets.values()) if cdr_sets else set()
        cdr_vals = [by_res[vpos[p]] for p in cdr_all
                    if p in vpos and vpos[p] in by_res]
        res["cdr_plddt"] = (round(sum(cdr_vals) / len(cdr_vals), 1)
                            if cdr_vals else None)
        res["variant_chain"] = vlabel
        results[os.path.basename(vpath)] = res

    if not results:
        print("[humanize] no variant structures evaluated", file=sys.stderr)
        return 1

    def _f(x):
        return f"{x:.2f}" if isinstance(x, (int, float)) else "-"

    print(f"[humanize] CDR-RMSD vs donor {os.path.basename(args.donor)} "
          f"(chain {dlabel}, {args.scheme}, framework-superposed)")
    print(f"{'variant':<28} {'CDR-RMSD':>9} {'FR-RMSD':>8} "
          f"{'CDR1':>6} {'CDR2':>6} {'CDR3':>6} {'CDR pLDDT':>10} "
          f"{'n_cdr':>6} {'clashes':>8} {'worst':>6}")
    for name, r in results.items():
        per = r.get("per_cdr") or {}
        print(f"{name:<28} {_f(r.get('cdr_rmsd')):>9} {_f(r.get('fr_rmsd')):>8} "
              f"{_f(per.get('CDR1')):>6} {_f(per.get('CDR2')):>6} "
              f"{_f(per.get('CDR3')):>6} {_f(r.get('cdr_plddt')):>10} "
              f"{r.get('n_cdr', 0):>6} "
              f"{r.get('n_clashes', 0):>8} {_f(r.get('worst_clash')):>6}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(results, fh, indent=2)
        print(f"\n[humanize] wrote {args.out}")
    return 0


def cmd_setup_germline(args):
    from humanize.germline import download_germline_db
    d = args.dir or os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "germline")
    print(f"[humanize] downloading NCBI IgBLAST germline into {d} ...")
    db = download_germline_db(d)
    print(f"[humanize] done: {len(db.v_for('H'))} IGHV, "
          f"{len(db.v_for('L'))} IG[KL]V, {len(db.j_for('H'))} IGHJ, "
          f"{len(db.j_for('L'))} IG[KL]J genes")
    return 0


def cmd_learn(args):
    from humanize.learning import (
        compute_position_effects,
        parse_experiments,
        write_calibration,
    )
    from humanize.germline import load_germline_db
    db = load_germline_db(args.germline_dir or os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "..", "data", "germline"))
    records = parse_experiments(args.experiments)
    effects, warnings = compute_position_effects(records, db)
    write_calibration(args.out, effects, meta={
        "n_experiments": len(records),
        "n_positions": len(effects),
    })
    for w in warnings:
        print(f"[humanize] warn: {w}")
    print(f"[humanize] learned effects for {len(effects)} framework positions "
          f"from {len(records)} experiments -> {args.out}")
    n_sig = sum(1 for e in effects.values() if abs(e.effect) >= 0.20)
    print(f"[humanize] positions with |ddG| >= 0.20 kcal/mol: {n_sig}")
    return 0


def cmd_stability(args):
    """Predict protein stability (ΔG/ΔΔG) for antibody structures."""
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools"))
    
    from stability_predictor.stability_analyzer import (
        predict_dg_esm3dg, predict_dg_saprodg, predict_ddg, analyze_batch,
        format_result_text, format_batch_text, format_result_json, format_batch_json,
        check_esm3dg_available, check_saprodg_available
    )
    
    # Check model availability
    print("[stability] Checking model availability...")
    esm3d_available = check_esm3dg_available()
    saprodg_available = check_saprodg_available()
    
    print(f"  ESM3dG: {'available' if esm3d_available else 'NOT INSTALLED'}")
    print(f"  SaProtΔG: {'available' if saprodg_available else 'NOT INSTALLED'}")
    
    if not esm3d_available and not saprodg_available:
        print("\n[stability] ERROR: No prediction models available!")
        print("  Install ESM3dG (recommended):")
        print("    pip install git+https://github.com/yehlincho/absolute-stability-predictor.git")
        return 1
    
    # Determine model
    model = args.model
    if model == "ESM3dG" and not esm3d_available:
        if saprodg_available:
            print("[stability] ESM3dG not available, using SaProtΔG")
            model = "SaProtΔG"
        else:
            print("[stability] ERROR: No models available")
            return 1
    elif model == "SaProtΔG" and not saprodg_available:
        if esm3d_available:
            print("[stability] SaProtΔG not available, using ESM3dG")
            model = "ESM3dG"
        else:
            print("[stability] ERROR: No models available")
            return 1
    
    # Run analysis
    if args.wt and args.mutant:
        print(f"\n[stability] Predicting ΔΔG: {os.path.basename(args.wt)} -> {os.path.basename(args.mutant)}")
        result = predict_ddg(args.wt, args.mutant, model)
        if args.json:
            print(json.dumps(format_result_json(result), indent=2))
        else:
            print(format_result_text(result))
    elif args.pdb:
        print(f"\n[stability] Analyzing: {os.path.basename(args.pdb)}")
        if model == "ESM3dG":
            result = predict_dg_esm3dg(args.pdb, args.chain)
        else:
            result = predict_dg_saprodg(args.pdb, args.chain)
        if args.json:
            print(json.dumps(format_result_json(result), indent=2))
        else:
            print(format_result_text(result))
    elif args.pdb_dir:
        print(f"\n[stability] Batch analyzing: {args.pdb_dir}")
        batch_result = analyze_batch(args.pdb_dir, model)
        if args.json:
            print(json.dumps(format_batch_json(batch_result), indent=2))
        else:
            print(format_batch_text(batch_result))
    
    # Save results
    if args.output:
        os.makedirs(args.output, exist_ok=True)
        if args.pdb:
            result = predict_dg_esm3dg(args.pdb, args.chain) if model == "ESM3dG" else predict_dg_saprodg(args.pdb, args.chain)
            data = format_result_json(result)
            output_file = os.path.join(args.output, "stability_result.json")
        elif args.wt and args.mutant:
            result = predict_ddg(args.wt, args.mutant, model)
            data = format_result_json(result)
            output_file = os.path.join(args.output, "stability_result.json")
        else:
            batch_result = analyze_batch(args.pdb_dir, model)
            data = format_batch_json(batch_result)
            output_file = os.path.join(args.output, "stability_batch_results.json")
        
        with open(output_file, "w") as f:
            json.dump(data, f, indent=2)
        print(f"\n[stability] Results saved to: {output_file}")
    
    return 0


def cmd_setup_check(args):
    checks = [
        ("python3", sys.executable),
        ("biopython", None),
        ("anarci", None),
        ("igblastn", shutil.which("igblastn")),
        ("protein_mpnn.py", shutil.which("protein_mpnn.py")),
        ("run_alphafold.py", shutil.which("run_alphafold.py")),
    ]
    print("[humanize] setup-check")
    for name, path in checks:
        if name == "biopython":
            try:
                import Bio  # type: ignore
                path = Bio.__version__
            except ImportError:
                path = None
        if name == "anarci":
            try:
                import anarci  # type: ignore
                path = "python module OK"
            except ImportError:
                path = None
        print(f"  {name:20s}: {path or 'NOT FOUND (optional)'}")
    print("\n  portable mode (no external tools) is fully functional.")
    return 0


def cmd_compare(args):
    """Lightweight germline comparison across 9 strategies (Step 1).
    
    This command only evaluates germline candidates without running the full
    humanization pipeline. It displays a comparison table showing the top 5
    candidates for each of the 9 strategies.
    """
    from humanize.sequences import parse_input
    from humanize.numbering import number_heavy, number_light
    from humanize.germline import load_germline_db
    from humanize.multi_strategy_germline import choose_germlines_multi_strategy, format_multi_strategy_report
    
    # Load germline database
    germline_dir = args.germline_dir or os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "data", "germline"
    )
    db = load_germline_db(germline_dir)
    
    # Parse input sequences
    fmt, chains = parse_input(args.input, args.format)
    
    print(f"\n{'='*100}")
    print(f"Germline Comparison Report (Step 1)")
    print(f"{'='*100}\n")
    
    for chain in chains:
        ctype = chain.chain_type
        if ctype == 'H':
            numbered = number_heavy(chain.sequence)
        else:
            numbered = number_light(chain.sequence)
        
        # Run multi-strategy germline selection (single pass, all 9 strategies)
        result = choose_germlines_multi_strategy(numbered, db, is_vhh=chain.is_vhh)
        
        # Display comparison table
        print(format_multi_strategy_report(result))
    
    print(f"\n{'='*100}")
    print("Next step: Run full pipeline with chosen germline strategy:")
    print("  python -m scripts.humanize.cli run --input <fasta> --germline-strategy <strategy>")
    print("  or force specific germline:")
    print("  python -m scripts.humanize.cli run --input <fasta> --force-germline VH=<gene> VL=<gene>")
    print(f"{'='*100}\n")
    
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="humanize", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="run the humanization pipeline")
    p_run.add_argument("--input", required=True, help="input FASTA (VH+VL or VHH)")
    p_run.add_argument("--outdir", default="outputs")
    p_run.add_argument("--format", default="auto", choices=["auto", "fab", "vhh"])
    p_run.add_argument("--scheme", default="kabat",
                       choices=["kabat", "chothia", "abm", "imgt"],
                       help="CDR definition used for the variant ladder")
    p_run.add_argument("--germline-dir", default="", help="NCBI germline FASTA dir")
    p_run.add_argument("--germline-strategy", default="auto",
                       choices=["fr_best", "cdr_best", "composite", "cvi_best", 
                                "min_backmutations", "current", "auto",
                                "adimab_frequency", "pioneer_frequency", "composite_3axis"],
                       help="germline selection strategy: "
                            "fr_best=FR identity, cdr_best=CDR identity, "
                            "composite=0.7*FR+0.3*CDR, cvi_best=CVI homology, "
                            "min_backmutations=fewest back-mutations, "
                            "current=current system (top 30%% FR with max CDR), "
                             "auto=adimab_frequency for VH, current for VL, "
                            "adimab_frequency=Adimab recommended + frequency, "
                            "pioneer_frequency=Pioneer library + frequency, "
                            "composite_3axis=0.5*CVI+0.3*freq+0.2*FR")
    p_run.add_argument("--force-germline", nargs="+", metavar="CHAIN=GENE",
                       help="force specific germline(s), e.g. --force-germline VH=IGHV1-3*01 VL=IGKV1-39*01")
    p_run.add_argument("--antigen", default=None, help="antigen sequence (AF3 complex)")
    p_run.add_argument("--calibration", default=None,
                       help="calibration.json from `humanize learn` (empirical scoring)")
    p_run.add_argument("--biophi-env", default=None,
                       help="conda env name containing biophi (Sapiens humanness cross-check)")
    p_run.add_argument("--oasis-db", default=None,
                       help="OASis 9-mer DB path (biophi oasis identity)")
    p_run.add_argument("--donor-structure", default=None, help="donor PDB/CIF")
    p_run.add_argument("--af3-mode", default="off", choices=["off", "local", "api"])
    p_run.add_argument("--af3-binary", default="", help="path to run_alphafold.py")
    p_run.add_argument("--af3-api", default="", help="AF3 API base URL")
    p_run.add_argument("--mpnn-mode", default="off", choices=["off", "local"])
    p_run.add_argument("--mpnn-script", default="", help="path to protein_mpnn.py")
    p_run.add_argument("--interactive-indel", action="store_true",
                       help="enable interactive FR indel selection (VH+VL)")
    p_run.add_argument("--af3-rmsd", action="store_true",
                       help="predict each variant with AF3 and compute CDR-RMSD "
                            "vs the donor model (requires --af3-mode)")
    p_run.add_argument("--design-panel", action="store_true",
                       help="emit a structure-guided V_opt panel of variants")
    p_run.add_argument("--design-panel-steps", type=int, default=3,
                       help="number of V_opt panel variants (default 3)")
    p_run.add_argument("--immunogenicity-json", default=None,
                       help="precomputed {kabat_position: score} MHC-II epitope map")
    p_run.add_argument("--netmhciipan", default=None,
                       help="NetMHCIIpan executable (per-position epitope scores)")
    p_run.add_argument("--netmhciipan-alleles", default="",
                       help="comma-separated HLA alleles for NetMHCIIpan")
    p_run.add_argument("--netmhciipan-env", default=None,
                       help="conda env holding NetMHCIIpan")
    p_run.set_defaults(func=cmd_run)

    # ---- compare: lightweight germline evaluation (Step 1) ----
    p_compare = sub.add_parser("compare", help="compare germline candidates across 9 strategies (Step 1)")
    p_compare.add_argument("--input", required=True, help="input FASTA (VH+VL or VHH)")
    p_compare.add_argument("--format", default="auto", choices=["auto", "fab", "vhh"])
    p_compare.add_argument("--germline-dir", default="", help="NCBI germline FASTA dir")
    p_compare.set_defaults(func=cmd_compare)

    # ---- rmsd: compare external variant structures to a donor structure ----
    p_rmsd = sub.add_parser(
        "rmsd", help="CDR-RMSD of variant structures vs a donor structure")
    p_rmsd.add_argument("--input", required=True,
                        help="donor FASTA (for sequence/numbering)")
    p_rmsd.add_argument("--donor", required=True,
                        help="donor reference structure (PDB/CIF)")
    p_rmsd.add_argument("--variants", nargs="+", required=True,
                        help="variant structures (PDB/CIF), any chain composition")
    p_rmsd.add_argument("--chain", default="H", choices=["H", "L"],
                        help="chain type being evaluated")
    p_rmsd.add_argument("--scheme", default="kabat",
                        choices=["kabat", "chothia", "abm", "imgt"])
    p_rmsd.add_argument("--format", default="auto", choices=["auto", "fab", "vhh"])
    p_rmsd.add_argument("--donor-chain", default=None,
                        help="force the donor PDB chain id")
    p_rmsd.add_argument("--variant-chain", default=None,
                        help="force the variant PDB chain id")
    p_rmsd.add_argument("--out", default=None, help="write results JSON here")
    p_rmsd.set_defaults(func=cmd_rmsd)

    p_g = sub.add_parser("setup-germline", help="download NCBI IgBLAST germline")
    p_g.add_argument("--dir", default="")
    p_g.set_defaults(func=cmd_setup_germline)

    p_c = sub.add_parser("setup-check", help="report available tools")
    p_c.set_defaults(func=cmd_setup_check)

    p_l = sub.add_parser(
        "learn", help="fit empirical position effects from experiment data")
    p_l.add_argument("--experiments", required=True,
                     help="experiments JSON (parent + variants + KD)")
    p_l.add_argument("--out", default="calibration.json",
                     help="output calibration file")
    p_l.add_argument("--germline-dir", default="")
    p_l.set_defaults(func=cmd_learn)

    # ---- stability: predict ΔG/ΔΔG for antibody structures ----
    p_stab = sub.add_parser(
        "stability", help="predict protein stability (ΔG/ΔΔG)")
    p_stab.add_argument("--pdb", help="single PDB file to analyze (ΔG)")
    p_stab.add_argument("--pdb-dir", help="directory of PDB files for batch analysis")
    p_stab.add_argument("--wt", help="wildtype PDB file (for ΔΔG)")
    p_stab.add_argument("--mutant", help="mutant PDB file (for ΔΔG)")
    p_stab.add_argument("--model", default="ESM3dG", choices=["ESM3dG", "SaProtΔG"],
                        help="prediction model (default: ESM3dG)")
    p_stab.add_argument("--chain", default="A", help="chain ID (default: A)")
    p_stab.add_argument("--output", help="output directory for results")
    p_stab.add_argument("--json", action="store_true", help="output JSON format")
    p_stab.set_defaults(func=cmd_stability)

    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except RuntimeError as e:
        print(f"[humanize] error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
