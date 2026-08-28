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
                            "auto=cvi_best for VH, cdr_best for VL, "
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
    p_run.set_defaults(func=cmd_run)

    # ---- compare: lightweight germline evaluation (Step 1) ----
    p_compare = sub.add_parser("compare", help="compare germline candidates across 9 strategies (Step 1)")
    p_compare.add_argument("--input", required=True, help="input FASTA (VH+VL or VHH)")
    p_compare.add_argument("--format", default="auto", choices=["auto", "fab", "vhh"])
    p_compare.add_argument("--germline-dir", default="", help="NCBI germline FASTA dir")
    p_compare.set_defaults(func=cmd_compare)

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

    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except RuntimeError as e:
        print(f"[humanize] error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
