#!/usr/bin/env python3
"""Tests for the stability_predictor ThermoMPNN analyzer.

Covers the fixes that matter for back-mutation decisions:
  * ThermoMPNN 0-based sequence positions -> true Kabat positions
  * gapped/missing residue handling
  * --positions means true Kabat positions
  * targeted human->donor / donor->human back-mutation evaluation
  * invocation helpers (direct script, conda-run splitting, checkpoint choice)
"""

import csv
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from tools.stability_predictor.thermompnn_analyzer import (  # noqa: E402
    BackMutationCandidate,
    MutationResult,
    StabilityReport,
    _python_parts,
    _result_for_pending,
    _summary_backmutation,
    _prepare_backmutation_targets,
    apply_annotations,
    build_position_table,
    find_checkpoint,
    format_text_report,
    load_backmutation_candidates,
    parse_positions,
    parse_thermompnn_csv,
    normalize_positions,
    resolve_kabat_positions,
    save_markdown_report,
)

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


THREE = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY",
    "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN",
    "P": "PRO", "Q": "GLN", "R": "ARG", "S": "SER", "T": "THR", "V": "VAL",
    "W": "TRP", "Y": "TYR",
}

VH_4D5 = (
    "EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTKYNQKFKG"
    "KATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS"
)
VH_AMG110 = (
    "QVQLLEQSGAEVVKPGSSVKISCKASGYAFTNYWLGWVKQAPGQGLEWIGDIFPGSGNIHYNEKFKG"
    "KATLTADKSTSTAYMELSSLRSEDTAVYFCARLRNWDEPMDYWGQGTTVTVSS"
)


def write_pdb(path, seq, resseqs=None, chain="H"):
    """Write a minimal CA-only PDB for mapping tests."""
    if resseqs is None:
        resseqs = list(range(1, len(seq) + 1))
    serial = 1
    with open(path, "w") as fh:
        for aa, resseq in zip(seq, resseqs):
            fh.write(
                f"ATOM  {serial:5d}  CA  {THREE[aa]:>3s} {chain}{resseq:4d}"
                f"    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n"
            )
            serial += 1


def temp_pdb(seq, resseqs=None, chain="H"):
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    write_pdb(path, seq, resseqs, chain)
    return path


def test_parse_positions():
    check("empty positions", parse_positions("") == [])
    check("numeric -> chain type",
          parse_positions("37,44,45", "H") == ["H37", "H44", "H45"])
    check("explicit prefix kept",
          parse_positions("H37,L24,37", "H") == ["H37", "L24", "H37"])
    check("light chain default",
          parse_positions("24,27", "L") == ["L24", "L27"])
    check("insertion code kept",
          parse_positions("H35A", "H") == ["H35A"])
    check("normalize per chain",
          normalize_positions(["37", "L24"], "L") == ["L37", "L24"])


def test_parse_csv_uses_sequence_positions():
    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "Model", "Dataset", "ddG_pred", "position", "wildtype",
            "mutation", "pdb", "chain",
        ])
        w.writeheader()
        w.writerow({"Model": "ThermoMPNN", "Dataset": "x", "ddG_pred": "-0.8",
                    "position": "36", "wildtype": "V", "mutation": "F",
                    "pdb": "x", "chain": "H"})
        w.writerow({"Model": "ThermoMPNN", "Dataset": "x", "ddG_pred": "0.0",
                    "position": "37", "wildtype": "A", "mutation": "A",
                    "pdb": "x", "chain": "H"})
        w.writerow({"Model": "ThermoMPNN", "Dataset": "x", "ddG_pred": "0.4",
                    "position": "68", "wildtype": "A", "mutation": "V",
                    "pdb": "x", "chain": "H"})
    try:
        rows = parse_thermompnn_csv(path, "x.pdb", "H")
        # Self-mutation (A->A) is dropped.
        check("parse csv count", len(rows) == 2, str(len(rows)))
        check("position is 0-based seq index", rows[0].position_seq == 36)
        check("raw parse has no Kabat label", rows[0].position_kabat == "")
        check("ddG parsed", rows[0].ddG == -0.8)
        check("stabilizing flag", rows[0].is_stabilizing is True)
    finally:
        os.unlink(path)


def test_4d5_mapping_off_by_one_fixed():
    path = temp_pdb(VH_4D5)
    try:
        cid, ctype, table, warnings, ungapped = build_position_table(path, "H", "H")
        by_pos = {cp.kabat_pos: cp for cp in table if cp is not None and cp.kabat_pos}
        check("chain resolved", cid == "H")
        check("chain type H", ctype == "H")
        check("ungapped sequence preserved", ungapped == VH_4D5)
        check("4D5 H37 -> seq_pos 36 (not 37)",
              by_pos["H37"].seq_pos == 36, str(by_pos["H37"].seq_pos))
        check("4D5 H37 PDB resseq 37", by_pos["H37"].resseq == 37)
        check("4D5 H37 region FR2", by_pos["H37"].kabat_region == "FR2")
        # A row at ThermoMPNN position 36 must become H37, not H36.
        rows = [MutationResult(
            pdb_file="x", chain="H", position_kabat="", wildtype="V",
            mutation="F", ddG=-0.5, position_seq=36,
        )]
        apply_annotations(rows, table, "H")
        check("position 36 annotated as H37", rows[0].position_kabat == "H37")
        check("position 36 annotated FR2", rows[0].region_kabat == "FR2")
        indices, missing = resolve_kabat_positions(table, ["H37"], "H")
        check("--positions H37 selects index 36", indices == {36}, str(indices))
        check("--positions H37 not missing", missing == [])
    finally:
        os.unlink(path)


def test_gapped_structure_mapping():
    # Renumber 4D5 with residue 10 missing: sequence is unchanged, but
    # ThermoMPNN inserts a gap at seq_pos 9 and all later seq_pos shift by 1.
    resseqs = list(range(1, 10)) + list(range(11, 11 + (len(VH_4D5) - 9)))
    assert len(resseqs) == len(VH_4D5), (len(resseqs), len(VH_4D5))
    path = temp_pdb(VH_4D5, resseqs=resseqs)
    try:
        _, _, table, _, _ = build_position_table(path, "H", "H")
        by_pos = {cp.kabat_pos: cp for cp in table if cp is not None and cp.kabat_pos}
        check("gap represented at seq_pos 9", table[9] is None)
        check("gap shifts H37 to seq_pos 37",
              by_pos["H37"].seq_pos == 37, str(by_pos["H37"].seq_pos))
        check("gap shifts H37 PDB resseq to 38",
              by_pos["H37"].resseq == 38, str(by_pos["H37"].resseq))
        check("H37 still FR2 after gap", by_pos["H37"].kabat_region == "FR2")
    finally:
        os.unlink(path)


def test_amg110_region_mapping_fixed():
    path = temp_pdb(VH_AMG110)
    try:
        _, _, table, _, _ = build_position_table(path, "H", "H")
        by_pos = {cp.kabat_pos: cp for cp in table if cp is not None and cp.kabat_pos}
        check("AMG110 H67 -> seq_pos 68",
              by_pos["H67"].seq_pos == 68, str(by_pos["H67"].seq_pos))
        check("AMG110 H67 residue A", by_pos["H67"].aa == "A")
        check("AMG110 H67 region FR3", by_pos["H67"].kabat_region == "FR3")
        check("AMG110 H66 -> seq_pos 67",
              by_pos["H66"].seq_pos == 67, str(by_pos["H66"].seq_pos))
        # The old bug selected index 67 and then reported region from index 66.
        rows = [MutationResult(
            pdb_file="x", chain="H", position_kabat="", wildtype="K",
            mutation="R", ddG=0.1, position_seq=67,
        )]
        apply_annotations(rows, table, "H")
        check("index 67 annotated H66 not H67", rows[0].position_kabat == "H66")
        check("index 67 annotated FR3 not CDR2", rows[0].region_kabat == "FR3")
        indices, _ = resolve_kabat_positions(table, ["H67"], "H")
        check("--positions H67 selects index 68", indices == {68}, str(indices))
    finally:
        os.unlink(path)


def test_targeted_backmutation_direction():
    # V0-like: position H67 carries the human residue V.
    v0 = list(VH_AMG110)
    v0[68] = "V"
    path = temp_pdb("".join(v0))
    try:
        _, _, table, _, _ = build_position_table(path, "H", "H")
        cands = [BackMutationCandidate(
            position="H67", donor_aa="A", human_aa="V", tier="T1",
        )]
        pending, skipped = _prepare_backmutation_targets(table, cands, "H")
        check("V0 candidate pending", len(pending) == 1 and not skipped)
        p = pending[0]
        check("V0 direction human_to_donor", p.direction == "human_to_donor")
        check("V0 target is donor A", p.target_aa == "A")
        row = MutationResult(
            pdb_file="x", chain="H", position_kabat="H67", wildtype="V",
            mutation="A", ddG=-0.8, position_seq=68,
        )
        res = _result_for_pending(p, row, "x.pdb", "H")
        check("V0 ddG_backmutation = ddG", res.ddG_backmutation == -0.8)
        check("V0 verdict supports", res.verdict == "supports back-mutation")
        check("V0 recommendation keep", res.recommendation == "keep")
    finally:
        os.unlink(path)

    # V2-like: position H67 already carries the donor residue A.
    path = temp_pdb(VH_AMG110)
    try:
        _, _, table, _, _ = build_position_table(path, "H", "H")
        cands = [BackMutationCandidate(
            position="H67", donor_aa="A", human_aa="V", tier="T1",
        )]
        pending, _ = _prepare_backmutation_targets(table, cands, "H")
        p = pending[0]
        check("V2 direction donor_to_human", p.direction == "donor_to_human")
        check("V2 target is human V", p.target_aa == "V")
        row = MutationResult(
            pdb_file="x", chain="H", position_kabat="H67", wildtype="A",
            mutation="V", ddG=0.5, position_seq=68,
        )
        res = _result_for_pending(p, row, "x.pdb", "H")
        # Removing a back-mutation by 0.5 kcal/mol means the back-mutation
        # itself was stabilizing by ~0.5 kcal/mol.
        check("V2 sign inverted for back-mutation effect",
              res.ddG_backmutation == -0.5, str(res.ddG_backmutation))
        check("V2 verdict supports", res.verdict == "supports back-mutation")
        summary = _summary_backmutation([res])
        check("targeted summary supports", summary["supports"] == 1)
    finally:
        os.unlink(path)


def test_targeted_skips_indel_and_cdr():
    path = temp_pdb(VH_AMG110)
    try:
        _, _, table, _, _ = build_position_table(path, "H", "H")
        cands = [
            BackMutationCandidate(position="H6", donor_aa="E", human_aa="-"),
            BackMutationCandidate(position="H999", donor_aa="A", human_aa="V"),
        ]
        pending, skipped = _prepare_backmutation_targets(table, cands, "H")
        check("no pending for invalid candidates", pending == [])
        reasons = {c.position: r for c, r in skipped}
        check("FR indel skipped", "FR indel" in reasons["H6"])
        check("missing position skipped", "not found" in reasons["H999"])
    finally:
        os.unlink(path)


def test_load_candidates_and_tiers():
    payload = {"chains": [{
        "name": "X.H", "chain_type": "H",
        "backmutations": [
            {"position": "H5", "donor_aa": "L", "human_aa": "V", "tier": "T1"},
            {"position": "H6", "donor_aa": "E", "human_aa": "-", "tier": "T3"},
        ],
    }]}
    fd, path = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    with open(path, "w") as fh:
        json.dump(payload, fh)
    try:
        all_c = load_backmutation_candidates(pipeline_json=path)
        t12 = load_backmutation_candidates(pipeline_json=path, tiers={"T1", "T2"})
        check("load all candidates", len(all_c) == 2)
        check("tier filter", len(t12) == 1 and t12[0].position == "H5")
    finally:
        os.unlink(path)

    fd, path = tempfile.mkstemp(suffix=".csv")
    os.close(fd)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "position", "donor_aa", "human_aa", "tier", "composite",
        ])
        w.writeheader()
        w.writerow({"position": "L1", "donor_aa": "E", "human_aa": "D",
                    "tier": "T2", "composite": "74.0"})
    try:
        c = load_backmutation_candidates(backmutations_csv=path)
        check("csv candidate loaded", len(c) == 1 and c[0].position == "L1")
        check("csv composite parsed", c[0].composite == 74.0)
    finally:
        os.unlink(path)


def test_report_formatting():
    rows = [
        MutationResult(
            pdb_file="x.pdb", chain="H", position_kabat="H67",
            wildtype="A", mutation="V", ddG=0.5, region_kabat="FR3",
            position_seq=68, position_pdb=68, tier="T1",
            donor_aa="A", human_aa="V", direction="donor_to_human",
            ddG_backmutation=-0.5, verdict="supports back-mutation",
            recommendation="keep",
        ),
    ]
    report = StabilityReport(
        pdb_file="x.pdb", chain="H", chain_type="H", mode="backmutation",
        context="back-mutated context (V1/V2/V3-like)", mutations=rows,
        summary=_summary_backmutation(rows),
    )
    text = format_text_report(report)
    check("targeted text title", "Back-mutation Stability Analysis" in text)
    check("targeted text has Kabat", "H67" in text)
    check("targeted text has dG_bm", "dG_bm" in text)
    check("targeted text has verdict", "supports back-mutation" in text)

    fd, path = tempfile.mkstemp(suffix=".md")
    os.close(fd)
    try:
        save_markdown_report(report, path)
        with open(path) as fh:
            md = fh.read()
        check("markdown title", "# Back-mutation Stability Analysis" in md)
        check("markdown table row", "H67" in md and "-0.50" in md)
    finally:
        os.unlink(path)


def test_invocation_helpers():
    parts = _python_parts("conda run -n thermoMPNN python")
    check("conda run split", parts[:2] == ["conda", "run"], str(parts))
    parts2 = _python_parts("~/miniconda3/bin/python")
    check("tilde expanded when path exists",
          parts2 == [os.path.expanduser("~/miniconda3/bin/python")], str(parts2))

    from pathlib import Path
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "models").mkdir()
        (root / "models" / "thermoMPNN_default.pt").write_text("x")
        (root / "vanilla_model_weights").mkdir()
        (root / "vanilla_model_weights" / "v_48_020.pt").write_text("x")
        ck = find_checkpoint(tmp)
        check("prefers thermoMPNN checkpoint",
              ck is not None and ck.name == "thermoMPNN_default.pt")
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "vanilla_model_weights").mkdir()
        (root / "vanilla_model_weights" / "v_48_020.pt").write_text("x")
        check("never uses vanilla weights as stability model",
              find_checkpoint(tmp) is None)


def test_end_to_end_with_fake_thermompnn():
    """Exercise the full analyze_structure path without a real ThermoMPNN env."""
    from pathlib import Path
    import stat

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        repo = root / "repo"
        (repo / "analysis").mkdir(parents=True)
        (repo / "models").mkdir()
        (repo / "analysis" / "custom_inference.py").write_text("")
        (repo / "local.yaml").write_text("")
        (repo / "models" / "thermoMPNN_default.pt").write_text("")

        fake_py = root / "fake_python"
        fake_py.write_text(
            "#!/usr/bin/env python3\n"
            "import sys, os, json, csv\n"
            "args = sys.argv[1:]\n"
            "if args and args[0] == '-c':\n"
            "    print('OK'); sys.exit(0)\n"
            "def val(flag):\n"
            "    return args[args.index(flag) + 1]\n"
            "out = val('--out-dir'); pdb = val('--pdb')\n"
            "targets = json.load(open(val('--targets'))) if '--targets' in args else None\n"
            "os.makedirs(out, exist_ok=True)\n"
            "rows = targets if targets else []\n"
            "path = os.path.join(out, 'ThermoMPNN_inference_%s.csv' % os.path.splitext(os.path.basename(pdb))[0])\n"
            "with open(path, 'w', newline='') as fh:\n"
            "    w = csv.writer(fh)\n"
            "    w.writerow(['Model','Dataset','ddG_pred','position','wildtype','mutation','pdb','chain'])\n"
            "    for wt, pos, mut in rows:\n"
            "        ddg = -0.8 if mut == 'A' else 0.5\n"
            "        w.writerow(['ThermoMPNN','x',ddg,pos,wt,mut,'x','H'])\n"
        )
        fake_py.chmod(fake_py.stat().st_mode | stat.S_IEXEC)

        v0 = list(VH_AMG110)
        v0[68] = "V"
        pdb = root / "v0.pdb"
        write_pdb(str(pdb), "".join(v0))

        from tools.stability_predictor.thermompnn_analyzer import analyze_structure
        cands = [BackMutationCandidate(
            position="H67", donor_aa="A", human_aa="V", tier="T1",
        )]
        report = analyze_structure(
            str(pdb), "H", str(repo), str(fake_py), str(root / "out"),
            positions=None, chain_type="H", candidates=cands, timeout=60,
        )
        check("e2e no error", report.error is None, str(report.error))
        check("e2e mode backmutation", report.mode == "backmutation")
        check("e2e one result", len(report.mutations) == 1, str(len(report.mutations)))
        if report.mutations:
            r = report.mutations[0]
            check("e2e kabat H67", r.position_kabat == "H67", r.position_kabat)
            check("e2e ddG_backmutation -0.8", r.ddG_backmutation == -0.8,
                  str(r.ddG_backmutation))
            check("e2e verdict supports", r.verdict == "supports back-mutation")
        target_file = root / "out" / "targeted_mutations.json"
        check("e2e targeted fast path used", target_file.exists())
        if target_file.exists():
            targets = json.loads(target_file.read_text())
            check("e2e target is V68A", targets == [["V", 68, "A"]], str(targets))

        # Generic --positions mode: 19 substitutions at true Kabat H67.
        report2 = analyze_structure(
            str(pdb), "H", str(repo), str(fake_py), str(root / "out2"),
            positions=["H67"], chain_type="H", candidates=None, timeout=60,
        )
        check("e2e ssm no error", report2.error is None, str(report2.error))
        check("e2e ssm 19 substitutions", len(report2.mutations) == 19,
              str(len(report2.mutations)))
        check("e2e ssm all H67",
              all(m.position_kabat == "H67" for m in report2.mutations))
        check("e2e ssm all FR3",
              all(m.region_kabat == "FR3" for m in report2.mutations))


def test_direct_script_entrypoint():
    script = os.path.join(
        ROOT, "tools", "stability_predictor", "stability_analyzer.py"
    )
    proc = subprocess.run(
        [sys.executable, script, "--help"],
        capture_output=True, text=True, timeout=60,
    )
    check("direct script --help exits 0", proc.returncode == 0, proc.stderr[-300:])
    check("direct script shows targeted options",
          "--pipeline-json" in proc.stdout and "--positions" in proc.stdout)


def run_all_tests():
    tests = [
        test_parse_positions,
        test_parse_csv_uses_sequence_positions,
        test_4d5_mapping_off_by_one_fixed,
        test_gapped_structure_mapping,
        test_amg110_region_mapping_fixed,
        test_targeted_backmutation_direction,
        test_targeted_skips_indel_and_cdr,
        test_load_candidates_and_tiers,
        test_report_formatting,
        test_invocation_helpers,
        test_end_to_end_with_fake_thermompnn,
        test_direct_script_entrypoint,
    ]
    for test in tests:
        print(f"\n{test.__name__}")
        try:
            test()
        except Exception as exc:  # noqa: BLE001
            import traceback
            traceback.print_exc()
            FAILURES.append(f"{test.__name__}: {exc}")

    print(f"\n{'='*60}")
    if FAILURES:
        print(f"FAILED ({len(FAILURES)}):")
        for f in FAILURES:
            print(f"  - {f}")
        return False
    print("All stability tests passed.")
    return True


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
