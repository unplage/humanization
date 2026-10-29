#!/usr/bin/env python3
"""Regression tests for tools/tap_profiler (plain script; exits nonzero on failure).

Covers the TAP rewrite:
  * published metric formulas (pairwise 1/r**2 PSH/PPC/PNC, SFvCSP as a
    *product* of chain charges, not a difference)
  * threshold flag logic for every bundled threshold set
  * neighbour distance / salt-bridge geometry on known coordinates
  * IMGT numbering consistency between structure and sequence modes
  * Fv scope excludes antigen chains
  * report renderers (text/JSON/Markdown/Word)
  * CLI --check and an end-to-end structure run
  * optional exact-parity check against the Exscientia reference numbers
    (set TAP_PARITY_PDB / TAP_PARITY_PSA to enable)

The end-to-end tests build a synthetic IMGT-ANARCI-numberable Fv from the
tracked 4D5 example sequences, so they need no external structure files.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "tap_profiler")
sys.path.insert(0, TOOL)

from tap_core import (  # noqa: E402
    METRIC_ORDER,
    NORMALISED_HYDROPHOBICITY,
    RESIDUE_CHARGES,
    THRESHOLD_SETS,
    Residue,
    analyze_sequences,
    analyze_structure,
    annotate_residues,
    build_pair_distances,
    compute_surface_scores,
    flag_value,
    load_residues,
    number_chain_imgt,
    risk_summary,
    write_work_pdb,
)

VH_4D5 = ("EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTK"
          "YNQKFKGKATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS")
VL_4D5 = ("DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGV"
          "PSRFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK")

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}

FAILURES = []


def check(condition, message):
    if not condition:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# Synthetic structure builder
# ---------------------------------------------------------------------------

def _write_synthetic_pdb(path, chains):
    """Write a deterministic toy structure for a list of (chain, sequence).

    Residues follow a serpentine path 3.8 A apart so every chain has real
    neighbours inside the 7.5 A TAP cutoff, and each residue carries a
    backbone plus CB around its CA. Coordinates depend only on (chain, index,
    residue) so repeated runs are bit-identical.
    """
    serial = 1
    with open(path, "w") as fh:
        for chain_offset, (chain_id, seq) in enumerate(chains):
            # distinct lateral offset per chain so different chains can never
            # generate coincident atoms (which would give zero distances)
            lateral = chain_offset * 5.1 + 0.37
            for i, aa in enumerate(seq, start=1):
                resname = ONE_TO_THREE[aa]
                x = (i % 20) * 3.8 + lateral
                y = (i // 20) * 3.8 + lateral
                z = ((i * 7 + chain_offset * 3 + 1) % 9) * 0.9 + chain_offset * 0.13
                atoms = (
                    (" CA ", 0.0, 0.0, 0.0, "C"),
                    (" N  ", -1.2, 0.5, 0.1, "N"),
                    (" C  ", 1.2, 0.5, -0.1, "C"),
                    (" O  ", 1.7, 1.6, -0.1, "O"),
                    (" CB ", 0.1, -1.5, 0.6, "C"),
                )
                for name, dx, dy, dz, element in atoms:
                    fh.write("ATOM  {:5d} {:4s} {:>3s} {:1s}{:4d}    "
                             "{:8.3f}{:8.3f}{:8.3f}  1.00 20.00          {:>2s}\n"
                             .format(serial, name, resname, chain_id, i,
                                     x + dx, y + dy, z + dz, element))
                    serial += 1
        fh.write("END\n")


def _residue(chain, resseq, resname, aa, coords):
    atoms = [(" CA ", "C", *coords)]
    res = Residue(chain=chain, resseq=resseq, icode=" ",
                   resname=resname, aa=aa, atoms=atoms)
    return res


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_threshold_flags():
    """Every threshold set must flag green / amber / red correctly."""
    # Published 2024 bands
    spec = THRESHOLD_SETS["2024"]
    check(flag_value(49, spec["Ltot"]) == "green", "Ltot 49 should be green")
    check(flag_value(40, spec["Ltot"]) == "amber", "Ltot 40 should be amber")
    check(flag_value(60, spec["Ltot"]) == "amber", "Ltot 60 should be amber")
    check(flag_value(30, spec["Ltot"]) == "red", "Ltot 30 should be red")
    check(flag_value(70, spec["Ltot"]) == "red", "Ltot 70 should be red")
    check(flag_value(140.0, spec["PSH"]) == "green", "PSH 140 green")
    check(flag_value(100.0, spec["PSH"]) == "amber", "PSH 100 amber")
    check(flag_value(90.0, spec["PSH"]) == "red", "PSH 90 red")
    check(flag_value(210.0, spec["PSH"]) == "red", "PSH 210 red")
    check(flag_value(1.0, spec["PPC"]) == "green", "PPC 1.0 green")
    check(flag_value(3.0, spec["PPC"]) == "amber", "PPC 3.0 amber")
    check(flag_value(5.0, spec["PPC"]) == "red", "PPC 5.0 red")
    check(flag_value(1.5, spec["PNC"]) == "green", "PNC 1.5 green")
    check(flag_value(3.0, spec["PNC"]) == "amber", "PNC 3.0 amber")
    check(flag_value(6.0, spec["PNC"]) == "red", "PNC 6.0 red")
    check(flag_value(2.0, spec["SFvCSP"]) == "green", "SFvCSP +2 green")
    check(flag_value(-10.0, spec["SFvCSP"]) == "amber", "SFvCSP -10 amber")
    check(flag_value(-40.0, spec["SFvCSP"]) == "red", "SFvCSP -40 red")
    # Every bundled set must assign a flag to representative values
    for name, set_spec in THRESHOLD_SETS.items():
        for metric in METRIC_ORDER:
            flag = flag_value(50.0 if metric in ("Ltot",) else 100.0
                              if metric == "PSH" else 1.0, set_spec[metric])
            check(flag in ("green", "amber", "red"),
                  f"{name}/{metric} produced invalid flag {flag!r}")
    # Risk escalation follows the published rule of thumb
    check(risk_summary({"Ltot": "green", "PSH": "green"}) == "LOW RISK",
          "all green should be LOW RISK")
    check(risk_summary({"Ltot": "amber", "PSH": "green"}) == "LOW-MEDIUM RISK",
          "one amber should be LOW-MEDIUM RISK")
    check(risk_summary({"Ltot": "amber", "PSH": "amber"}) == "MEDIUM RISK",
          "two amber should be MEDIUM RISK")
    check(risk_summary({"Ltot": "red", "PSH": "green"}) == "MEDIUM-HIGH RISK",
          "one red should be MEDIUM-HIGH RISK")
    check(risk_summary({"Ltot": "red", "PSH": "red"}) == "HIGH RISK",
          "two red should be HIGH RISK")


def test_build_pair_distances():
    """Closest heavy-atom distances between residues."""
    a = _residue("H", 1, "ALA", "A", (0.0, 0.0, 0.0))
    b = _residue("H", 2, "LEU", "L", (3.0, 4.0, 0.0))     # exactly 5.0 A
    c = _residue("H", 3, "TRP", "W", (0.0, 0.0, 50.0))    # far away
    dists = build_pair_distances([a, b, c])
    check((0, 1) in dists, "pair (0,1) must be found")
    check(abs(dists[(0, 1)] - 5.0) < 1e-9, f"expected 5.0, got {dists[(0, 1)]}")
    check((0, 2) not in dists, "50 A pair must be outside the 7.5 A cutoff")
    check((1, 2) not in dists, "far pair must be outside the cutoff")
    # Directional lookup must be symmetric for the caller
    dists2 = build_pair_distances([b, a, c])
    check((0, 1) in dists2 and abs(dists2[(0, 1)] - 5.0) < 1e-9,
          "pair distances must not depend on atom traversal order")


def test_surface_score_math():
    """PSH/PPC/PNC must equal the published pairwise x1*x2/r**2 sums."""
    # Two exposed LYS residues 3 A apart, both in the CDR vicinity.
    r1 = _residue("H", 27, "LYS", "K", (0.0, 0.0, 0.0))
    r2 = _residue("H", 28, "LYS", "K", (3.0, 0.0, 0.0))
    for r in (r1, r2):
        r.is_surface = True
        r.in_cdr_vicinity = True
        r.hydrophobicity = NORMALISED_HYDROPHOBICITY["LYS"]
        r.charge = RESIDUE_CHARGES["LYS"]
    dists = {(0, 1): 3.0}
    psh, ppc, pnc = compute_surface_scores([r1, r2], dists)
    h = NORMALISED_HYDROPHOBICITY["LYS"]
    check(abs(psh - 2 * h * h / 9.0) < 1e-9,
          f"PSH {psh} != {2 * h * h / 9.0}")
    check(abs(ppc - 2 * 1.0 * 1.0 / 9.0) < 1e-9,
          f"PPC {ppc} != {2 / 9.0}")
    check(abs(pnc) < 1e-12, f"PNC should be 0 for two positives, got {pnc}")

    # Mixed signs: no PPC/PNC, but PSH still counts the pair.
    r2.charge = RESIDUE_CHARGES["ASP"]
    r2.hydrophobicity = NORMALISED_HYDROPHOBICITY["ASP"]
    r2.resname, r2.aa = "ASP", "D"
    psh, ppc, pnc = compute_surface_scores([r1, r2], dists)
    h1 = NORMALISED_HYDROPHOBICITY["LYS"]
    h2 = NORMALISED_HYDROPHOBICITY["ASP"]
    check(abs(psh - 2 * h1 * h2 / 9.0) < 1e-9, "PSH must sum all vicinity pairs")
    check(abs(ppc) < 1e-12 and abs(pnc) < 1e-12,
          "opposite signs must contribute to neither PPC nor PNC")

    # Two ASP residues -> PNC only.
    r2.resname, r2.aa, r2.charge = "ASP", "D", RESIDUE_CHARGES["ASP"]
    r2.hydrophobicity = NORMALISED_HYDROPHOBICITY["ASP"]
    _psh, ppc, pnc = compute_surface_scores([r1, r2], dists)
    # r1 is still LYS -> still mixed; flip it too
    r1.resname, r1.aa, r1.charge = "ASP", "D", RESIDUE_CHARGES["ASP"]
    r1.hydrophobicity = NORMALISED_HYDROPHOBICITY["ASP"]
    _psh, ppc, pnc = compute_surface_scores([r1, r2], dists)
    check(abs(pnc - 2 * 1.0 / 9.0) < 1e-12,
          f"PNC {pnc} != {2 / 9.0} for two acidic residues")
    check(abs(ppc) < 1e-12, "PPC must stay 0 for two acidic residues")

    # A residue outside the CDR vicinity must not contribute.
    r2.in_cdr_vicinity = False
    r1.in_cdr_vicinity = True
    _psh, _ppc, pnc = compute_surface_scores([r1, r2], dists)
    check(abs(pnc) < 1e-12, "residues outside the CDR vicinity must be ignored")


def test_salt_bridge_zeroes_charge():
    """A K/R ... D/E pair within 3.2 A must zero both charges."""
    r1 = _residue("H", 100, "LYS", "K", (0.0, 0.0, 0.0))
    r2 = _residue("L", 50, "ASP", "D", (3.0, 0.0, 0.0))
    for r in (r1, r2):
        r.is_surface = True
        r.charge = RESIDUE_CHARGES[r.resname]
        r.hydrophobicity = NORMALISED_HYDROPHOBICITY[r.resname]
    # give both residues their defining side-chain atom names
    r1.atoms = [(" NZ ", "N", 0.0, 0.0, 0.0)]
    r2.atoms = [(" OD1", "O", 3.0, 0.0, 0.0)]
    from tap_core import assign_salt_bridges
    n = assign_salt_bridges([r1, r2])
    check(n == 1, f"expected one salt bridge, got {n}")
    check(r1.salt_bridge_partner == r2.key, "LYS must point at ASP")
    check(r2.salt_bridge_partner == r1.key, "ASP must point at LYS")
    # annotation re-applies charge/hydrophobicity after the bridge search
    annotate_residues([r1, r2])
    check(r1.charge == 0.0 and r2.charge == 0.0,
          f"bridged charges must be 0 (got {r1.charge}, {r2.charge})")
    check(abs(r1.hydrophobicity - NORMALISED_HYDROPHOBICITY["GLY"]) < 1e-12,
          "bridged residues must take glycine hydrophobicity")


def test_sequence_only_mode():
    """Sequence-only TAP must give an exact ANARCI Ltot."""
    result = analyze_sequences(VH_4D5, VL_4D5, name="4d5")
    check(result.metrics["Ltot"] > 0, "Ltot must be positive")
    for metric in METRIC_ORDER[1:]:
        check(result.metrics[metric] is None, f"{metric} must be n/a")
        check(result.flags[metric] == "unknown", f"{metric} flag must be unknown")
    total = sum(sum(rows.values()) for rows in result.details["cdr_lengths"].values())
    check(int(result.metrics["Ltot"]) == total,
          "Ltot must equal the sum of the per-CDR lengths")
    check("unknown" in result.risk_summary.lower(),
          "sequence-only risk must say it is unknown")


def test_synthetic_fv_end_to_end(tmpdir=None):
    """Full pipeline on a synthetic Fv: invariants, determinism, consistency."""
    workdir = tmpdir or tempfile.mkdtemp(prefix="tap_test_")
    pdb = os.path.join(workdir, "synthetic_fv.pdb")
    _write_synthetic_pdb(pdb, [("H", VH_4D5), ("L", VL_4D5)])

    r1 = analyze_structure(pdb, sasa_engine="freesasa", numbering="anarci")
    r2 = analyze_structure(pdb, sasa_engine="freesasa", numbering="anarci")

    for metric in METRIC_ORDER:
        check(metric in r1.metrics, f"missing metric {metric}")
    check(r1.metrics["Ltot"] > 0, "Ltot must be positive")
    for metric in ("PSH", "PPC", "PNC"):
        check(r1.metrics[metric] >= 0, f"{metric} must be non-negative")
    check(r1.metrics["PSH"] > 0, "PSH must be > 0 for a structured Fv")

    # Deterministic: identical second run
    for metric in METRIC_ORDER:
        check(r1.metrics[metric] == r2.metrics[metric],
              f"{metric} is not deterministic: {r1.metrics[metric]} vs {r2.metrics[metric]}")
    check(r1.flags == r2.flags, "flags must be deterministic")

    # Structure-mode Ltot must equal sequence-mode Ltot (same ANARCI numbering)
    seq = analyze_sequences(VH_4D5, VL_4D5)
    check(int(r1.metrics["Ltot"]) == int(seq.metrics["Ltot"]),
          f"structure Ltot {r1.metrics['Ltot']} != sequence Ltot {seq.metrics['Ltot']}")

    # SFvCSP is the product of the reported chain charges
    d = r1.details
    expected = d["net_charge_vh"] * d["net_charge_vl"]
    check(abs(r1.metrics["SFvCSP"] - expected) < 1e-9,
          f"SFvCSP {r1.metrics['SFvCSP']} != QH*QL {expected}")

    # Flags follow the threshold table
    spec = THRESHOLD_SETS[d["thresholds"]]
    for metric in METRIC_ORDER:
        check(r1.flags[metric] == flag_value(r1.metrics[metric], spec[metric]),
              f"{metric} flag {r1.flags[metric]} disagrees with thresholds")

    # ANARCI numbering must have been applied
    check("ANARCI" in d["numbering"], f"expected ANARCI numbering, got {d['numbering']}")
    check(d["vh_chain"] == "H" and d["vl_chain"] == "L",
          f"chain detection failed: VH={d['vh_chain']} VL={d['vl_chain']}")
    check(d["n_residues_analysed"] == len(VH_4D5) + len(VL_4D5),
          f"expected {len(VH_4D5) + len(VL_4D5)} analysed residues, "
          f"got {d['n_residues_analysed']}")


def test_antigen_chain_excluded():
    """Fv scope must ignore non-antibody chains (antigen, constant region)."""
    workdir = tempfile.mkdtemp(prefix="tap_test_ag_")
    pdb = os.path.join(workdir, "fv_plus_antigen.pdb")
    antigen = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPF"
    _write_synthetic_pdb(pdb, [("A", antigen), ("H", VH_4D5), ("L", VL_4D5)])

    r = analyze_structure(pdb, sasa_engine="freesasa", numbering="anarci",
                          scope="fv")
    d = r.details
    check(d["vh_chain"] == "H" and d["vl_chain"] == "L",
          f"antigen chain A was mistaken for an Fv chain: "
          f"VH={d['vh_chain']} VL={d['vl_chain']}")
    check(d["n_residues_analysed"] == len(VH_4D5) + len(VL_4D5),
          f"fv scope must analyse only the Fv, got {d['n_residues_analysed']} residues")

    # Same file in complex scope keeps the antigen for surface burial.
    r2 = analyze_structure(pdb, sasa_engine="freesasa", numbering="anarci",
                           scope="complex")
    check(r2.details["n_residues_analysed"] > d["n_residues_analysed"],
          "complex scope must analyse more residues than fv scope")

    # Ltot is identical in both scopes (CDR count only spans the Fv).
    check(int(r2.metrics["Ltot"]) == int(r.metrics["Ltot"]),
          "Ltot must not depend on scope")


def test_reports_render():
    """Text / JSON / Markdown / Word renderers must not crash."""
    sys.path.insert(0, os.path.join(ROOT, "tools", "tap_profiler"))
    from tap_analyzer import (format_text_report, save_json_report,
                              save_markdown_report, save_word_report)
    workdir = tempfile.mkdtemp(prefix="tap_test_rep_")
    pdb = os.path.join(workdir, "fv.pdb")
    _write_synthetic_pdb(pdb, [("H", VH_4D5), ("L", VL_4D5)])
    result = analyze_structure(pdb, sasa_engine="freesasa")

    text = format_text_report(result, "TAP Analysis: test")
    check("TAP METRICS" in text, "text report must contain the metrics block")
    check("SFvCSP" in text, "text report must list SFvCSP")
    check("RISK ASSESSMENT" in text, "text report must contain the risk block")

    json_path = os.path.join(workdir, "out.json")
    save_json_report(result, json_path)
    with open(json_path) as fh:
        payload = json.load(fh)
    for metric in METRIC_ORDER:
        check(metric in payload["metrics"], f"JSON must contain {metric}")
    check(payload["method"]["framework"].startswith("TAP"), "JSON must cite TAP")

    md_path = os.path.join(workdir, "out.md")
    save_markdown_report(result, md_path, "TAP Analysis: test")
    with open(md_path) as fh:
        md = fh.read()
    check("| Ltot |" in md, "Markdown must contain the metric table")

    word_path = os.path.join(workdir, "out.docx")
    save_word_report(result, word_path, "TAP Analysis: test")
    check(os.path.getsize(word_path) > 1000, "Word report must be non-trivial")

    # Sequence-only result must render too (n/a metrics)
    seq = analyze_sequences(VH_4D5, VL_4D5, name="seq")
    text = format_text_report(seq, "TAP Analysis: seq")
    check("n/a" in text, "sequence-only text report must show n/a metrics")


def test_cli_end_to_end():
    """The CLI must run --check and a full structure analysis."""
    env = dict(os.environ)
    env["PYTHONPATH"] = ROOT
    proc = subprocess.run(
        [sys.executable, os.path.join(TOOL, "tap_analyzer.py"), "--check"],
        capture_output=True, text=True, env=env, timeout=300)
    check(proc.returncode == 0, f"--check failed: {proc.stderr}")
    check("ANARCI" in proc.stdout, "--check must report ANARCI status")

    workdir = tempfile.mkdtemp(prefix="tap_test_cli_")
    pdb = os.path.join(workdir, "fv.pdb")
    _write_synthetic_pdb(pdb, [("H", VH_4D5), ("L", VL_4D5)])
    outdir = os.path.join(workdir, "out")
    proc = subprocess.run(
        [sys.executable, os.path.join(TOOL, "tap_analyzer.py"),
         "--pdb", pdb, "--all-formats", "--output", outdir],
        capture_output=True, text=True, env=env, timeout=600)
    check(proc.returncode == 0, f"CLI run failed: {proc.stderr}")
    check("TAP METRICS" in proc.stdout, "CLI must print the metrics table")
    files = os.listdir(outdir)
    for suffix in ("_tap.json", "_tap.md", "_tap.docx"):
        check(any(f.endswith(suffix) for f in files),
              f"missing {suffix} in {files}")


def test_reference_parity():
    """Optional: exact match against the Exscientia reference implementation.

    Enable with:
      TAP_PARITY_PDB=/path/test_complex_reference.pdb
      TAP_PARITY_PSA=/path/to/psa
    Expected values (Exscientia/ab-characterisation, MIT, run on that PDB):
      Ltot 64, PSH 199.813, PPC 0.323, PNC 0.466, SFvCSP -2.880
    """
    pdb = os.environ.get("TAP_PARITY_PDB")
    psa = os.environ.get("TAP_PARITY_PSA")
    if not pdb or not psa:
        print("SKIP: test_reference_parity (set TAP_PARITY_PDB / TAP_PARITY_PSA)")
        return
    if not (os.path.isfile(pdb) and os.path.isfile(psa)):
        print("SKIP: test_reference_parity (paths not found)")
        return
    r = analyze_structure(pdb, sasa_engine="psa", psa_binary=psa,
                          numbering="native", scope="complex",
                          thresholds="exscientia")
    expected = {"Ltot": 64.0, "PSH": 199.813, "PPC": 0.323,
                "PNC": 0.466, "SFvCSP": -2.880}
    for metric, want in expected.items():
        got = r.metrics[metric]
        check(abs(got - want) < 5e-4,
              f"{metric}: got {got}, reference {want}")


def run_all_tests():
    tests = [
        test_threshold_flags,
        test_build_pair_distances,
        test_surface_score_math,
        test_salt_bridge_zeroes_charge,
        test_sequence_only_mode,
        test_synthetic_fv_end_to_end,
        test_antigen_chain_excluded,
        test_reports_render,
        test_cli_end_to_end,
        test_reference_parity,
    ]
    passed = failed = 0
    for test in tests:
        try:
            test()
            print(f"PASS: {test.__name__}")
            passed += 1
        except Exception as exc:
            print(f"FAIL: {test.__name__}: {exc}")
            import traceback
            traceback.print_exc()
            failed += 1
    print(f"\n{passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    ok = run_all_tests()
    sys.exit(0 if ok else 1)
