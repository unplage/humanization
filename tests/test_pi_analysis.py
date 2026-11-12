#!/usr/bin/env python3
"""Regression tests for tools/pi_analysis (plain script; exits nonzero on failure).

Covers:
  * pKa-set values / aliases and sequence pI against Biopython
  * free vs capped termini, disulfide handling, VH/VL charge descriptors
  * FASTA chain-order inference (VL-first files)
  * structure-aware pI, TAP metrics, interface flags and patch candidates
  * mutation filters (CDR/FR4/VHH hallmark/interface, patch sign)
  * text/JSON/Markdown/CSV renderers and the CLI
  * PROPKA .pka parsing (no PROPKA binary required)

Structure tests build a deterministic synthetic Fv from the tracked 4D5
sequences, so they do not depend on af3/ or inputs/ (both gitignored).
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "pi_analysis")
sys.path.insert(0, TOOL)
sys.path.insert(0, os.path.join(ROOT, "tools", "tap_profiler"))

import pi_core  # noqa: E402
import pi_analyzer  # noqa: E402
import tap_core  # noqa: E402
from pi_core import (  # noqa: E402
    ChainInfo,
    Patch,
    analyze_sequence,
    analyze_structure,
    isoelectric_point,
    mutation_candidates,
    net_charge,
    parse_propka_pka,
    parse_propka_summary,
    sequence_pi_table,
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


def check(condition, message):
    if not condition:
        raise AssertionError(message)


def _write_synthetic_fv(path):
    serial = 1
    with open(path, "w") as fh:
        for chain_offset, (chain_id, seq) in enumerate((("H", VH_4D5), ("L", VL_4D5))):
            lateral = chain_offset * 5.0 + 0.37
            for i, aa in enumerate(seq, start=1):
                resname = ONE_TO_THREE[aa]
                x = (i % 20) * 3.8 + lateral
                y = (i // 20) * 3.8 + lateral
                z = ((i * 7 + chain_offset * 3 + 1) % 9) * 0.9 + chain_offset * 0.13
                atoms = [
                    (" CA ", 0.0, 0.0, 0.0, "C"),
                    (" N  ", -1.2, 0.5, 0.1, "N"),
                    (" C  ", 1.2, 0.5, -0.1, "C"),
                    (" O  ", 1.7, 1.6, -0.1, "O"),
                ]
                if aa != "G":
                    atoms.append((" CB ", 0.1, -1.5, 0.6, "C"))
                for name, dx, dy, dz, element in atoms:
                    fh.write(
                        "ATOM  {:5d} {:4s} {:>3s} {:1s}{:4d}    "
                        "{:8.3f}{:8.3f}{:8.3f}  1.00 20.00          {:>2s}\n"
                        .format(serial, name, resname, chain_id, i,
                                x + dx, y + dy, z + dz, element))
                    serial += 1
        fh.write("END\n")


def _residue(chain, resseq, aa, imgt_pos, is_cdr=False, is_anchor=False,
             in_cdr_vicinity=False, icode=" "):
    from tap_core import Residue
    return Residue(
        chain=chain, resseq=resseq, icode=icode,
        resname=ONE_TO_THREE[aa], aa=aa, atoms=[],
        imgt_pos=imgt_pos, is_cdr=is_cdr, is_anchor=is_anchor,
        in_cdr_vicinity=in_cdr_vicinity, rel_sasa=50.0,
    )


def test_pka_tables_and_aliases():
    expected = {
        "Bjellqvist": (7.5, 3.55, 4.05, 4.45, 5.98, 10.0, 12.0),
        "IPC_protein": (9.094, 2.869, 3.872, 4.412, 5.637, 9.052, 11.84),
        "EMBOSS": (8.6, 3.6, 3.9, 4.1, 6.5, 10.8, 12.5),
        "Sillero": (8.2, 3.2, 4.0, 4.5, 6.4, 10.4, 12.0),
        "Nozaki_Tanford": (7.5, 3.8, 4.0, 4.4, 6.3, 10.4, 12.0),
    }
    for name, values in expected.items():
        table = pi_core.PKA_SETS[name]
        got = tuple(table[k] for k in ("Nterm", "Cterm", "D", "E", "H", "K", "R"))
        check(got == values, "{} pKa values changed: {}".format(name, got))
    check(pi_core.normalize_pka_set("expasy") == "Bjellqvist", "ExPASy alias")
    check(pi_core.normalize_pka_set("IPC") == "IPC_protein", "IPC alias")
    check(pi_core.normalize_pka_set("nozaki-tanford") == "Nozaki_Tanford",
          "Nozaki alias")


def test_sequence_pi_matches_biopython():
    try:
        from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
    except ImportError:
        print("SKIP: Biopython not available")
        return
    for name, seq in (("VH", VH_4D5), ("VL", VL_4D5)):
        table = sequence_pi_table([ChainInfo(name=name, chain_type="H", sequence=seq)],
                                  "Bjellqvist")
        want = IsoelectricPoint(seq).pi()
        got = table["pI_free"]
        check(abs(got - want) < 0.02,
              "{} Bjellqvist pI {:.4f} != Biopython {:.4f}".format(name, got, want))

    chains = [
        ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
        ChainInfo(name="VL", chain_type="L", sequence=VL_4D5),
    ]
    combined = sequence_pi_table(chains, "Bjellqvist")["pI_free"]
    vh_pi = sequence_pi_table([chains[0]], "Bjellqvist")["pI_free"]
    vl_pi = sequence_pi_table([chains[1]], "Bjellqvist")["pI_free"]
    check(min(vh_pi, vl_pi) <= combined <= max(vh_pi, vl_pi),
          "combined Fv pI should lie between chain pIs")

    groups, termini = pi_core._make_groups_from_sequence(chains, "Bjellqvist")
    q = net_charge(groups, termini, combined)
    check(abs(q) < 1e-5, "net charge at pI should be ~0, got {}".format(q))


def test_sequence_charge_descriptors():
    chains = [
        ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
        ChainInfo(name="VL", chain_type="L", sequence=VL_4D5),
    ]
    result = analyze_sequence(chains, pka_set="IPC_protein", all_pka_sets=True)
    fv = result["fv"]
    check(fv["pI_free"] is not None and 0.0 < fv["pI_free"] < 14.0,
          "invalid Fv pI")
    qh = result["chains"]["VH"]["net_charge_ph74"]
    ql = result["chains"]["VL"]["net_charge_ph74"]
    check(abs(fv["Fv_chml"] - (qh - ql)) < 1e-9, "Fv_chml must be VH-VL")
    check(abs(fv["SFvCSP"] - (qh * ql)) < 1e-9, "SFvCSP must be VH*VL")
    check(len(result["sequence_pI"]) == len(pi_core.PKA_SETS),
          "all_pka_sets should populate every set")
    for name, entry in result["sequence_pI"].items():
        check(entry["pI_free"] is not None and math.isfinite(entry["pI_free"]),
              "{} pI not finite".format(name))


def test_vhh_sequence_mode():
    result = analyze_sequence(
        [ChainInfo(name="VHH", chain_type="H", sequence=VH_4D5)],
        pka_set="Bjellqvist")
    check(result["fv"]["Fv_chml"] is None, "VHH has no VH/VL charge asymmetry")
    check(result["fv"]["pI_free"] is not None, "VHH pI missing")
    check(result["patches"] == [], "sequence mode should not invent patches")


def test_disulfide_assignment():
    chains = [
        ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
        ChainInfo(name="VL", chain_type="L", sequence=VL_4D5),
    ]
    found, warnings = pi_core.find_disulfide_cys_sequence(chains)
    check(len(found) == 4, "expected two disulfides (4 Cys), got {}".format(found))
    auto = analyze_sequence(chains, pka_set="Bjellqvist", disulfide_mode="auto")
    free = analyze_sequence(chains, pka_set="Bjellqvist", disulfide_mode="free")
    check(auto["fv"]["pI_free"] != free["fv"]["pI_free"],
          "auto vs free Cys should change pI")
    check(auto["disulfide_mode"] == "auto", "auto disulfide mode")
    check(free["disulfide_mode"] == "free", "free disulfide mode")
    check(not any("no conserved disulfide" in w.lower() for w in auto["warnings"]),
          "auto mode should find the conserved disulfide")


def test_fasta_chain_order():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "vl_first.fasta")
        with open(path, "w") as fh:
            fh.write(">x_VL\n{}\n>x_VH\n{}\n".format(VL_4D5, VH_4D5))
        chains = pi_analyzer._load_fasta_chains(path)
        check(chains[0].name == "VH" and chains[0].chain_type == "H",
              "VL-first FASTA must still yield VH first")
        check(chains[1].name == "VL" and chains[1].chain_type == "L",
              "VL-first FASTA must still yield VL second")
        check(chains[0].sequence == VH_4D5, "VH sequence mapping")
        check(chains[1].sequence == VL_4D5, "VL sequence mapping")

        path = os.path.join(tmp, "dot_headers.fasta")
        with open(path, "w") as fh:
            fh.write(">antibody.H\n{}\n>antibody.L\n{}\n".format(VH_4D5, VL_4D5))
        chains = pi_analyzer._load_fasta_chains(path)
        check(chains[0].chain_type == "H" and chains[1].chain_type == "L",
              "dot-separated H/L headers")

        path = os.path.join(tmp, "single.fasta")
        with open(path, "w") as fh:
            fh.write(">nanobody\n{}\n".format(VH_4D5))
        chains = pi_analyzer._load_fasta_chains(path)
        check(len(chains) == 1 and chains[0].name == "VHH",
              "single FASTA record should be VHH")


def test_propka_parser():
    text = """
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   q    RMS  (SD)
  ASP  25 A   3.80    3.80
  GLU  31 A   4.90    4.50
  HIS  40 A   6.10    6.50
  CYS  23 A   9.00    9.00
  TYR  37 A  10.20   10.00
  LYS  43 A  10.40   10.50
  ARG  45 A  12.30   12.00
"""
    got = parse_propka_pka(text)
    check(len(got) == 7, "expected 7 PROPKA residues, got {}".format(got))
    check(got[("A", 25, " ")] == 3.80, "ASP pKa parse")
    check(got[("A", 45, " ")] == 12.30, "ARG pKa parse")


def test_structure_end_to_end():
    try:
        import anarci  # noqa: F401
        import freesasa  # noqa: F401
    except ImportError:
        print("SKIP: ANARCI/FreeSASA not available")
        return
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "synthetic_fv.pdb")
        _write_synthetic_fv(pdb)
        result = analyze_structure(
            pdb, pka_set="IPC_protein", all_pka_sets=True,
            structure_model="heuristic")
    check(result["input"]["kind"] == "fab", "synthetic Fv should be Fab")
    check(result["input"]["vh_chain"] == "H", "VH chain detection")
    check(result["input"]["vl_chain"] == "L", "VL chain detection")
    check(result["structure_pI"], "structure pI table missing")
    check(result["fv"]["pI_free"] is not None, "Fv pI missing")
    check(result["fv"]["pI_selected"] is not None, "selected pI missing")
    check(result["termini_mode"] == "capped", "Fab should default to capped termini")
    check(result["fv"]["Fv_chml"] is not None, "Fv_chml missing")
    check(result["fv"]["Fv_chml_total"] is not None, "total Fv_chml missing")
    check(result["fv"]["SFvCSP"] is not None, "SFvCSP missing")
    check(result["fv"]["SFvCSP_total"] is not None, "total SFvCSP missing")
    check(result["sensitivity"]["pI_selected_range"] is not None,
          "pKa sensitivity range missing")
    check(len(result["sensitivity"]["pka_sets"]) == len(pi_core.PKA_SETS),
          "sensitivity should cover every bundled pKa set")
    curve = result["charge_curve"]
    check(len(curve) >= 2 and all("net_charge" in p for p in curve),
          "charge curve missing")
    check(all(curve[i]["net_charge"] >= curve[i + 1]["net_charge"] - 1e-6
              for i in range(len(curve) - 1)),
          "net charge should decrease with pH")
    check(result["structure_quality"]["n_residues"] > 0,
          "structure quality missing")
    check(result["tap"]["n_salt_bridges_tap"] is not None
          and result["tap"]["n_salt_bridges_pka"] is not None,
          "both salt-bridge definitions should be reported")
    check(result["tap"]["metrics"]["Ltot"] > 0, "TAP Ltot missing")
    check(result["patches"], "expected charge patches")
    for patch in result["patches"]:
        check(patch["sign"] in ("positive", "negative", "mixed"), "patch sign")
        check(patch["residues"], "empty patch")
        for row in patch["residues"]:
            check("is_interface" in row, "patch residue lacks is_interface")
    check(any(r.get("is_interface") for p in result["patches"]
              for r in p.get("residues", [])),
          "Fab should have at least one VH/VL interface patch residue")
    check(any("pI" in r["category"] for r in result["recommendations"]),
          "missing pI recommendation")


def test_patch_candidate_filters():
    residues = [
        _residue("H", 10, "D", 10),                       # framework acid
        _residue("H", 11, "K", 11),                       # framework base
        _residue("H", 12, "D", 12, is_cdr=True),           # CDR acid
        _residue("H", 13, "D", 118),                       # FR4 acid
        _residue("H", 14, "K", 42),                        # VHH hallmark
        _residue("H", 15, "D", 15, in_cdr_vicinity=True),  # vicinity acid
        _residue("H", 16, "D", 1),                        # N-terminal acid
    ]
    patch = Patch(
        patch_id=1, sign="negative", ph=6.0, total_charge=-3.0,
        residues=[
            {"label": "H10", "aa": "D", "region": "FR1", "rel_sasa": 50.0,
             "charge": -1.0, "in_cdr_vicinity": False},
            {"label": "H11", "aa": "K", "region": "FR1", "rel_sasa": 50.0,
             "charge": 1.0, "in_cdr_vicinity": False},
            {"label": "H12", "aa": "D", "region": "CDR1", "rel_sasa": 50.0,
             "charge": -1.0, "in_cdr_vicinity": False},
            {"label": "H13", "aa": "D", "region": "FR4", "rel_sasa": 50.0,
             "charge": -1.0, "in_cdr_vicinity": False},
            {"label": "H14", "aa": "K", "region": "FR2", "rel_sasa": 50.0,
             "charge": 1.0, "in_cdr_vicinity": False},
            {"label": "H15", "aa": "D", "region": "FR1", "rel_sasa": 50.0,
             "charge": -1.0, "in_cdr_vicinity": True},
            {"label": "H16", "aa": "D", "region": "FR1", "rel_sasa": 50.0,
             "charge": -1.0, "in_cdr_vicinity": False},
        ],
    )
    cands = mutation_candidates(
        patch, residues,
        interface_keys={("H", 15, " ")},
        vhh_hallmark_keys={("H", 14, " ")},
        terminal_keys={("H", 16, " ")})
    labels = {c["label"] for c in cands}
    check("H10" in labels, "framework acid should be a candidate")
    check("H11" not in labels, "positive residue must not be mutated in a negative patch")
    check("H12" not in labels, "CDR residue must be excluded")
    check("H13" not in labels, "FR4 residue must be excluded")
    check("H14" not in labels, "VHH hallmark must be excluded")
    check("H15" in labels, "vicinity acid should remain a candidate")
    check("H16" not in labels, "N-terminal residue must be excluded")
    h15 = [c for c in cands if c["label"] == "H15"][0]
    check(h15["risk"] == "high", "interface candidate should be high risk")
    check(isinstance(h15.get("delta_pI_estimate"), (int, float))
          and isinstance(h15.get("delta_charge_formulation_estimate"), (int, float)),
          "candidate should carry what-if estimates")
    effects = h15.get("substitution_effects", {})
    check(effects and set(effects) <= {"Q", "N", "S", "T", "A"},
          "candidate should carry per-substitution effects")

    positive = Patch(
        patch_id=2, sign="positive", ph=6.0, total_charge=2.0,
        residues=[
            {"label": "H10", "aa": "D", "region": "FR1", "rel_sasa": 50.0,
             "charge": -1.0},
            {"label": "H11", "aa": "K", "region": "FR1", "rel_sasa": 50.0,
             "charge": 1.0},
        ],
    )
    labels = {c["label"] for c in mutation_candidates(positive, residues)}
    check("H11" in labels, "positive residue should be a candidate in a positive patch")
    check("H10" not in labels, "negative residue must not be mutated in a positive patch")


def test_reports_and_cli():
    result = analyze_sequence(
        [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
         ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)],
        pka_set="IPC_protein", all_pka_sets=True, name="4D5")
    text = pi_analyzer.format_text_report(result, max_patches=2)
    check("pI BY PKA SET" in text and "RECOMMENDATIONS" in text,
          "text report sections missing")
    with tempfile.TemporaryDirectory() as tmp:
        json_path = os.path.join(tmp, "out.json")
        md_path = os.path.join(tmp, "out.md")
        csv_path = os.path.join(tmp, "out.csv")
        patch_csv = os.path.join(tmp, "patches.csv")
        pi_analyzer.save_json_report(result, json_path)
        pi_analyzer.save_markdown_report(result, md_path)
        pi_analyzer.save_csv_reports(result, csv_path, patch_csv)
        with open(json_path) as fh:
            payload = json.load(fh)
        check(payload["name"] == "4D5", "JSON name")
        check(payload["fv"]["pI_free"] is not None, "JSON pI")
        with open(md_path) as fh:
            check("# pI / Fv Charge Analysis" in fh.read(), "Markdown header")
        with open(csv_path) as fh:
            header = fh.readline()
        check("pI_selected" in header and "Fv_chml_total" in header,
              "CSV header")

        proc = subprocess.run(
            [sys.executable, os.path.join(TOOL, "pi_analyzer.py"), "--check"],
            capture_output=True, text=True)
        check(proc.returncode == 0, "CLI --check failed: {}".format(proc.stderr))
        check("pKa sets" in proc.stdout, "CLI --check output")

        outdir = os.path.join(tmp, "cli")
        proc = subprocess.run(
            [sys.executable, os.path.join(TOOL, "pi_analyzer.py"),
             "--vh", VH_4D5, "--vl", VL_4D5, "--json", "--output", outdir,
             "--name", "cli_4d5"],
            capture_output=True, text=True)
        check(proc.returncode == 0, "CLI sequence run failed: {}".format(proc.stderr))
        check(os.path.exists(os.path.join(outdir, "cli_4d5_pi.json")),
              "CLI JSON output missing")


def test_structure_model_none():
    try:
        import anarci  # noqa: F401
        import freesasa  # noqa: F401
    except ImportError:
        print("SKIP: ANARCI/FreeSASA not available")
        return
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "synthetic_fv.pdb")
        _write_synthetic_fv(pdb)
        heuristic = analyze_structure(pdb, pka_set="Bjellqvist",
                                      structure_model="heuristic")
        none = analyze_structure(pdb, pka_set="Bjellqvist",
                                 structure_model="none")
    check(heuristic["structure_pI"]["Bjellqvist"]["model"] == "heuristic",
          "heuristic model tag")
    check(none["structure_pI"]["Bjellqvist"]["model"] == "none",
          "none model tag")
    check(heuristic["fv"]["pI_free"] != none["fv"]["pI_free"],
          "microenvironment correction should change pI")


def test_termini_and_charge_curve():
    fab = analyze_sequence(
        [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
         ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)],
        pka_set="IPC_protein", all_pka_sets=True, termini_mode="auto")
    check(fab["termini_mode"] == "capped", "Fab auto termini should be capped")
    check(abs(fab["fv"]["pI_selected"] - fab["fv"]["pI_capped"]) < 1e-9,
          "Fab selected pI should equal capped pI")

    vhh = analyze_sequence(
        [ChainInfo(name="VHH", chain_type="H", sequence=VH_4D5)],
        pka_set="IPC_protein", all_pka_sets=True, termini_mode="auto")
    check(vhh["termini_mode"] == "free", "VHH auto termini should be free")
    check(abs(vhh["fv"]["pI_selected"] - vhh["fv"]["pI_free"]) < 1e-9,
          "VHH selected pI should equal free pI")

    free_fab = analyze_sequence(
        [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
         ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)],
        pka_set="IPC_protein", termini_mode="free")
    check(fab["fv"]["pI_selected"] > free_fab["fv"]["pI_selected"],
          "capped Fab pI should exceed free Fab pI")

    curve = fab["charge_curve"]
    check(len(curve) >= 10, "default pH grid should have >= 10 points")
    check(all(curve[i]["net_charge"] >= curve[i + 1]["net_charge"] - 1e-6
              for i in range(len(curve) - 1)),
          "net charge should be monotonic decreasing in pH")
    rng = fab["sensitivity"]["pI_selected_range"]
    check(rng[0] - 0.01 <= fab["fv"]["pI_selected"] <= rng[1] + 0.01,
          "selected pI should lie inside the sensitivity range")


def test_propka_termini_parser():
    text = """
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   q    RMS  (SD)
  ASP  25 A   3.80    3.80
  N+    1 A   8.10    8.00
  C-  127 A   3.30    3.20
"""
    side, termini = parse_propka_summary(text)
    check(side[("A", 25, " ")] == 3.80, "side-chain pKa parse")
    check(termini[("A", "N")] == 8.10, "N-terminal PROPKA pKa parse")
    check(termini[("A", "C")] == 3.30, "C-terminal PROPKA pKa parse")
    check(parse_propka_pka(text) == side, "parse_propka_pka wrapper")


def test_patch_risk_calibration():
    low, basis = pi_core._patch_risk(0.40, "positive", False, 3.0)
    elevated, _ = pi_core._patch_risk(0.60, "positive", False, 3.0)
    elevated_cdr, basis_cdr = pi_core._patch_risk(0.60, "positive", True, 3.0)
    high, _ = pi_core._patch_risk(2.50, "positive", False, 3.0)
    check(low == "low" and elevated == "elevated" and high == "high",
          "positive patch risk bands")
    check(elevated_cdr in ("elevated", "high"),
          "CDR vicinity should not lower patch risk")
    check("cdr" in basis_cdr.lower(), "CDR penalty should be documented in basis")
    balanced, balanced_basis = pi_core._patch_risk(0.10, "mixed", False, 0.0)
    check(balanced == "low" and "balanced" in balanced_basis.lower(),
          "balanced mixed patch should be low risk")


# ---------------------------------------------------------------------------
# Regression tests for the defects found in the scientific review
# ---------------------------------------------------------------------------

def test_propka_sensitivity_is_not_degenerate():
    """PROPKA mode must still vary with the bundled pKa set.

    Regression: the same PROPKA pKa dict was passed for every pKa set, so
    all five sets returned an identical pI and the reported uncertainty was
    a meaningless zero-width range.  PROPKA's *model* pKa is now used to
    separate the environmental contribution, which is transferred onto each
    bundled set.
    """
    import tools.pi_analysis.pi_core as core

    def fake_propka(residues, timeout=600):
        pka, model = {}, {}
        for r in residues:
            if r.aa in ("D", "E"):
                v, m = 3.9, 3.8
            elif r.aa in ("K", "R"):
                v, m = 10.5, 10.5
            elif r.aa == "H":
                v, m = 6.0, 6.5
            elif r.aa in ("C", "Y"):
                v, m = 9.0, 9.0
            else:
                continue
            pka[r.key] = v
            model[r.key] = m
        termini, tmodel = {}, {}
        for ch in {r.chain for r in residues}:
            termini[(ch, "N")] = 8.5
            tmodel[(ch, "N")] = 8.0
            termini[(ch, "C")] = 3.6
            tmodel[(ch, "C")] = 3.2
        return pka, model, termini, tmodel, None

    original_run = core._run_propka_on_fv
    original_find = core.find_propka
    core._run_propka_on_fv = fake_propka
    core.find_propka = lambda: "/usr/bin/propka3"
    try:
        with tempfile.TemporaryDirectory() as tmp:
            pdb = os.path.join(tmp, "fv.pdb")
            _write_synthetic_fv(pdb)
            result = analyze_structure(
                pdb, pka_set="IPC_protein", all_pka_sets=True,
                structure_model="propka")
    finally:
        core._run_propka_on_fv = original_run
        core.find_propka = original_find

    values = result["sensitivity"]["pka_sets"]
    check(len(values) == len(pi_core.PKA_SETS),
          "every pKa set should be reported")
    check(len(set(round(v, 6) for v in values.values())) > 1,
          "PROPKA mode must not collapse all pKa sets onto one value: {}".format(
              values))
    rng = result["sensitivity"]["pI_selected_range"]
    check(rng[1] > rng[0], "pKa-set sensitivity range must be non-degenerate")


def test_patch_graph_uses_documented_8A_cutoff():
    """Pairs at 7.5-8.0 A must join a patch (the documented cutoff is 8 A).

    Regression: the patch graph reused the TAP 7.5 A neighbour map, so a
    pair at 7.7 A was split into two patches.
    """
    from tap_core import Residue

    def res(resseq, atoms):
        r = Residue(chain="H", resseq=resseq, icode=" ", resname="ASP",
                    aa="D", atoms=atoms)
        r.imgt_pos = 30
        r.rel_sasa = 50.0
        return r

    d1 = [(" OD1", "O", 0.0, 0.0, 0.0), (" OD2", "O", 0.0, 1.0, 0.0),
          (" CG ", "C", -1.5, 0.5, 0.0), (" CB ", "C", -2.5, 1.5, 0.0)]
    d2 = [(" OD1", "O", 7.7, 0.0, 0.0), (" OD2", "O", 7.7, 1.0, 0.0),
          (" CG ", "C", 9.2, 0.5, 0.0), (" CB ", "C", 10.2, 1.5, 0.0)]
    residues = [res(10, d1), res(20, d2)]

    def groups_for(rs):
        return {r.key: pi_core.IonizableGroup(
            chain=r.chain, resseq=r.resseq, icode=r.icode, aa=r.aa,
            kind="acid", pka_model=4.0, pka_eff=4.0, active=True,
            structure_chain=r.chain) for r in rs}

    near = tap_core.build_pair_distances(residues, cutoff=7.5)
    wide = tap_core.build_pair_distances(residues, cutoff=pi_core.PATCH_CUTOFF)
    patches_near = pi_core.cluster_patches(
        residues, near, 6.0, groups_by_key=groups_for(residues))
    patches_wide = pi_core.cluster_patches(
        residues, wide, 6.0, groups_by_key=groups_for(residues))
    check(len(patches_near) == 2,
          "7.5 A map must not connect a 7.7 A pair (control)")
    check(len(patches_wide) == 1,
          "8.0 A map must connect a 7.7 A pair into one patch")


def test_same_chain_for_vh_and_vl_is_rejected():
    """--vh-chain == --vl-chain must fail loudly, not double-count."""
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "fv.pdb")
        _write_synthetic_fv(pdb)
        try:
            analyze_structure(pdb, vh_chain="H", vl_chain="H",
                              structure_model="heuristic")
        except ValueError as exc:
            check("different chains" in str(exc),
                  "error should explain that the chains must differ")
        else:
            raise AssertionError("same chain for VH and VL must be rejected")


def test_text_report_header_matches_columns():
    """The non-structure pKa table header must list every data column."""
    from tools.pi_analysis import pi_analyzer
    result = analyze_sequence(
        [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
         ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)],
        pka_set="IPC_protein")
    text = pi_analyzer.format_text_report(result)
    header = [ln for ln in text.splitlines() if "seq pI" in ln and "struct pI" in ln]
    check(len(header) == 1, "expected exactly one pKa-set header")
    for column in ("pKa set", "seq pI", "struct pI", "capped", "q pH74"):
        check(column in header[0],
              "header must list the {!r} column: {!r}".format(column, header[0]))
    data = [ln for ln in text.splitlines()
            if ln.strip().startswith("IPC_protein")]
    check(len(data) == 1, "expected one data row")
    check(len(data[0].split()) == 5,
          "the data row must have five values: {!r}".format(data[0]))


def test_ph_curve_reaches_saved_text_report():
    """--ph-curve must also appear in the saved text report and in batch mode."""
    out = tempfile.mkdtemp()
    try:
        pdb = os.path.join(out, "fv.pdb")
        _write_synthetic_fv(pdb)
        proc = subprocess.run(
            [sys.executable, os.path.join(TOOL, "pi_analyzer.py"),
             "--pdb", pdb, "--structure-model", "heuristic",
             "--text", "--ph-curve", "--output", out],
            capture_output=True, text=True)
        check(proc.returncode == 0,
              "CLI failed: {}".format(proc.stderr[-400:]))
        with open(os.path.join(out, "fv_pi.txt")) as fh:
            text = fh.read()
        check("NET CHARGE VS pH" in text,
              "saved text report must contain the charge curve")
    finally:
        for name in os.listdir(out):
            os.remove(os.path.join(out, name))
        os.rmdir(out)


def test_clean_sequence_warnings_are_reported():
    """U/O substitution and dropped residues must be reported."""
    seq = "EVQLUESGGGXQVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRYDFWSGYPPGEEYFQHWGQGTLVTVSS"
    result = analyze_sequence(
        [ChainInfo(name="VHH", chain_type="H", sequence=seq)],
        pka_set="Bjellqvist")
    joined = " ".join(result["warnings"]).lower()
    check("selenocysteine" in joined,
          "U -> C substitution must be warned about: {}".format(result["warnings"]))
    check("non-standard" in joined,
          "dropped residue must be warned about: {}".format(result["warnings"]))
    check("U" not in result["chains"]["VHH"]["sequence"],
          "selenocysteine must be substituted, not kept")


def test_bfactor_provenance_semantics():
    """A crystal B-factor must not be read as a low-confidence pLDDT."""
    from tap_core import Residue

    def write(path, header, bfactor):
        with open(path, "w") as fh:
            fh.write(header)
            serial = 1
            for cid, seq in (("H", VH_4D5), ("L", VL_4D5)):
                for i, aa in enumerate(seq, 1):
                    x = (i % 20) * 3.8
                    y = (i // 20) * 3.8
                    z = ((i * 7) % 9) * 0.9
                    for name, dx, dy, dz, el in (
                            (" CA ", 0, 0, 0, "C"), (" N  ", -1.2, 0.5, 0.1, "N"),
                            (" C  ", 1.2, 0.5, -0.1, "C"), (" O  ", 1.7, 1.6, -0.1, "O")):
                        fh.write(
                            "ATOM  {:5d} {:4s} {:>3s} {:1s}{:4d}    "
                            "{:8.3f}{:8.3f}{:8.3f}  1.00 {:6.2f}          {:>2s}\n".format(
                                serial, name, ONE_TO_THREE[aa], cid, i,
                                x + dx, y + dy, z + dz, bfactor, el))
                        serial += 1
            fh.write("END\n")

    with tempfile.TemporaryDirectory() as tmp:
        crystal = os.path.join(tmp, "crystal.pdb")
        write(crystal, "HEADER    X\nEXPDTA    X-RAY DIFFRACTION\n"
                       "REMARK   2 RESOLUTION.    1.90 ANGSTROMS\n", 20.0)
        predicted = os.path.join(tmp, "predicted.pdb")
        write(predicted, "HEADER    ALPHAFOLD MONOMER\n"
                         "REMARK   1 pLDDT per-residue confidence\n", 55.0)
        unknown = os.path.join(tmp, "unknown.pdb")
        write(unknown, "HEADER    SOMETHING\n", 30.0)

        r_cryst = analyze_structure(crystal, structure_model="heuristic")
        r_pred = analyze_structure(predicted, structure_model="heuristic")
        r_unk = analyze_structure(unknown, structure_model="heuristic")

    check(r_cryst["structure_quality"]["provenance"] == "experimental",
          "EXPDTA must mark the structure as experimental")
    check(r_pred["structure_quality"]["provenance"] == "predicted",
          "AlphaFold marker must mark the structure as predicted")
    check(r_unk["structure_quality"]["provenance"] == "unknown",
          "an unrecognised header must stay 'unknown'")
    check("crystallographic" in r_cryst["structure_quality"]["bfactor_meaning"],
          "crystal B-factor meaning")
    check("pLDDT" in r_pred["structure_quality"]["bfactor_meaning"],
          "predicted-model B-factor meaning")

    crystal_warnings = " ".join(r_cryst["warnings"]).lower()
    check("low mean ca b-factor" not in crystal_warnings,
          "a well-ordered crystal structure must not be called low-confidence")
    check("plddt" not in crystal_warnings,
          "crystal B-factors must not be reported as pLDDT")
    pred_warnings = " ".join(r_pred["warnings"]).lower()
    check("plddt" in pred_warnings,
          "a low-pLDDT model must be reported as low confidence")


def test_terminus_completeness_is_chain_type_specific():
    """VH ends at IMGT 128, VL at 127; both must be checked correctly."""
    check(pi_core.TRUE_TERMINUS_C_IMGT["H"] == 128,
          "heavy-chain FR4 ends at IMGT 128")
    check(pi_core.TRUE_TERMINUS_C_IMGT["L"] == 127,
          "light-chain FR4 ends at IMGT 127")
    check(pi_core._true_c_terminus("H") == 128, "H terminus threshold")
    check(pi_core._true_c_terminus("L") == 127, "L terminus threshold")

    from tap_core import Residue

    def chain_ending_at(imgt_pos, chain_id="H"):
        """A minimal two-residue chain: IMGT 1 at the N end, imgt_pos at C."""
        residues = []
        for resseq, pos in ((1, 1), (2, imgt_pos)):
            r = Residue(chain=chain_id, resseq=resseq, icode=" ",
                        resname="ALA", aa="A", atoms=[])
            r.imgt_pos = pos
            residues.append(r)
        return residues

    # A heavy chain truncated at IMGT 127 is missing its last FR4 residue.
    quality = pi_core.assess_structure_quality(
        chain_ending_at(127), "fab", chain_types={"H": "H"},
        provenance="experimental")
    check(any("C(IMGT 127" in t for t in quality["truncated_termini"]),
          "a VH truncated at IMGT 127 must be reported: {}".format(
              quality["truncated_termini"]))

    # A complete heavy chain at IMGT 128 must not be reported.
    quality = pi_core.assess_structure_quality(
        chain_ending_at(128), "fab", chain_types={"H": "H"},
        provenance="experimental")
    check(not quality["truncated_termini"],
          "a complete VH at IMGT 128 must not be reported as truncated: "
          "{}".format(quality["truncated_termini"]))

    # A light chain is complete at IMGT 127.
    quality = pi_core.assess_structure_quality(
        chain_ending_at(127, "L"), "fab", chain_types={"L": "L"},
        provenance="experimental")
    check(not quality["truncated_termini"],
          "a complete VL at IMGT 127 must not be reported as truncated: "
          "{}".format(quality["truncated_termini"]))

    # A light chain truncated at IMGT 126 must be reported.
    quality = pi_core.assess_structure_quality(
        chain_ending_at(126, "L"), "fab", chain_types={"L": "L"},
        provenance="experimental")
    check(any("C(IMGT 126" in t for t in quality["truncated_termini"]),
          "a VL truncated at IMGT 126 must be reported: {}".format(
              quality["truncated_termini"]))


def test_hbond_requires_donor_acceptor_polarity():
    """D...D and K...K contacts are not hydrogen bonds."""
    from tap_core import Residue

    def res(resseq, aa, atoms, imgt):
        r = Residue(chain="H", resseq=resseq, icode=" ", resname="XXX",
                    aa=aa, atoms=atoms)
        r.imgt_pos = imgt
        return r

    acid = [(" OD1", "O", 0.0, 0.0, 0.0), (" OD2", "O", 0.0, 1.0, 0.0),
            (" CG ", "C", -1.5, 0.5, 0.0), (" CB ", "C", -2.5, 1.5, 0.0)]
    acid2 = [(" OD1", "O", 3.0, 0.0, 0.0), (" OD2", "O", 3.0, 1.0, 0.0),
             (" CG ", "C", 4.5, 0.5, 0.0), (" CB ", "C", 5.5, 1.5, 0.0)]
    base = [(" NZ ", "N", 0.0, 0.0, 0.0), (" CE ", "C", -1.5, 0.5, 0.0),
            (" CB ", "C", -2.5, 1.5, 0.0)]
    base2 = [(" NZ ", "N", 3.0, 0.0, 0.0), (" CE ", "C", 4.5, 0.5, 0.0),
             (" CB ", "C", 5.5, 1.5, 0.0)]
    tyr = [(" OH ", "O", 0.0, 0.0, 0.0), (" CZ ", "C", -1.2, 0.6, 0.0),
           (" CB ", "C", -2.4, 0.2, 0.6)]

    dd = [res(10, "D", acid, 33), res(20, "D", acid2, 50)]
    kk = [res(11, "K", base, 34), res(21, "K", base2, 51)]
    ty = [res(12, "Y", tyr, 35), res(22, "D", acid2, 52)]

    for label, pair in (("D...D", dd), ("K...K", kk)):
        salt, _pairs, _d = pi_core.detect_salt_bridges(pair)
        flagged, _dist = pi_core.detect_hbonds(pair, salt)
        check(not flagged,
              "{} is not a hydrogen bond (no donor/acceptor complement)".format(label))
    salt, _pairs, _d = pi_core.detect_salt_bridges(ty)
    flagged, _dist = pi_core.detect_hbonds(ty, salt)
    check(flagged, "Tyr-OH donor -> Asp acceptor must be detected")


def test_conserved_cys_requires_both_partners():
    """A lone conserved Cys is a free thiol, not a disulfide partner."""
    full = VL_4D5
    broken = full.replace("YFCQQHYTTPPTFGQGTKVEIK", "YFAQQHYTTPPTFGQGTKVEIK")
    found, warnings = pi_core.find_disulfide_cys_sequence(
        [ChainInfo(name="VL", chain_type="L", sequence=full)])
    check(len(found) == 2, "both conserved VL Cys should be excluded")
    found_broken, warnings_broken = pi_core.find_disulfide_cys_sequence(
        [ChainInfo(name="VL", chain_type="L", sequence=broken)])
    check(len(found_broken) == 0,
          "a lone conserved Cys must stay ionizable, got {}".format(found_broken))
    check(warnings_broken, "an unpaired conserved Cys must produce a warning")


def test_propka_parser_handles_blank_chain_column():
    """PROPKA may leave the chain column empty; the pKa must not be misread."""
    text = """
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa   ligand atom-type
   ASP  25      5.07       3.80
   N+    1      8.96       8.00
"""
    side, model, termini, termini_model = pi_core._parse_propka_rows(text)
    check(side[("", 25, " ")] == 5.07,
          "the pKa must be read as 5.07, not the model-pKa: {}".format(side))
    check(model[("", 25, " ")] == 3.80,
          "the model-pKa must be captured separately: {}".format(model))
    check(termini[("", "N")] == 8.96, "terminal pKa with a blank chain column")


def test_region_and_patch_id_robustness():
    """Out-of-range IMGT positions and duplicate patch IDs are handled."""
    check(pi_core._region_imgt(200) == "unknown",
          "an out-of-range position must not be called FR4")
    check(pi_core._region_imgt(128) == "FR4", "FR4 range")
    check(pi_core._region_imgt(None) == "", "unknown position")


def test_patch_ids_are_unique_across_ph():
    """Patch IDs must be unique in the JSON/CSV output."""
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "fv.pdb")
        _write_synthetic_fv(pdb)
        result = analyze_structure(pdb, structure_model="heuristic")
    ids = [p["patch_id"] for p in result["patches"]]
    check(len(ids) == len(set(ids)),
          "patch IDs must be unique across pH values: {}".format(ids))


def test_anchors_are_never_proposed_for_mutation():
    """IMGT anchors are excluded from mutation candidates, not just flagged."""
    from tap_core import Residue
    anchor = Residue(chain="H", resseq=26, icode=" ", resname="LYS", aa="K",
                      atoms=[])
    anchor.imgt_pos = 26
    anchor.is_anchor = True
    anchor.rel_sasa = 50.0
    patch = pi_core.Patch(patch_id=1, sign="positive", ph=6.0)
    patch.residues = [{"label": "H26", "aa": "K", "region": "FR1",
                       "rel_sasa": 50.0, "charge": 1.0,
                       "in_cdr_vicinity": True}]
    candidates = pi_core.mutation_candidates(patch, [anchor], kind="fab")
    check(not candidates,
          "an IMGT anchor must not be proposed for mutation")


def test_patch_calibration_file_overrides_bands():
    """A calibration file replaces the default patch risk bands."""
    cal = pi_core.PatchCalibration(
        pos={"elevated": 1.0, "high": 4.0},
        neg={"elevated": 1.2, "high": 5.0},
        cdr_factor=1.0, source="test")
    risk, basis = pi_core._patch_risk(2.5, "positive", False, 3.0,
                                      calibration=cal)
    check(risk == "elevated",
          "calibrated bands must be applied, got {}".format(risk))
    default_risk, _ = pi_core._patch_risk(2.5, "positive", False, 3.0)
    check(default_risk == "high",
          "the default bands must still say high for 2.5")


def test_charge_descriptor_semantics():
    """TAP thresholds must only be applied to the surface descriptors."""
    chains = [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
              ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)]
    seq_result = analyze_sequence(chains, pka_set="IPC_protein")
    fv = seq_result["fv"]
    check(fv["SFvCSP_surface"] is None,
          "sequence mode has no surface SFvCSP")
    check(fv["SFvCSP_basis"] == "total",
          "sequence mode must record that the alias holds the total value")
    check(fv["SFvCSP_total"] is not None, "total SFvCSP must be reported")


def test_sequence_mode_declares_what_it_cannot_compute():
    """A sequence-only run must state which outputs need a structure.

    Regression: the sequence-only report simply omitted the TAP flags, the
    charge patches and the structure-aware pI, which reads as "no flags
    raised" rather than "not computed".  The result now carries an explicit
    inventory and the reports say so.
    """
    from tools.pi_analysis import pi_analyzer
    chains = [ChainInfo(name="VH", chain_type="H", sequence=VH_4D5),
              ChainInfo(name="VL", chain_type="L", sequence=VL_4D5)]
    result = analyze_sequence(chains, pka_set="IPC_protein")

    check(result["analysis_type"] == "sequence_only",
          "sequence mode must declare its analysis type")
    missing = result["unavailable_metrics"]
    check(isinstance(missing, list) and missing,
          "sequence mode must list the metrics it cannot compute")
    keys = {m["key"] for m in missing}
    for expected in ("tap", "patches", "structure_pI"):
        check(expected in keys,
              "{} must be listed as structure-dependent".format(expected))
    for item in missing:
        check(item.get("label") and item.get("reason"),
              "every entry needs a label and a reason: {}".format(item))

    # The absence must also be stated in the human-readable reports.
    text = pi_analyzer.format_text_report(result)
    check("SCOPE OF THIS ANALYSIS" in text,
          "the text report must contain a scope block")
    check("NOT computed" in text,
          "the scope block must say the metrics were not computed")
    check("not a clean flag" in text.lower(),
          "the scope block must warn that an absent flag is not a clean flag")
    check("TAP flags" in text,
          "the scope block must name the TAP metrics")

    # A machine-readable warning must exist too.
    scope_recs = [r for r in result["recommendations"]
                   if r["category"] == "scope"]
    check(any("NOT computed" in r["message"] for r in scope_recs),
          "a sequence-mode recommendation must state the omission")

    # Structure mode must not carry the warning.
    with tempfile.TemporaryDirectory() as tmp:
        pdb = os.path.join(tmp, "fv.pdb")
        _write_synthetic_fv(pdb)
        sres = analyze_structure(pdb, structure_model="heuristic")
    check(sres["analysis_type"] == "structure",
          "structure mode must declare its analysis type")
    check(sres["unavailable_metrics"] == [],
          "structure mode has no structure-dependent omission")
    check("SCOPE OF THIS ANALYSIS" not in pi_analyzer.format_text_report(sres),
          "the structure report must not carry the sequence-mode scope block")


def test_multichar_cif_chain_ids_are_handled():
    """mmCIF chain IDs longer than one character must not corrupt the PDB.

    Regression: the SASA work PDB was written with the raw chain ID, so a
    two-character ID pushed the resseq and coordinate columns out of
    alignment and the structure could not be analysed at all.
    """
    cols = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
            "label_comp_id", "label_asym_id", "label_entity_id",
            "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y",
            "Cartn_z", "occupancy", "B_iso_or_equiv", "auth_seq_id",
            "auth_comp_id", "auth_asym_id", "pdbx_PDB_model_num"]
    header = "data_test\nloop_\n" + "".join(
        "_atom_site.{}\n".format(c) for c in cols)
    rows = []
    serial = 1
    for chain_id, seq in (("HHL", VH_4D5), ("LGT", VL_4D5)):
        for i, aa in enumerate(seq, 1):
            for name, element, dx, dy, dz in (
                    ("CA", "C", 0, 0, 0), ("N", "N", -1.2, 0.5, 0.1),
                    ("C", "C", 1.2, 0.5, -0.1), ("O", "O", 1.7, 1.6, -0.1)):
                rows.append(" ".join([
                    "ATOM", str(serial), element, name, ".", ONE_TO_THREE[aa],
                    chain_id, "1", str(i), "?",
                    "{:.3f}".format((i % 20) * 3.8 + dx),
                    "{:.3f}".format((i // 20) * 3.8 + dy),
                    "{:.3f}".format(((i * 7) % 9) * 0.9 + dz),
                    "1.00", "20.00", str(i), ONE_TO_THREE[aa], chain_id, "1"]))
                serial += 1
    with tempfile.TemporaryDirectory() as tmp:
        cif = os.path.join(tmp, "multi.cif")
        with open(cif, "w") as fh:
            fh.write(header + "\n".join(rows) + "\n#\n")
        result = analyze_structure(cif, pka_set="IPC_protein",
                                   structure_model="heuristic")
    check(result["input"]["chains"] == ["HHL", "LGT"],
          "both multi-character chains must be found: {}".format(
              result["input"]["chains"]))
    check(result["structure_quality"]["n_residues"] > 0,
          "residues must be parsed from a multi-character chain file")
    check(result["fv"]["pI_selected"] is not None,
          "pI must be computed for a multi-character chain file")


def run_all_tests():
    tests = [
        test_pka_tables_and_aliases,
        test_sequence_pi_matches_biopython,
        test_sequence_charge_descriptors,
        test_vhh_sequence_mode,
        test_disulfide_assignment,
        test_fasta_chain_order,
        test_propka_parser,
        test_structure_end_to_end,
        test_patch_candidate_filters,
        test_reports_and_cli,
        test_structure_model_none,
        test_termini_and_charge_curve,
        test_propka_termini_parser,
        test_patch_risk_calibration,
        # regression tests for the scientific review
        test_propka_sensitivity_is_not_degenerate,
        test_patch_graph_uses_documented_8A_cutoff,
        test_same_chain_for_vh_and_vl_is_rejected,
        test_text_report_header_matches_columns,
        test_ph_curve_reaches_saved_text_report,
        test_clean_sequence_warnings_are_reported,
        test_bfactor_provenance_semantics,
        test_terminus_completeness_is_chain_type_specific,
        test_hbond_requires_donor_acceptor_polarity,
        test_conserved_cys_requires_both_partners,
        test_propka_parser_handles_blank_chain_column,
        test_region_and_patch_id_robustness,
        test_patch_ids_are_unique_across_ph,
        test_anchors_are_never_proposed_for_mutation,
        test_patch_calibration_file_overrides_bands,
        test_charge_descriptor_semantics,
        test_multichar_cif_chain_ids_are_handled,
        test_sequence_mode_declares_what_it_cannot_compute,
    ]
    passed = failed = 0
    for test in tests:
        try:
            test()
            print("PASS: {}".format(test.__name__))
            passed += 1
        except Exception as exc:
            print("FAIL: {}: {}".format(test.__name__, exc))
            import traceback
            traceback.print_exc()
            failed += 1
    print("\n{} passed, {} failed".format(passed, failed))
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all_tests() else 1)
