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

import pi_core  # noqa: E402
import pi_analyzer  # noqa: E402
from pi_core import (  # noqa: E402
    ChainInfo,
    Patch,
    analyze_sequence,
    analyze_structure,
    isoelectric_point,
    mutation_candidates,
    net_charge,
    parse_propka_pka,
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
    check(result["fv"]["Fv_chml"] is not None, "Fv_chml missing")
    check(result["fv"]["SFvCSP"] is not None, "SFvCSP missing")
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
        ],
    )
    cands = mutation_candidates(
        patch, residues,
        interface_keys={("H", 15, " ")},
        vhh_hallmark_keys={("H", 14, " ")})
    labels = {c["label"] for c in cands}
    check("H10" in labels, "framework acid should be a candidate")
    check("H11" not in labels, "positive residue must not be mutated in a negative patch")
    check("H12" not in labels, "CDR residue must be excluded")
    check("H13" not in labels, "FR4 residue must be excluded")
    check("H14" not in labels, "VHH hallmark must be excluded")
    check("H15" in labels, "vicinity acid should remain a candidate")
    h15 = [c for c in cands if c["label"] == "H15"][0]
    check(h15["risk"] == "high", "interface candidate should be high risk")

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
        check("pI_free" in header and "Fv_chml" in header, "CSV header")

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
