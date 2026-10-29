#!/usr/bin/env python3
"""Regression tests for tools/cdr_homology (plain script; exits nonzero on failure).

Covers the bugs fixed in the 2.0 rewrite:
  * CDR3 uses local alignment; identity is not diluted by length mismatch
  * ANARCI return-shape probing (legacy + modern)
  * bispecific chains are searched and tagged, no loop truncation
  * --top validation, deterministic tie-break
  * collision-free CSV names, Markdown pipe escaping
  * VHH detection and query naming
  * provenance in the report / JSON
  * ANARCI failure surfaces as a non-zero exit code
  * create_paired_fasta index isolation
"""

from __future__ import annotations

import importlib.util
import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOL = os.path.join(ROOT, "tools", "cdr_homology")
sys.path.insert(0, TOOL)

FAILURES = []


def check(group, name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail and not condition else ""))
    if not condition:
        FAILURES.append(f"{group}: {name}")


def load(name, filename):
    path = os.path.join(TOOL, filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


core = load("cdr_core", "cdr_core.py")
query = load("cdr_query", "query.py")
paired = load("cdr_paired", "create_paired_fasta.py")
preprocess = load("cdr_preprocess", "preprocess.py")


VH = ("EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTK"
      "YNQKFKGKATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS")
VL = ("DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPS"
      "RFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK")
VHH = ("QVQLVESGGGLVQPGGSLRLSCAASGRTFSQYAMGWFRQAPGKEREFVAAISWSGGSTY"
       "YADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA")


# ---------------------------------------------------------------------------
def test_alignment_modes():
    print("alignment modes (CDR3 local, CDR1/2 global)")
    # A query CDR3 found intact inside a longer subject CDR3.
    local = core.align_sequences("YYDILTNY", "YYDILTNYGFV", "local")
    glob = core.align_sequences("YYDILTNY", "YYDILTNYGFV", "global")
    check("align", "local identity is 100% for an intact 8-residue core",
          abs(local["identity"] - 1.0) < 1e-9, str(local["identity"]))
    check("align", "global identity is 8/11 for the same pair",
          abs(glob["identity"] - round(8 / 11, 4)) < 1e-9, str(glob["identity"]))
    check("align", "local coverage is 1.0 when the whole query matches",
          abs(local["coverage"] - 1.0) < 1e-9, str(local["coverage"]))
    check("align", "local mode is recorded in the result", local["mode"] == "local")
    check("align", "global mode is recorded in the result", glob["mode"] == "global")

    # Subject shorter than query: coverage < 1, identity still over the core.
    short = core.align_sequences("YYDILTNYGFV", "YYDILTNY", "local")
    check("align", "coverage reflects the matched fraction of the query",
          abs(short["coverage"] - round(8 / 11, 4)) < 1e-9, str(short["coverage"]))
    check("align", "identity_over_query is matches/len(query)",
          abs(short["identity_over_query"] - round(8 / 11, 4)) < 1e-9,
          str(short["identity_over_query"]))

    check("align", "CDR3 defaults to local", core.cdr_mode("CDR3") == "local")
    check("align", "CDR1 defaults to global", core.cdr_mode("CDR1") == "global")
    check("align", "an explicit mode overrides the default",
          core.cdr_mode("CDR3", "global") == "global")
    check("align", "empty input yields zero metrics",
          core.align_sequences("", "AAAA")["identity"] == 0.0)

    # The aligner cache must be reused, not rebuilt per comparison.
    first = core._get_aligner("global")
    check("align", "aligner is cached per mode", core._get_aligner("global") is first)


def test_anarci_shapes():
    print("ANARCI return-shape probing")
    legacy = [([((1, " "), "Q"), ((2, " "), "V")], 0, 1)]
    found = core.extract_numbering(legacy)
    check("anarci", "legacy shape (anarci 0.1.x) is recognised",
          found == [((1, " "), "Q"), ((2, " "), "V")], repr(found))

    modern = ([[([((1, " "), "Q")], 0, 0)]], [], [])
    found = core.extract_numbering(modern)
    check("anarci", "modern shape (anarci >= 2024) is recognised",
          found == [((1, " "), "Q")], repr(found))

    check("anarci", "garbage input returns None",
          core.extract_numbering([None, [], []]) is None)

    if core.anarci_available():
        posmap = core.number_sequence_imgt(VH, "H")
        check("anarci", "a real VH numbers to a non-empty posmap",
              bool(posmap) and len(posmap) > 100, str(len(posmap or {})))
        check("anarci", "chain type is detected as heavy",
              core.detect_chain_type(VH) == "H", str(core.detect_chain_type(VH)))
        check("anarci", "chain type is detected as light",
              core.detect_chain_type(VL) == "L", str(core.detect_chain_type(VL)))


def test_cdr_extraction():
    print("IMGT CDR extraction")
    cdrs = core.extract_cdrs(VH, "H")
    check("cdr", "VH CDR3 matches the 4D5 reference",
          cdrs["CDR3"] == "SRWGGDGFYAMDY", cdrs["CDR3"])
    check("cdr", "VH CDR1 is the IMGT 27-38 window",
          cdrs["CDR1"] == "GYTFTDTY", cdrs["CDR1"])
    check("cdr", "VH CDR2 is the IMGT 56-65 window",
          cdrs["CDR2"] == "INPYNGVT", cdrs["CDR2"])
    cdrs_l = core.extract_cdrs(VL, "L")
    check("cdr", "VL CDR3 extracted", cdrs_l["CDR3"] == "QQHYTTPPT", cdrs_l["CDR3"])
    check("cdr", "CDR3 insertion order is preserved (111/112 walk)",
          core.extract_cdrs_from_posmap(
              {"H111": "A", "H111A": "B", "H112C": "C", "H112B": "D",
               "H112A": "E", "H112": "F"}, "H")["CDR3"] == "ABCDEF")


def test_vhh_detection():
    print("VHH detection and query kind")
    posmap = core.number_sequence_imgt(VHH, "H")
    is_vhh, count, matched = core.detect_vhh(posmap, "H")
    check("vhh", "camelid hallmark detected (>=3 of 4)", is_vhh, str(matched))
    check("vhh", "human VH is not called a VHH",
          core.detect_vhh(core.number_sequence_imgt(VH, "H"), "H")[0] is False)
    kind, _ = core.chain_kind(VHH, "")
    check("vhh", "single camelid chain classified as vhh", kind == "vhh", kind)
    kind, _ = core.chain_kind(VH, "")
    check("vhh", "single non-camelid chain classified as vhh_suspect",
          kind == "vhh_suspect", kind)
    kind, _ = core.chain_kind(VH, VL)
    check("vhh", "VH+VL classified as fab", kind == "fab", kind)


def test_names():
    print("naming and escaping")
    check("name", "chain suffix and description stripped",
          core.base_name("4D5_VH mouse anti-HER2 heavy chain") == "4D5",
          core.base_name("4D5_VH mouse anti-HER2 heavy chain"))
    check("name", "VHH suffix stripped", core.base_name("1BZQ_VHH") == "1BZQ",
          core.base_name("1BZQ_VHH"))
    check("name", "ladder index kept",
          core.base_name("AMG110_V0_VH") == "AMG110_V0",
          core.base_name("AMG110_V0_VH"))

    used = set()
    a = core.safe_filename("A/B", used)
    b = core.safe_filename("A_B", used)
    check("name", "colliding sanitised names are de-duplicated", a != b, f"{a} vs {b}")
    long_a = "X" * 60 + "AAA"
    long_b = "X" * 60 + "BBB"
    used2 = set()
    check("name", "truncated long names stay distinct",
          core.safe_filename(long_a, used2) != core.safe_filename(long_b, used2))

    check("name", "pipe is escaped for Markdown",
          core.markdown_cell("A|B") == "A\\|B", core.markdown_cell("A|B"))
    check("name", "newlines are flattened for Markdown",
          core.markdown_cell("A\nB") == "A B")


def _tiny_db():
    return [
        {"name": "Alpha", "format": "Fab", "genetics": "", "target": "",
         "highest_trial": "",
         "chains": {"VH": {"CDR1": "GYTFTDTY", "CDR2": "INPYNGVT",
                           "CDR3": "SRWGGDGFYAMDY"}},
         "bispecific": {"VH": {"CDR1": "GYTFTDTY"}}},
        {"name": "Beta", "format": "Fab", "genetics": "", "target": "",
         "highest_trial": "",
         "chains": {"VH": {"CDR1": "AAAAAAAA", "CDR2": "BBBBBBBB",
                           "CDR3": "CCCCCCCC"}}},
        {"name": "Gamma|Pipe", "format": "Fab", "genetics": "", "target": "",
         "highest_trial": "",
         "chains": {"VH": {"CDR1": "GGGGGGGG", "CDR2": "HHHHHHHH",
                           "CDR3": "IIIIIIII"}}},
    ]


def test_search():
    print("search behaviour")
    db = _tiny_db()
    res = query.search_cdr("SRWGGDGFYAMDY", db, "VH", "CDR3", top_n=5)
    names = [h["therapeutic"] for h in res["top_matches"]]
    check("search", "every entry is scanned (no loop truncation)",
          {"Alpha", "Beta", "Gamma|Pipe"} <= set(names), str(names))
    check("search", "the exact match ranks first", names[0] == "Alpha", str(names))
    res_cdr1 = query.search_cdr("GYTFTDTY", db, "VH", "CDR1", top_n=5)
    check("search", "candidate count includes bispecific chains",
          res_cdr1["n_candidates"] == 4, str(res_cdr1["n_candidates"]))
    check("search", "bispecific hits are tagged",
          any(h["chain_source"] == "bispecific"
              for h in res_cdr1["top_matches"]))
    check("search", "background median is reported",
          res["background_median_identity"] is not None)
    check("search", "top-hit percentile is reported",
          res["top_match_percentile"] is not None)

    # Deterministic ordering for ties.
    tie_db = [{"name": n, "chains": {"VH": {"CDR3": "AAAA"}}} for n in ("Z", "A", "M")]
    first = [h["therapeutic"] for h in
             query.search_cdr("AAAA", tie_db, "VH", "CDR3")["top_matches"]]
    check("search", "ties are broken deterministically by name",
          first == sorted(first), str(first))

    # Thresholds filter and are counted.
    filtered = query.search_cdr("SRWGGDGFYAMDY", db, "VH", "CDR3",
                                min_identity=0.9)
    check("search", "min_identity filters weak hits",
          all(h["identity"] >= 0.9 for h in filtered["top_matches"]))
    check("search", "filtered candidates are counted",
          filtered["n_filtered"] > 0, str(filtered["n_filtered"]))

    # --top must be a positive integer.
    for bad in ("0", "-1"):
        try:
            query._positive_int(bad)
            ok = False
        except Exception:
            ok = True
        check("search", f"--top {bad} is rejected", ok)
    check("search", "--top 3 is accepted", query._positive_int("3") == 3)


def test_reports():
    print("report generation")
    db = _tiny_db()
    results = query.run_query([("4D5", VH, VL, {})], db, top_n=3)
    tmp = tempfile.mkdtemp()

    class Args:
        database = os.path.join(TOOL, "..", "..", "data", "therasabdab_cdrs.json")
        top = 3
        align = "auto"
        min_score = 0.0
        min_identity = 0.0

    meta = query.build_metadata(Args(), db, {"source_file": "x.xlsx",
                                             "source_sha256": "abc"})
    md = query.generate_report(results, db, {"total_entries": len(db),
                                             "format_counts": {"Fab": 3}},
                               meta, tmp, 3)
    text = open(md, encoding="utf-8").read()
    check("report", "provenance section present", "## 1. Provenance" in text)
    check("report", "database SHA-256 recorded", "SHA-256" in text)
    check("report", "method section present", "## 2. Method" in text)
    check("report", "CDR3 mean column present", "CDR3 mean id" in text)
    check("report", "background median column present", "Background median" in text)
    check("report", "appendix G is populated",
          "## G. Top 30 Drugs by Entry Count" in text and "| Alpha | 1 |" in text)
    check("report", "pipe in a therapeutic name is escaped",
          "Gamma\\|Pipe" in text and "| Gamma|Pipe |" not in text)

    # CSV names must not collide.
    csv_results = [
        {"name": "A/B", "matches": {}},
        {"name": "A_B", "matches": {}},
    ]
    query.write_csv_details(csv_results, tmp)
    files = os.listdir(os.path.join(tmp, "detail"))
    check("report", "CSV filenames are collision-free", len(files) == 2, str(files))

    # JSON round-trip with provenance.
    payload = {"metadata": meta, "queries": results}
    path = core.write_json(os.path.join(tmp, "r.json"), payload)
    loaded = json.load(open(path, encoding="utf-8"))
    check("report", "JSON carries metadata", "metadata" in loaded)
    check("report", "JSON metadata has parameters",
          "parameters" in loaded["metadata"])


def test_cli_exit_codes():
    print("CLI exit codes")
    tmp = tempfile.mkdtemp()
    db = os.path.join(TOOL, "..", "..", "data", "therasabdab_cdrs.json")
    if not os.path.exists(db):
        print("  [SKIP] database not present")
        return

    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        code = query.main(["--vh", VH, "--vl", VL, "--name", "ok",
                           "--top", "2", "--output", tmp, "--json"])
    check("cli", "successful run exits 0", code == 0, str(code))

    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        code = query.main(["--vh", "XXXXXNOTASEQUENCEXXXXX", "--name", "bad",
                           "--top", "2", "--output", tmp, "--json"])
    check("cli", "un-numberable query exits 2", code == 2, str(code))
    check("cli", "JSON is still written on partial failure",
          os.path.exists(os.path.join(tmp, "cdr_homology_results.json")))

    buf = io.StringIO()
    with redirect_stdout(buf), redirect_stderr(buf):
        code = query.main(["--vh", VH, "--database", "/nonexistent.json",
                           "--output", tmp])
    check("cli", "missing database exits 1", code == 1, str(code))


def test_paired_fasta():
    print("create_paired_fasta grouping")
    records = [
        ("AMG110.H|H_V0|desc", "AAAA"),
        ("OTHER.H|H_V0|desc", "BBBB"),
        ("AMG110.L|L_V0|desc", "CCCC"),
        ("OTHER.L|L_V0|desc", "DDDD"),
    ]
    errors = []
    groups = paired.build_pairs(records, errors)
    check("paired", "two antibodies sharing a ladder index stay separate",
          groups.get("AMG110_V0", {}).get("H") == "AAAA"
          and groups.get("OTHER_V0", {}).get("H") == "BBBB", str(groups))
    check("paired", "no spurious errors for a clean file", errors == [], str(errors))

    errors = []
    groups = paired.build_pairs([("A_VH", "AAAA"), ("A_VH", "BBBB")], errors)
    check("paired", "duplicate chain is reported, not silently dropped",
          len(errors) == 1 and groups["A"]["H"] == "AAAA", str(errors))

    tmp = tempfile.mkdtemp()
    path = paired.write_paired_fasta(
        paired.build_pairs(records), os.path.join(tmp, "p.fa"), "AMG110")
    headers = [line[1:].strip() for line in open(path, encoding="utf-8")
               if line.startswith(">")]
    check("paired", "output headers follow the documented form",
          "AMG110_V0_VH" in headers and "AMG110_V0_VL" in headers, str(headers))


def test_statistics():
    print("database statistics")
    db = [{"name": "a", "format": "Fab", "genetics": "Humanized", "vd_lc": "Kappa",
           "highest_trial": "Approved", "target": "HER2",
           "chains": {"VH": {"CDR1": "A" * n, "CDR2": "B", "CDR3": "C"}}}
          for n in range(1, 11)]
    stats = preprocess.build_statistics(db)
    check("stats", "median uses statistics.median (5.5 for 1..10)",
          stats["cdr_length_summary"]["VH_CDR1"]["median"] == 5.5,
          str(stats["cdr_length_summary"]["VH_CDR1"]["median"]))

    db_with_bispec = [dict(entry) for entry in db]
    db_with_bispec[0] = dict(db_with_bispec[0])
    db_with_bispec[0]["bispecific"] = {"VH": {"CDR1": "A", "CDR2": "B", "CDR3": "C"}}
    stats2 = preprocess.build_statistics(db_with_bispec)
    check("stats", "bispecific chains are included in the summary",
          stats2["cdr_length_summary"]["VH_CDR1"]["count"] == 11,
          str(stats2["cdr_length_summary"]["VH_CDR1"]["count"]))
    check("stats", "primary-only summary excludes bispecific chains",
          stats2["cdr_length_summary_primary"]["VH_CDR1"]["count"] == 10)
    check("stats", "schema version recorded", stats["schema_version"] == 2)
    check("stats", "NaN-like values are cleaned",
          preprocess._clean(float("nan")) == "" and preprocess._clean("NA") == "")


def main() -> int:
    for test in (test_alignment_modes, test_anarci_shapes, test_cdr_extraction,
                 test_vhh_detection, test_names, test_search, test_reports,
                 test_cli_exit_codes, test_paired_fasta, test_statistics):
        test()
        print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for failure in FAILURES:
            print("  -", failure)
        return 1
    print("ALL CDR-HOMOLOGY TESTS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
