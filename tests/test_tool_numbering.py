#!/usr/bin/env python3
"""Tool numbering/partitioning parity tests.

Run with `python3 tests/test_tool_numbering.py` (plain script, no pytest).

The analysis tools must use the same Kabat/IMGT numbering and region
partitioning as the main pipeline. These tests pin the regressions found in
the tool audit:

  * IMGT CDR3 insertion order (111, 111A.., 112C..112A, 112) was transposed
    by numeric iteration in cdr_homology / pataa_homology / humanness_score
  * humanness_score report CDR extraction now reads pipeline region labels
  * ptm_exposure region maps match the pipeline engines residue-for-residue
  * humanness_score IMGT scoring requires ANARCI instead of silently
    converting approximate Kabat numbering
"""
import importlib.util
import io
import json
import os
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from humanize.anarci_adapter import number_with_anarci_imgt
from humanize.numbering import number_heavy, number_light

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


def load_module(name, rel_path):
    path = os.path.join(ROOT, rel_path)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


VH_PREFIX = ("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQGLEWMGWINAGNGNTKYSQKFQGRV"
             "TITRDTSASTAYMELSSLRSEDTAVYYC")
VH_FR4 = "WGQGTLVTVSS"
VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG"
      "TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
LAMBDA = ("QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKS"
          "GTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL")
# Amostomig VH: IMGT CDR2 has 60A/60B and 61C/61B/61A insertions, which ANARCI
# emits in sequence order (60, 60A, 60B, 61C, 61B, 61A, 61) -- not numeric.
VH_CDR2_INS = ("QVQLVESGGGLVQPGGSLRLSCAASGFTFSRYWMYWVRQAPGKGLEWVSSIEDSSINTGGGTETTY"
               "YTDSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCAKGDYYCTTYECQHYYGMDYWGQGTTVTVSS")


def _vh_with_cdr3(n):
    body = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(n))
    return VH_PREFIX + body + VH_FR4


def _cdr_from_chain(chain, region):
    return "".join(r.aa for r in chain.residues if r.region == region)


def _imgt_posmap(seq, chain_type="H"):
    chain = number_with_anarci_imgt(seq, chain_type)
    return chain, chain.posmap()


def test_cdr_homology_cdr3_order():
    print("cdr_homology CDR3 follows ANARCI sequence order (111/112 insertions)")
    q = load_module("cdr_query", "tools/cdr_homology/query.py")
    p = load_module("cdr_preprocess", "tools/cdr_homology/preprocess.py")
    for n in (13, 14, 17, 20, 23):
        seq = _vh_with_cdr3(n)
        chain, posmap = _imgt_posmap(seq)
        expected = _cdr_from_chain(chain, "CDR3")
        got_q = q._extract_cdrs_from_posmap(posmap, "H")["CDR3"]
        got_p = p.extract_cdrs_from_posmap(posmap, "H")["CDR3"]
        check(f"query n={n}", got_q == expected, f"{got_q} != {expected}")
        check(f"preprocess n={n}", got_p == expected, f"{got_p} != {expected}")
        for region in ("CDR1", "CDR2"):
            want = _cdr_from_chain(chain, region)
            check(f"query {region} n={n}",
                  q._extract_cdrs_from_posmap(posmap, "H")[region] == want)


def test_cdr2_insertion_order():
    print("IMGT CDR2 insertion order (60A/60B, 61C..61A) is preserved")
    q = load_module("cdr_query_c2", "tools/cdr_homology/query.py")
    chain, posmap = _imgt_posmap(VH_CDR2_INS)
    expected = _cdr_from_chain(chain, "CDR2")
    got = q._extract_cdrs_from_posmap(posmap, "H")["CDR2"]
    check("cdr_homology", got == expected, f"{got} != {expected}")

    g = load_module("humanness_reports_c2",
                    "tools/humanness_score/generate_all_reports.py")
    got_g = g.get_cdr_sequences(VH_CDR2_INS, "H")["CDR2"]
    check("humanness report", got_g == expected, f"{got_g} != {expected}")


def test_pataa_cdr3_order():
    print("pataa_homology CDR3 follows ANARCI sequence order")
    q = load_module("pataa_query", "tools/pataa_homology/query_pataa.py")
    pre = load_module("pataa_preprocess", "tools/pataa_homology/preprocess_pataa.py")
    for n in (13, 17, 23):
        seq = _vh_with_cdr3(n)
        chain, posmap = _imgt_posmap(seq)
        expected = _cdr_from_chain(chain, "CDR3")
        got_q = q._extract_cdrs_from_posmap(posmap, "H")["CDR3"]
        check(f"query n={n}", got_q == expected, f"{got_q} != {expected}")
        records = pre._number_batch([("q", seq)])
        got_pre = records[0]["cdrs"]["CDR3"] if records else None
        check(f"preprocess n={n}", got_pre == expected, f"{got_pre} != {expected}")


def test_pataa_repair_stage():
    print("pataa repair re-extracts CDRs from stored numbering, then builds")
    pre = load_module("pataa_preprocess_repair",
                      "tools/pataa_homology/preprocess_pataa.py")
    seq = _vh_with_cdr3(17)
    chain = number_with_anarci_imgt(seq, "H")
    numbered = chain.posmap()
    correct = pre._extract_cdrs_from_posmap(numbered, "H")

    def old_extract(posmap, prefix):
        """Pre-fix behavior: numeric walk with ascending insertion letters."""
        out = {}
        for cdr, (start, end) in pre._IMGT_CDR_RANGES.items():
            frag = ""
            for pos in range(start, end + 1):
                aa = posmap.get(f"{prefix}{pos}")
                if aa:
                    frag += aa
                for letter in "ABCDEFGHIJK":
                    iaa = posmap.get(f"{prefix}{pos}{letter}")
                    if iaa:
                        frag += iaa
            out[cdr] = frag
        return out

    stale = old_extract(numbered, "H")
    check("test data actually stale", stale["CDR3"] != correct["CDR3"],
          f"{stale['CDR3']} == {correct['CDR3']}")

    tmp = tempfile.mkdtemp(prefix="pataa_repair_")
    paths = {
        "CDR_RECORDS": Path(tmp) / "cdr_records.jsonl",
        "POSITIVES_TSV": Path(tmp) / "positives.tsv",
        "INDEX_TSV": Path(tmp) / "index.tsv",
        "OUT_JSON": Path(tmp) / "pataa_cdrs.json",
        "OUT_STATS": Path(tmp) / "pataa_stats.json",
    }
    long_numbered = dict(numbered)
    for letter in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        long_numbered[f"H111{letter}"] = "A"
        long_numbered[f"H112{letter}"] = "A"
    with open(paths["CDR_RECORDS"], "w") as fh:
        fh.write(json.dumps({
            "qid": "p000000001", "chain": "VH", "sequence": seq,
            "numbered": numbered, "cdrs": stale,
        }) + "\n")
        fh.write(json.dumps({
            "qid": "p000000002", "chain": "VH", "sequence": seq,
            "numbered": long_numbered, "cdrs": stale,
        }) + "\n")
    with open(paths["POSITIVES_TSV"], "w") as fh:
        fh.write("p000000001\t200.0\tIGHV1-69*01\t85.0\n")
        fh.write("p000000002\t150.0\tIGHV1-69*01\t85.0\n")
    with open(paths["INDEX_TSV"], "w") as fh:
        fh.write("p000000001\tSequence 1 from patent US 1234567\n")
        fh.write("p000000002\tSequence 2 from patent US 1234567\n")

    saved = {k: getattr(pre, k) for k in paths}
    try:
        for k, v in paths.items():
            setattr(pre, k, v)
        pre.stage_repair()
        with open(paths["CDR_RECORDS"]) as fh:
            repaired = json.loads(fh.readline())
        check("JSONL repaired",
              repaired["cdrs"]["CDR3"] == correct["CDR3"],
              f"{repaired['cdrs']['CDR3']} != {correct['CDR3']}")
        with open(paths["OUT_JSON"]) as fh:
            data = json.load(fh)
        check("built DB count", len(data) == 1)
        got = data[0]["chains"]["VH"]["CDR3"]
        check("built DB CDR3", got == correct["CDR3"], f"{got} != {correct['CDR3']}")
        with open(paths["OUT_STATS"]) as fh:
            stats = json.load(fh)
        check("long CDR3 filtered",
              stats["skipped_long_cdr3"] == 1)
        check("dedup input count",
              stats["total_positive_records"] == 1)
    finally:
        for k, v in saved.items():
            setattr(pre, k, v)
        shutil.rmtree(tmp, ignore_errors=True)


def test_pataa_fresh_setup():
    print("pataa fresh-machine setup: tool fallback + archive auto-extract")
    pre = load_module("pataa_preprocess_fresh",
                      "tools/pataa_homology/preprocess_pataa.py")
    tmp = tempfile.mkdtemp(prefix="pataa_fresh_")
    try:
        found = pre._tool_path(Path(tmp) / "does-not-exist", "python3")
        check("tool falls back to PATH", found.exists(), str(found))

        db = Path(tmp) / "database" / "pataa"
        db.parent.mkdir(parents=True)
        with tarfile.open(db.parent / "pataa.tar.gz", "w:gz") as tf:
            for member, data in (("pataa.pin", b"pin"), ("pataa.psq", b"psq"),
                                 ("taxdb.btd", b"tax")):
                info = tarfile.TarInfo(member)
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
        check("db extracted", pre.ensure_pataa_db(db) == db)
        check("pin present", db.with_suffix(".pin").exists())
        check("taxdb present", (db.parent / "taxdb.btd").exists())
        pre.ensure_pataa_db(db)  # idempotent

        missing = Path(tmp) / "missing" / "pataa"
        try:
            pre.ensure_pataa_db(missing)
        except RuntimeError as exc:
            check("missing archive message", "pataa.tar.gz" in str(exc))
        else:
            check("missing archive raises", False)

        fake = Path(tmp) / "blastdbcmd"
        fake.write_text(
            '#!/usr/bin/env python3\nprint(">seq1\\nEVQLVESGGG")\n')
        fake.chmod(0o755)
        saved = (pre.FASTA, pre.WORK, pre.BLASTDBCMD)
        pre.FASTA = Path(tmp) / "pataa.fasta"
        pre.WORK = Path(tmp)
        pre.BLASTDBCMD = fake
        try:
            pre.stage_dump(db)
            check("dump wrote FASTA",
                  pre.FASTA.read_text() == ">seq1\nEVQLVESGGG\n")
            check("dump left no temp",
                  not (Path(tmp) / "pataa.fasta.tmp").exists())
        finally:
            pre.FASTA, pre.WORK, pre.BLASTDBCMD = saved
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_humanness_report_cdrs():
    print("humanness_score report CDRs match pipeline regions (IMGT + Kabat)")
    g = load_module("humanness_reports", "tools/humanness_score/generate_all_reports.py")

    for seq, chain_type in ((_vh_with_cdr3(17), "H"), (VL, "L")):
        imgt = number_with_anarci_imgt(seq, chain_type)
        got = g.get_cdr_sequences(seq, chain_type)
        for region in ("CDR1", "CDR2", "CDR3"):
            check(f"IMGT {chain_type} {region}",
                  got[region] == _cdr_from_chain(imgt, region),
                  f"{got[region]} != {_cdr_from_chain(imgt, region)}")

        kabat = number_heavy(seq) if chain_type == "H" else number_light(seq)
        got_k = g.get_kabat_cdr_sequences(seq, chain_type)
        for region in ("CDR1", "CDR2", "CDR3"):
            check(f"Kabat {chain_type} {region}",
                  got_k[region] == _cdr_from_chain(kabat, region),
                  f"{got_k[region]} != {_cdr_from_chain(kabat, region)}")


def test_ptm_exposure_regions():
    print("ptm_exposure Kabat/IMGT region maps match pipeline engines")
    p = load_module("ptm_analyze", "tools/ptm_exposure/analyze_ptm.py")

    for seq, chain_type in ((_vh_with_cdr3(17), "H"), (VL, "L"), (LAMBDA, "L")):
        rk, ri, detected, _ = p.get_chain_regions(seq, chain_type)
        kabat = number_heavy(seq) if chain_type == "H" else number_light(seq)
        imgt = number_with_anarci_imgt(seq, chain_type)
        check(f"{chain_type} detected", detected == kabat.chain_type,
              f"{detected} != {kabat.chain_type}")
        ok_k = all(rk[r.index + 1] == r.region for r in kabat.residues)
        ok_i = all(ri[r.index + 1] == r.region for r in imgt.residues)
        check(f"{chain_type} Kabat map", ok_k)
        check(f"{chain_type} IMGT map", ok_i)


def test_humanness_imgt_requires_anarci():
    print("humanness_score IMGT path refuses approximate fallback")
    e = load_module("humanness_eval", "tools/humanness_score/evaluate_humanness.py")
    original = e.number_with_anarci_imgt
    try:
        e.number_with_anarci_imgt = lambda seq, chain_type: None
        scored, numbered = e.evaluate_sequence_imgt("H", _vh_with_cdr3(13))
        check("returns no scores", scored is None and numbered is None,
              f"{scored}, {numbered}")
    finally:
        e.number_with_anarci_imgt = original


def test_humanness_match_counts_and_coverage():
    print("humanness_score reports exact match counts and coverage")
    e = load_module("humanness_eval_counts", "tools/humanness_score/evaluate_humanness.py")
    from humanize.sequences import parse_fasta

    records = parse_fasta(os.path.join(ROOT, "data/examples/mouse_4d5_fab.fasta"))
    vh = next(seq for name, seq in records if "VH" in name)
    scored, numbered = e.evaluate_sequence("H", vh)
    gene, scores = scored[0]
    q_map = numbered.posmap()
    g_map = gene.numbered.posmap()

    for region, key in (("FR1", "fr1"), ("FR2", "fr2"), ("FR3", "fr3"),
                        ("CDR1", "cdr1"), ("CDR2", "cdr2")):
        positions = [p for p in q_map if p in g_map and numbered.region_of(p) == region]
        expected = sum(1 for p in positions if q_map[p] == g_map[p])
        check(f"Kabat {region} match",
              scores[f"{key}_match"] == expected and scores[f"n_{key}"] == len(positions),
              f"{scores[f'{key}_match']}/{scores[f'n_{key}']} != {expected}/{len(positions)}")
    check("4D5 FR2 match is not truncated",
          scores["fr2_match"] == 8 and scores["n_fr2"] == 14,
          f"{scores['fr2_match']}/{scores['n_fr2']}")
    check("Kabat coverage fields",
          scores["n_query"] == 98 and scores["n_compared"] == 98 and scores["coverage"] == 1.0,
          f"{scores['n_query']}, {scores['n_compared']}, {scores['coverage']}")

    formatted = e.format_results("H", vh, scored, numbered, top_n=1)
    fr2_lines = [line for line in formatted.splitlines() if line.strip().startswith("FR2")]
    check("Kabat FR2 table shows match/total",
          fr2_lines and "8" in fr2_lines[0] and "14" in fr2_lines[0],
          fr2_lines[0] if fr2_lines else "missing FR2 line")

    scored_imgt, numbered_imgt = e.evaluate_sequence_imgt("H", vh)
    if scored_imgt:
        stats = scored_imgt[0][1].get("imgt_region_stats", {})
        check("IMGT stats have match fields",
              all("match" in s and 0 <= s["match"] <= s["count"] for s in stats.values()))
        check("IMGT coverage fields",
              "n_query" in scored_imgt[0][1] and "coverage" in scored_imgt[0][1])
    else:
        check("IMGT path available", False, "ANARCI/IMGT unavailable")


def test_humanness_true_imgt_alignment():
    print("humanness_score true IMGT alignment preserves insertion order")
    e = load_module("humanness_eval_imgt_align", "tools/humanness_score/evaluate_humanness.py")
    seq = VH_CDR2_INS
    scored, numbered = e.evaluate_sequence_imgt("H", seq)
    check("IMGT scores available", scored is not None)
    if not scored:
        return

    gene = scored[0][0]
    alignment = e._format_imgt_alignment(numbered, gene._imgt_posmap, "H")
    cdr2_expected = _cdr_from_chain(numbered, "CDR2")
    cdr2_lines = [line for line in alignment.splitlines() if line.strip().startswith("CDR2")]
    check("alignment has CDR2", bool(cdr2_lines), alignment)
    if cdr2_lines:
        query_part = cdr2_lines[0].split(":")[1].split("|")[0].strip()
        check("IMGT CDR2 is complete", query_part == cdr2_expected,
              f"{query_part} != {cdr2_expected}")
        pieces = [
            part.split(":", 1)[-1].strip() if ":" in part else part.strip()
            for part in cdr2_lines[0].split("|")
        ]
        check("alignment columns are equal length",
              len({len(piece) for piece in pieces}) == 1, cdr2_lines[0])

    kabat_scored, kabat_numbered = e.evaluate_sequence("H", seq)
    table = e.format_sequence_comparison_table(
        kabat_numbered, kabat_scored[0][0].numbered, "H")
    check("approximate IMGT table removed", "IMGT region boundaries" not in table)


def test_humanness_word_fallback_and_names():
    print("humanness_score Word fallback shows N/A and sanitizes names")
    g = load_module("humanness_reports_na", "tools/humanness_score/generate_all_reports.py")
    seq = _vh_with_cdr3(13)
    original = g.evaluate_sequence_imgt
    tmp = tempfile.mkdtemp(prefix="humanness_na_")
    try:
        g.evaluate_sequence_imgt = lambda *args, **kwargs: (None, None)
        results = g.evaluate_antibody("test", seq, "H")
        check("IMGT unavailable flag", results["imgt"]["available"] is False)
        check("IMGT best gene is N/A", results["imgt"]["best_gene"] == "N/A")

        out_path = os.path.join(tmp, "combined.docx")
        g.create_combined_report([("test", seq, "H", results)], out_path)
        doc = g.Document(out_path)
        imgt_row = doc.tables[1].rows[1]
        check("IMGT summary cells are N/A",
              all(cell.text == "N/A" for cell in imgt_row.cells[2:]),
              " | ".join(cell.text for cell in imgt_row.cells))
    finally:
        g.evaluate_sequence_imgt = original
        shutil.rmtree(tmp, ignore_errors=True)

    used = set()
    first = g._safe_report_stem("a/b c", used)
    second = g._safe_report_stem("a/b c", used)
    check("filename sanitized", first == "a_b_c", first)
    check("duplicate names unique", first != second, f"{first}, {second}")


def test_humanness_input_names_and_errors():
    print("humanness_score input loader deduplicates names and skips bad records")
    e = load_module("humanness_inputs", "tools/humanness_score/evaluate_humanness.py")
    tmp = tempfile.mkdtemp(prefix="humanness_input_")
    fasta_path = os.path.join(tmp, "mixed.fasta")
    with open(fasta_path, "w") as fh:
        fh.write(f">dup\n{_vh_with_cdr3(17)}\n")
        fh.write(">bad\nACDE\n")
        fh.write(f">dup\n{VL}\n")
    try:
        sequences, errors = e.load_sequences(fasta_path, skip_invalid=True)
        names = [name for name, _, _, _ in sequences]
        check("duplicate names made unique", names == ["dup", "dup_1"], str(names))
        check("bad record reported", len(errors) == 1 and errors[0][0] == "bad",
              str(errors))
        check("valid records classified",
              [chain for _, chain, _, _ in sequences] == ["H", "L"],
              str([(c, n) for n, c, _, _ in sequences]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    test_cdr_homology_cdr3_order()
    test_cdr2_insertion_order()
    test_pataa_cdr3_order()
    test_pataa_repair_stage()
    test_pataa_fresh_setup()
    test_humanness_report_cdrs()
    test_ptm_exposure_regions()
    test_humanness_imgt_requires_anarci()
    test_humanness_match_counts_and_coverage()
    test_humanness_true_imgt_alignment()
    test_humanness_word_fallback_and_names()
    test_humanness_input_names_and_errors()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES: {FAILURES}")
        return 1
    print("ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
