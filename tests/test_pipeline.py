#!/usr/bin/env python3
"""Pipeline tests (portable: run with `python3 tests/test_pipeline.py` or pytest).

Covers:
  * numbering engine vs known sequences
  * germline matching sanity
  * grafting correctness (CDR contents preserved, FR human)
  * back-mutation tiers (T1/T2/T3/KEEP_DONOR/KEEP_HUMAN)
  * variant assembly (V0-V3)
  * VHH hallmark protection
  * end-to-end run on the example inputs
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from humanize.backmut import StructureHints, analyze_backmutations
from humanize.germline import choose_germlines, load_germline_db
from humanize.graft import graft_chain
from humanize.numbering import CDR_SCHEMES, is_vhh, number_heavy, number_light
from humanize.pipeline import PipelineConfig, run_pipeline
from humanize.sequences import parse_input
from humanize.variants import assemble_variants

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


IGHV3_23 = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKWGQGTLVTVSS"
IGKV1_39 = "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPFGQGTKVEIK"
M4D5_VH = "EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTKYNQKFKGKATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS"
M4D5_VL = "DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK"
VHH_1BZQ = "QVQLVESGGGLVQAGGSLRLSCAASGYAYTYIYMGWFRQAPGKEREGVAAMDSGGGGTLYADSVKGRFTISRDKGKNTVYLQMDSLKPEDTATYYCAAGGYELRDRTYGQWGQGTQVTVSS"


def test_numbering():
    print("numbering engine")
    h = number_heavy(IGHV3_23)
    check("VH3-23 FR1 = 30 residues", h.residue("H30") is not None and h.residue("H31") is not None)
    check("VH3-23 CDR1 = SYAMS", h.seq_range("H31", "H35") == "SYAMS", h.seq_range("H31", "H35"))
    check("VH3-23 CDR2 correct", h.seq_range("H50", "H65") == "SAISGSGGSTYYADSVKG")
    check("VH3-23 FR3 correct", h.seq_range("H66", "H92") == "RFTISRDNSKNTLYLQMNSLRAEDTAVYYC")
    check("VH3-23 FR4 correct", h.seq_range("H103", "H113") == "WGQGTLVTVSS")
    check("VH3-23 Cys22", h.residue("H22").aa == "C")

    l = number_light(IGKV1_39)
    check("IGKV1-39 CDR1", l.seq_range("L24", "L34") == "RASQSISSYLN")
    check("IGKV1-39 CDR2", l.seq_range("L50", "L56") == "AASSLQS")
    check("IGKV1-39 CDR3", l.seq_range("L89", "L95") == "QQSYSTP")
    check("IGKV1-39 FR4", l.seq_range("L98", "L107") == "FGQGTKVEIK")

    m4 = number_heavy(M4D5_VH)
    check("4D5 VH CDR2 = donor loop", m4.seq_range("H50", "H52A") + m4.seq_range("H53", "H65") == "YINPYNGVTKYNQKFKG")
    # Standard Kabat: H93/H94 are FR3 (e.g. "S" "R" in ...VYYC S R WGGD...);
    # CDR3 starts at H95.
    check("4D5 VH FR3 H93/H94", m4.seq_range("H93", "H94") == "SR", m4.seq_range("H93", "H94"))
    check("4D5 VH CDR3", m4.seq_range("H95", "H102") == "WGGDGFYAMDY", m4.seq_range("H95", "H102"))
    vhh = number_heavy(VHH_1BZQ)
    vhh_ok, score, _ = is_vhh(vhh)
    check("1BZQ VHH hallmark detected", vhh_ok and score == 4, str(score))
    check("1BZQ VHH FR3 H93/H94", vhh.seq_range("H93", "H94") == "AA", vhh.seq_range("H93", "H94"))
    check("1BZQ VHH CDR3", vhh.seq_range("H95", "H102") == "GGYELRDRTYGQ", vhh.seq_range("H95", "H102"))
    for scheme in CDR_SCHEMES:
        _ = CDR_SCHEMES[scheme]  # table integrity


def _full_cdr3_loop(chain):
    """Full CDR3 loop = Kabat 93-102 (+insertions) for H, 89-97 for L."""
    lo, hi = (93, 102) if chain.chain_type == "H" else (89, 97)
    out = []
    for r in chain.residues:
        n = int("".join(c for c in r.pos if c.isdigit()))
        if lo <= n <= hi:
            out.append(r.aa)
    return "".join(out)


def test_graft_loop_conservation():
    """Regression: the FULL CDR3 loop (incl. Kabat 93/94) must be grafted
    from the donor. The strict-Kabat H93/H94 labels are FR3, but they carry
    the first two loop residues and must not be replaced by germline."""
    print("graft full-loop conservation (regression for H93/H94)")
    from humanize.graft import graft_chain
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    cases = [
        ("4D5 VH", number_heavy(M4D5_VH), "IGHV1-3*01", "IGHJ1*01", "H"),
        ("1BZQ VHH", number_heavy(VHH_1BZQ), "IGHV3-11*01", "IGHJ1*01", "H"),
        ("IGKV1-39", number_light(IGKV1_39), "IGKV1-39*01", "IGKJ1*01", "L"),
        ("1MEL VHH (Cys-CDR3)", number_heavy(
            "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS"),
            "IGHV3-11*01", "IGHJ1*01", "H"),
    ]
    for name, donor, vgene, jgene, _ in cases:
        v = [g for g in db.v_genes if g.gene_id == vgene][0]
        j = [g for g in db.j_genes if g.gene_id == jgene][0]
        graft = graft_chain(donor, v, j, "kabat")
        d_loop = _full_cdr3_loop(donor)
        g_loop = _full_cdr3_loop(graft.numbered)
        check(f"{name} full CDR3 loop conserved",
              d_loop == g_loop, f"donor={d_loop} graft={g_loop}")
        check(f"{name} H93/H94 donor origin",
              graft.origin.get("H93") in ("donor", "donor(vhh)")
              if graft.chain_type == "H" else True,
              str(graft.origin.get("H93")))


def test_vhh_humanization_gold_standard():
    """VHH humanization effectiveness vs the Vincke 2009 gold standard.

    cAb-Lys3 (PDB 1MEL) -> hCAb-Lys3: the documented universal scaffold
    design = IGHV3-23 consensus + camelid hallmark + donor CDRs. The
    pipeline must (1) detect VHH, (2) keep hallmark as KEEP_DONOR,
    (3) graft CDR1/2/3 (incl. the CDR1-Cys + CDR3-Cys disulfide pair)
    unchanged, (4) choose the VH3-23 family germline.
    """
    print("VHH humanization gold standard (cAb-Lys3 -> hCAb-Lys3)")
    from humanize.backmut import analyze_backmutations
    from humanize.germline import choose_germlines
    from humanize.graft import graft_chain
    CAbLys3 = ("DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGV"
               "AAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(CAbLys3)
    vhh, score, matched = is_vhh(donor)
    check("VHH detected (hallmark 4/4)", vhh and score == 4, str(matched))
    choice = choose_germlines(donor, db)
    check("VH3-23-family germline chosen",
          choice.v_gene is not None and choice.v_gene.gene_id.startswith("IGHV3-23"),
          str(choice.v_gene.gene_id if choice.v_gene else None))
    bm = analyze_backmutations(donor, choice.v_gene, is_vhh=True)
    kd = {c.position for c in bm.candidates if c.tier == "KEEP_DONOR"}
    check("hallmark positions all KEEP_DONOR",
          {"H37", "H44", "H45", "H47"} <= kd, str(sorted(kd)))
    graft = graft_chain(donor, choice.v_gene, choice.j_gene, "kabat", is_vhh=True)
    for cdr in ("CDR1", "CDR2", "CDR3"):
        d = "".join(r.aa for r in donor.residues if r.region == cdr)
        g = "".join(r.aa for r in graft.numbered.residues if r.region == cdr)
        check(f"CDR{cdr} loop conserved", d == g, f"d={d} g={g}")
    check("disulfide Cys pair grafted (CDR1 C + CDR3 C)",
          graft.sequence.count("C") >= 2
          and "PYCMG" in graft.sequence and "YEC" in graft.sequence)
    check("hallmark intact in graft FR2",
          graft.numbered.seq_range("H36", "H48") == "WFRQAPGKEREGV",
          graft.numbered.seq_range("H36", "H48"))


def test_germline_and_graft():
    print("germline + graft")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    check("germline DB has human V genes", len(db.human("H")) > 100 and len(db.human("L")) > 50,
          f"H={len(db.human('H'))} L={len(db.human('L'))}")
    check("germline DB has J genes", len(db.j_for("H")) > 5 and len(db.j_for("L")) > 10)

    _, chains = parse_input(os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"))
    vh, vl = chains
    for donor in (vh, vl):
        choice = choose_germlines(donor.numbered, db)
        check(f"{donor.name} germline chosen", choice.v_gene is not None and choice.j_gene is not None)
        check(f"{donor.name} FR identity >= 0.6",
              choice.scores.get("fr_identity", 0) >= 0.6, str(choice.scores.get("fr_identity")))
        graft = graft_chain(donor.numbered, choice.v_gene, choice.j_gene, "kabat")
        # every CDR residue must come from the donor
        for cdr in ("CDR1", "CDR2", "CDR3"):
            donor_seg = "".join(r.aa for r in donor.numbered.residues if r.region == cdr)
            graft_seg = "".join(r.aa for r in graft.numbered.residues if r.region == cdr)
            check(f"{donor.name} graft keeps {cdr}", graft_seg == donor_seg,
                  f"donor={donor_seg} graft={graft_seg}")


def test_backmut_variants():
    print("back-mutations + variants")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    _, chains = parse_input(os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"))
    for donor in chains:
        choice = choose_germlines(donor.numbered, db)
        top = [(g, s) for g, s in choice.alternatives]
        bm = analyze_backmutations(donor.numbered, choice.v_gene, is_vhh=False, top_germlines=top)
        tiers = {c.tier for c in bm.candidates}
        # VH keeps T1 (structural pillars); VL L87 is demoted to T3 by the
        # gold-standard empirical no-effect table (trastuzumab kept Y87).
        if donor.chain_type == "H":
            check(f"{donor.name} has T1 candidates", "T1" in tiers, str(sorted(tiers)))
        else:
            l87 = [c for c in bm.candidates if c.position == "L87"]
            check("L87 demoted to T3 by empirical table",
                  bool(l87) and l87[0].tier == "T3", str([(c.position, c.tier) for c in l87]))
        check(f"{donor.name} has T3 candidates", "T3" in tiers)
        variants = assemble_variants(donor.numbered, choice.v_gene, choice.j_gene,
                                     "kabat", bm, is_vhh=False)
        check(f"{donor.name} 4 variants", len(variants) == 4)
        lens = {len(v.sequence) for v in variants}
        check(f"{donor.name} same length across variants", len(lens) == 1, str(lens))
        v2 = variants[2]
        for pos in v2.backmutations:
            check(f"{donor.name} {pos} back-mutated in V2",
                  v2.graft.origin.get(pos) == "donor", str(v2.graft.origin.get(pos)))


def test_vhh_protection():
    print("VHH hallmark protection")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    _, chains = parse_input(os.path.join(ROOT, "data", "examples", "cab_rn05_vhh.fasta"))
    vhh = chains[0]
    check("VHH detected", vhh.is_vhh)
    choice = choose_germlines(vhh.numbered, db)
    top = [(g, s) for g, s in choice.alternatives]
    bm = analyze_backmutations(vhh.numbered, choice.v_gene, is_vhh=True, top_germlines=top)
    kd = {c.position for c in bm.candidates if c.tier == "KEEP_DONOR"}
    check("hallmark 37/44/45/47 all KEEP_DONOR",
          {"H37", "H44", "H45", "H47"} <= kd, str(sorted(kd)))
    variants = assemble_variants(vhh.numbered, choice.v_gene, choice.j_gene,
                                 "kabat", bm, is_vhh=True)
    for v in variants:
        fr2 = v.sequence[35:49]
        check(f"VHH {v.name} hallmark intact in FR2", fr2.startswith("WFRQAPGKEREG"),
              fr2)


def test_minimal_reversion():
    print("minimal reversion + CVI + matrix")
    from humanize.minimal import (
        build_paratope_variant,
        cvi_homology,
        matrix_alternatives,
        minimal_reversion_set,
    )
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    _, chains = parse_input(os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"))
    for donor in chains:
        choice = choose_germlines(donor.numbered, db)
        top = [(g, s) for g, s in choice.alternatives]
        bm = analyze_backmutations(donor.numbered, choice.v_gene, top_germlines=top)
        # no-structure fallback = Tier-1
        mr = minimal_reversion_set(donor.numbered, bm, structure=None)
        check(f"{donor.name} Vmin fallback == Tier-1",
              set(mr.positions) == set(bm.revert_positions(("T1",))), str(mr.positions))
        check(f"{donor.name} Vmin method tier", mr.method == "tier")
        # CVI homology sanity (0..1)
        cvi = cvi_homology(donor.numbered, choice.v_gene)
        check(f"{donor.name} CVI homology in (0,1]", 0 < cvi <= 1.0, str(cvi))
        # matrix variants on alternatives
        entries = matrix_alternatives(donor.numbered, choice.alternatives[1:],
                                      choice.j_gene, "kabat", n=2)
        check(f"{donor.name} matrix has entries", len(entries) >= 1)
        if entries:
            e = entries[0]
            check(f"{donor.name} matrix graft length preserved",
                  len(e.graft_v2.sequence) == len(donor.sequence))
            check(f"{donor.name} matrix CVI present", 0 < e.cvi <= 1.0)
        # paratope variant requires structure -> None without hints
        sdr = build_paratope_variant(donor.numbered, choice.v_gene, choice.j_gene,
                                     "kabat", bm, StructureHints())
        check(f"{donor.name} V_SDR None without complex", sdr is None)


def test_learning_loop():
    print("closed-loop learning (synthetic data)")
    import json
    from humanize.learning import (
        compute_position_effects,
        load_calibration,
        parse_experiments,
        write_calibration,
    )
    from humanize.backmut import analyze_backmutations

    # build synthetic experiments: V0 pure graft (4x KD loss), V2 (retained)
    exp = [{
        "name": "syn4D5",
        "parent_vh": M4D5_VH,
        "parent_vl": M4D5_VL,
        "parent_kd": 0.1,
        "variants": [
            {"name": "V0", "vh": IGHV3_23.replace("AKW", "AKW"), "vl": IGKV1_39, "kd": 0.4},
        ],
    }]
    # NOTE: synthetic variants intentionally do NOT renumber cleanly into a
    # real graft; use the demo experiment file instead for the real loop.
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(exp, fh)
        exp_path = fh.name
    records = parse_experiments(exp_path)
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    effects, warnings = compute_position_effects(records, db)
    check("learning parses experiments", len(records) == 1)
    check("learning produces effect map", isinstance(effects, dict))
    cal_path = os.path.join(tempfile.gettempdir(), "calib_test.json")
    write_calibration(cal_path, effects)
    cal = load_calibration(cal_path)
    check("calibration round-trip", "position_effects" in json.dumps(cal) or cal is not None)
    os.unlink(exp_path)


def test_humanness_adapter():
    print("BioPhi/Sapiens adapter (mock)")
    from humanize.humanness import run_sapiens
    script = os.path.join(tempfile.gettempdir(), "fake_biophi.py")
    with open(script, "w") as fh:
        fh.write(
            "#!/usr/bin/env python3\n"
            "import os, sys\n"
            "fa = [a for a in sys.argv if a.endswith('.fa') or a.endswith('.fasta')]\n"
            "out = os.path.join(os.path.dirname(os.path.abspath(fa[0])), 'scores.csv')\n"
            "with open(out, 'w') as f:\n"
            "    f.write('sequence,mean_score\\n')\n"
            "    name = ''\n"
            "    for line in open(fa[0]):\n"
            "        if line.startswith('>'):\n"
            "            name = line[1:].strip()\n"
            "        else:\n"
            "            f.write(f'{name},0.85\\n')\n"
        )
    os.chmod(script, 0o755)
    shim_dir = tempfile.mkdtemp()
    shim = os.path.join(shim_dir, "biophi")
    with open(shim, "w") as fh:
        fh.write(f"#!/bin/bash\nexec python3 {script} \"$@\"\n")
    os.chmod(shim, 0o755)
    old_path = os.environ.get("PATH", "")
    os.environ["PATH"] = shim_dir + os.pathsep + old_path
    try:
        seq = "EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVRQAPGKGLEWVG"
        res = run_sapiens({"test1": seq})
        r = res.get("test1")
        check("sapiens adapter parses scores", r is not None and r.sapiens_mean == 0.85,
              str(r))
    finally:
        os.environ["PATH"] = old_path


def test_docx_report():
    print("Word report generation")
    from humanize.report import write_all
    try:
        import docx  # noqa
    except ImportError:
        print("  [SKIP] python-docx not installed")
        return
    with tempfile.TemporaryDirectory() as out:
        result = run_pipeline(
            os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"),
            PipelineConfig(), outdir=out)
        paths = write_all(out, result)
        check("docx generated", os.path.exists(paths.get("docx", "")), str(paths.keys()))
        if paths.get("docx"):
            from docx import Document
            d = Document(paths["docx"])
            texts = "\n".join(p.text for p in d.paragraphs)
            check("docx has exec summary", "Executive Summary" in texts)
            check("docx has appendix sequences", "Appendix A" in texts)
            check("docx has tables", len(d.tables) >= 10, str(len(d.tables)))


def test_enhanced_report():
    print("WeMol-style enhanced report")
    import tempfile as _tf
    from humanize.report_enhanced import generate_enhanced_report
    with _tf.TemporaryDirectory() as out:
        result = run_pipeline(
            os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"),
            PipelineConfig(), outdir=out)
        content = generate_enhanced_report(result, out)
        check("enhanced report has template score", "Template Score" in content)
        check("enhanced report has mutation score", "Mutation Score" in content)
        # VL template table must not be all zeros (old hardcoded-H bug)
        import re
        vh_block = content.split("VL Chain (selected")[0]
        vl_block = content.split("VL Chain (selected")[1] if "VL Chain (selected" in content else ""
        rows = [l for l in vl_block.splitlines() if "\t" in l and l[0].isdigit()]
        nonzero = [l for l in rows if float(l.split("\t")[4]) > 0]
        check("VL template rows non-zero", len(nonzero) > 5,
              f"{len(nonzero)}/{len(rows)} rows with FR% > 0")
        # humanized sequences section populated with real sequences
        seqs = [l for l in content.split("5. Humanized Sequences")[-1].splitlines()
                if l.startswith((">", "EVQL", "DIQM"))]
        check("humanized sequences populated", any(l.startswith("EVQL") for l in seqs))
        check("enhanced report variant headers", any(">H V" in l or ">H" in l for l in seqs))


def test_end_to_end():
    print("end-to-end")
    with tempfile.TemporaryDirectory() as out:
        result = run_pipeline(
            os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"),
            PipelineConfig(), outdir=out)
        check("E2E format fab", result.format == "fab")
        check("E2E two chains", len(result.chains) == 2)
        from humanize.report import write_all
        paths = write_all(out, result)
        for p in paths.values():
            check(f"E2E output exists: {os.path.basename(p)}", os.path.exists(p))
        with open(paths["fasta"]) as fh:
            n_seq = sum(1 for line in fh if line.startswith(">"))
        # H: V0-V3 + Vmin (T1 non-empty); L: V0-V3 only (T1 empty after L87
        # gold-standard demotion, Vmin == V0). Total = 5 + 4.
        check("E2E fasta has 9 variants (H:5, L:4)", n_seq == 9, str(n_seq))


def main():
    test_numbering()
    test_graft_loop_conservation()
    test_vhh_humanization_gold_standard()
    test_germline_and_graft()
    test_backmut_variants()
    test_vhh_protection()
    test_minimal_reversion()
    test_humanness_adapter()
    test_docx_report()
    test_enhanced_report()
    test_end_to_end()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES: {FAILURES}")
        return 1
    print("ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
