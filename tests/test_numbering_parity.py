#!/usr/bin/env python3
"""Numbering parity/regression tests.

Run with `python3 tests/test_numbering_parity.py` (plain script, no pytest).

Covers the anchor-engine fixes verified against ANARCI / the `number`
reference implementation:
  * VL CDR3 >= 14 no longer swallows the J region (FR4 preserved)
  * lambda FR1 deletion numbered at L10, conserved Cys stays on L23
  * VL FR1 capped when the Cys anchor shifts (no duplicated labels)
  * V-region-only VL warns instead of silently folding the tail into CDR3
  * built-in VH CDR3-H short-loop gap placement matches ANARCI
  * Chothia CDR-H1 = Kabat 26-32 (AbM stays 26-35)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from humanize import anarci_adapter
from humanize.germline import load_germline_db
from humanize.graft import CDR_POS_SETS, cdr_positions_for
from humanize.numbering import CDR_SCHEMES, number_heavy, number_light

FAILURES = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    print(f"  [{status}] {name}" + (f" -- {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(name)


VL = ("DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG"
      "TDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK")
LAMBDA = ("QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKS"
          "GTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL")
VH_PREFIX = ("QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQGLEWMGWINAGNGNTKYSQKFQGRV"
             "TITRDTSASTAYMELSSLRSEDTAVYYC")
VH_FR4 = "WGQGTLVTVSS"


def test_vl_cdr3_ladder():
    print("VL CDR3 insertion ladder keeps L96/L97 + FR4")
    for extra in range(0, 9):
        cdr3 = "QQHYTTP" + "A" * extra + "PT"
        seq = VL.replace("QQHYTTPPT", cdr3)
        ch = number_light(seq)
        cdr3_r = [r.pos for r in ch.residues if r.region == "CDR3"]
        fr4_r = [r.pos for r in ch.residues if r.region == "FR4"]
        check(f"ladder +{extra}: CDR3 has {len(cdr3)} residues",
              len(cdr3_r) == len(cdr3), str(len(cdr3_r)))
        check(f"ladder +{extra}: CDR3 ends at L97", cdr3_r[-1] == "L97", cdr3_r[-1])
        check(f"ladder +{extra}: FR4 is L98-L107",
              fr4_r == [f"L{i}" for i in range(98, 108)], str(fr4_r))

    # The regression: at 14 residues the old 16-residue J window lost "FGQG"
    # and turned the whole tail (incl. human FR4) into CDR3.
    cdr3 = "QQHYTTP" + "A" * 5 + "PT"
    ch = number_light(VL.replace("QQHYTTPPT", cdr3))
    cdr3_r = [r.pos for r in ch.residues if r.region == "CDR3"]
    check("14-residue CDR3 keeps 14 residues", len(cdr3_r) == 14, str(len(cdr3_r)))
    check("14-residue CDR3 keeps insertion codes",
          cdr3_r[7:12] == ["L95A", "L95B", "L95C", "L95D", "L95E"], str(cdr3_r))


def test_lambda_fr1():
    print("lambda FR1: Kabat deletion at L10, Cys on L23")
    ch = number_light(LAMBDA)
    pm = ch.posmap()
    check("no residue on L10", pm.get("L10") is None, str(pm.get("L10")))
    check("L11 is V", pm.get("L11") == "V", str(pm.get("L11")))
    check("L22 is S", pm.get("L22") == "S", str(pm.get("L22")))
    check("conserved Cys on L23", pm.get("L23") == "C", str(pm.get("L23")))
    check("CDR1 starts on L24", pm.get("L24") == "S", str(pm.get("L24")))
    labels = [r.pos for r in ch.residues]
    check("no duplicated position labels", len(labels) == len(set(labels)))
    check("L10 deletion is warned about",
          any("L10" in w for w in ch.warnings), str(ch.warnings))


def test_vl_fr1_cap():
    print("VL FR1 cap when Cys anchor shifts")
    ch = number_light("X" + VL)
    labels = [r.pos for r in ch.residues]
    check("no duplicated position labels", len(labels) == len(set(labels)))
    check("exactly one L24", labels.count("L24") == 1, str(labels.count("L24")))


def test_vl_v_region_only_warns():
    print("V-region-only VL does not silently lose FR4")
    ch = number_light(VL[:VL.index("FGQGTKVEIK")])
    check("FR4 is empty", not any(r.region == "FR4" for r in ch.residues))
    check("J-anchor warning emitted",
          any("Phe98" in w for w in ch.warnings), str(ch.warnings))


def _with_anarci_disabled(fn, *args, **kwargs):
    saved = anarci_adapter.is_anarci_available
    try:
        anarci_adapter.is_anarci_available = lambda: False
        return fn(*args, **kwargs)
    finally:
        anarci_adapter.is_anarci_available = saved


def test_vh_cdr3_gap_matches_anarci():
    print("built-in VH CDR3-H gap placement matches ANARCI")
    if not anarci_adapter.is_anarci_available():
        print("  [SKIP] ANARCI not installed")
        return
    from anarci import anarci as _anarci
    for loop_len in range(3, 11):
        seq = VH_PREFIX + "AR" + "G" * loop_len + VH_FR4
        ch = _with_anarci_disabled(number_heavy, seq)
        mine = {r.pos: r.aa for r in ch.residues
                if r.region == "CDR3" and r.pos.startswith("H")}
        res = _anarci([("q", seq)], scheme="kabat", output=False)
        nl = res[0][0][0][0]
        ref = {f"H{pn}{fl.strip()}": aa for (pn, fl), aa in nl
               if aa != "-" and 95 <= pn <= 102}
        check(f"loop {loop_len + 2}: labels/AA identical", mine == ref,
              f"{sorted(mine)} vs {sorted(ref)}")


def test_chothia_h1_window():
    print("Chothia CDR-H1 = Kabat 26-32")
    check("CDR_SCHEMES chothia H1 is H26-H32",
          CDR_SCHEMES["chothia"]["H"]["CDR1"] == ("H26", "H32"),
          str(CDR_SCHEMES["chothia"]["H"]["CDR1"]))
    check("graft CDR_POS_SETS chothia H1 is 26-32",
          CDR_POS_SETS["chothia"]["H"]["CDR1"] == (26, 32),
          str(CDR_POS_SETS["chothia"]["H"]["CDR1"]))
    pos = cdr_positions_for("chothia", "H")
    check("chothia H1 covers 26-32",
          set(range(26, 33)) <= pos, str(sorted(pos)))
    check("chothia H1 excludes 33-35",
          not (set(range(33, 36)) & pos), str(sorted(pos)))
    check("AbM H1 remains 26-35 (unchanged)",
          CDR_POS_SETS["abm"]["H"]["CDR1"] == (26, 35))


def test_germline_lambda_fr1_consistent():
    """The bundled V genes must use the same lambda FR1 convention as donors.

    Guards the germline numbering cache version: stale caches built before the
    L10 fix store L1-L22 with Cys on L22, which would silently misalign every
    lambda FR1 comparison.
    """
    print("bundled lambda germline FR1 matches the donor convention")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    ilgv = [g for g in db.v_genes if g.gene_id.startswith("IGLV") and g.numbered]
    check("bundled DB contains numbered IGLV genes", len(ilgv) > 0, str(len(ilgv)))
    bad = []
    for g in ilgv:
        pm = g.numbered.posmap()
        if pm.get("L10") is not None or pm.get("L23") != "C":
            bad.append(g.gene_id)
    check("every IGLV gene has an L10 gap and Cys on L23", not bad, str(bad[:5]))


def main():
    test_vl_cdr3_ladder()
    test_lambda_fr1()
    test_vl_fr1_cap()
    test_vl_v_region_only_warns()
    test_vh_cdr3_gap_matches_anarci()
    test_chothia_h1_window()
    test_germline_lambda_fr1_consistent()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES: {FAILURES}")
        return 1
    print("ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
