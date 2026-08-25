#!/usr/bin/env python3
"""Retrospective backtest: does the pipeline recover the actual (clinically
approved) humanization decisions?

Cases (parent -> actual humanized, sequences from PDB structures):
  Case 1   mouse anti-HER2 4D5       -> trastuzumab (1N8Z / 1FVC)
  Case 2   mouse anti-VEGF A4.6.1    -> humanized anti-VEGF Fab (1BJ1)

For every framework position we compare what the ACTUAL humanized sequence
carries (donor residue / human germline residue / neither) against the
pipeline recommendation (T1/T2 revert, T3 optional, keep-human) and compute
precision / recall / minimality.

Two modes:
  A - full pipeline (germline chosen by the pipeline itself)
  B - forced germline (the framework family used historically)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "scripts"))

from humanize.backmut import analyze_backmutations
from humanize.germline import GermlineDB, load_germline_db
from humanize.graft import graft_chain
from humanize.numbering import number_heavy, number_light
from humanize.pipeline import human_likeness_percent
from humanize.sequences import InputChain, classify_sequences, detect_format

# ---------------------------------------------------------------------------
# Gold-standard sequences (verified from RCSB structures)
# ---------------------------------------------------------------------------

CASE1 = {
    "name": "4D5 -> trastuzumab (anti-HER2)",
    "parent": {
        "VH": "EVQLQQSGPELVKPGASVKMSCKASGYTFTDTYIHWVKQSHGKSLEWIGYINPYNGVTKYNQKFKGKATLTSDKSSSTAYMELSSLTSEDSAVYYCSRWGGDGFYAMDYWGQGTSVTVSS",
        "VL": "DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK",
    },
    "actual": {   # 1N8Z chains B/A (V region + CH1/CL tail, ignored by engine)
        "VH": "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEP",
        "VL": "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
    },
    "forced_germlines": {"VH": "IGHV3-66*01", "VL": "IGKV1-39*01"},
    "cdr_expected": {"VH": True, "VL": True},
}

CASE_VHH = {
    "name": "cAb-Lys3 -> hCAb-Lys3 (VHH humanization, Vincke 2009)",
    "is_vhh": True,
    "parent": {   # PDB 1MEL (camelid anti-lysozyme VHH)
        "VH": "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS",
    },
    # Reconstructed from the documented Vincke 2009 universal scaffold design:
    #   IGHV3-23-consensus framework + camelid hallmark (37F/44E/45R/47G) + S49
    #   + cAb-Lys3 CDRs grafted unchanged (incl. the disulfide CDR3).
    # Sequence should be verified against the paper before production use.
    "actual": {
        "VH": "EVQLVESGGGLVQPGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVSAINMGGGITYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSS",
    },
    "forced_germlines": {"VH": "IGHV3-23*01"},
    "cdr_expected": {"VH": True},
}


CASE2 = {
    "name": "A4.6.1 -> humanized anti-VEGF (bevacizumab lineage)",
    "parent": {   # 1CZ8 chains D/C (mouse A4.6.1)
        "VH": "EVQLVESGGGLVQPGGSLRLSCAASGYDFTHYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPYYYGTSHWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHL",
        "VL": "DIQLTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKVLIYFTSSLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
    },
    "actual": {   # 1BJ1 chains B/A (humanized anti-VEGF Fab)
        "VH": "EVQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHT",
        "VL": "DIQMTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKVLIYFTSSLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
    },
    "forced_germlines": {"VH": "IGHV3-23*01", "VL": "IGKV1-16*01"},
    "cdr_expected": {"VH": False, "VL": False},  # CDRs were affinity-matured
}


def get_gene(db: GermlineDB, gene_id: str):
    for g in db.v_genes + db.j_genes:
        if g.gene_id == gene_id:
            return g
    return None


def number_pair(vh_seq, vl_seq):
    vh = number_heavy(vh_seq)
    vl = number_light(vl_seq)
    return vh, vl


def analyze_case(case, db, force_germline=False):
    from humanize.backmut import StructureHints
    from humanize.minimal import cvi_homology
    from humanize.numbering import is_vhh

    is_vhh_case = case.get("is_vhh", False)
    chains = ("VH",) if is_vhh_case else ("VH", "VL")
    results = {}
    for ctype in chains:
        parent = number_heavy(case["parent"][ctype]) if ctype == "VH" \
            else number_light(case["parent"][ctype])
        actual = number_heavy(case["actual"][ctype]) if ctype == "VH" \
            else number_light(case["actual"][ctype])

        # germline choice
        if force_germline:
            from humanize.germline import GermlineChoice, compare_to_germline, score_j_match
            v_gene = get_gene(db, case["forced_germlines"][ctype])
            j_genes = db.j_for("H" if ctype == "VH" else "L")
            j_best = max(j_genes, key=lambda jg: score_j_match(actual, jg)[0])
            choice = GermlineChoice(v_gene=v_gene, j_gene=j_best)
            scores = compare_to_germline(parent, v_gene)
        else:
            from humanize.germline import choose_germlines
            choice = choose_germlines(parent, db)
            scores = choice.scores
        v_gene = choice.v_gene
        if v_gene is None or v_gene.numbered is None:
            results[ctype] = {"error": "no germline"}
            continue

        # pipeline recommendation on the parent (T1+T2 = V2 class)
        bm = analyze_backmutations(parent, v_gene, is_vhh=is_vhh_case)
        t12 = set(bm.revert_positions(("T1", "T2")))

        # what does the ACTUAL humanized sequence carry at each FR position?
        dmap, amap = parent.posmap(), actual.posmap()
        gmap = v_gene.numbered.posmap()
        tp = fp = fn = tn = 0
        actual_bm, compromise = set(), set()
        details = []
        for pos in sorted(set(dmap) & set(gmap), key=lambda p: (p[0], int("".join(c for c in p if c.isdigit())))):
            num = int("".join(c for c in pos if c.isdigit()))
            if num >= (103 if ctype == "VH" else 98):
                continue
            if parent.region_of(pos) not in ("FR1", "FR2", "FR3"):
                continue
            # H93/H94 carry the first two CDR3-loop residues; they are
            # grafted from the donor (P0 fix) and never back-mutation targets.
            if ctype == "VH" and num in (93, 94):
                continue
            if dmap[pos] == gmap[pos] or pos not in amap:
                continue
            actual_aa = amap[pos]
            if actual_aa == dmap[pos]:
                is_bm = True
                kind = "donor"                 # actual kept the donor residue
            elif actual_aa == gmap[pos]:
                is_bm = False
                kind = "germline"              # actual kept the human residue
            else:
                is_bm = False
                kind = "compromise"            # actual = neither (engineered)
            rec = pos in t12                   # pipeline recommends revert
            # VHH hallmark positions are kept via KEEP_DONOR graft protection
            # (not T1/T2 reversion); retaining the donor there is satisfied.
            kd_pos = False
            if is_vhh_case and ctype == "VH":
                c_obj = next((c for c in bm.candidates if c.position == pos), None)
                kd_pos = c_obj is not None and c_obj.tier == "KEEP_DONOR"
            if is_bm and kd_pos:
                rec = True
            if is_bm:
                actual_bm.add(pos)
            elif kind == "compromise":
                compromise.add(pos)
            if is_bm and rec:
                tp += 1
            elif is_bm and not rec:
                fn += 1
            elif not is_bm and rec:
                fp += 1
            else:
                tn += 1
            details.append((pos, dmap[pos], gmap[pos], actual_aa, kind, rec))

        # CDR identity between graft and actual
        graft = graft_chain(parent, v_gene, choice.j_gene, "kabat")
        gseq = graft.numbered.posmap()
        cdr_diff = []
        for pos, aa in gseq.items():
            if pos in amap and amap[pos] != aa and \
                    actual.region_of(pos).startswith("CDR"):
                cdr_diff.append((pos, aa, amap[pos]))

        results[ctype] = {
            "germline": v_gene.gene_id,
            "detail_bm_candidates": bm.candidates,
            "fr_identity_parent": round(scores.get("fr_identity", 0), 3),
            "cvi": cvi_homology(parent, v_gene),
            "actual_backmutations": sorted(actual_bm, key=lambda p: (p[0], int("".join(c for c in p if c.isdigit())))),
            "compromise_positions": sorted(compromise, key=lambda p: (p[0], int("".join(c for c in p if c.isdigit())))),
            "our_T1T2": sorted(t12, key=lambda p: (p[0], int("".join(c for c in p if c.isdigit())))),
            "tp": tp, "fp": fp, "fn": fn, "tn": tn,
            "precision": round(tp / (tp + fp), 3) if tp + fp else None,
            "recall": round(tp / (tp + fn), 3) if tp + fn else None,
            "minimality": f"{len(t12)} recommended vs {len(actual_bm)} actual",
            "cdr_differences": cdr_diff,
            "hl_vs_actual": round(human_likeness_percent(actual, v_gene), 1),
            "details": details,
        }
    return results


def print_case(case, db, force=False):
    label = "MODE B (forced historical germline)" if force else "MODE A (pipeline germline choice)"
    print(f"\n{'='*78}\n{case['name']}  [{label}]\n{'='*78}")
    res = analyze_case(case, db, force_germline=force)
    if case.get("is_vhh"):
        from humanize.numbering import is_vhh
        parent = number_heavy(case["parent"]["VH"])
        vhh, score, matched = is_vhh(parent)
        print(f"\n  VHH hallmark detected: {vhh} ({score}/4: {matched})")
        kd = [c.position for c in res["VH"]["detail_bm_candidates"]
              if c.tier == "KEEP_DONOR"] if "detail_bm_candidates" in res["VH"] else None
        if kd:
            print(f"  hallmark positions kept KEEP_DONOR: {sorted(kd)}")
    for ctype, r in res.items():
        if "error" in r:
            print(f"  {ctype}: {r['error']}")
            continue
        print(f"\n  {ctype}:")
        print(f"    germline: {r['germline']}  (parent FR id {r['fr_identity_parent']}, CVI {r['cvi']})")
        print(f"    actual back-mutations: {', '.join(r['actual_backmutations']) or '(none)'}")
        print(f"    actual compromise positions: "
              f"{', '.join(r['compromise_positions']) or '(none)'}")
        print(f"    our T1+T2 reversion : {', '.join(r['our_T1T2']) or '(none)'}")
        print(f"    precision {r['precision']} | recall {r['recall']} | {r['minimality']}")
        print(f"    human-likeness of ACTUAL vs chosen germline: {r['hl_vs_actual']}%")
        if r["cdr_differences"]:
            print(f"    CDR diffs graft-vs-actual: "
                  f"{', '.join(f'{p} {a}->{b}' for p, a, b in r['cdr_differences'][:8])}")
        misses = [d for d in r["details"] if d[4] == "donor" and not d[5]]
        over = [d for d in r["details"] if d[4] == "germline" and d[5]]
        comp = [d for d in r["details"] if d[4] == "compromise" and d[5]]
        if misses:
            print(f"    missed (actual reverted, we didn't): "
                  f"{', '.join(f'{d[0]}({d[1]}->{d[3]})' for d in misses)}")
        if over:
            print(f"    over-reverted (we did, actual=germline): "
                  f"{', '.join(f'{d[0]}({d[1]}->{d[3]})' for d in over)}")
        if comp:
            print(f"    compromise zones (actual=neither, we recommend donor): "
                  f"{', '.join(f'{d[0]}({d[1]}->{d[3]})' for d in comp)}")


def main():
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    print_case(CASE1, db, force=False)
    print_case(CASE1, db, force=True)
    print_case(CASE2, db, force=False)
    print_case(CASE2, db, force=True)
    print_case(CASE_VHH, db, force=False)
    print_case(CASE_VHH, db, force=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
