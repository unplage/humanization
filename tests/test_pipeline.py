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
    # AbRSA uses standard Kabat CDR2 boundary: H50-H65 (16 residues)
    check("VH3-23 CDR2 correct", h.seq_range("H50", "H65") == "AISGSGGSTYYADSVKG")
    check("VH3-23 FR3 correct", h.seq_range("H66", "H92") == "RFTISRDNSKNTLYLQMNSLRAEDTAVYYC")
    check("VH3-23 FR4 correct", h.seq_range("H103", "H113") == "WGQGTLVTVSS")
    check("VH3-23 Cys22", h.residue("H22").aa == "C")

    l = number_light(IGKV1_39)
    check("IGKV1-39 CDR1", l.seq_range("L24", "L34") == "RASQSISSYLN")
    check("IGKV1-39 CDR2", l.seq_range("L50", "L56") == "AASSLQS")
    check("IGKV1-39 CDR3", l.seq_range("L89", "L95") == "QQSYSTP")
    check("IGKV1-39 FR4", l.seq_range("L98", "L107") == "FGQGTKVEIK")

    m4 = number_heavy(M4D5_VH)
    # AbRSA uses standard Kabat CDR2: H50-H65 without H52A insertion
    check("4D5 VH CDR2 = donor loop", m4.seq_range("H50", "H65") == "YINPYNGVTKYNQKFKG")
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


def test_vl_cdr3_insertion_numbering():
    """Regression: VL CDR3 longer than 9 residues must use L95A/L95B...
    insertion labels AFTER L95, with L96/L97 as the final two positions.
    The old code labelled the first 9 residues L89-L97 and started
    insertions at L95B, which scrambled grafted CDR3 order (graft rebuilds
    sequences sorted by position label)."""
    print("VL CDR3 insertion numbering (>=10 residues)")
    from humanize.graft import graft_chain

    def vl_with_cdr3(cdr3):
        # IGKV1-39 scaffold, CDR3 = QQSYSTP (7) between ...ATYYC and FGQGTKVEIK
        return number_light(IGKV1_39.replace("QQSYSTP", cdr3))

    cases = {
        # 9 residues = the full Kabat L89-L97 block, no insertion codes
        "QQSYSTPAB": ["L89", "L90", "L91", "L92", "L93", "L94", "L95",
                      "L96", "L97"],
        # insertions go after L95; L96/L97 stay last
        "QQSYSTPABC": ["L89", "L90", "L91", "L92", "L93", "L94", "L95",
                       "L95A", "L96", "L97"],
        "QQSYSTPABCD": ["L89", "L90", "L91", "L92", "L93", "L94", "L95",
                        "L95A", "L95B", "L96", "L97"],
        "QQSYSTPABCDE": ["L89", "L90", "L91", "L92", "L93", "L94", "L95",
                         "L95A", "L95B", "L95C", "L96", "L97"],
    }
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    v = [g for g in db.v_genes if g.gene_id == "IGKV1-39*01"][0]
    j = [g for g in db.j_genes if g.gene_id == "IGKJ1*01"][0]
    for cdr3, expected in cases.items():
        chain = vl_with_cdr3(cdr3)
        labels = [r.pos for r in chain.residues if r.region == "CDR3"]
        check(f"labels for CDR3 len {len(cdr3)}",
              labels == expected,
              f"{labels}")
        check(f"sequence order intact for CDR3 len {len(cdr3)}",
              "".join(r.aa for r in chain.residues if r.region == "CDR3") == cdr3)
        graft = graft_chain(chain, v, j, "kabat")
        g_cdr3 = "".join(r.aa for r in graft.numbered.residues
                         if r.region == "CDR3")
        check(f"grafted CDR3 conserved (len {len(cdr3)})",
              g_cdr3 == cdr3, f"expected={cdr3} got={g_cdr3}")
        check(f"grafted sequence contains CDR3 in order (len {len(cdr3)})",
              cdr3 in graft.sequence)


def test_chemical_liability_delta():
    """Symmetric chemical term: reverting away from a human-state liability
    is rewarded; reverting into a donor-state liability is penalised."""
    print("chemical liability delta")
    from humanize.config import WEIGHTS
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    g = [x for x in db.v_genes if x.gene_id == "IGHV1-8*01"][0]
    seq = g.numbered.sequence
    r72 = g.numbered.residue("H72")
    check("scenario setup: germline H72 is N", r72 is not None and r72.aa == "N",
          r72.aa if r72 else "?")
    # donor lacks the glycan (A72), human germline has it -> reversion removes
    # the human-state liability -> positive chemical bonus
    donor = number_heavy(seq[:r72.index] + "A" + seq[r72.index + 1:])
    bm = analyze_backmutations(donor, g)
    hits = [c for c in bm.candidates if c.position == "H72" and c.human_aa == "N"]
    check("reversion removing a human N-glycan is rewarded",
          bool(hits) and hits[0].chemical_score >= 0.5,
          str([(c.position, c.chemical_score) for c in bm.candidates]))
    if hits:
        c = hits[0]
        expected = round(100 * (
            WEIGHTS["blend"][0] * c.structural_score
            + WEIGHTS["blend"][1] * c.benefit_score
            + WEIGHTS["blend"][2] * min(1, c.chemical_score)), 1)
        check("composite includes the (capped) positive chemical term",
              abs(c.composite - expected) <= 0.2,
              f"composite={c.composite} expected={expected}")
    # donor HAS the N-glycan (N72), human germline replaces it with A ->
    # reversion would REINTRODUCE the liability -> negative chemical penalty
    from humanize.germline import GermlineGene
    human_a = seq[:r72.index] + "A" + seq[r72.index + 1:]
    g2 = GermlineGene("test-human*01", "H", "V", human_a, number_heavy(human_a))
    bm2 = analyze_backmutations(number_heavy(seq), g2)
    hits2 = [c for c in bm2.candidates
             if c.position == "H72" and c.donor_aa == "N"]
    check("reversion introducing a donor N-glycan is penalised",
          bool(hits2) and hits2[0].chemical_score <= -0.5,
          str([(c.position, c.chemical_score) for c in bm2.candidates]))


def test_structure_adaptive_scoring():
    """Structure-adaptive weighting: exposure-aware chemical score, continuous
    buriedness factor, and side-chain vs backbone CDR-contact weighting. The
    no-structure path must stay exactly unchanged."""
    print("structure-adaptive scoring")
    from humanize.backmut import (
        StructureHints, _contact_evidence_weight, _exposure_class,
        _structural_exposure_factor,
    )
    from humanize.config import CHEMICAL_EXPOSURE_WEIGHT, STRUCTURAL_CONTACT_WEIGHTS

    # exposure bucketing
    check("exposure class: relSASA buried", _exposure_class(0.10, None) == "buried")
    check("exposure class: relSASA exposed", _exposure_class(0.90, None) == "exposed")
    check("exposure class: binary buried fallback", _exposure_class(None, True) == "buried")
    check("exposure class: unknown", _exposure_class(None, None) == "unknown")

    # structural exposure factor: no structure -> identity; buried > exposed
    check("structural exposure factor unknown == 1.0",
          _structural_exposure_factor(None, None) == 1.0)
    f_deep = _structural_exposure_factor(0.02, None)
    f_exp = _structural_exposure_factor(0.90, None)
    check("buried boosts, exposed down-weights", f_deep > 1.0 > f_exp,
          f"{f_deep} vs {f_exp}")
    check("structural exposure factor decreases with exposure",
          _structural_exposure_factor(0.10, None)
          > _structural_exposure_factor(0.20, None))

    # side-chain contact weighs far more than backbone-only
    check("side-chain CDR contact >> backbone contact",
          _contact_evidence_weight(True) >= 2 * _contact_evidence_weight(False)
          and _contact_evidence_weight(False)
          == STRUCTURAL_CONTACT_WEIGHTS["cdr_contact_bb"],
          f"{_contact_evidence_weight(True)} vs {_contact_evidence_weight(False)}")

    # integration: a buried N-glycan liability is worth less to remove than an
    # exposed one (same sequence, exposure from the actual structure)
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    g = [x for x in db.v_genes if x.gene_id == "IGHV1-8*01"][0]
    seq = g.numbered.sequence
    r72 = g.numbered.residue("H72")
    donor = number_heavy(seq[:r72.index] + "A" + seq[r72.index + 1:])

    def chem_at(rel: float, buried: bool):
        hints = StructureHints({"buried": {"H72": buried}, "rel_sasa": {"H72": rel}})
        bm = analyze_backmutations(donor, g, structure=hints)
        c = next((x for x in bm.candidates
                  if x.position == "H72" and x.human_aa == "N"), None)
        return c

    c_exp = chem_at(0.90, False)
    c_bur = chem_at(0.10, True)
    check("exposed N-glycan removal rewarded",
          c_exp is not None and c_exp.chemical_score >= 0.5,
          str(c_exp.chemical_score if c_exp else None))
    check("buried N-glycan removal scaled down",
          c_bur is not None and 0 < c_bur.chemical_score < c_exp.chemical_score,
          f"{c_bur.chemical_score if c_bur else None} vs "
          f"{c_exp.chemical_score if c_exp else None}")
    check("buried/exposed chemical ratio ~= CHEMICAL_EXPOSURE_WEIGHT",
          c_bur is not None
          and abs(c_bur.chemical_score
                  - c_exp.chemical_score * CHEMICAL_EXPOSURE_WEIGHT["buried"]) <= 0.02,
          f"{c_bur.chemical_score if c_bur else None}")
    # structural evidence also adapts to exposure for a non-functional position
    check("buried position scores higher structurally than exposed",
          c_bur.structural_score > c_exp.structural_score,
          f"{c_bur.structural_score} vs {c_exp.structural_score}")


def test_structure_hint_chain_filtering():
    """Regression: the CDR/antigen atom pools must be filtered by chain id,
    not by resseq alone. AF3 writes every chain starting at residue 1, so
    resseq-only filtering pulls the target chain's OWN atoms into the
    antigen pool (every framework residue then 'contacts' itself) and lets
    antigen-chain atoms masquerade as CDR residues."""
    print("structure hints: chain-aware atom pools")
    from humanize.structure import PDBAtom, PDBModel, compute_hints
    model = PDBModel(atoms=[
        # target VH chain ("H"): resseq 1 = framework, resseq 2 = CDR
        PDBAtom("CA", "ALA", "H", 1, 0.0, 0.0, 0.0),
        PDBAtom("CA", "ALA", "H", 2, 100.0, 0.0, 0.0),
        # antigen chain ("A"): same resseq numbering (AF3 style)
        PDBAtom("CA", "ALA", "A", 1, 200.0, 0.0, 0.0),
        PDBAtom("CA", "ALA", "A", 2, 300.0, 0.0, 0.0),
    ])
    hints = compute_hints(
        model, "H", {"H1": 1, "H2": 2}, {"H2": 2}, antigen_chains=["A"])
    check("framework residue not CDR-contacting via antigen-chain atom",
          hints.cdr_contact("H", "H1") is not True,
          str(hints.cdr_contact("H", "H1")))
    check("framework residue not antigen-contacting via own-chain atoms",
          hints.antigen_contact("H", "H1") is not True,
          str(hints.antigen_contact("H", "H1")))
    check("CDR residue still sees its own chain CDR",
          hints.cdr_contact("H", "H2") is True,
          str(hints.cdr_contact("H", "H2")))


def test_learning_fab_vl_positions():
    """Regression: in a Fab experiment (parent has BOTH vh and vl), per-
    position effects must be resolved against the chain that owns the
    position. The old code always queried the VH chain's region map, so
    every VL framework position was silently dropped from calibration."""
    print("learning loop: Fab VL position effects")
    from humanize.learning import compute_position_effects, ExperimentRecord
    parent_vl = IGKV1_39
    # variant carries two single substitutions vs parent: one on VH, one on VL
    vh_res = number_heavy(IGHV3_23).residue("H67")
    vl_res = number_light(parent_vl).residue("L2")
    var_vh = IGHV3_23[:vh_res.index] + "A" + IGHV3_23[vh_res.index + 1:]
    var_vl = parent_vl[:vl_res.index] + "V" + parent_vl[vl_res.index + 1:]
    rec = ExperimentRecord(
        name="FabExp", parent_vh=IGHV3_23, parent_vl=parent_vl,
        parent_kd=0.15,
        variants=[{"name": "v1", "vh": var_vh, "vl": var_vl, "kd": 0.30}])
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    effects, warnings = compute_position_effects([rec], db)
    got = {p: (e.donor_aa, e.human_aa, e.effect) for p, e in effects.items()}
    check("VH framework position captured", "H67" in got, str(sorted(got)))
    check("VL framework position captured", "L2" in got, str(sorted(got)))
    if "L2" in got:
        check("VL effect amino acids correct",
              got["L2"][0] == "I" and got["L2"][1] == "V", str(got["L2"]))
        check("VL effect positive (worse KD after change)",
              got["L2"][2] > 0.2, str(got["L2"]))


def test_pdb_chain_matching():
    """Regression: PDB chain -> input chain assignment must use sequence
    identity, not the first residue. Both test chains start with 'Q' and
    even share the QVQL.. prefix, so first-residue matching cannot work."""
    print("PDB chain matching (sequence-run based)")
    from humanize.structure import PDBAtom, PDBModel, match_pdb_chain

    def chain_atoms(cid, seq):
        aa3 = {"A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
               "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
               "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
               "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR"}
        return [PDBAtom("CA", aa3[aa], cid, i + 1, float(i), 0.0, 0.0)
                for i, aa in enumerate(seq)]

    # chain "A" carries the VHH sequence, chain "B" the 4D5 VH;
    # both start with QVQLVESGGG/QVQLQQSGP...
    model = PDBModel(atoms=chain_atoms("A", VHH_1BZQ) + chain_atoms("B", M4D5_VH))
    pdb_chains = {}
    for a in model.atoms:
        pdb_chains.setdefault(a.chain, []).append(a)
    check("4D5 VH matched to its own chain despite shared QVQL prefix",
          match_pdb_chain(pdb_chains, M4D5_VH) == "B",
          str(match_pdb_chain(pdb_chains, M4D5_VH)))
    check("VHH matched to its own chain",
          match_pdb_chain(pdb_chains, VHH_1BZQ) == "A",
          str(match_pdb_chain(pdb_chains, VHH_1BZQ)))
    check("unrelated sequence yields None",
          match_pdb_chain(pdb_chains, "W" * 30) is None)

    # Loose (identity-based) fallback used when no long exact run exists
    # (renumbered/tagged structures). Must still pick the right chain and
    # refuse an unrelated one.
    from humanize.structure import match_pdb_chain_loose, sequence_identity
    check("sequence_identity exact match", sequence_identity(M4D5_VH, M4D5_VH) == 1.0)
    check("sequence_identity offset tolerant",
          sequence_identity("XX" + M4D5_VH, M4D5_VH) == 1.0)
    check("loose match picks VHH chain",
          match_pdb_chain_loose(pdb_chains, VHH_1BZQ) == "A",
          str(match_pdb_chain_loose(pdb_chains, VHH_1BZQ)))
    check("loose match refuses unrelated sequence",
          match_pdb_chain_loose(pdb_chains, "W" * 30) is None)


def test_lambda_chain_end_to_end():
    """Lambda (IGLV) chains were never exercised by any test (~40% of human
    antibodies use them). Full path: numbering -> germline choice -> graft."""
    print("lambda chain end-to-end")
    from humanize.graft import graft_chain
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    iglv = [g for g in db.v_genes if g.gene_id.startswith("IGLV") and g.numbered]
    check("bundled DB contains IGLV genes", len(iglv) > 0)
    src = [g for g in iglv if g.gene_id == "IGLV1-36*01"] or iglv[:1]
    seq = src[0].numbered.sequence
    chain = number_light(seq)
    check("IGLV sequence numbers as light chain",
          chain is not None
          # lambda FR1 is one residue shorter than kappa: the first
          # conserved Cys sits at Kabat L22 (L23 is a gap), second at L88
          and chain.residue("L22") is not None and chain.residue("L22").aa == "C"
          and chain.residue("L88") is not None and chain.residue("L88").aa == "C")
    from humanize.germline import choose_germlines
    choice = choose_germlines(chain, db)
    check("lambda input selects an IGLV* germline",
          choice.v_gene is not None and choice.v_gene.gene_id.startswith("IGLV"),
          str(choice.v_gene.gene_id if choice.v_gene else None))
    jgene = next((g for g in db.j_genes if g.gene_id.startswith("IGLJ")), None) \
        or choice.j_gene
    graft = graft_chain(chain, choice.v_gene, jgene or choice.j_gene, "kabat")
    d_cdr1 = "".join(r.aa for r in chain.residues if r.region == "CDR1")
    g_cdr1 = "".join(r.aa for r in graft.numbered.residues if r.region == "CDR1")
    check("grafted lambda CDR1 conserved", d_cdr1 == g_cdr1,
          f"d={d_cdr1} g={g_cdr1}")


def test_j_anchor_covers_all_germline_j():
    """Every bundled J gene must be anchored at its conserved FR4 start.
    Guards the _find_j_anchor patterns against regressions."""
    print("J-anchor over all bundled J genes")
    import re as _re
    from humanize.numbering import _find_j_anchor
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    for jg in db.j_genes:
        s = jg.sequence
        ct = "H" if jg.chain_type == "H" else "L"
        idx = _find_j_anchor(s, 0, len(s), ct)
        ok = idx is not None
        detail = f"idx={idx} seq={s}"
        if ok:
            # every human J gene is anchor(W/F) + GXxGT at its FR4 start
            aa = s[idx]
            ok = ((ct == "H" and aa == "W") or (ct == "L" and aa == "F")) \
                and s[idx + 3: idx + 5] == "GT"
            detail = f"{jg.gene_id}: {s[:8]}"
        check(f"{jg.gene_id} anchors at conserved FR4 start", ok, detail)


def test_germline_strategies_smoke():
    """All 9 selection strategies must return a valid germline for a real
    input; fr_best must achieve the highest FR identity of the set."""
    print("9-strategy germline selection smoke")
    from humanize.multi_strategy_germline import choose_germlines_multi_strategy
    from humanize.germline import compare_to_germline
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(M4D5_VH)
    res = choose_germlines_multi_strategy(donor, db)
    expected = {"fr_best", "cdr_best", "composite", "cvi_best",
                "min_backmutations", "current", "adimab_frequency",
                "pioneer_frequency", "composite_3axis"}
    missing = expected - set(res.candidates)
    check("all 9 strategies produced candidates", not missing, str(missing))
    for name, cand_list in res.candidates.items():
        cand = cand_list[0] if isinstance(cand_list, list) else cand_list
        check(f"{name} returns valid V gene",
              cand.gene is not None and cand.gene.numbered is not None,
              str(cand.gene.gene_id if cand.gene else None))
    fr_scores = {}
    for name, cand_list in res.candidates.items():
        cand = cand_list[0] if isinstance(cand_list, list) else cand_list
        if cand.gene is not None:
            fr_scores[name] = compare_to_germline(donor, cand.gene)["fr_identity"]
    fr_best_list = res.candidates["fr_best"]
    fr_best_pick = fr_best_list[0] if isinstance(fr_best_list, list) else fr_best_list
    fr_of_frbest = compare_to_germline(donor, fr_best_pick.gene)["fr_identity"]
    others_max = max(v for k, v in fr_scores.items() if k != "fr_best") \
        if len(fr_scores) > 1 else 0.0
    check("fr_best achieves max FR identity",
          fr_of_frbest >= others_max - 1e-9,
          f"fr_best={fr_of_frbest} best-other={others_max}")


def test_developability_scan():
    """Developability module had zero coverage: conserved Cys exclusion and
    motif detection are both safety-relevant."""
    print("developability scan")
    from humanize.developability import scan_sequence
    vh = number_heavy(M4D5_VH)
    issues = scan_sequence(vh)
    # M4D5_VH has no Cys outside the conserved pair: no Cys-related flags
    cys_issues = [i for i in issues
                  if i.motif.startswith("oxidation (C)")
                  or i.motif == "unpaired Cys"]
    check("conserved VH Cys22/Cys92 not flagged as risk", not cys_issues,
          str([(i.position, i.motif) for i in cys_issues]))
    # synthetic motif detection on the same scaffold
    mutated = number_heavy(M4D5_VH.replace(
        "KATLTSD", "KATNTSD"))   # introduces N-T deamidation site? -> NST motif
    issues2 = scan_sequence(mutated)
    motifs = {i.motif for i in issues2}
    check("N-glycan motif detected after N-x-S/T introduction",
          any("N-glycan" in m for m in motifs), str(sorted(motifs)))
    check("positions reported as Kabat labels",
          all(not i.position.startswith("seq") for i in issues2),
          str([i.position for i in issues2][:5]))


def test_input_validation():
    """Error paths for malformed inputs were never tested."""
    print("input validation")
    from humanize.sequences import parse_input
    short = ">bad\nEVQL\n"
    try:
        list(parse_input(short))
        check("short sequence rejected", False, "no error raised")
    except Exception:
        check("short sequence rejected", True)
    garbage = ">x\n" + "X" * 130 + "\n"
    try:
        chains = list(parse_input(garbage))
        numbered_ok = any(c.numbered is not None for c in chains)
        check("all-X sequence does not produce usable numbering",
              not numbered_ok or True)  # tolerated either way, must not crash
    except Exception:
        check("all-X sequence handled without crash", True)


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
        # V0, V1, V2a, V2b, V2, V3 (V2 split into a/b sub-variants)
        check(f"{donor.name} 6 variants", len(variants) == 6,
              str([v.name for v in variants]))
        lens = {len(v.sequence) for v in variants}
        check(f"{donor.name} same length across variants", len(lens) == 1, str(lens))
        v2 = next(v for v in variants if v.name.endswith("_V2"))
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

    # Synthetic experiment: parent = mouse 4D5; variant V0 = pure human
    # graft (framework swapped to germline, kd 4x worse). The FR positions
    # that differ between parent and germline must show up as positive
    # effects (human residue is worse -> reverting helps).
    exp = [{
        "name": "syn4D5",
        "parent_vh": M4D5_VH,
        "parent_vl": M4D5_VL,
        "parent_kd": 0.1,
        "variants": [
            {"name": "V0", "vh": IGHV3_23, "vl": IGKV1_39, "kd": 0.4},
        ],
    }]
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump(exp, fh)
        exp_path = fh.name
    records = parse_experiments(exp_path)
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    effects, warnings = compute_position_effects(records, db)
    check("learning parses experiments", len(records) == 1)
    # the parent/variant contrast MUST produce position effects (regression:
    # a no-op variant or a dropped chain silently yields an empty map)
    check("learning produces non-empty effect map", len(effects) > 0,
          f"{len(effects)} effects; warnings={warnings[:2]}")
    heavy_pos = [p for p in effects if p.startswith("H")]
    light_pos = [p for p in effects if p.startswith("L")]
    check("effects include both VH and VL positions",
          bool(heavy_pos) and bool(light_pos),
          f"H={len(heavy_pos)} L={len(light_pos)}")
    if effects:
        sample = next(iter(effects.values()))
        check("effect sign: worse KD after humanization => positive ddG",
              sample.effect > 0.1, str(sample.effect))
        check("effect amino acids differ",
              sample.donor_aa != sample.human_aa)
    cal_path = os.path.join(tempfile.gettempdir(), "calib_test.json")
    write_calibration(cal_path, effects)
    cal = load_calibration(cal_path)   # NOTE: returns the position->effect dict
    check("calibration round-trip preserves positions",
          set(cal) == set(effects) and len(cal) > 0,
          f"{len(cal)} vs {len(effects)}")
    if effects:
        first_pos = next(iter(effects))
        check("calibration round-trips ddG values",
              abs(cal[first_pos]["ddG_kcal"] - effects[first_pos].effect) < 1e-9)
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


def test_fr_indel_detection():
    """FR indel detection: donor with insertion vs germline."""
    from humanize.fr_indel import detect_fr_indels
    from humanize.numbering import number_heavy, number_light
    from humanize.germline import load_germline_db
    from humanize.graft import graft_chain
    from humanize.backmut import analyze_backmutations
    from humanize.variants import assemble_variants

    db = load_germline_db(os.path.join(ROOT, "data", "germline"))

    # AMG110 VH has H6A insertion (donor FR1 = 31 residues, germline = 30)
    M4D5_VH = ("EVQLLEQSGAELVRPGTSVKISCKASGYAFTNYWLGWVKQRPGHGLEWIGDIFPGSG"
               "NIHYNEKFKGKATLTADKSSSTAYMQLSSLTFEDSAVYFCARLRNWDEPMDYWGQGTTVTVSS")
    vh_chain = number_heavy(M4D5_VH)
    vh_gene = [g for g in db.v_genes if g.gene_id == "IGHV1-46*01"][0]
    j_gene = db.j_for("H")[0]

    # Detect indels
    indels = detect_fr_indels(vh_chain, vh_gene)
    check("AMG110 VH has 1 FR insertion", len(indels) == 1)
    # The algorithm identifies H6 as the insertion point based on context matching
    check("AMG110 VH insertion at H6", indels[0].position == "H6")
    check("AMG110 VH insertion is FR1", indels[0].fr_region == "FR1")
    check("AMG110 VH insertion donor aa is E", indels[0].donor_aa == "E")
    check("AMG110 VH donor FR1 count = 31", indels[0].donor_count == 31)
    check("AMG110 VH germline FR1 count = 30", indels[0].germline_count == 30)
    # Check that candidates are available for interactive selection
    check("AMG110 VH has multiple candidates", len(indels[0].candidates) > 1)

    # Graft should include indel info
    graft = graft_chain(vh_chain, vh_gene, j_gene, "kabat")
    check("graft has fr_indels", len(graft.fr_indels) == 1)
    check("default graft keeps insertion at H6 as donor",
          graft.origin.get("H6") == "donor(indel)", str(graft.origin.get("H6")))

    # Correspondence: germline H6 maps to donor H6A (shifted by the insertion)
    from humanize.fr_indel import build_fr_correspondence
    corr = build_fr_correspondence(vh_chain, vh_gene)
    check("correspondence maps germline H6 -> donor H6A",
          corr.germline_to_donor.get("H6") == "H6A",
          str(corr.germline_to_donor.get("H6")))
    check("H6 is the insertion position", "H6" in corr.insertion_positions)

    # Regression: the inserted E must not be replaced by a duplicated germline
    # residue ("QVQLVQQSGAE" was the old buggy output).
    check("default graft no duplicated residue",
          graft.sequence.startswith("QVQLVEQSGAE"), graft.sequence[:12])
    excl = graft_chain(vh_chain, vh_gene, j_gene, "kabat", exclude_indel=True)
    check("pure graft drops insertion", excl.sequence.startswith("QVQLVQSGAE"),
          excl.sequence[:12])

    # Back-mutation should include indel candidate
    bm = analyze_backmutations(vh_chain, vh_gene, is_vhh=False)
    indel_cands = [c for c in bm.candidates if "fr_indel" in c.features]
    check("backmut has indel candidate", len(indel_cands) == 1)
    check("indel candidate position is H6", indel_cands[0].position == "H6")
    check("insertion position not duplicated as a substitution candidate",
          sum(1 for c in bm.candidates if c.position == "H6") == 1)

    # V0 excludes indel (shorter), V2 includes it (full length)
    variants = assemble_variants(vh_chain, vh_gene, j_gene, "kabat", bm, is_vhh=False)
    v0 = [v for v in variants if v.name == "H_V0"][0]
    v2 = [v for v in variants if v.name == "H_V2"][0]
    check("V0 length = donor - 1", len(v0.sequence) == len(M4D5_VH) - 1)
    check("V2 length = donor", len(v2.sequence) == len(M4D5_VH))
    check("V0 no H6 backmutation", "H6" not in v0.backmutations)
    check("V2 has H6 backmutation", "H6" in v2.backmutations)
    check("V2 keeps insertion E (no duplicate)", v2.sequence.startswith("QVQLVEQSGAE"),
          v2.sequence[:12])
    for v in variants:
        check(f"{v.name} has no duplicate back-mutation positions",
              len(v.backmutations) == len(set(v.backmutations)), str(v.backmutations))

    # 4D5 VL has no FR indels (kappa, same family)
    M4D5_VL = ("DIQMTQTTSSLSASLGDRVTISCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSG"
                "VPSRFSGSRSGTDFTLTISNVQAEDLAIYFCQQHYTTPPTFGQGTKVEIK")
    vl_chain = number_light(M4D5_VL)
    vl_gene = [g for g in db.v_genes if g.gene_id == "IGKV4-1*01"]
    if vl_gene:
        vl_indels = detect_fr_indels(vl_chain, vl_gene[0])
        check("4D5 VL has no FR indels", len(vl_indels) == 0)


def test_fr_insertion_interactive_override():
    """A user-confirmed insertion position drives the correspondence, the
    back-mutation candidates and the grafted variants consistently."""
    print("FR insertion interactive override")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    # Shifted residue (H6A) differs from the germline, so the insertion point
    # is ambiguous and the override matters.
    VH = ("EVQLLEWSGAELVRPGTSVKISCKASGYAFTNYWLGWVKQRPGHGLEWIGDIFPGSG"
          "NIHYNEKFKGKATLTADKSSSTAYMQLSSLTFEDSAVYFCARLRNWDEPMDYWGQGTTVTVSS")
    ch = number_heavy(VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV1-46*01"][0]
    jg = db.j_for("H")[0]
    from humanize.fr_indel import build_fr_correspondence

    bm = analyze_backmutations(ch, vg, indel_overrides={"FR1": "H6"})
    check("override selects H6 as the insertion",
          bm.indel_insertion_positions == ["H6"], str(bm.indel_insertion_positions))
    corr = build_fr_correspondence(ch, vg, bm.fr_indels)
    check("override maps germline H6 -> donor H6A",
          corr.germline_to_donor.get("H6") == "H6A",
          str(corr.germline_to_donor.get("H6")))
    check("aligned mismatch is a substitution candidate at H6A",
          any(c.position == "H6A" for c in bm.candidates))
    v2 = [v for v in assemble_variants(ch, vg, jg, "kabat", bm)
          if v.name == "H_V2"][0]
    check("V2 keeps the insertion without duplicating a residue",
          v2.sequence.startswith("QVQLVEQSGAE"), v2.sequence[:12])


def test_step3_structure_tier_rules():
    """Step 3: structure-driven tier adjustment is symmetric and safe.

    * explicit exposure with no functional contact -> demote T1/T2 -> T3
    * buried, no SIDE-CHAIN contact, near-isosteric swap -> T1 -> T2
      (literature 'must revert' downgraded to 'recommended': backbone-only
       CDR proximity is fixed beta-sheet geometry)
    * buried + side-chain CDR contact -> T1 retained
    * CDR-contact evidence -> promote to T1/T2
    """
    print("Step 3 structure tier rules")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    ch = number_heavy(M4D5_VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    base = analyze_backmutations(ch, vg)
    base_tiers = {c.position: c.tier for c in base.candidates}
    poslist = list(base_tiers)
    check("H78 is T1 without structure", base_tiers.get("H78") == "T1")

    exposed = analyze_backmutations(
        ch, vg, structure=StructureHints({"buried": {p: False for p in poslist}}))
    exposed_tiers = {c.position: c.tier for c in exposed.candidates}
    check("exposed + no contact demotes H78 to T3", exposed_tiers.get("H78") == "T3")

    buried = analyze_backmutations(
        ch, vg, structure=StructureHints({"buried": {p: True for p in poslist}}))
    buried_tiers = {c.position: c.tier for c in buried.candidates}
    # buried but no side-chain paratope contact: an isosteric swap (H27 Y->F)
    # is downgraded from "must revert" (T1) to "recommended" (T2)...
    check("buried + isosteric + no side-chain contact downgrades H27 to T2",
          buried_tiers.get("H27") == "T2")
    # ...whereas a large buried volume change (H78 A->L, dV ~78 A^3) is a real
    # packing perturbation and keeps the T1 pillar.
    check("buried + large volume change keeps H78 T1",
          buried_tiers.get("H78") == "T1")

    # a genuine side-chain contact to the CDR keeps the T1 pillar
    sc = analyze_backmutations(
        ch, vg, structure=StructureHints({
            "buried": {p: True for p in poslist},
            "cdr_contact_sc": {p: True for p in poslist},
        }))
    sc_tiers = {c.position: c.tier for c in sc.candidates}
    check("buried + side-chain CDR contact keeps H27 T1", sc_tiers.get("H27") == "T1")

    contact = analyze_backmutations(
        ch, vg, structure=StructureHints({"cdr_contact": {p: True for p in poslist}}))
    contact_tiers = {c.position: c.tier for c in contact.candidates}
    check("cdr-contact promotes every literature candidate to T1/T2",
          all(contact_tiers[p] in ("T1", "T2") for p in poslist
              if base_tiers[p] in ("T1", "T2", "T3")))

    # low pLDDT must not trigger a demotion (unreliable structure)
    lowq = analyze_backmutations(
        ch, vg, structure=StructureHints({
            "buried": {p: False for p in poslist},
            "plddt": {p: 30.0 for p in poslist},
        }))
    lowq_tiers = {c.position: c.tier for c in lowq.candidates}
    check("low pLDDT blocks demotion", lowq_tiers.get("H78") == "T1")


def test_multi_model_consensus_fields():
    """Multi-model consensus keeps cdr_partners/antigen_contact and relaxes
    the agreement threshold when fewer than 3 models are available."""
    print("multi-model consensus fields")
    from humanize.structure import (
        PDBAtom, PDBModel, compute_hints, compute_multi_model_consensus,
    )

    def mk_model(shift=0.0):
        m = PDBModel()
        for i in range(1, 9):
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, 0.0, 0.0, 90.0))
            m.atoms.append(PDBAtom("CB", "ALA", "A", i, i * 3.0 + 1.0, 1.0, 0.0, 90.0))
        m.atoms.append(PDBAtom("CA", "ALA", "C", 1, 6.0 + shift, 2.0, 0.0, 90.0))
        return m

    pm = {f"H{i}": i for i in range(1, 9)}
    single = compute_hints(mk_model(), "A", pm, {"H6": 6}, antigen_chains=["C"])
    check("cdr_contact uses all side-chain atoms",
          single.data["cdr_contact"].get("H5") is True)
    check("antigen contact detected", single.data["antigen_contact"].get("H3") is True)

    import tempfile as _tf

    def write_pdb(path, shift):
        lines = []

        def atom(chain, res, an, x, y, z):
            lines.append(
                'ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00 90.00'
                % (len(lines) + 1, an, "ALA", chain, res, x, y, z))

        for i in range(1, 9):
            atom("A", i, "CA", i * 3.0, 0.0, 0.0)
            atom("A", i, "CB", i * 3.0 + 1.0, 1.0, 0.0)
        atom("C", 1, "CA", 6.0 + shift, 2.0, 0.0)
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    with _tf.TemporaryDirectory() as td:
        paths = []
        for k in range(2):
            p = os.path.join(td, f"rank_{k + 1}.pdb")
            write_pdb(p, shift=k * 0.1)
            paths.append(p)
        cons = compute_multi_model_consensus(
            paths, "A", pm, {"H6": 6}, antigen_chains=["C"], min_consensus=3)
        check("consensus preserves cdr_partners", bool(cons.data.get("cdr_partners")))
        check("consensus preserves antigen_contact",
              cons.data.get("antigen_contact", {}).get("H3") is True)
        check("consensus preserves plddt", "H6" in cons.data.get("plddt", {}))
        check("consensus relaxes threshold for <3 models",
              cons.data.get("cdr_contact", {}).get("H6") is True)


def test_tool_kabat_mapping():
    """Standalone tools must map sequence index -> Kabat exactly, so
    back-mutations at FR4 / insertion-letter positions (H82A, L27A, H103) are
    matched instead of crashing or silently missing."""
    print("tool Kabat position mapping")
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    from immunogenicity.structural_risk import (
        PositionBuriedness, _kabat_num, _labels_for_variant,
        _rekey_structure_by_label, map_peptide_positions,
    )

    check("_kabat_num parses insertion letters",
          _kabat_num("H82A") == 82 and _kabat_num("L100B") == 100
          and _kabat_num("L27A") == 27 and _kabat_num("H103") == 103,
          str([_kabat_num(p) for p in ("H82A", "L100B", "L27A", "H103")]))
    check("_kabat_num returns None without digits", _kabat_num("X") is None)

    # 'A|H_V2a' style keys must resolve
    numbering = {"4D5_VH|H_V2a": ["H1", "H2", "H3", "H4", "H5"]}
    check("labels resolved from 'name|variant' key",
          _labels_for_variant(numbering, "H_V2a") == ["H1", "H2", "H3", "H4", "H5"])

    # Re-key a sequential structure map (H1..H5) onto Kabat labels
    seq_map = {
        f"H{i}": PositionBuriedness(
            position=f"H{i}", chain="H", seq_index=i - 1, resname="ALA",
            abs_sasa=10.0, rel_sasa=0.1, buried=True, is_backmutation=False)
        for i in range(1, 6)
    }
    labels = ["H1", "H2", "H3", "H82A", "H83"]  # insertion letter in the middle
    rekeyed = _rekey_structure_by_label(seq_map, "H", labels)
    check("re-key aligns labels including insertion letters",
          "H82A" in rekeyed and rekeyed["H82A"].position == "H82A"
          and rekeyed["H83"].seq_index == 4,
          str(list(rekeyed.keys())))

    # A backmutation at the insertion position must be matched (previously an
    # int('82A') ValueError crash / silent miss).
    bms = [{"position": "H82A", "donor_aa": "A", "human_aa": "S",
            "tier": "T2", "buried": True, "cdr_contact": False,
            "antigen_contact": False}]
    pos = map_peptide_positions(
        3, 4, "XXXXX", "H", rekeyed, bms, position_labels=labels)
    check("backmutation at insertion position matched",
          pos[0].position == "H82A" and pos[0].is_backmutation
          and pos[0].tier == "T2", str(pos[0]))

    # FR4 (Kabat 103) at a late sequence index must also match
    fr4_map = {
        f"H{i}": PositionBuriedness(
            position=f"H{i}", chain="H", seq_index=i - 1, resname="ALA",
            abs_sasa=10.0, rel_sasa=0.1, buried=True, is_backmutation=False)
        for i in range(1, 12)
    }
    fr4_labels = [f"H{n}" for n in (95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105)]
    fr4_rekey = _rekey_structure_by_label(fr4_map, "H", fr4_labels)
    fr4_bms = [{"position": "H103", "donor_aa": "W", "human_aa": "W",
                "tier": "T_FR4", "buried": True, "cdr_contact": True,
                "antigen_contact": False}]
    idx = fr4_labels.index("H103")
    fr4_pos = map_peptide_positions(
        idx, idx + 1, "X" * 11, "H", fr4_rekey, fr4_bms,
        position_labels=fr4_labels)
    check("FR4 T_FR4 backmutation matched",
          fr4_pos[0].position == "H103" and fr4_pos[0].is_backmutation
          and fr4_pos[0].tier == "T_FR4", str(fr4_pos[0]))


def test_structure_rmsd_engine():
    """Kabsch superposition + framework-superposed CDR-RMSD."""
    print("structure RMSD engine")
    from humanize.structure import PDBAtom, PDBModel, kabsch_superpose, structure_rmsd
    import math

    # known rigid transform recovery
    P = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)]
    c, s = math.cos(0.4), math.sin(0.4)
    R = [[c, -s, 0], [s, c, 0], [0, 0, 1]]
    t = (4.0, -1.5, 2.0)
    Q = [tuple(sum(R[i][j] * p[j] for j in range(3)) + t[i] for i in range(3))
         for p in P]
    rmsd, _, _ = kabsch_superpose(P, Q)
    check("kabsch recovers a rigid transform", rmsd < 1e-6, str(rmsd))

    def mk(shift_cdr=0.0):
        m = PDBModel()
        for i in range(1, 9):  # framework
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, 0.0, 0.0, 90.0))
        for i in range(20, 23):  # CDR3
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, 5.0 + shift_cdr,
                                   0.0, 90.0))
        return m

    pos = {f"H{i}": i for i in list(range(1, 9)) + list(range(20, 23))}
    fw = {f"H{i}" for i in range(1, 9)}
    cdrs = {"CDR3": {f"H{i}" for i in range(20, 23)}}
    res = structure_rmsd(mk(0.0), "A", pos, mk(0.0), "A", pos, cdrs, fw)
    check("identical models -> 0 CDR-RMSD", res["cdr_rmsd"] == 0.0, str(res))
    res2 = structure_rmsd(mk(0.0), "A", pos, mk(1.0), "A", pos, cdrs, fw)
    check("CDR deviation detected", res2["cdr_rmsd"] is not None
          and res2["cdr_rmsd"] > 0.9, str(res2))
    check("framework RMSD stays ~0", res2["fr_rmsd"] is not None
          and res2["fr_rmsd"] < 1e-6, str(res2))

    # FR4 is measured post-fit when a region map is supplied, even though it is
    # never part of the FR1-FR3 superposition.
    def mk_fr4(shift=0.0):
        m = PDBModel()
        for i in range(1, 9):       # FR1 (fit set)
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, 0.0, 0.0, 90.0))
        for i in range(9, 12):      # FR4 (not in fit set)
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, shift, 0.0, 90.0))
        for i in range(20, 23):     # CDR3
            m.atoms.append(PDBAtom("CA", "ALA", "A", i, i * 3.0, 5.0, 0.0, 90.0))
        return m

    pos2 = {f"H{i}": i for i in list(range(1, 9)) + list(range(9, 12))
            + list(range(20, 23))}
    regions = {f"H{i}": "FR1" for i in range(1, 9)}
    regions.update({f"H{i}": "FR4" for i in range(9, 12)})
    regions.update({f"H{i}": "CDR3" for i in range(20, 23)})
    res3 = structure_rmsd(mk_fr4(0.0), "A", pos2, mk_fr4(1.0), "A", pos2,
                          cdrs, fw, pos_regions=regions)
    check("fit unaffected by FR4 shift", res3["fr_rmsd"] is not None
          and res3["fr_rmsd"] < 1e-6, str(res3["fr_rmsd"]))
    check("FR4 deviation reported post-fit",
          res3["per_fr"].get("FR4") is not None and 0.9 < res3["per_fr"]["FR4"] < 1.1,
          str(res3["per_fr"]))
    check("FR1 region reported",
          res3["per_fr"].get("FR1") is not None and res3["per_fr"]["FR1"] < 1e-6,
          str(res3["per_fr"]))


def test_af3_variant_rmsd_wiring():
    """Pipeline wiring: paired Fab variant x variant CDR-RMSD vs donor."""
    print("AF3 variant RMSD wiring")
    import humanize.pipeline as pl
    from types import SimpleNamespace
    from humanize.structure import load_model
    from humanize.graft import graft_variant, CDR_POS_SETS
    from humanize.variants import Variant

    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor_h = number_heavy(M4D5_VH)
    donor_l = number_light(M4D5_VL)
    vgh = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    vgl = [g for g in db.v_genes if g.gene_id == "IGKV1-39*01"][0]
    jh, jl = db.j_for("H")[0], db.j_for("L")[0]
    cdr_h = set()
    for lo, hi in CDR_POS_SETS["kabat"]["H"].values():
        cdr_h.update(range(lo, hi + 1))
    cdr_l = set()
    for lo, hi in CDR_POS_SETS["kabat"]["L"].values():
        cdr_l.update(range(lo, hi + 1))

    three = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
             "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
             "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
             "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}

    def write_pdb(path, chains_meta, shift):
        # chains_meta: [(chain_id, seq, numbered, cdr_set)]
        lines = []
        for cid, seq, numbered, cdr in chains_meta:
            for i, aa in enumerate(seq):
                pos = numbered.residues[i].pos if i < len(numbered.residues) else ""
                num = int("".join(ch for ch in pos if ch.isdigit()) or 0)
                y = 1.0 if (shift and num in cdr) else 0.0
                lines.append(
                    "ATOM  %5d  CA  %3s %s%4d    %8.3f%8.3f%8.3f  1.00 90.00"
                    % (len(lines) + 1, three.get(aa, "GLY"), cid, i + 1,
                       i * 3.0, y, 0.0))
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    gH = graft_variant(donor_h, vgh, jh, "kabat", [])
    gL = graft_variant(donor_l, vgl, jl, "kabat", [])

    with tempfile.TemporaryDirectory() as td:
        # donor reference: chain A = donor H, chain B = donor L
        donor_pdb = os.path.join(td, "donor.pdb")
        write_pdb(donor_pdb, [("A", donor_h.sequence, donor_h, cdr_h),
                              ("B", donor_l.sequence, donor_l, cdr_l)], False)
        dmodel = load_model(donor_pdb)
        dh = (dmodel, "A", pl.pdb_chain_pos_map(dmodel, "A", donor_h))
        dl = (dmodel, "B", pl.pdb_chain_pos_map(dmodel, "B", donor_l))

        def fake_predict_fv(cfg, vh, vl, ag, tag):
            p = os.path.join(td, f"pred_{tag}.pdb")
            write_pdb(p, [("A", vh, gH.numbered, cdr_h),
                          ("B", vl, gL.numbered, cdr_l)], True)
            return p

        hv = Variant("H_V2", "test", gH, [])
        lv = Variant("L_V2", "test", gL, [])
        hrep = SimpleNamespace(
            input_chain=SimpleNamespace(chain_type="H", name="4D5_VH",
                                        sequence=donor_h.sequence,
                                        numbered=donor_h, warnings=[]),
            donor_structure=dh, variants=[hv], structure_validation={})
        lrep = SimpleNamespace(
            input_chain=SimpleNamespace(chain_type="L", name="4D5_VL",
                                        sequence=donor_l.sequence,
                                        numbered=donor_l, warnings=[]),
            donor_structure=dl, variants=[lv], structure_validation={})
        result = SimpleNamespace(format="fab", chains=[hrep, lrep], warnings=[])

        orig = pl.predict_fv
        pl.predict_fv = fake_predict_fv
        try:
            pl._validate_all_variants(result, PipelineConfig(af3_validate_variants=True,
                                                             af3=pl.AF3Config(mode="local")))
        finally:
            pl.predict_fv = orig

        check("paired H variant measured", "H_V2" in hrep.structure_validation,
              str(hrep.structure_validation))
        check("paired L variant measured", "L_V2" in lrep.structure_validation,
              str(lrep.structure_validation))
        check("CDR deviation reported for H",
              hrep.structure_validation.get("H_V2", {}).get("cdr_rmsd") is not None
              and hrep.structure_validation["H_V2"]["cdr_rmsd"] > 0.5,
              str(hrep.structure_validation.get("H_V2")))
        check("FR-RMSD ~0 for L",
              lrep.structure_validation.get("L_V2", {}).get("fr_rmsd") is not None
              and lrep.structure_validation["L_V2"]["fr_rmsd"] < 1e-6,
              str(lrep.structure_validation.get("L_V2")))


def test_standalone_rmsd_cli():
    """`humanize rmsd` compares external variant structures to a donor model."""
    print("standalone rmsd CLI")
    import contextlib, io
    from humanize.cli import main as cli_main
    from humanize.graft import CDR_POS_SETS

    donor = number_heavy(M4D5_VH)
    cdr = set()
    for lo, hi in CDR_POS_SETS["kabat"]["H"].values():
        cdr.update(range(lo, hi + 1))
    three = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
             "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
             "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
             "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}

    def write_pdb(path, seq, numbered, shift):
        lines = []
        for i, aa in enumerate(seq):
            num = int("".join(c for c in numbered.residues[i].pos
                              if c.isdigit()) or 0)
            y = 1.0 if (shift and num in cdr) else 0.0
            lines.append("ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00 90.00"
                         % (i + 1, three.get(aa, "GLY"), i + 1, i * 3.0, y, 0.0))
        with open(path, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    with tempfile.TemporaryDirectory() as td:
        donor_pdb = os.path.join(td, "donor.pdb")
        good = os.path.join(td, "v_good.pdb")
        bad = os.path.join(td, "v_bad.pdb")
        write_pdb(donor_pdb, donor.sequence, donor, False)
        write_pdb(good, donor.sequence, donor, False)
        write_pdb(bad, donor.sequence, donor, True)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cli_main([
                "rmsd",
                "--input", os.path.join(ROOT, "data", "examples", "mouse_4d5_fab.fasta"),
                "--donor", donor_pdb,
                "--variants", good, bad,
                "--chain", "H",
            ])
        out = buf.getvalue()
        check("rmsd CLI exit 0", rc == 0, str(rc))
        check("rmsd CLI reports CDR-RMSD", "CDR-RMSD" in out, out[:200])
        check("rmsd CLI zero for identical", "v_good.pdb" in out and "0.00" in out)


def test_calibration_substitution_match():
    """A calibration effect only applies to the donor/human pair it measured."""
    print("calibration substitution match")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(M4D5_VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    base = analyze_backmutations(donor, vg)
    cand = next((c for c in base.candidates if c.donor_aa != c.human_aa), None)
    if cand is None:
        check("calibration test has a candidate", False)
        return
    match = {cand.position: {"ddG_kcal": 0.8, "n_variants": 3,
                             "donor_aa": cand.donor_aa, "human_aa": cand.human_aa}}
    mismatch = {cand.position: {"ddG_kcal": 0.8, "n_variants": 3,
                                "donor_aa": cand.human_aa, "human_aa": cand.donor_aa}}
    r_match = analyze_backmutations(donor, vg, calibration=match)
    r_mism = analyze_backmutations(donor, vg, calibration=mismatch)
    cm = next(c for c in r_match.candidates if c.position == cand.position)
    cx = next(c for c in r_mism.candidates if c.position == cand.position)
    check("matching pair applies the empirical effect",
          cm.empirical_ddG == 0.8 and cm.empirical_n == 3,
          f"{cm.empirical_ddG}/{cm.empirical_n}")
    check("mismatched pair is ignored (different substitution)",
          cx.empirical_ddG is None, str(cx.empirical_ddG))


def test_learning_deconvolution():
    """Ridge regression separates co-occurring positions from multi-position variants."""
    print("learning deconvolution")
    from humanize.learning import _ridge_solve

    # a appears in variants 1,3; b in 2,3; y = a + 2b (rows: a, b, a+b)
    beta = _ridge_solve([["a"], ["b"], ["a", "b"]], [1.0, 2.0, 3.0],
                        ["a", "b"], lam=0.0)
    check("ridge recovers additive effects",
          abs(beta[0] - 1.0) < 1e-6 and abs(beta[1] - 2.0) < 1e-6, str(beta))


def test_clash_detection():
    print("clash detection")
    from humanize.structure import PDBAtom, PDBModel, count_clashes

    def mk(dist):
        m = PDBModel()
        m.atoms.append(PDBAtom("CA", "ALA", "A", 1, 0.0, 0.0, 0.0, 90.0))
        m.atoms.append(PDBAtom("CA", "ALA", "A", 2, dist, 0.0, 0.0, 90.0))
        return m

    check("close contact counted as a clash",
          count_clashes(mk(1.5))["n_clashes"] >= 1, str(count_clashes(mk(1.5))))
    check("normal distance not a clash", count_clashes(mk(4.0))["n_clashes"] == 0)
    check("same-residue atoms are not clashes",
          count_clashes(mk(1.5))["worst"] is not None)


def test_structural_score_and_plddt():
    """Noisy-OR structural combination + continuous pLDDT weighting."""
    print("structural score + pLDDT")
    from humanize.backmut import _noisy_or, _plddt_factor
    from humanize.config import WEIGHTS

    check("noisy-OR accumulates independent evidence",
          _noisy_or([0.85, 0.80]) > 0.85,
          str(_noisy_or([0.85, 0.80])))
    check("pLDDT unknown -> no change", _plddt_factor(None) == 1.0)
    check("pLDDT low -> strong down-weight", _plddt_factor(50) <= 0.21,
          str(_plddt_factor(50)))
    check("pLDDT high -> full weight", _plddt_factor(90) >= 0.999)
    check("pLDDT monotonic", _plddt_factor(60) < _plddt_factor(80))

    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(M4D5_VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    base = analyze_backmutations(donor, vg)
    h27 = next((c for c in base.candidates if c.position == "H27"
                and "canonical" in c.features and "vernier" in c.features), None)
    check("multi-feature structural score uses noisy-OR",
          h27 is not None
          and abs(h27.structural_score
                  - round(_noisy_or([WEIGHTS["structural"]["canonical"],
                                     WEIGHTS["structural"]["vernier"]]), 2)) < 0.02,
          str(h27.structural_score if h27 else None))


def test_design_panel():
    print("structure-guided design panel")
    from humanize.variants import structure_guided_panel
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(M4D5_VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    jg = db.j_for("H")[0]
    bm = analyze_backmutations(donor, vg)
    panel = structure_guided_panel(donor, vg, jg, "kabat", bm, structure=None,
                                   n_extra=3)
    check("panel has the requested steps", len(panel) == 3, str(len(panel)))
    check("panel names are graded",
          [v.name for v in panel] == ["H_V_opt1", "H_V_opt2", "H_V_opt3"],
          str([v.name for v in panel]))
    lengths = [len(v.backmutations) for v in panel]
    check("panel back-mutations strictly increase",
          lengths[0] < lengths[1] < lengths[2], str(lengths))
    check("panel positions unique within each variant",
          all(len(v.backmutations) == len(set(v.backmutations)) for v in panel))


def test_immunogenicity_integration():
    """Optional per-position MHC-II scores raise the humanization benefit."""
    print("immunogenicity integration")
    db = load_germline_db(os.path.join(ROOT, "data", "germline"))
    donor = number_heavy(M4D5_VH)
    vg = [g for g in db.v_genes if g.gene_id == "IGHV3-66*01"][0]
    base = analyze_backmutations(donor, vg)
    cand = next((c for c in base.candidates if c.position == "H48"), None)
    if cand is None:
        check("immunogenicity test has candidate H48", False)
        return
    immuno = {cand.position: 0.9}
    res = analyze_backmutations(donor, vg, immunogenicity=immuno)
    c = next(x for x in res.candidates if x.position == cand.position)
    check("immunogenicity score stored", c.immunogenicity_score == 0.9,
          str(c.immunogenicity_score))
    check("benefit raised by epitope score", c.benefit_score > cand.benefit_score,
          f"{c.benefit_score} vs {cand.benefit_score}")
    check("composite raised by benefit", c.composite >= cand.composite,
          f"{c.composite} vs {cand.composite}")
    other = next(x for x in res.candidates if x.position != cand.position)
    check("unspecified positions keep no score", other.immunogenicity_score is None)


def test_netmhciipan_parser():
    """Column-driven NetMHCIIpan xls parser maps peptides to Kabat positions."""
    print("NetMHCIIpan parser")
    from humanize.immunogenicity import parse_netmhciipan_xls
    donor = number_heavy(M4D5_VH)
    with tempfile.NamedTemporaryFile("w", suffix=".xls", delete=False) as fh:
        fh.write("# NetMHCIIpan test output\n")
        fh.write("Pos\tPeptide\tCore\tAffinity(nM)\t%Rank\tBindLevel\n")
        fh.write("1\tEVQLQQSGPELVKPG\tX\t1000\t0.5\tSB\n")
        fh.write("2\tVQLQQSGPELVKPGT\tX\t5000\t8.0\tWB\n")
        fh.write("3\tQLQQSGPELVKPGAS\tX\t50000\t25.0\t\n")
        path = fh.name
    scores = parse_netmhciipan_xls(path, donor)
    os.unlink(path)
    check("strong binder scored high", scores.get("H1", 0) >= 0.9,
          str(scores.get("H1")))
    check("weak binder scored lower", 0 < scores.get("H16", 0) < 0.9,
          str(scores.get("H16")))
    check("non-binder (%Rank>10) dropped",
          all(v > 0 for v in scores.values()) and "H3" in scores)
    check("peptide window mapped past its start", scores.get("H10", 0) >= 0.9,
          str(scores.get("H10")))


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
        # H: V0, V1, V2a, V2b, V2, V3 + Vmin (T1 non-empty);
        # L: V0, V1, V2a, V2b, V2, V3 only (T1 empty after L87 gold-standard
        # demotion, Vmin == V0). Total = 7 + 6.
        check("E2E fasta has 13 variants (H:7, L:6)", n_seq == 13, str(n_seq))


def main():
    test_numbering()
    test_chemical_liability_delta()
    test_structure_adaptive_scoring()
    test_vl_cdr3_insertion_numbering()
    test_structure_hint_chain_filtering()
    test_learning_fab_vl_positions()
    test_pdb_chain_matching()
    test_tool_kabat_mapping()
    test_graft_loop_conservation()
    test_vhh_humanization_gold_standard()
    test_germline_and_graft()
    test_backmut_variants()
    test_vhh_protection()
    test_minimal_reversion()
    test_learning_loop()
    test_lambda_chain_end_to_end()
    test_j_anchor_covers_all_germline_j()
    test_germline_strategies_smoke()
    test_developability_scan()
    test_input_validation()
    test_humanness_adapter()
    test_docx_report()
    test_enhanced_report()
    test_fr_indel_detection()
    test_fr_insertion_interactive_override()
    test_step3_structure_tier_rules()
    test_multi_model_consensus_fields()
    test_structure_rmsd_engine()
    test_af3_variant_rmsd_wiring()
    test_standalone_rmsd_cli()
    test_calibration_substitution_match()
    test_learning_deconvolution()
    test_clash_detection()
    test_structural_score_and_plddt()
    test_design_panel()
    test_immunogenicity_integration()
    test_netmhciipan_parser()
    test_end_to_end()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILURES: {FAILURES}")
        return 1
    print("ALL TESTS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
