#!/usr/bin/env python3
"""Structure-aware pI and Fv charge-patch analysis for antibodies.

The tool accepts a PDB/mmCIF structure (preferred) or VH/VL/VHH sequences.
When a structure is available it combines TAP/ANARCI/FreeSASA annotations with
an empirical pKa microenvironment model (or optional PROPKA) to estimate pI,
net charge, charge patches and developability risks. Without a structure it
falls back to sequence-based pI with a selectable pKa set.

Usage:
    python3 tools/pi_analysis/pi_analyzer.py --check
    python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --all-formats
    python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb \
        --pka-set Bjellqvist --all-pka-sets --structure-model heuristic
    python3 tools/pi_analysis/pi_analyzer.py --fasta inputs/antibody.fasta
    python3 tools/pi_analysis/pi_analyzer.py --vh <VH_SEQ> --vl <VL_SEQ>
    python3 tools/pi_analysis/pi_analyzer.py --pdb-dir af3/ --all-formats

pKa sets: Bjellqvist (ExPASy), IPC_protein, EMBOSS, Sillero, Nozaki_Tanford.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import pi_core  # noqa: E402
from pi_core import (  # noqa: E402
    PKA_SETS,
    ChainInfo,
    analyze_sequence,
    analyze_structure,
    clean_sequence,
    find_propka,
    load_patch_calibration,
    normalize_pka_set,
)

_STRUCTURE_EXTS = (".pdb", ".ent", ".cif", ".mmcif")
_METHOD = {
    "framework": "structure-aware antibody Fv pI / charge-patch analysis",
    "pka_sets": {
        "Bjellqvist": "ExPASy/Bjellqvist (1993) pKa set",
        "IPC_protein": "IPC_protein (Kozlowski 2016/2021)",
        "EMBOSS": "EMBOSS iep pKa set",
        "Sillero": "Sillero & Ribeiro (1989) pKa set",
        "Nozaki_Tanford": "Nozaki & Tanford (1967) pKa set",
    },
    "structure_model": (
        "heuristic: SASA burial + salt bridge + H-bond pKa shifts "
        "(clamped +/-2 pH units); propka: PROPKA3 pKa values"
    ),
    "references": [
        "Thorsteinson et al. mAbs 2021;13(1):1863582 (structure-based Fv charge/pI)",
        "Gupta et al. mAbs 2022;14(1):2020203 (Fv pI 7.5-9.0 developability window)",
        "Raybould et al. PNAS 2019;116(10):4025-4030 (TAP)",
        "Raybould et al. Commun Biol 2024;7:62 (TAP 2024 thresholds)",
        "Olsson et al. J Chem Theory Comput 2011;7:525-537 (PROPKA3)",
    ],
}


def _fmt(value: object, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    try:
        return "{:.{}f}".format(float(str(value)), digits)
    except (TypeError, ValueError):
        return str(value)


def _signed(value: object, digits: int = 2) -> str:
    text = _fmt(value, digits)
    if text == "n/a" or text.startswith("-"):
        return text
    return "+" + text


def _truncate(text: str, width: int = 120) -> str:
    return text if len(text) <= width else text[: width - 3] + "..."


def _input_description(result: Dict[str, object]) -> str:
    inp = result.get("input", {})
    if not isinstance(inp, dict):
        return str(result.get("name", ""))
    if inp.get("type") == "structure":
        kind = inp.get("kind", "structure")
        vh = inp.get("vh_chain") or "-"
        vl = inp.get("vl_chain") or "-"
        return "{} ({}; VH={}, VL={})".format(inp.get("path"), kind, vh, vl)
    chains = inp.get("chains", [])
    if isinstance(chains, list):
        return "sequence: " + ", ".join(
            "{} ({} aa)".format(c.get("name"), c.get("length"))
            for c in chains if isinstance(c, dict))
    return "sequence"


def _pka_sets(result: Dict[str, object]) -> List[str]:
    seq = result.get("sequence_pI", {})
    struct = result.get("structure_pI", {})
    computed = set()
    if isinstance(seq, dict):
        computed.update(seq.keys())
    if isinstance(struct, dict):
        computed.update(struct.keys())
    available = result.get("pka_sets_available", [])
    if isinstance(available, list) and available:
        shown = [str(s) for s in available if str(s) in computed]
        if shown:
            return shown
    return [str(result.get("pka_set", "IPC_protein"))]


def _pI_row(result: Dict[str, object], pka_set: str) -> Dict[str, object]:
    seq = result.get("sequence_pI", {})
    struct = result.get("structure_pI", {})
    seq_entry = seq.get(pka_set, {}) if isinstance(seq, dict) else {}
    struct_entry = struct.get(pka_set, {}) if isinstance(struct, dict) else {}
    if not isinstance(seq_entry, dict):
        seq_entry = {}
    if not isinstance(struct_entry, dict):
        struct_entry = {}

    use_struct = bool(struct_entry.get("pI_selected") is not None
                      or struct_entry.get("pI_free") is not None)
    entry = struct_entry if use_struct else seq_entry
    return {
        "pka_set": pka_set,
        "seq_pI": seq_entry.get("pI_selected", seq_entry.get("pI_free")),
        "struct_pI": struct_entry.get("pI_selected",
                                      struct_entry.get("pI_free")),
        "pI": entry.get("pI_selected", entry.get("pI_free")),
        "pI_free": entry.get("pI_free"),
        "pI_capped": entry.get("pI_capped"),
        "delta": struct_entry.get("delta_vs_sequence"),
        "q_ph74": (struct_entry.get("net_charge_ph74")
                   if use_struct else seq_entry.get("net_charge_selected_ph",
                                                    seq_entry.get("net_charge_ph74_free"))),
        "q_form": (struct_entry.get("net_charge_formulation_ph")
                   if use_struct else seq_entry.get("net_charge_selected_formulation_ph")),
        "n_groups": struct_entry.get("n_groups"),
    }


def check_tools() -> None:
    print("\n" + "=" * 68)
    print("  pI / Fv Charge Analyzer - Tool Availability")
    print("=" * 68)
    try:
        import anarci  # noqa: F401
        anarci_ok = True
    except ImportError:
        try:
            from scripts.humanize.anarci_adapter import is_anarci_available
            anarci_ok = bool(is_anarci_available())
        except Exception:
            anarci_ok = False
    try:
        import freesasa  # noqa: F401
        freesasa_ok = True
    except ImportError:
        freesasa_ok = False
    propka = find_propka()

    print("\n  Numbering:")
    print("    ANARCI IMGT:   {}".format("available" if anarci_ok else "not installed"))
    print("\n  Surface calculation:")
    print("    FreeSASA:      {}".format("available" if freesasa_ok else "not installed"))
    print("\n  Optional structure pKa:")
    print("    PROPKA3:       {}".format(propka if propka else "not installed (heuristic fallback)"))
    print("\n  pKa sets:")
    for name in PKA_SETS:
        print("    - {}".format(name))
    print("\n  Notes:")
    print("    - Structure mode needs ANARCI + FreeSASA.")
    print("    - Sequence mode is fully portable (stdlib only).")
    print("    - Default structure pKa model is PROPKA3; install propka3")
    print("      (missing -> automatic heuristic fallback with a warning).")
    print("")


def format_text_report(
    result: Dict[str, object],
    title: str = "pI / Fv Charge Analysis",
    max_patches: int = 12,
    show_ionizable: bool = False,
    show_charge_curve: bool = False,
) -> str:
    lines: List[str] = []
    lines.append("=" * 88)
    lines.append("  " + title)
    lines.append("=" * 88)
    lines.append("")
    lines.append("  Name:          {}".format(result.get("name", "")))
    lines.append("  Input:         {}".format(_input_description(result)))
    lines.append("  pKa set:       {}".format(result.get("pka_set", "")))
    lines.append("  pH:            {} (formulation pH {})".format(
        _fmt(result.get("ph"), 1), _fmt(result.get("formulation_ph"), 1)))
    if result.get("structure_pI"):
        lines.append("  Structure pKa: {} (SASA engine: {})".format(
            result.get("structure_model", "heuristic"),
            result.get("sasa_engine", "n/a")))
    lines.append("  Disulfides:    {}".format(result.get("disulfide_mode", "auto")))
    lines.append("  Termini mode:  {}   (ionic strength {})".format(
        result.get("termini_mode", "free"),
        _fmt(result.get("ionic_strength"), 3)))
    lines.append("")

    # Scope block: state explicitly which outputs a sequence-only run could
    # not compute, so that an absent developability flag is not read as a
    # clean one.
    if result.get("analysis_type") == "sequence_only":
        missing = result.get("unavailable_metrics")
        lines.append("-" * 88)
        lines.append("  SCOPE OF THIS ANALYSIS")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  Sequence-only analysis: no three-dimensional structure was used.")
        lines.append("  The following outputs require a structure and were NOT computed:")
        lines.append("")
        if isinstance(missing, list):
            for item in missing:
                if not isinstance(item, dict):
                    continue
                lines.append("    - {}".format(item.get("label", "")))
                lines.append("      {}".format(item.get("reason", "")))
        lines.append("")
        lines.append("  An absent flag is NOT a clean flag: an antibody with a real")
        lines.append("  hydrophobic-patch or charge-asymmetry flag is indistinguishable")
        lines.append("  from a clean one in this mode.")
        lines.append("")

    fv = result.get("fv", {})
    if isinstance(fv, dict):
        lines.append("-" * 88)
        lines.append("  Fv SUMMARY")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  pI (selected, {}):     {}".format(
            fv.get("termini_mode", result.get("termini_mode", "free")),
            _fmt(fv.get("pI_selected"))))
        rng = fv.get("pI_range")
        if isinstance(rng, list) and len(rng) == 2:
            lines.append("  pI range across sets:     {}-{}".format(
                _fmt(rng[0]), _fmt(rng[1])))
        lines.append("  pI (free termini):        {}".format(_fmt(fv.get("pI_free"))))
        lines.append("  pI (capped C-terminus):   {}".format(_fmt(fv.get("pI_capped"))))
        lines.append("  Net charge @ pH {}:      {}".format(
            _fmt(result.get("ph"), 1), _signed(fv.get("net_charge_ph74"))))
        lines.append("  Net charge @ formulation: {}   (pH {})".format(
            _signed(fv.get("net_charge_formulation_ph")),
            _fmt(result.get("formulation_ph"), 1)))
        if fv.get("Fv_chml_surface") is not None:
            lines.append("  Fv_chml total / surface:  {} / {}".format(
                _signed(fv.get("Fv_chml_total")),
                _signed(fv.get("Fv_chml_surface"))))
        elif fv.get("Fv_chml_total") is not None:
            lines.append("  Fv_chml (total charge):   {}".format(
                _signed(fv.get("Fv_chml_total"))))
        if fv.get("SFvCSP_surface") is not None:
            lines.append("  SFvCSP total / surface:   {} / {}".format(
                _fmt(fv.get("SFvCSP_total")), _fmt(fv.get("SFvCSP_surface"))))
        elif fv.get("SFvCSP_total") is not None:
            lines.append("  SFvCSP (total charge):    {}".format(
                _fmt(fv.get("SFvCSP_total"))))
        if fv.get("surface_charge_vh") is not None:
            lines.append("  Surface charge VH / VL:   {} / {}".format(
                _signed(fv.get("surface_charge_vh")),
                _signed(fv.get("surface_charge_vl"))))
        if fv.get("total_charge_vh") is not None:
            lines.append("  Total charge VH / VL:     {} / {}".format(
                _signed(fv.get("total_charge_vh")),
                _signed(fv.get("total_charge_vl"))))
        lines.append("")

    lines.append("-" * 88)
    lines.append("  pI BY PKA SET")
    lines.append("-" * 88)
    lines.append("")
    if result.get("structure_pI"):
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "pKa set", "seq pI", "struct pI", "delta", "struct cap", "q pH74", "q form"))
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "-" * 14, "-" * 9, "-" * 9, "-" * 9, "-" * 9, "-" * 9, "-" * 9))
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
                s, _fmt(row["seq_pI"]), _fmt(row["struct_pI"]),
                _signed(row["delta"]), _fmt(row["pI_capped"]),
                _signed(row["q_ph74"]), _signed(row["q_form"])))
    else:
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "pKa set", "seq pI", "struct pI", "capped", "q pH74"))
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "-" * 14, "-" * 9, "-" * 9, "-" * 9, "-" * 9))
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("  {:<14s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
                s, _fmt(row["seq_pI"]), _fmt(row["struct_pI"]),
                _fmt(row["pI_capped"]), _signed(row["q_ph74"])))
    lines.append("")

    chains = result.get("chains", {})
    if isinstance(chains, dict) and chains:
        lines.append("-" * 88)
        lines.append("  CHAINS")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "Chain", "Type", "Length", "pI sel", "pI free", "pI capped", "q pH74", "q form"))
        lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "-" * 8, "-" * 5, "-" * 7, "-" * 9, "-" * 9, "-" * 9, "-" * 9, "-" * 9))
        for name, c in chains.items():
            if not isinstance(c, dict):
                continue
            lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
                str(name), str(c.get("chain_type", "")), str(c.get("length", "")),
                _fmt(c.get("pI_selected", c.get("pI_free"))),
                _fmt(c.get("pI_free")), _fmt(c.get("pI_capped")),
                _signed(c.get("net_charge_ph74")),
                _signed(c.get("net_charge_formulation_ph"))))
        lines.append("")

    tap = result.get("tap", {})
    if isinstance(tap, dict) and tap.get("metrics"):
        lines.append("-" * 88)
        lines.append("  TAP METRICS")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  {:<10s} {:>12s} {:>10s}".format("Metric", "Value", "Flag"))
        lines.append("  {:<10s} {:>12s} {:>10s}".format("-" * 10, "-" * 12, "-" * 10))
        metrics = tap.get("metrics", {})
        flags = tap.get("flags", {})
        for metric in ("Ltot", "PSH", "PPC", "PNC", "SFvCSP"):
            if not isinstance(metrics, dict):
                continue
            value = metrics.get(metric)
            shown = "n/a" if value is None else (
                str(int(float(value))) if metric == "Ltot" else _fmt(value, 3))
            flag = flags.get(metric, "unknown") if isinstance(flags, dict) else "unknown"
            lines.append("  {:<10s} {:>12s} {:>10s}".format(metric, shown, str(flag).upper()))
        lines.append("")
        lines.append("  Overall TAP risk: {}".format(tap.get("risk_summary", "UNKNOWN")))
        lines.append("  Surface residues: {}; CDR vicinity: {}".format(
            tap.get("n_surface", "n/a"), tap.get("n_cdr_vicinity", "n/a")))
        lines.append("  Salt bridges: TAP {} (3.2 A exposed) / pKa microenvironment {} "
                     "(4.0 A, {} partner residues)".format(
                         tap.get("n_salt_bridges_tap", tap.get("n_salt_bridges", "n/a")),
                         tap.get("n_salt_bridges_pka", "n/a"),
                         tap.get("n_salt_bridge_residues_pka", "n/a")))
        lines.append("")

    quality = result.get("structure_quality")
    if isinstance(quality, dict):
        lines.append("-" * 88)
        lines.append("  STRUCTURE QUALITY")
        lines.append("-" * 88)
        lines.append("")
        bf = quality.get("bfactor_mean")
        lines.append("  CA B-factor / pLDDT proxy: {} (range {}-{}, pLDDT-like: {})".format(
            _fmt(bf, 1), _fmt(quality.get("bfactor_min"), 1),
            _fmt(quality.get("bfactor_max"), 1),
            quality.get("bfactor_plddt_like", "n/a")))
        lines.append("  Missing side chains: {}/{} ({} ionizable atom(s) missing)".format(
            quality.get("n_missing_sidechain", "n/a"),
            quality.get("n_residues", "n/a"),
            quality.get("n_missing_ionizable_atoms", "n/a")))
        lines.append("  Chain breaks (>{} A CA-CA): {}".format(
            quality.get("ca_break_cutoff", 4.5),
            quality.get("n_chain_breaks", "n/a")))
        lines.append("  Truncated termini: {}".format(
            ", ".join(str(t) for t in quality.get("truncated_termini", [])) or "none"))
        lines.append("")

    if show_charge_curve:
        curve = result.get("charge_curve")
        if isinstance(curve, list) and curve:
            lines.append("-" * 88)
            lines.append("  NET CHARGE VS pH")
            lines.append("-" * 88)
            lines.append("")
            for point in curve:
                if isinstance(point, dict):
                    lines.append("  pH {:>5s}  q = {:>8s}".format(
                        _fmt(point.get("ph"), 2),
                        _signed(point.get("net_charge"))))
            lines.append("")

    patches = result.get("patches", [])
    if isinstance(patches, list) and patches:
        lines.append("-" * 88)
        lines.append("  CHARGE PATCHES")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  {:<5s} {:<9s} {:>8s} {:>8s} {:<9s} {:<8s} {}".format(
            "pH", "Sign", "Charge", "Score", "Risk", "CDR/FR", "Residues / candidates"))
        lines.append("  {:<5s} {:<9s} {:>8s} {:>8s} {:<9s} {:<8s} {}".format(
            "-" * 5, "-" * 9, "-" * 8, "-" * 8, "-" * 9, "-" * 8, "-" * 24))
        ordered = sorted(
            [p for p in patches if isinstance(p, dict)],
            key=lambda p: (-abs(float(p.get("total_charge") or 0.0)),
                           -float(p.get("score") or 0.0), str(p.get("ph"))))
        shown = ordered if max_patches <= 0 else ordered[:max_patches]
        for p in shown:
            residues = p.get("residues", [])
            res_txt = " ".join(
                "{}{}".format(r.get("label"), r.get("aa"))
                for r in residues if isinstance(r, dict))
            cands = p.get("mutation_candidates", [])
            cand_txt = " ".join(
                "{}->{}".format(
                    c.get("label"),
                    "/".join((c.get("preferred_substitutions") or [])[:3]))
                for c in cands if isinstance(c, dict))
            detail = res_txt + (" | " + cand_txt if cand_txt else "")
            lines.append("  {:<5s} {:<9s} {:>8s} {:>8s} {:<9s} {:<8s} {}".format(
                _fmt(p.get("ph"), 1), str(p.get("sign", "")),
                _signed(p.get("total_charge")), _fmt(p.get("score")),
                str(p.get("risk", "")),
                "{}/{}".format(p.get("n_cdr", 0), p.get("n_framework", 0)),
                _truncate(detail, 100)))
        if max_patches > 0 and len(ordered) > max_patches:
            lines.append("  ... {} more patches in JSON/CSV".format(len(ordered) - max_patches))
        lines.append("")

    if show_ionizable:
        ions = result.get("ionizable_residues", [])
        if isinstance(ions, list) and ions:
            lines.append("-" * 88)
            lines.append("  IONIZABLE RESIDUES (structure-aware pKa)")
            lines.append("-" * 88)
            lines.append("")
            lines.append("  {:<10s} {:>4s} {:<8s} {:>7s} {:>7s} {:>7s} {:<12s} {}".format(
                "Residue", "AA", "Region", "pKa", "pKa eff", "d pKa", "Exposure", "Reason"))
            lines.append("  {:<10s} {:>4s} {:<8s} {:>7s} {:>7s} {:>7s} {:<12s} {}".format(
                "-" * 10, "-" * 4, "-" * 8, "-" * 7, "-" * 7, "-" * 7, "-" * 12, "-" * 20))
            for r in ions:
                if not isinstance(r, dict):
                    continue
                lines.append("  {:<10s} {:>4s} {:<8s} {:>7s} {:>7s} {:>7s} {:<12s} {}".format(
                    str(r.get("label", "")), str(r.get("aa", "")),
                    str(r.get("region", "")), _fmt(r.get("pka_model")),
                    _fmt(r.get("pka_eff")), _signed(r.get("delta_pka")),
                    str(r.get("exposure", "")), _truncate(str(r.get("reason", "")), 40)))
            lines.append("")

    recs = result.get("recommendations", [])
    if isinstance(recs, list) and recs:
        lines.append("-" * 88)
        lines.append("  RECOMMENDATIONS")
        lines.append("-" * 88)
        lines.append("")
        for rec in recs:
            if not isinstance(rec, dict):
                continue
            lines.append("  [{}] {}: {}".format(
                str(rec.get("severity", "info")).upper(),
                rec.get("category", ""), rec.get("message", "")))
            details = str(rec.get("details", "") or "")
            if details:
                for chunk in _wrap(details, 84):
                    lines.append("      " + chunk)
        lines.append("")

    warnings = result.get("warnings", [])
    if isinstance(warnings, list) and warnings:
        lines.append("-" * 88)
        lines.append("  WARNINGS")
        lines.append("-" * 88)
        lines.append("")
        for w in warnings:
            for chunk in _wrap(str(w), 84):
                lines.append("  NOTE: " + chunk)
        lines.append("")

    lines.append("-" * 88)
    lines.append("  METHOD")
    lines.append("-" * 88)
    lines.append("")
    lines.append("  pI is solved by bisection on the net-charge curve.")
    if result.get("analysis_type") == "structure":
        lines.append("  Structure mode uses per-residue SASA burial, salt bridges and H-bonds;")
        lines.append("  the heuristic pKa shift is clamped to +/-2 pH units and screened by")
        lines.append("  Debye-Huckel ionic strength. PROPKA is optional. Patch clustering uses")
        lines.append("  exposed charged residues (rSASA >= 7.5 %) within 8 A.")
    else:
        lines.append("  Sequence mode applies the selected pKa set to every ionizable group")
        lines.append("  without a microenvironment correction (no solvent exposure, salt bridge")
        lines.append("  or hydrogen-bond terms), so the reported pI is the free-solution value")
        lines.append("  of the modelled variable domains.")
    lines.append("  Computational pI is a screening estimate; confirm by cIEF/IEX.")
    lines.append("")
    return "\n".join(lines)


def _wrap(text: str, width: int) -> List[str]:
    words = text.split()
    out: List[str] = []
    line = ""
    for word in words:
        if not line:
            line = word
        elif len(line) + 1 + len(word) <= width:
            line += " " + word
        else:
            out.append(line)
            line = word
    if line:
        out.append(line)
    return out or [""]


def save_json_report(result: Dict[str, object], output_path: str) -> None:
    payload = dict(result)
    payload["method"] = _METHOD
    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)


def save_markdown_report(
    result: Dict[str, object],
    output_path: str,
    title: str = "pI / Fv Charge Analysis",
    max_patches: int = 12,
    show_ionizable: bool = False,
) -> None:
    lines: List[str] = ["# " + title, ""]
    lines.append("- **Name**: {}".format(result.get("name", "")))
    lines.append("- **Input**: {}".format(_input_description(result)))
    lines.append("- **pKa set**: {}".format(result.get("pka_set", "")))
    lines.append("- **pH**: {} (formulation {})".format(
        _fmt(result.get("ph"), 1), _fmt(result.get("formulation_ph"), 1)))
    if result.get("structure_pI"):
        lines.append("- **Structure pKa**: {} (SASA: {})".format(
            result.get("structure_model", "heuristic"),
            result.get("sasa_engine", "n/a")))
    lines.append("- **Termini mode**: {} (ionic strength {})".format(
        result.get("termini_mode", "free"),
        _fmt(result.get("ionic_strength"), 3)))
    lines.append("")

    if result.get("analysis_type") == "sequence_only":
        missing = result.get("unavailable_metrics")
        lines.append("> **Sequence-only analysis.** No three-dimensional "
                     "structure was used, so the outputs below require a "
                     "structure and were **not computed**:")
        lines.append(">")
        if isinstance(missing, list):
            for item in missing:
                if not isinstance(item, dict):
                    continue
                lines.append("> - **{}** — {}".format(
                    item.get("label", ""), item.get("reason", "")))
        lines.append(">")
        lines.append("> An absent flag is **not** a clean flag: an antibody "
                     "with a real hydrophobic-patch or charge-asymmetry flag "
                     "is indistinguishable from a clean one in this mode.")
        lines.append("")

    fv = result.get("fv", {})
    if isinstance(fv, dict):
        lines.append("## Fv summary")
        lines.append("")
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        lines.append("| pI (selected, {}) | {} |".format(
            fv.get("termini_mode", result.get("termini_mode", "free")),
            _fmt(fv.get("pI_selected"))))
        rng = fv.get("pI_range")
        if isinstance(rng, list) and len(rng) == 2:
            lines.append("| pI range across pKa sets | {}-{} |".format(
                _fmt(rng[0]), _fmt(rng[1])))
        lines.append("| pI (free termini) | {} |".format(_fmt(fv.get("pI_free"))))
        lines.append("| pI (capped C-terminus) | {} |".format(_fmt(fv.get("pI_capped"))))
        lines.append("| Net charge @ pH {} | {} |".format(
            _fmt(result.get("ph"), 1), _signed(fv.get("net_charge_ph74"))))
        lines.append("| Net charge @ formulation | {} |".format(
            _signed(fv.get("net_charge_formulation_ph"))))
        if fv.get("Fv_chml_total") is not None:
            lines.append("| Fv_chml total / surface | {} / {} |".format(
                _signed(fv.get("Fv_chml_total")),
                _signed(fv.get("Fv_chml_surface", "n/a"))))
        if fv.get("SFvCSP_total") is not None:
            lines.append("| SFvCSP total / surface | {} / {} |".format(
                _fmt(fv.get("SFvCSP_total")),
                _fmt(fv.get("SFvCSP_surface", "n/a"))))
        lines.append("")

    lines.append("## pI by pKa set")
    lines.append("")
    if result.get("structure_pI"):
        lines.append("| pKa set | seq pI | struct pI | delta | struct capped | q pH74 | q form |")
        lines.append("|---------|--------|-----------|-------|---------------|--------|--------|")
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                s, _fmt(row["seq_pI"]), _fmt(row["struct_pI"]),
                _signed(row["delta"]), _fmt(row["pI_capped"]),
                _signed(row["q_ph74"]), _signed(row["q_form"])))
    else:
        lines.append("| pKa set | seq pI | struct pI | capped | q pH74 |")
        lines.append("|---------|--------|-----------|--------|--------|")
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("| {} | {} | {} | {} | {} |".format(
                s, _fmt(row["seq_pI"]), _fmt(row["struct_pI"]),
                _fmt(row["pI_capped"]), _signed(row["q_ph74"])))
    lines.append("")

    chains = result.get("chains", {})
    if isinstance(chains, dict) and chains:
        lines.append("## Chains")
        lines.append("")
        lines.append("| Chain | Type | Length | pI sel | pI free | pI capped | q pH74 | q form |")
        lines.append("|-------|------|--------|--------|---------|-----------|--------|--------|")
        for name, c in chains.items():
            if not isinstance(c, dict):
                continue
            lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
                name, c.get("chain_type", ""), c.get("length", ""),
                _fmt(c.get("pI_selected", c.get("pI_free"))),
                _fmt(c.get("pI_free")), _fmt(c.get("pI_capped")),
                _signed(c.get("net_charge_ph74")),
                _signed(c.get("net_charge_formulation_ph"))))
        lines.append("")

    tap = result.get("tap", {})
    if isinstance(tap, dict) and tap.get("metrics"):
        lines.append("## TAP metrics")
        lines.append("")
        lines.append("| Metric | Value | Flag |")
        lines.append("|--------|-------|------|")
        metrics = tap.get("metrics", {})
        flags = tap.get("flags", {})
        for metric in ("Ltot", "PSH", "PPC", "PNC", "SFvCSP"):
            if not isinstance(metrics, dict):
                continue
            value = metrics.get(metric)
            shown = "n/a" if value is None else (
                str(int(float(value))) if metric == "Ltot" else _fmt(value, 3))
            flag = flags.get(metric, "unknown") if isinstance(flags, dict) else "unknown"
            lines.append("| {} | {} | {} |".format(metric, shown, str(flag).upper()))
        lines.append("")
        lines.append("**Overall TAP risk**: {}".format(tap.get("risk_summary", "UNKNOWN")))
        lines.append("")
        lines.append("Salt bridges: TAP {} (3.2 A exposed) / pKa microenvironment {} "
                     "(4.0 A, {} partner residues).".format(
                         tap.get("n_salt_bridges_tap", tap.get("n_salt_bridges", "n/a")),
                         tap.get("n_salt_bridges_pka", "n/a"),
                         tap.get("n_salt_bridge_residues_pka", "n/a")))
        lines.append("")

    quality = result.get("structure_quality")
    if isinstance(quality, dict):
        lines.append("## Structure quality")
        lines.append("")
        lines.append("- **CA B-factor / pLDDT proxy**: {} (range {}-{}, pLDDT-like: {})".format(
            _fmt(quality.get("bfactor_mean"), 1),
            _fmt(quality.get("bfactor_min"), 1),
            _fmt(quality.get("bfactor_max"), 1),
            quality.get("bfactor_plddt_like", "n/a")))
        lines.append("- **Missing side chains**: {}/{} ({} ionizable atom(s) missing)".format(
            quality.get("n_missing_sidechain", "n/a"),
            quality.get("n_residues", "n/a"),
            quality.get("n_missing_ionizable_atoms", "n/a")))
        lines.append("- **Chain breaks**: {} (>{} A CA-CA)".format(
            quality.get("n_chain_breaks", "n/a"),
            quality.get("ca_break_cutoff", 4.5)))
        lines.append("- **Truncated termini**: {}".format(
            ", ".join(str(t) for t in quality.get("truncated_termini", [])) or "none"))
        lines.append("")

    patches = result.get("patches", [])
    if isinstance(patches, list) and patches:
        lines.append("## Charge patches")
        lines.append("")
        lines.append("| pH | Sign | Charge | Score | Risk | CDR/FR | Residues | Candidates |")
        lines.append("|----|------|--------|-------|------|--------|----------|------------|")
        ordered = sorted(
            [p for p in patches if isinstance(p, dict)],
            key=lambda p: (-abs(float(p.get("total_charge") or 0.0)),
                           -float(p.get("score") or 0.0), str(p.get("ph"))))
        shown = ordered if max_patches <= 0 else ordered[:max_patches]
        for p in shown:
            residues = p.get("residues", [])
            res_txt = " ".join(
                "{}{}".format(r.get("label"), r.get("aa"))
                for r in residues if isinstance(r, dict))
            cands = p.get("mutation_candidates", [])
            cand_txt = " ".join(
                "{}->{}".format(
                    c.get("label"),
                    "/".join((c.get("preferred_substitutions") or [])[:3]))
                for c in cands if isinstance(c, dict))
            lines.append("| {} | {} | {} | {} | {} | {}/{} | {} | {} |".format(
                _fmt(p.get("ph"), 1), p.get("sign", ""),
                _signed(p.get("total_charge")), _fmt(p.get("score")),
                p.get("risk", ""), p.get("n_cdr", 0), p.get("n_framework", 0),
                res_txt or "-", cand_txt or "-"))
        if max_patches > 0 and len(ordered) > max_patches:
            lines.append("")
            lines.append("_{} more patches in JSON/CSV._".format(len(ordered) - max_patches))
        lines.append("")

    if show_ionizable:
        ions = result.get("ionizable_residues", [])
        if isinstance(ions, list) and ions:
            lines.append("## Ionizable residues")
            lines.append("")
            lines.append("| Residue | AA | Region | pKa | pKa eff | d pKa | Exposure | Reason |")
            lines.append("|---------|----|--------|-----|---------|-------|----------|--------|")
            for r in ions:
                if not isinstance(r, dict):
                    continue
                lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
                    r.get("label", ""), r.get("aa", ""), r.get("region", ""),
                    _fmt(r.get("pka_model")), _fmt(r.get("pka_eff")),
                    _signed(r.get("delta_pka")), r.get("exposure", ""),
                    r.get("reason", "")))
            lines.append("")

    recs = result.get("recommendations", [])
    if isinstance(recs, list) and recs:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recs:
            if not isinstance(rec, dict):
                continue
            lines.append("- **{} / {}**: {}".format(
                str(rec.get("severity", "info")).upper(),
                rec.get("category", ""), rec.get("message", "")))
            details = str(rec.get("details", "") or "")
            if details:
                lines.append("  - {}".format(details))
        lines.append("")

    warnings = result.get("warnings", [])
    if isinstance(warnings, list) and warnings:
        lines.append("## Warnings")
        lines.append("")
        for w in warnings:
            lines.append("> NOTE: {}".format(w))
        lines.append("")

    lines.append("## Method")
    lines.append("")
    lines.append("pI is solved by bisection on the net-charge curve. Structure mode uses")
    lines.append("per-residue SASA burial, salt bridges and H-bonds; the heuristic pKa shift")
    lines.append("is clamped to +/-2 pH units and screened by Debye-Huckel ionic strength.")
    lines.append("PROPKA is optional. Patch clustering uses exposed charged residues")
    lines.append("(rSASA >= 7.5 %) within 8 A. Computational pI is a screening estimate;")
    lines.append("confirm by cIEF/IEX.")
    lines.append("")
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))


def save_csv_reports(
    result: Dict[str, object],
    summary_path: str,
    patch_path: Optional[str] = None,
) -> None:
    fv = result.get("fv", {})
    tap = result.get("tap", {})
    if not isinstance(fv, dict):
        fv = {}
    if not isinstance(tap, dict):
        tap = {}
    inp = result.get("input", {})
    if not isinstance(inp, dict):
        inp = {}
    row = {
        "name": result.get("name", ""),
        "input_type": inp.get("type", ""),
        "input": inp.get("path", ""),
        "pka_set": result.get("pka_set", ""),
        "structure_model": result.get("structure_model", ""),
        "ph": result.get("ph", ""),
        "formulation_ph": result.get("formulation_ph", ""),
        "termini_mode": result.get("termini_mode", ""),
        "ionic_strength": result.get("ionic_strength", ""),
        "pI_selected": fv.get("pI_selected", ""),
        "pI_free": fv.get("pI_free", ""),
        "pI_capped": fv.get("pI_capped", ""),
        "pI_range_low": (fv.get("pI_range") or [""])[0]
                        if isinstance(fv.get("pI_range"), list) else "",
        "pI_range_high": (fv.get("pI_range") or ["", ""])[1]
                         if isinstance(fv.get("pI_range"), list) else "",
        "net_charge_ph74": fv.get("net_charge_ph74", ""),
        "net_charge_formulation_ph": fv.get("net_charge_formulation_ph", ""),
        "Fv_chml_total": fv.get("Fv_chml_total", ""),
        "Fv_chml_surface": fv.get("Fv_chml_surface", ""),
        "SFvCSP_total": fv.get("SFvCSP_total", ""),
        "SFvCSP_surface": fv.get("SFvCSP_surface", ""),
        "tap_risk": tap.get("risk_summary", ""),
        "n_salt_bridges_tap": tap.get("n_salt_bridges_tap", ""),
        "n_salt_bridges_pka": tap.get("n_salt_bridges_pka", ""),
        "n_patches": len(result.get("patches", []) or []),
        "n_warnings": len(result.get("warnings", []) or []),
    }
    with open(summary_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)

    if patch_path:
        patches = result.get("patches", [])
        with open(patch_path, "w", newline="", encoding="utf-8") as fh:
            fields = [
                "name", "patch_id", "ph", "sign", "total_charge", "score",
                "score_tap", "risk", "risk_basis", "n_cdr", "n_framework",
                "n_anchor", "in_cdr_vicinity", "terminal_residues",
                "chains", "residues", "candidates",
            ]
            writer = csv.DictWriter(fh, fieldnames=fields)
            writer.writeheader()
            if isinstance(patches, list):
                for p in patches:
                    if not isinstance(p, dict):
                        continue
                    residues = p.get("residues", [])
                    res_txt = " ".join(
                        "{}{}".format(r.get("label"), r.get("aa"))
                        for r in residues if isinstance(r, dict))
                    cands = p.get("mutation_candidates", [])
                    cand_txt = " ".join(
                        "{}->{}".format(
                            c.get("label"),
                            "/".join((c.get("preferred_substitutions") or [])[:3]))
                        for c in cands if isinstance(c, dict))
                    writer.writerow({
                        "name": result.get("name", ""),
                        "patch_id": p.get("patch_id", ""),
                        "ph": p.get("ph", ""),
                        "sign": p.get("sign", ""),
                        "total_charge": p.get("total_charge", ""),
                        "score": p.get("score", ""),
                        "score_tap": p.get("score_tap", ""),
                        "risk": p.get("risk", ""),
                        "risk_basis": p.get("risk_basis", ""),
                        "n_cdr": p.get("n_cdr", ""),
                        "n_framework": p.get("n_framework", ""),
                        "n_anchor": p.get("n_anchor", ""),
                        "in_cdr_vicinity": p.get("in_cdr_vicinity", ""),
                        "terminal_residues": " ".join(
                            str(r) for r in p.get("terminal_residues", [])),
                        "chains": " ".join(str(c) for c in p.get("chains", [])),
                        "residues": res_txt,
                        "candidates": cand_txt,
                    })


def _parse_ph_grid(text: str) -> List[float]:
    """Parse START:STOP:STEP (inclusive of START, not of STOP)."""
    parts = [p.strip() for p in text.split(":")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "--ph-grid must be START:STOP:STEP, e.g. 4:10:0.5")
    try:
        start, stop, step = (float(p) for p in parts)
    except ValueError:
        raise argparse.ArgumentTypeError("--ph-grid values must be numeric")
    if step <= 0 or stop <= start:
        raise argparse.ArgumentTypeError(
            "--ph-grid requires STOP > START and STEP > 0")
    values: List[float] = []
    x = start
    while x <= stop + 1e-9:
        values.append(round(x, 6))
        x += step
    return values


def _parse_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    seq: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(";"):
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name = line[1:].strip() or "sequence_{}".format(len(records) + 1)
                seq = []
            else:
                seq.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        records.append((name, "".join(seq)))
    if not records:
        raise ValueError("No FASTA records found in {}".format(path))
    return records


def _infer_chain_type(name: str, index: int, total: int) -> Tuple[str, str]:
    """Return (chain type, source) for one FASTA record.

    ``source`` is ``"header"`` when the record name states the chain type
    explicitly and ``"order"`` when the file order was used as a fallback.
    """
    upper = name.upper()
    if re.search(r"(^|[_\s.-])VH([_\s.-]|$)|HEAVY|(^|[_\s.-])H([_\s.-]|$)", upper):
        return "H", "header"
    if re.search(r"(^|[_\s.-])VL([_\s.-]|$)|LIGHT|KAPPA|LAMBDA|(^|[_\s.-])L([_\s.-]|$)", upper):
        return "L", "header"
    if total == 1:
        return "H", "order"
    return ("H" if index == 0 else "L"), "order"


def _infer_chain_type_from_sequence(seq: str) -> Optional[str]:
    """Classify a sequence as a heavy or a light chain from its content.

    A record name is only a hint; the sequence itself decides.  Two
    classifiers are tried, in order of confidence:

    1. ANARCI's HMM assignment (the same engine the structure path uses),
       which returns the best domain hit together with its chain type;
    2. the repo's portable numbering engine, whose conserved FR3 Cys anchor
       separates the two: a heavy V domain carries Cys92 and a light one
       Cys88, so for a well-formed domain exactly one of ``number_heavy`` /
       ``number_light`` succeeds.

    Returns ``"H"``, ``"L"`` or ``None`` when the sequence cannot be
    classified.  ``None`` means "no opinion", not "not an antibody".
    """
    if not seq:
        return None
    try:
        import tap_core
        ctype = tap_core._anarci_chain_type(seq)
    except Exception:
        ctype = None
    if ctype in ("H", "K", "L"):
        return "H" if ctype == "H" else "L"
    try:
        from scripts.humanize import numbering as _numbering
    except Exception:
        return None
    heavy_ok = light_ok = False
    try:
        _numbering.number_heavy(seq)
        heavy_ok = True
    except Exception:
        heavy_ok = False
    try:
        _numbering.number_light(seq)
        light_ok = True
    except Exception:
        light_ok = False
    if heavy_ok and not light_ok:
        return "H"
    if light_ok and not heavy_ok:
        return "L"
    if heavy_ok and light_ok:
        # A heavy-like domain that also numbers as a light chain is VHH-like.
        return "H"
    return None


def _load_fasta_chains(path: str) -> List[ChainInfo]:
    records = _parse_fasta(path)
    if len(records) > 2:
        raise ValueError(
            "FASTA has {} records; provide one VHH or two VH/VL records, "
            "or use --vh/--vl.".format(len(records)))
    cleaned = []
    for name, seq in records:
        clean, _warnings = clean_sequence(seq)
        if not clean:
            raise ValueError("Empty sequence in FASTA record {!r}".format(name))
        # The raw sequence is kept so that analyze_sequence can report the
        # cleaning warnings (dropped/substituted non-standard residues).
        cleaned.append((name, seq))
    if len(cleaned) == 1:
        return [ChainInfo(name="VHH", chain_type="H", sequence=cleaned[0][1])]

    inferred = [_infer_chain_type(name, i, 2) for i, (name, _seq) in enumerate(cleaned)]
    types = [t for t, _src in inferred]
    sources = [src for _t, src in inferred]

    # An explicit header wins, but the sequence content is used whenever the
    # header does not state the chain type: file order alone is a guess, and
    # a two-record FASTA whose light chain comes first would otherwise be
    # assigned backwards.
    if sources[0] == "order" or sources[1] == "order":
        content = [_infer_chain_type_from_sequence(seq)
                   for _name, seq in cleaned]
        if content[0] and content[1] and content[0] != content[1]:
            types = [content[0], content[1]]

    if types[0] == types[1]:
        # Still ambiguous: preserve file order as VH, VL.
        return [
            ChainInfo(name="VH", chain_type="H", sequence=cleaned[0][1]),
            ChainInfo(name="VL", chain_type="L", sequence=cleaned[1][1]),
        ]
    by_type = {t: seq for t, (_name, seq) in zip(types, cleaned)}
    return [
        ChainInfo(name="VH", chain_type="H", sequence=by_type["H"]),
        ChainInfo(name="VL", chain_type="L", sequence=by_type["L"]),
    ]


def _result_stem(result: Dict[str, object]) -> str:
    inp = result.get("input", {})
    path = inp.get("path") if isinstance(inp, dict) else ""
    if path:
        stem = Path(str(path)).stem
    else:
        stem = str(result.get("name", "pi"))
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem) or "pi"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Structure-aware pI and Fv charge-patch analysis for antibodies",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 tools/pi_analysis/pi_analyzer.py --check
  python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --all-formats
  python3 tools/pi_analysis/pi_analyzer.py --pdb af3/antibody.pdb --all-pka-sets
  python3 tools/pi_analysis/pi_analyzer.py --fasta inputs/antibody.fasta
  python3 tools/pi_analysis/pi_analyzer.py --vh <VH_SEQ> --vl <VL_SEQ>
  python3 tools/pi_analysis/pi_analyzer.py --pdb-dir af3/ --all-formats

pKa sets: Bjellqvist (ExPASy), IPC_protein, EMBOSS, Sillero, Nozaki_Tanford.
""")
    parser.add_argument("--check", action="store_true", help="Check available tools")
    parser.add_argument("--pdb", help="Path to a PDB/mmCIF structure")
    parser.add_argument("--pdb-dir", help="Directory of structures (recursive)")
    parser.add_argument("--fasta", help="FASTA with one VHH or two VH/VL records")
    parser.add_argument("--vh", help="VH amino acid sequence")
    parser.add_argument("--vl", help="VL amino acid sequence (omit for VHH)")
    parser.add_argument("--name", default="", help="Name for the analysis")
    parser.add_argument("--output", "-o", default="outputs/pi_analysis",
                        help="Output directory (default: outputs/pi_analysis)")
    parser.add_argument("--pka-set", default="IPC_protein", metavar="SET",
                        help="pKa set: Bjellqvist/ExPASy, IPC_protein, EMBOSS, "
                             "Sillero, Nozaki_Tanford (default: IPC_protein)")
    parser.add_argument("--all-pka-sets", action="store_true",
                        help="Report every bundled pKa set")
    parser.add_argument("--ph", type=float, default=7.4,
                        help="Physiological pH for net-charge reporting (default: 7.4)")
    parser.add_argument("--formulation-ph", type=float, default=6.0,
                        help="Formulation pH for patch/charge reporting (default: 6.0)")
    parser.add_argument("--disulfide", choices=("auto", "exclude", "free"),
                        default="auto",
                        help="Cys handling: auto = detect conserved disulfide "
                             "(default), exclude = all Cys non-ionizable, "
                             "free = all Cys ionizable")
    parser.add_argument("--structure-model", choices=("heuristic", "propka", "none"),
                        default="propka",
                        help="Structure pKa model: propka (default; falls back to "
                             "heuristic with a warning if PROPKA3 is missing), "
                             "heuristic, or none (no microenvironment correction)")
    parser.add_argument("--termini", choices=("auto", "free", "capped", "blocked"),
                        default="auto",
                        help="Terminus model: auto (Fab capped, VHH free), free, "
                             "capped C-terminus, or blocked N- and C-termini "
                             "(default: auto)")
    parser.add_argument("--ionic-strength", type=float, default=0.15, metavar="M",
                        help="Ionic strength for Debye-Huckel screening (default: 0.15 M)")
    parser.add_argument("--ph-grid", default=None, metavar="START:STOP:STEP",
                        help="pH grid for the charge-vs-pH curve "
                             "(default: 4:10:0.5)")
    parser.add_argument("--ph-curve", action="store_true",
                        help="Print the net-charge-vs-pH curve in the text report")
    parser.add_argument("--vh-chain", default=None,
                        help="Explicit structure chain id for VH")
    parser.add_argument("--vl-chain", default=None,
                        help="Explicit structure chain id for VL")
    parser.add_argument("--sasa", choices=("auto", "psa", "freesasa"),
                        default="auto",
                        help="Surface engine for structure mode (default: auto)")
    parser.add_argument("--psa-binary", default=None,
                        help="Path to the JOY psa executable")
    parser.add_argument("--tap-thresholds", choices=("2024", "2019", "exscientia", "legacy"),
                        default="2024", help="TAP flag thresholds (default: 2024)")
    parser.add_argument("--provenance", choices=("auto", "predicted", "experimental"),
                        default="auto",
                        help="How to read the B-factor column: auto (header, "
                             "then file content), predicted (pLDDT), or "
                             "experimental (crystallographic B-factor)")
    parser.add_argument("--patch-cutoff", type=float, default=pi_core.PATCH_CUTOFF,
                        metavar="A",
                        help="3D proximity cutoff for the charge-patch graph "
                             "(default: {:.1f} A)".format(pi_core.PATCH_CUTOFF))
    parser.add_argument("--temperature", type=float, default=298.15, metavar="K",
                        help="Temperature for the Debye-Huckel screening "
                             "(default: 298.15 K)")
    parser.add_argument("--propka-inactive-threshold", type=float,
                        default=pi_core.PROPKA_INACTIVE_PKA, metavar="PKA",
                        help="PROPKA pKa at or above which a group is treated "
                             "as non-ionizable (default: {:.0f})".format(
                                 pi_core.PROPKA_INACTIVE_PKA))
    parser.add_argument("--patch-calibration", default=None, metavar="FILE",
                        help="JSON file with local patch risk bands "
                             "({'positive': {...}, 'negative': {...}, "
                             "'cdr_factor': 1.25}); replaces the documented "
                             "default bands")
    parser.add_argument("--max-patches", type=int, default=12,
                        help="Patches shown in text/Markdown; 0 = all (default: 12)")
    parser.add_argument("--ionizable", action="store_true",
                        help="Include the ionizable-residue table in text/Markdown")
    parser.add_argument("--json", action="store_true", help="Save JSON report")
    parser.add_argument("--markdown", action="store_true", help="Save Markdown report")
    parser.add_argument("--csv", action="store_true", help="Save CSV summary + patches")
    parser.add_argument("--text", action="store_true", help="Save the text report")
    parser.add_argument("--all-formats", action="store_true",
                        help="Save text + JSON + Markdown + CSV")

    args = parser.parse_args()

    if args.check:
        check_tools()
        return
    if not any((args.pdb, args.pdb_dir, args.fasta, args.vh)):
        parser.error("Specify --check, --pdb, --pdb-dir, --fasta, or --vh/--vl")
    if args.fasta and (args.vh or args.vl):
        parser.error("--fasta cannot be combined with --vh/--vl")
    if args.vl_chain and not args.vh_chain:
        parser.error("--vl-chain requires --vh-chain")
    if (args.vh_chain or args.vl_chain) and not (args.pdb or args.pdb_dir):
        parser.error("--vh-chain/--vl-chain require --pdb or --pdb-dir")
    ph_grid = _parse_ph_grid(args.ph_grid) if args.ph_grid else None

    try:
        args.pka_set = normalize_pka_set(args.pka_set)
    except ValueError as exc:
        parser.error(str(exc))

    try:
        patch_calibration = load_patch_calibration(args.patch_calibration)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error("cannot read --patch-calibration: {}".format(exc))

    os.makedirs(args.output, exist_ok=True)
    results: List[Dict[str, object]] = []

    if args.pdb:
        print("\nAnalysing structure: {}".format(args.pdb))
        try:
            result = analyze_structure(
                args.pdb,
                pka_set=args.pka_set,
                all_pka_sets=args.all_pka_sets,
                ph=args.ph,
                formulation_ph=args.formulation_ph,
                vh_chain=args.vh_chain,
                vl_chain=args.vl_chain,
                sasa_engine=args.sasa,
                psa_binary=args.psa_binary,
                disulfide_mode=args.disulfide,
                structure_model=args.structure_model,
                tap_thresholds=args.tap_thresholds,
                termini_mode=args.termini,
                ionic_strength=args.ionic_strength,
                ph_grid=ph_grid,
                name=args.name,
                patch_calibration=patch_calibration,
                provenance_override=args.provenance,
                patch_cutoff=args.patch_cutoff,
                temperature=args.temperature,
                propka_inactive_threshold=args.propka_inactive_threshold,
            )
        except Exception as exc:
            print("  ERROR: {}".format(exc), file=sys.stderr)
            sys.exit(1)
        results.append(result)
        print(format_text_report(
            result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
            max_patches=args.max_patches, show_ionizable=args.ionizable,
            show_charge_curve=args.ph_curve))

    elif args.pdb_dir:
        files = [p for p in Path(args.pdb_dir).rglob("*")
                 if p.suffix.lower() in _STRUCTURE_EXTS]
        if not files:
            print("No PDB/CIF files found in {}".format(args.pdb_dir), file=sys.stderr)
            sys.exit(1)
        print("\nFound {} structure files".format(len(files)))
        for pdb_path in files:
            print("\nAnalysing: {}".format(pdb_path))
            try:
                result = analyze_structure(
                    str(pdb_path),
                    pka_set=args.pka_set,
                    all_pka_sets=args.all_pka_sets,
                    ph=args.ph,
                    formulation_ph=args.formulation_ph,
                    vh_chain=args.vh_chain,
                    vl_chain=args.vl_chain,
                    sasa_engine=args.sasa,
                    psa_binary=args.psa_binary,
                    disulfide_mode=args.disulfide,
                    structure_model=args.structure_model,
                    tap_thresholds=args.tap_thresholds,
                    termini_mode=args.termini,
                    ionic_strength=args.ionic_strength,
                    ph_grid=ph_grid,
                    name=args.name or pdb_path.stem,
                    patch_calibration=patch_calibration,
                    provenance_override=args.provenance,
                patch_cutoff=args.patch_cutoff,
                temperature=args.temperature,
                propka_inactive_threshold=args.propka_inactive_threshold,
                )
            except Exception as exc:
                print("  ERROR: {}".format(exc), file=sys.stderr)
                continue
            results.append(result)
            print(format_text_report(
                result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
                max_patches=args.max_patches, show_ionizable=args.ionizable,
                show_charge_curve=args.ph_curve))

    else:
        if args.fasta:
            try:
                chains = _load_fasta_chains(args.fasta)
            except Exception as exc:
                print("  ERROR: {}".format(exc), file=sys.stderr)
                sys.exit(1)
            input_path = args.fasta
            default_name = Path(args.fasta).stem
        else:
            clean_vh, _ = clean_sequence(args.vh)
            if not clean_vh:
                print("  ERROR: empty --vh sequence", file=sys.stderr)
                sys.exit(1)
            if args.vl:
                clean_vl, _ = clean_sequence(args.vl)
                if not clean_vl:
                    print("  ERROR: empty --vl sequence", file=sys.stderr)
                    sys.exit(1)
                # The raw sequences are passed on so that analyze_sequence
                # can report the cleaning warnings (dropped/substituted
                # non-standard residues).
                chains = [
                    ChainInfo(name="VH", chain_type="H", sequence=args.vh),
                    ChainInfo(name="VL", chain_type="L", sequence=args.vl),
                ]
            else:
                chains = [ChainInfo(name="VHH", chain_type="H", sequence=args.vh)]
            input_path = ""
            default_name = "sequence"

        print("\nAnalysing sequences (structure-aware fields unavailable)")
        try:
            result = analyze_sequence(
                chains,
                pka_set=args.pka_set,
                all_pka_sets=args.all_pka_sets,
                ph=args.ph,
                formulation_ph=args.formulation_ph,
                disulfide_mode=args.disulfide,
                termini_mode=args.termini,
                name=args.name or default_name,
                input_path=input_path,
                patch_calibration=patch_calibration,
            )
        except Exception as exc:
            print("  ERROR: {}".format(exc), file=sys.stderr)
            sys.exit(1)
        results.append(result)
        print(format_text_report(
            result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
            max_patches=args.max_patches, show_ionizable=args.ionizable,
            show_charge_curve=args.ph_curve))

    for result in results:
        stem = _result_stem(result)
        if args.text or args.all_formats:
            path = os.path.join(args.output, "{}_pi.txt".format(stem))
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(format_text_report(
                    result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
                    max_patches=args.max_patches, show_ionizable=args.ionizable,
                    show_charge_curve=args.ph_curve))
            print("  Text report saved: {}".format(path))
        if args.json or args.all_formats:
            path = os.path.join(args.output, "{}_pi.json".format(stem))
            save_json_report(result, path)
            print("  JSON report saved: {}".format(path))
        if args.markdown or args.all_formats:
            path = os.path.join(args.output, "{}_pi.md".format(stem))
            save_markdown_report(
                result, path, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
                max_patches=args.max_patches, show_ionizable=args.ionizable)
            print("  Markdown report saved: {}".format(path))
        if args.csv or args.all_formats:
            summary = os.path.join(args.output, "{}_pi.csv".format(stem))
            patches = os.path.join(args.output, "{}_pi_patches.csv".format(stem))
            save_csv_reports(result, summary, patches)
            print("  CSV reports saved: {}, {}".format(summary, patches))

    if len(results) > 1 and (args.json or args.all_formats):
        combined = os.path.join(args.output, "combined_pi_results.json")
        # The combined file is a cross-structure summary, but it must carry the
        # same contract as the individual reports: a consumer has to be able
        # to tell which pKa model produced a number and whether a metric was
        # computed at all.  An omitted field reads like a clean flag, so
        # `analysis_type` / `unavailable_metrics` / `structure_model` are
        # included here as well.
        with open(combined, "w", encoding="utf-8") as fh:
            json.dump([{
                "name": r.get("name", ""),
                "input": r.get("input", {}),
                "pka_set": r.get("pka_set", ""),
                "structure_model": r.get("structure_model", ""),
                "analysis_type": r.get("analysis_type", ""),
                "unavailable_metrics": r.get("unavailable_metrics", []),
                "termini_mode": r.get("termini_mode", ""),
                "ionic_strength": r.get("ionic_strength", ""),
                "fv": r.get("fv", {}),
                "chains": r.get("chains", {}),
                "sensitivity": r.get("sensitivity", {}),
                "charge_curve": r.get("charge_curve", []),
                "structure_quality": r.get("structure_quality", {}),
                "tap": r.get("tap", {}),
                "patch_calibration": r.get("patch_calibration", {}),
                "n_patches": len(r.get("patches", []) or []),
                "patches": r.get("patches", []),
                "recommendations": r.get("recommendations", []),
                "warnings": r.get("warnings", []),
            } for r in results], fh, indent=2, ensure_ascii=False)
        print("\nCombined JSON report saved: {}".format(combined))


if __name__ == "__main__":
    main()
