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

from pi_core import (  # noqa: E402
    PKA_SETS,
    ChainInfo,
    analyze_sequence,
    analyze_structure,
    clean_sequence,
    find_propka,
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

    use_struct = bool(struct_entry.get("pI_free") is not None)
    entry = struct_entry if use_struct else seq_entry
    return {
        "pka_set": pka_set,
        "seq_pI": seq_entry.get("pI_free"),
        "struct_pI": struct_entry.get("pI_free"),
        "pI": entry.get("pI_free"),
        "pI_capped": entry.get("pI_capped"),
        "delta": struct_entry.get("delta_vs_sequence"),
        "q_ph74": (struct_entry.get("net_charge_ph74")
                   if use_struct else seq_entry.get("net_charge_ph74_free")),
        "q_form": (struct_entry.get("net_charge_formulation_ph")
                   if use_struct else None),
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
    print("    - Install PROPKA3 for --structure-model propka.")
    print("")


def format_text_report(
    result: Dict[str, object],
    title: str = "pI / Fv Charge Analysis",
    max_patches: int = 12,
    show_ionizable: bool = False,
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
    lines.append("")

    fv = result.get("fv", {})
    if isinstance(fv, dict):
        lines.append("-" * 88)
        lines.append("  Fv SUMMARY")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  pI (free termini):        {}".format(_fmt(fv.get("pI_free"))))
        lines.append("  pI (capped C-terminus):   {}".format(_fmt(fv.get("pI_capped"))))
        lines.append("  Net charge @ pH {}:      {}".format(
            _fmt(result.get("ph"), 1), _signed(fv.get("net_charge_ph74"))))
        lines.append("  Net charge @ formulation: {}   (pH {})".format(
            _signed(fv.get("net_charge_formulation_ph")),
            _fmt(result.get("formulation_ph"), 1)))
        if fv.get("Fv_chml") is not None:
            lines.append("  Fv_chml (VH-VL charge):   {}".format(_signed(fv.get("Fv_chml"))))
        if fv.get("SFvCSP") is not None:
            lines.append("  SFvCSP (VH*VL charge):    {}".format(_fmt(fv.get("SFvCSP"))))
        if fv.get("surface_charge_vh") is not None:
            lines.append("  Surface charge VH / VL:   {} / {}".format(
                _signed(fv.get("surface_charge_vh")),
                _signed(fv.get("surface_charge_vl"))))
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
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s}".format(
            "pKa set", "seq pI", "capped", "q pH74"))
        lines.append("  {:<14s} {:>9s} {:>9s} {:>9s}".format(
            "-" * 14, "-" * 9, "-" * 9, "-" * 9))
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("  {:<14s} {:>9s} {:>9s} {:>9s}".format(
                s, _fmt(row["seq_pI"]), _fmt(row["pI_capped"]),
                _signed(row["q_ph74"])))
    lines.append("")

    chains = result.get("chains", {})
    if isinstance(chains, dict) and chains:
        lines.append("-" * 88)
        lines.append("  CHAINS")
        lines.append("-" * 88)
        lines.append("")
        lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "Chain", "Type", "Length", "pI free", "pI capped", "q pH74", "q form"))
        lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
            "-" * 8, "-" * 5, "-" * 7, "-" * 9, "-" * 9, "-" * 9, "-" * 9))
        for name, c in chains.items():
            if not isinstance(c, dict):
                continue
            lines.append("  {:<8s} {:>5s} {:>7s} {:>9s} {:>9s} {:>9s} {:>9s}".format(
                str(name), str(c.get("chain_type", "")), str(c.get("length", "")),
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
        lines.append("  Surface residues: {}; CDR vicinity: {}; salt bridges: {}".format(
            tap.get("n_surface", "n/a"), tap.get("n_cdr_vicinity", "n/a"),
            tap.get("n_salt_bridges", "n/a")))
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
    lines.append("  Structure mode uses per-residue SASA burial, salt bridges and H-bonds;")
    lines.append("  the heuristic pKa shift is clamped to +/-2 pH units. PROPKA is optional.")
    lines.append("  Patch clustering uses exposed charged residues (rSASA >= 7.5 %) within 8 A.")
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
    lines.append("")

    fv = result.get("fv", {})
    if isinstance(fv, dict):
        lines.append("## Fv summary")
        lines.append("")
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        lines.append("| pI (free termini) | {} |".format(_fmt(fv.get("pI_free"))))
        lines.append("| pI (capped C-terminus) | {} |".format(_fmt(fv.get("pI_capped"))))
        lines.append("| Net charge @ pH {} | {} |".format(
            _fmt(result.get("ph"), 1), _signed(fv.get("net_charge_ph74"))))
        lines.append("| Net charge @ formulation | {} |".format(
            _signed(fv.get("net_charge_formulation_ph"))))
        if fv.get("Fv_chml") is not None:
            lines.append("| Fv_chml | {} |".format(_signed(fv.get("Fv_chml"))))
        if fv.get("SFvCSP") is not None:
            lines.append("| SFvCSP | {} |".format(_fmt(fv.get("SFvCSP"))))
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
        lines.append("| pKa set | seq pI | capped | q pH74 |")
        lines.append("|---------|--------|--------|--------|")
        for s in _pka_sets(result):
            row = _pI_row(result, s)
            lines.append("| {} | {} | {} | {} |".format(
                s, _fmt(row["seq_pI"]), _fmt(row["pI_capped"]),
                _signed(row["q_ph74"])))
    lines.append("")

    chains = result.get("chains", {})
    if isinstance(chains, dict) and chains:
        lines.append("## Chains")
        lines.append("")
        lines.append("| Chain | Type | Length | pI free | pI capped | q pH74 | q form |")
        lines.append("|-------|------|--------|---------|-----------|--------|--------|")
        for name, c in chains.items():
            if not isinstance(c, dict):
                continue
            lines.append("| {} | {} | {} | {} | {} | {} | {} |".format(
                name, c.get("chain_type", ""), c.get("length", ""),
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
    lines.append("is clamped to +/-2 pH units. PROPKA is optional. Patch clustering uses")
    lines.append("exposed charged residues (rSASA >= 7.5 %) within 8 A. Computational pI is a")
    lines.append("screening estimate; confirm by cIEF/IEX.")
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
        "pI_free": fv.get("pI_free", ""),
        "pI_capped": fv.get("pI_capped", ""),
        "net_charge_ph74": fv.get("net_charge_ph74", ""),
        "net_charge_formulation_ph": fv.get("net_charge_formulation_ph", ""),
        "Fv_chml": fv.get("Fv_chml", ""),
        "SFvCSP": fv.get("SFvCSP", ""),
        "tap_risk": tap.get("risk_summary", ""),
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
                "risk", "n_cdr", "n_framework", "n_anchor", "in_cdr_vicinity",
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
                        "risk": p.get("risk", ""),
                        "n_cdr": p.get("n_cdr", ""),
                        "n_framework": p.get("n_framework", ""),
                        "n_anchor": p.get("n_anchor", ""),
                        "in_cdr_vicinity": p.get("in_cdr_vicinity", ""),
                        "chains": " ".join(str(c) for c in p.get("chains", [])),
                        "residues": res_txt,
                        "candidates": cand_txt,
                    })


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


def _infer_chain_type(name: str, index: int, total: int) -> str:
    upper = name.upper()
    if re.search(r"(^|[_\s.-])VH([_\s.-]|$)|HEAVY|(^|[_\s.-])H([_\s.-]|$)", upper):
        return "H"
    if re.search(r"(^|[_\s.-])VL([_\s.-]|$)|LIGHT|KAPPA|LAMBDA|(^|[_\s.-])L([_\s.-]|$)", upper):
        return "L"
    if total == 1:
        return "H"
    return "H" if index == 0 else "L"


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
        cleaned.append((name, clean))

    if len(cleaned) == 1:
        return [ChainInfo(name="VHH", chain_type="H", sequence=cleaned[0][1])]

    types = [_infer_chain_type(name, i, 2) for i, (name, _seq) in enumerate(cleaned)]
    if types[0] == types[1]:
        # Ambiguous headers: preserve file order as VH, VL.
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
                        default="heuristic",
                        help="Structure pKa model: heuristic (default), propka, "
                             "or none (no microenvironment correction)")
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
    if args.vh_chain and not args.vl_chain or args.vl_chain and not args.vh_chain:
        parser.error("--vh-chain and --vl-chain must be given together")
    if (args.vh_chain or args.vl_chain) and not (args.pdb or args.pdb_dir):
        parser.error("--vh-chain/--vl-chain require --pdb or --pdb-dir")

    try:
        args.pka_set = normalize_pka_set(args.pka_set)
    except ValueError as exc:
        parser.error(str(exc))

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
                name=args.name,
            )
        except Exception as exc:
            print("  ERROR: {}".format(exc), file=sys.stderr)
            sys.exit(1)
        results.append(result)
        print(format_text_report(
            result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
            max_patches=args.max_patches, show_ionizable=args.ionizable))

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
                    name=args.name or pdb_path.stem,
                )
            except Exception as exc:
                print("  ERROR: {}".format(exc), file=sys.stderr)
                continue
            results.append(result)
            print(format_text_report(
                result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
                max_patches=args.max_patches, show_ionizable=args.ionizable))

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
                chains = [
                    ChainInfo(name="VH", chain_type="H", sequence=clean_vh),
                    ChainInfo(name="VL", chain_type="L", sequence=clean_vl),
                ]
            else:
                chains = [ChainInfo(name="VHH", chain_type="H", sequence=clean_vh)]
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
                name=args.name or default_name,
                input_path=input_path,
            )
        except Exception as exc:
            print("  ERROR: {}".format(exc), file=sys.stderr)
            sys.exit(1)
        results.append(result)
        print(format_text_report(
            result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
            max_patches=args.max_patches, show_ionizable=args.ionizable))

    for result in results:
        stem = _result_stem(result)
        if args.text or args.all_formats:
            path = os.path.join(args.output, "{}_pi.txt".format(stem))
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(format_text_report(
                    result, "pI / Fv Charge Analysis: {}".format(result.get("name", "")),
                    max_patches=args.max_patches, show_ionizable=args.ionizable))
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
        with open(combined, "w", encoding="utf-8") as fh:
            json.dump([{
                "name": r.get("name", ""),
                "input": r.get("input", {}),
                "pka_set": r.get("pka_set", ""),
                "fv": r.get("fv", {}),
                "tap": r.get("tap", {}),
                "n_patches": len(r.get("patches", []) or []),
                "recommendations": r.get("recommendations", []),
                "warnings": r.get("warnings", []),
            } for r in results], fh, indent=2, ensure_ascii=False)
        print("\nCombined JSON report saved: {}".format(combined))


if __name__ == "__main__":
    main()
