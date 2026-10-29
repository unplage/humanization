#!/usr/bin/env python3
"""Post-humanization immunogenicity analysis tool.

Predicts MHC-II T-cell epitope burden for humanized antibody variants,
comparing them to the donor and V0 graft.

Backends: HLAIIPred (Apache-2.0, recommended) / heuristic (fallback).

Input:
  - variants.fasta  (from pipeline output or any multi-sequence FASTA)
  - output dir      (auto-discovers variants.fasta)
  - donor FASTA     (for donor-vs-humanized comparison)

Usage:
  python3 tools/immunogenicity/immunogenicity_analyzer.py --input variants.fasta
  python3 tools/immunogenicity/immunogenicity_analyzer.py --outdir outputs/amg110_step3/step3
  python3 tools/immunogenicity/immunogenicity_analyzer.py --input variants.fasta --backend heuristic --all-formats
  python3 tools/immunogenicity/immunogenicity_analyzer.py --check
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import textwrap
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Ensure tools dir is on path for backends import
_tools_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if _tools_dir not in sys.path:
    sys.path.insert(0, _tools_dir)

from immunogenicity.backends import (
    ALLELE_PANELS, DEFAULT_ALLELES, DEFAULT_STRONG_THRESHOLD,
    DEFAULT_WEAK_THRESHOLD, MAX_PEP_LEN, MIN_PEP_LEN,
    auto_select_backend, check_backends, get_backend,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRONG_THRESHOLD = DEFAULT_STRONG_THRESHOLD
WEAK_THRESHOLD = DEFAULT_WEAK_THRESHOLD


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> List[Tuple[str, str, str]]:
    """Parse FASTA into list of (header, sequence, name).

    Header format: '>name|variant|description' or '>name description'
    Returns: [(full_header_without_>, sequence, short_name)]
    """
    entries = []
    current_header = None
    current_seq = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    seq = "".join(current_seq)
                    # Parse name from header
                    hdr = current_header[1:].strip()
                    parts = hdr.split("|")
                    name = parts[0].strip() if parts else hdr
                    entries.append((hdr, seq, name))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_header is not None:
        seq = "".join(current_seq)
        hdr = current_header[1:].strip()
        parts = hdr.split("|")
        name = parts[0].strip() if parts else hdr
        entries.append((hdr, seq, name))
    return entries


def discover_variants(outdir: str) -> Optional[str]:
    """Look for variants.fasta in an output directory tree."""
    candidates = [
        os.path.join(outdir, "variants.fasta"),
    ]
    for sub in ["step2", "step3"]:
        candidates.append(os.path.join(outdir, sub, "variants.fasta"))
    for c in candidates:
        if os.path.isfile(c):
            return c
    # Recursive search (shallow)
    for root, dirs, files in os.walk(outdir):
        if "variants.fasta" in files:
            return os.path.join(root, "variants.fasta")
        # Don't go too deep
        depth = root.replace(outdir, "").count(os.sep)
        if depth > 3:
            dirs.clear()
    return None


# ---------------------------------------------------------------------------
# Sequence grouping
# ---------------------------------------------------------------------------

_H_RE = re.compile(r"(?:^|[|_.\s-])(VHH|VH|HEAVY|H)(?:[|_.\s-]|$)", re.I)
_L_RE = re.compile(r"(?:^|[|_.\s-])(VL|VK|LIGHT|KAPPA|LAMBDA|L)(?:[|_.\s-]|$)", re.I)
_VH_FIRST = {
    "EVQL", "QVQL", "EVKL", "QVTL", "QMQL", "QITL", "DVQL",
    "QLQL", "QSVE", "QSLE", "VQLV", "QSVQ", "EVQM",
}


def _infer_chain_type(hdr: str, seq: str, name: str) -> str:
    hdr_upper = hdr.upper()
    if _H_RE.search(hdr_upper):
        return "H"
    if _L_RE.search(hdr_upper):
        return "L"
    seq_upper = (seq or "").upper()
    if seq_upper[:4] in _VH_FIRST:
        return "H"
    if "VHH" in hdr_upper or "NANOBODY" in hdr_upper:
        return "H"
    if seq_upper[:4] in {"DIQM", "DIVM", "EIVL", "EIQM", "DIVL"}:
        return "L"
    return "H" if len(seq_upper) > 200 else "L"


def group_by_chain(entries: List[Tuple[str, str, str]]) -> Dict[str, Dict[str, str]]:
    """Group sequences by chain type (H/L) and variant name.

    Returns: {chain_type: {variant_name: sequence}}
    """
    chains: Dict[str, Dict[str, str]] = defaultdict(dict)
    for hdr, seq, name in entries:
        chain_type = _infer_chain_type(hdr, seq, name)
        parts = hdr.split("|")
        variant_name = parts[1].strip() if len(parts) > 1 else name
        chains[chain_type][variant_name] = seq
    return dict(chains)


def _is_reference_variant(vname: str, chain_type: str,
                          reference_variant: Optional[str]) -> bool:
    upper = (vname or "").upper()
    if reference_variant:
        ref = reference_variant.strip()
        return (vname == ref or upper == ref.upper()
                or vname.split("|")[-1].strip() == ref
                or upper.split("|")[-1].strip() == ref.upper())
    if upper in {"V0", "DONOR", "GRAFT", "REFERENCE"}:
        return True
    if upper == f"{chain_type.upper()}_V0":
        return True
    if upper in {"COMBO_V0", "COMBO_V0_V0"}:
        return True
    if upper.startswith("COMBO_"):
        return False
    return upper.endswith("_V0") and upper.count("_") == 1


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

@dataclass
class VariantAnalysis:
    name: str
    chain_type: str
    sequence: str
    epitope_scores: Dict[int, float] = field(default_factory=dict)
    n_strong: int = 0
    n_weak: int = 0
    n_total_windows: int = 0
    epitope_density: float = 0.0
    mean_epitope_load: float = 0.0
    raw_mean_epitope_load: float = 0.0
    max_epitope: float = 0.0
    immunogenicity_score: Optional[float] = 0.0
    per_peptide_hits: list = field(default_factory=list)
    positions_of_interest: List[Dict] = field(default_factory=list)
    error: Optional[str] = None
    evaluated: bool = True
    structural_adjustment: Optional[float] = None


def compute_epitope_metrics(per_residue: Dict[int, float],
                            seq_len: int,
                            per_peptide: Optional[list] = None,
                            strong_threshold: float = STRONG_THRESHOLD,
                            weak_threshold: float = WEAK_THRESHOLD) -> Dict:
    """Compute aggregate metrics.

    Metrics are peptide-window based whenever per-peptide predictions are
    available. ``density`` is the fraction of peptide windows predicted as
    strong or weak binders; it is not a residue-coverage fraction.
    """
    scores: List[float] = []
    if per_peptide:
        for hit in per_peptide:
            if isinstance(hit, dict):
                scores.append(float(hit.get("score", 0.0)))
            else:
                scores.append(float(getattr(hit, "score", 0.0)))
        n_total = len(scores)
    elif per_residue:
        scores = [float(v) for v in per_residue.values()]
        n_total = len(scores)
    else:
        return {"n_strong": 0, "n_weak": 0, "density": 0.0, "mean": 0.0,
                "max": 0.0, "n_total_windows": 0}

    n_strong = sum(1 for s in scores if s >= strong_threshold)
    n_weak = sum(1 for s in scores if weak_threshold <= s < strong_threshold)
    density = (n_strong + n_weak) / max(1, n_total)
    mean_load = sum(scores) / max(1, len(scores))
    return {
        "n_strong": n_strong,
        "n_weak": n_weak,
        "density": round(density, 4),
        "mean": round(mean_load, 4),
        "max": round(max(scores), 4),
        "n_total_windows": n_total,
    }


def compute_immunogenicity_score(metrics: Dict,
                                 donor_metrics: Optional[Dict] = None) -> float:
    """Relative epitope load; 100 = same as reference.

    Values above 100 are allowed and mean the variant has a higher predicted
    epitope load than the reference. Without a reference, the absolute mean
    window score is returned on a 0-100 scale.
    """
    variant_load = metrics.get("mean", 0.0)
    donor_load = (donor_metrics or {}).get("mean", 0.0)
    if donor_load and donor_load > 1e-12:
        return round(variant_load / donor_load * 100.0, 1)
    return round(variant_load * 100.0, 1)


def _align_index_pairs(a: str, b: str) -> List[Tuple[int, int]]:
    """Global alignment; return (index_in_a, index_in_b) pairs for aligned columns."""
    if a == b:
        return [(i, i) for i in range(len(a))]
    match, mismatch, gap = 2, -1, -2
    n, m = len(a), len(b)
    if n * m > 4_000_000:
        return [(i, i) for i in range(min(n, m))]
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap
        trace[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap
        trace[0][j] = 2
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = score[i - 1][j - 1] + (match if ai == b[j - 1] else mismatch)
            up = score[i - 1][j] + gap
            left = score[i][j - 1] + gap
            best = max(diag, up, left)
            score[i][j] = best
            trace[i][j] = 0 if best == diag else (1 if best == up else 2)
    pairs: List[Tuple[int, int]] = []
    i, j = n, m
    while i > 0 or j > 0:
        t = trace[i][j]
        if i > 0 and j > 0 and t == 0:
            pairs.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or t == 1):
            i -= 1
        else:
            j -= 1
    pairs.reverse()
    return pairs


def find_positions_of_interest(donor_seq: str, variant_seq: str,
                               donor_epi: Dict,
                               variant_epi: Dict) -> List[Dict]:
    """Find aligned positions where the epitope score changes.

    ``position`` is a 1-based index into the variant chain sequence (MHC
    backends score residue offsets). Use the structural-risk report for
    Kabat/IMGT labels.
    """
    if donor_epi and isinstance(next(iter(donor_epi)), str):
        donor_epi = {int(k): v for k, v in donor_epi.items()}
    if variant_epi and isinstance(next(iter(variant_epi)), str):
        variant_epi = {int(k): v for k, v in variant_epi.items()}

    positions = []
    for d_idx, v_idx in _align_index_pairs(donor_seq, variant_seq):
        d_aa = donor_seq[d_idx]
        v_aa = variant_seq[v_idx]
        if d_aa == v_aa:
            continue
        d_score = float(donor_epi.get(d_idx, 0.0))
        v_score = float(variant_epi.get(v_idx, 0.0))
        delta = v_score - d_score
        if abs(delta) <= 0.05:
            continue
        positions.append({
            "position": v_idx + 1,
            "donor_position": d_idx + 1,
            "donor_aa": d_aa,
            "variant_aa": v_aa,
            "donor_epitope": round(d_score, 3),
            "variant_epitope": round(v_score, 3),
            "delta": round(delta, 3),
            "direction": "reduced" if delta < 0 else "increased",
        })
    positions.sort(key=lambda x: abs(x["delta"]), reverse=True)
    return positions[:30]


# ---------------------------------------------------------------------------
# Main analysis function
# ---------------------------------------------------------------------------

def analyze_immunogenicity(
    input_path: Optional[str] = None,
    outdir: Optional[str] = None,
    donor_path: Optional[str] = None,
    backend_name: str = "auto",
    alleles: Optional[List[str]] = None,
    pep_len: int = 15,
    conda_envs: Optional[Dict[str, str]] = None,
    verbose: bool = False,
    reference_variant: Optional[str] = None,
    strong_threshold: float = STRONG_THRESHOLD,
    weak_threshold: float = WEAK_THRESHOLD,
) -> Dict:
    """Run immunogenicity analysis on variants.

    ``donor_path`` is always scored as the reference when supplied. Without it,
    an exact reference variant (default: V0 / H_V0 / L_V0) is used. Textual
    substring matching is deliberately not used because combined variants such
    as ``combo_V0_V1`` contain V0 in only one chain.
    """
    if alleles is None:
        alleles = DEFAULT_ALLELES
    if not (MIN_PEP_LEN <= int(pep_len) <= MAX_PEP_LEN):
        raise ValueError(
            f"--pep-len must be between {MIN_PEP_LEN} and {MAX_PEP_LEN}")
    if not (0.0 <= weak_threshold < strong_threshold <= 1.0):
        raise ValueError("require 0 <= weak_threshold < strong_threshold <= 1")

    fasta_path = input_path
    if not fasta_path and outdir:
        fasta_path = discover_variants(outdir)
    if not fasta_path or not os.path.isfile(fasta_path):
        raise FileNotFoundError(
            f"No variants.fasta found. Provide --input or --outdir containing variants.fasta")

    entries = parse_fasta(fasta_path)
    if not entries:
        raise ValueError(f"No sequences found in {fasta_path}")

    chains = group_by_chain(entries)
    if not chains:
        raise ValueError(f"No chain sequences found in {fasta_path}")

    donor_chains = {}
    if donor_path and os.path.isfile(donor_path):
        donor_chains = group_by_chain(parse_fasta(donor_path))
    elif donor_path:
        raise FileNotFoundError(f"Donor FASTA not found: {donor_path}")

    if backend_name == "auto":
        mhc_backend = auto_select_backend()
        if not getattr(mhc_backend, "is_neural", False):
            print("  WARNING: HLAIIPred is unavailable; using the allele-agnostic "
                  "heuristic backend (not a trained MHC-II model).", file=sys.stderr)
    else:
        mhc_backend = get_backend(backend_name)
        if not mhc_backend.is_available():
            raise RuntimeError(f"Backend '{backend_name}' is not available")

    if verbose:
        print(f"  Backend: {mhc_backend.describe()}")
        print(f"  Alleles: {', '.join(alleles)}")
        print(f"  Peptide length: {pep_len}")

    all_results = {
        "input_fasta": fasta_path,
        "alleles": alleles,
        "pep_len": pep_len,
        "mhc_backend": mhc_backend.name,
        "strong_threshold": strong_threshold,
        "weak_threshold": weak_threshold,
        "score_semantics": "percent of reference mean peptide-window score; >100 means higher than reference",
        "chains": {},
    }

    for chain_type in sorted(chains.keys()):
        chain_variants = chains[chain_type]
        if verbose:
            print(f"\n  Chain {chain_type}: {len(chain_variants)} variants")

        chain_types = {vname: chain_type for vname in chain_variants}
        mhc_results = mhc_backend.score(
            chain_variants, alleles, pep_len=pep_len, chain_types=chain_types)

        reference_name = next(
            (v for v in chain_variants
             if _is_reference_variant(v, chain_type, reference_variant)),
            None,
        )

        donor_seq = ""
        donor_source = "none"
        donor_metrics = None
        donor_epi: Dict[int, float] = {}

        if chain_type in donor_chains and donor_chains[chain_type]:
            donor_seq = next(iter(donor_chains[chain_type].values()))
            donor_source = "donor_file"
            donor_res = mhc_backend.score(
                {"__donor__": donor_seq}, alleles, pep_len=pep_len,
                chain_types={"__donor__": chain_type}).get("__donor__")
            if donor_res and not donor_res.error:
                donor_epi = donor_res.per_residue
                donor_metrics = compute_epitope_metrics(
                    donor_epi, len(donor_seq), donor_res.per_peptide,
                    strong_threshold, weak_threshold)
            elif verbose and donor_res:
                print(f"    donor scoring failed: {donor_res.error}")
        elif reference_name:
            donor_seq = chain_variants.get(reference_name, "")
            donor_source = "reference_variant"
            donor_res = mhc_results.get(reference_name)
            if donor_res and not donor_res.error:
                donor_epi = donor_res.per_residue
                donor_metrics = compute_epitope_metrics(
                    donor_epi, len(donor_seq), donor_res.per_peptide,
                    strong_threshold, weak_threshold)

        variant_analyses = {}
        for vname, seq in chain_variants.items():
            backend_res = mhc_results.get(vname)
            if backend_res is None:
                variant_analyses[vname] = VariantAnalysis(
                    name=vname, chain_type=chain_type, sequence=seq,
                    error="backend returned no result", evaluated=False,
                    immunogenicity_score=None)
                continue
            if backend_res.error:
                if verbose:
                    print(f"    {vname}: ERROR - {backend_res.error}")
                variant_analyses[vname] = VariantAnalysis(
                    name=vname, chain_type=chain_type, sequence=seq,
                    error=backend_res.error, evaluated=False,
                    immunogenicity_score=None)
                continue
            if not backend_res.per_peptide and not backend_res.per_residue:
                variant_analyses[vname] = VariantAnalysis(
                    name=vname, chain_type=chain_type, sequence=seq,
                    error="no peptide predictions", evaluated=False,
                    immunogenicity_score=None)
                continue

            epi = backend_res.per_residue
            metrics = compute_epitope_metrics(
                epi, len(seq), backend_res.per_peptide,
                strong_threshold, weak_threshold)
            imm_score = compute_immunogenicity_score(metrics, donor_metrics)
            poi = (find_positions_of_interest(donor_seq, seq, donor_epi, epi)
                   if donor_seq else [])

            va = VariantAnalysis(
                name=vname,
                chain_type=chain_type,
                sequence=seq,
                epitope_scores=epi,
                n_strong=metrics["n_strong"],
                n_weak=metrics["n_weak"],
                n_total_windows=metrics["n_total_windows"],
                epitope_density=metrics["density"],
                mean_epitope_load=metrics["mean"],
                raw_mean_epitope_load=metrics["mean"],
                max_epitope=metrics["max"],
                immunogenicity_score=imm_score,
                per_peptide_hits=[asdict(h) for h in backend_res.per_peptide],
                positions_of_interest=poi,
            )
            variant_analyses[vname] = va
            if verbose:
                print(f"    {vname}: score={imm_score}, strong={metrics['n_strong']}, "
                      f"density={metrics['density']:.3f}")

        all_results["chains"][chain_type] = {
            "variants": {k: asdict(v) for k, v in variant_analyses.items()},
            "donor_metrics": donor_metrics,
            "reference_name": reference_name,
            "reference_source": donor_source,
            "reference_length": len(donor_seq),
        }

    return all_results


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def _variant_key(name: str):
    return [int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", name)]


def _sorted_variants(variants: Dict) -> List[str]:
    return sorted(variants.keys(), key=_variant_key)


def _score_text(v: Dict) -> str:
    if v.get("error") or v.get("immunogenicity_score") is None:
        return "N/A"
    suffix = "*" if v.get("structural_score_applied") else ""
    return f"{v['immunogenicity_score']:.1f}%{suffix}"


def _adj_text(v: Dict) -> str:
    adj = v.get("structural_adjustment")
    return "-" if adj is None else f"{adj:.3f}"


def format_text_report(results: Dict, title: str = "Immunogenicity Analysis") -> str:
    """Format results as text report (stdout)."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"  {title}")
    lines.append("=" * 80)
    lines.append(f"  Input:    {results['input_fasta']}")
    lines.append(f"  Backend:  {results['mhc_backend']}")
    alleles = results["alleles"]
    lines.append(f"  Alleles:  {', '.join(alleles[:5])}"
                 + (" ..." if len(alleles) > 5 else ""))
    lines.append(f"  Pep len:  {results['pep_len']}")
    lines.append(f"  Strong/weak thresholds: {results.get('strong_threshold', 0.5):.2f} / "
                 f"{results.get('weak_threshold', 0.3):.2f}")
    lines.append("  Strong/Weak = peptide windows; Density = (strong+weak)/windows")
    lines.append("")

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        donor_m = cdata.get("donor_metrics")
        lines.append(f"  Chain {chain_type} ({len(variants)} variants)")
        lines.append("  " + "-" * 88)
        ref_name = cdata.get("reference_name") or "-"
        ref_src = cdata.get("reference_source") or "none"
        lines.append(f"  Reference: {ref_src} ({ref_name})")
        if cdata.get("structural_reference_mean") is not None:
            sref = cdata.get("structural_reference_name") or "raw donor mean"
            lines.append(f"  Structural reference: {sref}")
            lines.append("  (* = final score rewritten by structural adjustment)")
        lines.append("")

        lines.append(f"  {'Variant':<22} {'Score':>8} {'Strong':>7} {'Weak':>7} "
                     f"{'Density':>9} {'Mean EPI':>10} {'StructAdj':>10}")
        lines.append("  " + "-" * 84)
        for vname in _sorted_variants(variants):
            v = variants[vname]
            mean = v.get("mean_epitope_load")
            mean_txt = "N/A" if mean is None else f"{mean:.3f}"
            lines.append(
                f"  {vname:<22} {_score_text(v):>8} "
                f"{v.get('n_strong', 0):>7} {v.get('n_weak', 0):>7} "
                f"{v.get('epitope_density', 0.0):>8.3f} {mean_txt:>10} "
                f"{_adj_text(v):>10}"
            )
        lines.append("")

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("error"):
                lines.append(f"  {vname}: ERROR - {v['error']}")
            if v.get("structural_error"):
                lines.append(f"  {vname}: STRUCTURAL WARNING - {v['structural_error']}")
        if any(v.get("error") or v.get("structural_error")
               for v in variants.values()):
            lines.append("")

        if donor_m:
            lines.append(f"  Reference metrics: mean_epitope={donor_m['mean']:.3f}, "
                         f"strong={donor_m['n_strong']}, weak={donor_m['n_weak']}, "
                         f"density={donor_m['density']:.3f}")
            lines.append("")

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("positions_of_interest"):
                lines.append(f"  Top epitope-changing positions ({vname}):")
                lines.append("  (Pos = 1-based sequence index in the variant, not Kabat/IMGT)")
                lines.append(f"  {'VPos':>5} {'DPos':>5} {'Donor':>6} {'Variant':>8} "
                             f"{'DeltaEPI':>9} {'Direction':>10}")
                for p in v["positions_of_interest"][:10]:
                    lines.append(
                        f"  {p['position']:>5} {p.get('donor_position', 0):>5} "
                        f"{p['donor_aa']:>6} {p['variant_aa']:>8} "
                        f"{p['delta']:>+9.3f} {p['direction']:>10}"
                    )
                lines.append("")

        strong_thr = results.get("strong_threshold", 0.5)
        for vname in _sorted_variants(variants):
            peptides = variants[vname].get("per_peptide_hits", [])
            if not peptides:
                continue
            top = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top if p.get("score", 0) >= strong_thr]
            if not strong:
                continue
            lines.append(f"  Top predicted epitopes ({vname}):")
            lines.append(f"  {'Peptide':<18} {'Score':>7} {'Core':<12} {'Allele':<24}")
            lines.append("  " + "-" * 66)
            for p in strong:
                lines.append(
                    f"  {p.get('peptide', ''):<18} {p.get('score', 0):>7.3f} "
                    f"{p.get('core', ''):<12} {p.get('allele', ''):<24}"
                )
            lines.append("")

    return "\n".join(lines)


def save_json_report(results: Dict, output_path: str):
    """Save results as JSON."""
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False, default=str)


def save_markdown_report(results: Dict, output_path: str,
                         title: str = "Immunogenicity Analysis"):
    """Save results as Markdown."""
    lines = []
    lines.append(f"# {title}\n")
    lines.append(f"**Input:** `{results['input_fasta']}`")
    lines.append(f"**Backend:** {results['mhc_backend']}")
    alleles = results["alleles"]
    lines.append(f"**Alleles:** {', '.join(alleles[:5])}"
                 + (" ..." if len(alleles) > 5 else ""))
    lines.append(f"**Peptide length:** {results['pep_len']}")
    lines.append(f"**Strong/weak thresholds:** {results.get('strong_threshold', 0.5):.2f} / "
                 f"{results.get('weak_threshold', 0.3):.2f}")
    lines.append("**Score semantics:** "
                 f"{results.get('score_semantics', 'percent of reference mean peptide-window score')}.")
    lines.append("A `*` marks variants whose final score was rewritten after "
                 "structural adjustment.")
    lines.append("")

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        lines.append(f"## Chain {chain_type} ({len(variants)} variants)\n")
        lines.append(f"**Reference:** {cdata.get('reference_source', 'none')} "
                     f"(`{cdata.get('reference_name') or '-'}`)\n")
        if cdata.get("structural_reference_mean") is not None:
            sref = cdata.get("structural_reference_name") or "raw donor mean"
            lines.append(f"**Structural reference:** `{sref}`\n")

        lines.append("| Variant | Score | Strong | Weak | Density | Mean EPI | StructAdj |")
        lines.append("|---------|-------|--------|------|---------|----------|-----------|")
        for vname in _sorted_variants(variants):
            v = variants[vname]
            mean = v.get("mean_epitope_load")
            mean_txt = "N/A" if mean is None else f"{mean:.3f}"
            lines.append(
                f"| {vname} | {_score_text(v)} | "
                f"{v.get('n_strong', 0)} | {v.get('n_weak', 0)} | "
                f"{v.get('epitope_density', 0.0):.3f} | {mean_txt} | {_adj_text(v)} |"
            )
        lines.append("")

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("error"):
                lines.append(f"- **{vname}**: {v['error']}")
            if v.get("structural_error"):
                lines.append(f"- **{vname}** (structural): {v['structural_error']}")
        if any(v.get("error") or v.get("structural_error")
               for v in variants.values()):
            lines.append("")

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("positions_of_interest"):
                lines.append(f"### Epitope-changing positions ({vname})\n")
                lines.append("| VPos | DPos | Donor | Variant | DeltaEPI | Direction |")
                lines.append("|------|------|-------|---------|----------|-----------|")
                for p in v["positions_of_interest"][:15]:
                    lines.append(
                        f"| {p['position']} | {p.get('donor_position', 0)} | "
                        f"{p['donor_aa']} | {p['variant_aa']} | "
                        f"{p['delta']:+.3f} | {p['direction']} |"
                    )
                lines.append("")

        strong_thr = results.get("strong_threshold", 0.5)
        for vname in _sorted_variants(variants):
            peptides = variants[vname].get("per_peptide_hits", [])
            if not peptides:
                continue
            top = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top if p.get("score", 0) >= strong_thr]
            if not strong:
                continue
            lines.append(f"### Top predicted epitopes ({vname})\n")
            lines.append("| Peptide | Score | Core | Allele |")
            lines.append("|---------|-------|------|--------|")
            for p in strong:
                lines.append(
                    f"| `{p.get('peptide', '')}` | {p.get('score', 0):.3f} | "
                    f"`{p.get('core', '')}` | {p.get('allele', '')} |")
            lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def save_word_report(results: Dict, output_path: str,
                     title: str = "Immunogenicity Analysis"):
    """Save results as Word document."""
    try:
        from docx import Document
        from docx.shared import RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        print("  WARNING: python-docx not installed. Skipping Word report.",
              file=sys.stderr)
        return

    doc = Document()
    heading = doc.add_heading(title, 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_paragraph(f"Input: {results['input_fasta']}", style="Intense Quote")
    doc.add_paragraph(f"Backend: {results['mhc_backend']}")
    alleles = results["alleles"]
    doc.add_paragraph(f"Alleles: {', '.join(alleles[:5])}"
                      + (" ..." if len(alleles) > 5 else ""))
    doc.add_paragraph(f"Strong/weak thresholds: {results.get('strong_threshold', 0.5):.2f} / "
                      f"{results.get('weak_threshold', 0.3):.2f}")
    doc.add_paragraph(f"Score semantics: {results.get('score_semantics', 'percent of reference mean')}")
    doc.add_paragraph("* marks scores rewritten by structural adjustment.")

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        doc.add_heading(f"Chain {chain_type} ({len(variants)} variants)", level=1)
        doc.add_paragraph(f"Reference: {cdata.get('reference_source', 'none')} "
                          f"({cdata.get('reference_name') or '-'})")
        if cdata.get("structural_reference_mean") is not None:
            doc.add_paragraph(
                "Structural reference: "
                f"{cdata.get('structural_reference_name') or 'raw donor mean'}")

        table = doc.add_table(rows=1 + len(variants), cols=7)
        table.style = "Table Grid"
        headers = ["Variant", "Score", "Strong", "Weak", "Density", "Mean EPI", "StructAdj"]
        for i, h in enumerate(headers):
            cell = table.rows[0].cells[i]
            cell.text = h
            for p in cell.paragraphs:
                for run in p.runs:
                    run.bold = True
        for ri, vname in enumerate(_sorted_variants(variants), 1):
            v = variants[vname]
            mean = v.get("mean_epitope_load")
            vals = [vname, _score_text(v),
                    str(v.get("n_strong", 0)), str(v.get("n_weak", 0)),
                    f"{v.get('epitope_density', 0.0):.3f}",
                    "N/A" if mean is None else f"{mean:.3f}",
                    _adj_text(v)]
            for ci, val in enumerate(vals):
                table.rows[ri].cells[ci].text = val

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("error"):
                doc.add_paragraph(f"{vname}: ERROR - {v['error']}")
            if v.get("structural_error"):
                doc.add_paragraph(f"{vname} (structural): {v['structural_error']}")

        for vname in _sorted_variants(variants):
            v = variants[vname]
            if v.get("positions_of_interest"):
                doc.add_heading(f"Epitope-changing positions ({vname})", level=2)
                table2 = doc.add_table(rows=1 + min(15, len(v["positions_of_interest"])), cols=6)
                table2.style = "Table Grid"
                for i, h in enumerate(["VPos", "DPos", "Donor", "Variant", "DeltaEPI", "Direction"]):
                    cell = table2.rows[0].cells[i]
                    cell.text = h
                    for p in cell.paragraphs:
                        for run in p.runs:
                            run.bold = True
                for ri, p in enumerate(v["positions_of_interest"][:15], 1):
                    row = table2.rows[ri]
                    vals = [str(p['position']), str(p.get('donor_position', 0)),
                            p['donor_aa'], p['variant_aa'],
                            f"{p['delta']:+.3f}", p['direction']]
                    for ci, val in enumerate(vals):
                        row.cells[ci].text = val
                    delta_cell = row.cells[4]
                    for para in delta_cell.paragraphs:
                        for run in para.runs:
                            if p['delta'] < 0:
                                run.font.color.rgb = RGBColor(0, 128, 0)
                            elif p['delta'] > 0:
                                run.font.color.rgb = RGBColor(255, 0, 0)

        strong_thr = results.get("strong_threshold", 0.5)
        for vname in _sorted_variants(variants):
            peptides = variants[vname].get("per_peptide_hits", [])
            if not peptides:
                continue
            top = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top if p.get("score", 0) >= strong_thr]
            if not strong:
                continue
            doc.add_heading(f"Top predicted epitopes ({vname})", level=2)
            table3 = doc.add_table(rows=1 + len(strong), cols=4)
            table3.style = "Table Grid"
            for i, h in enumerate(["Peptide", "Score", "Core", "Allele"]):
                cell = table3.rows[0].cells[i]
                cell.text = h
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.bold = True
            for ri, p in enumerate(strong, 1):
                row = table3.rows[ri]
                vals = [p.get("peptide", ""), f"{p.get('score', 0):.3f}",
                        p.get("core", ""), p.get("allele", "")]
                for ci, val in enumerate(vals):
                    row.cells[ci].text = val

    doc.save(output_path)


def save_epitope_heatmap_csv(results: Dict, output_path: str):
    """Save per-position epitope scores as a chain-aware CSV heatmap."""
    columns = []
    for chain_type in sorted(results["chains"]):
        positions = set()
        for vdata in results["chains"][chain_type]["variants"].values():
            for key in vdata.get("epitope_scores", {}):
                positions.add(int(key))
        for p in sorted(positions):
            columns.append((chain_type, p))
    if not columns:
        return

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Chain", "Variant"] + [f"{c}:{p + 1}" for c, p in columns])
        for chain_type in sorted(results["chains"]):
            for vname in _sorted_variants(results["chains"][chain_type]["variants"]):
                vdata = results["chains"][chain_type]["variants"][vname]
                values = {int(k): v for k, v in vdata.get("epitope_scores", {}).items()}
                row = [chain_type, vname]
                for c, p in columns:
                    row.append(f"{values.get(p, 0.0):.3f}" if c == chain_type else "")
                writer.writerow(row)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Post-humanization immunogenicity analysis (MHC-II epitope + nativeness)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
        Examples:
          # Analyze pipeline output
          python3 tools/immunogenicity/immunogenicity_analyzer.py --outdir outputs/amg110_step3/step3

          # Analyze specific FASTA with all outputs
          python3 tools/immunogenicity/immunogenicity_analyzer.py --input variants.fasta --all-formats

          # Heuristic fallback (no external deps)
          python3 tools/immunogenicity/immunogenicity_analyzer.py --input variants.fasta --backend heuristic

          # With donor reference
          python3 tools/immunogenicity/immunogenicity_analyzer.py --input variants.fasta --donor donor.fasta

          # Check backend availability
          python3 tools/immunogenicity/immunogenicity_analyzer.py --check

        License:
          HLAIIPred: Apache-2.0 (commercial use OK)
        """)
    )
    parser.add_argument("--input", "-i", help="Path to variants FASTA file")
    parser.add_argument("--outdir", "-d",
                        help="Pipeline output directory (auto-discovers variants.fasta)")
    parser.add_argument("--donor", help="Donor antibody FASTA (for comparison)")
    parser.add_argument("--backend", "-b", default="auto",
                        choices=["auto", "hlaiipred", "heuristic"],
                        help="MHC-II prediction backend (default: auto)")
    parser.add_argument("--alleles", nargs="+", default=None,
                        help="HLA-II alleles; overrides --allele-panel")
    parser.add_argument("--allele-panel", default="9common",
                        choices=sorted(ALLELE_PANELS),
                        help="Named allele panel (default: 9common)")
    parser.add_argument("--pep-len", type=int, default=15,
                        help=f"Peptide window length, {MIN_PEP_LEN}-{MAX_PEP_LEN} (default: 15)")
    parser.add_argument("--strong-threshold", type=float,
                        default=DEFAULT_STRONG_THRESHOLD,
                        help="Strong presentation-probability cutoff (default: 0.5)")
    parser.add_argument("--weak-threshold", type=float,
                        default=DEFAULT_WEAK_THRESHOLD,
                        help="Weak presentation-probability cutoff (default: 0.3)")
    parser.add_argument("--reference-variant",
                        help="Explicit reference variant name (e.g. H_V0 or combo_V0_V0)")
    parser.add_argument("--output", "-o", default="outputs/immunogenicity",
                        help="Output directory (default: outputs/immunogenicity)")
    parser.add_argument("--json", action="store_true", help="Save JSON report")
    parser.add_argument("--markdown", action="store_true", help="Save Markdown report")
    parser.add_argument("--word", action="store_true", help="Save Word report")
    parser.add_argument("--csv", action="store_true", help="Save epitope heatmap CSV")
    parser.add_argument("--all-formats", action="store_true",
                        help="Save all formats (text + JSON + Markdown + Word + CSV)")
    parser.add_argument("--check", action="store_true",
                        help="Check backend availability and exit")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")

    # Structural risk adjustment options
    parser.add_argument("--structure-dir",
                        help="Directory containing variant PDB files for structural risk assessment")
    parser.add_argument("--backmutation-dir",
                        help="Directory containing backmutation CSVs (from pipeline output)")
    parser.add_argument("--structural-risk", action="store_true",
                        help="Enable structural risk adjustment (requires --structure-dir)")
    parser.add_argument("--numbering-json",
                        help="Per-variant Kabat labels for structural risk "
                             "(auto-discovered in --backmutation-dir when omitted)")
    parser.add_argument("--structural-method", default="relSASA",
                        choices=["relSASA", "binary"],
                        help="Structural risk adjustment method (default: relSASA)")

    args = parser.parse_args()

    # Check mode
    if args.check:
        print("Backend availability:")
        avail = check_backends()
        for name, ok in avail.items():
            status = "available" if ok else "NOT installed"
            print(f"  {name:<15} {status}")
        try:
            from immunogenicity.structural_risk import HAS_FREESASA
        except Exception:
            HAS_FREESASA = False
        print(f"  {'freesasa':<15} {'available' if HAS_FREESASA else 'NOT installed (structural step degraded)'}")
        try:
            import docx  # noqa: F401
            has_docx = True
        except Exception:
            has_docx = False
        print(f"  {'python-docx':<15} {'available' if has_docx else 'NOT installed (Word report skipped)'}")
        return 0

    # Validate input
    if not args.input and not args.outdir:
        parser.error("Either --input or --outdir must be specified")
    if args.structural_risk and not args.structure_dir:
        parser.error("--structural-risk requires --structure-dir")
    if not (MIN_PEP_LEN <= args.pep_len <= MAX_PEP_LEN):
        parser.error(f"--pep-len must be between {MIN_PEP_LEN} and {MAX_PEP_LEN}")
    if not (0.0 <= args.weak_threshold < args.strong_threshold <= 1.0):
        parser.error("require 0 <= weak-threshold < strong-threshold <= 1")

    alleles = args.alleles or ALLELE_PANELS[args.allele_panel]

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Run analysis
    try:
        results = analyze_immunogenicity(
            input_path=args.input,
            outdir=args.outdir,
            donor_path=args.donor,
            backend_name=args.backend,
            alleles=alleles,
            pep_len=args.pep_len,
            verbose=args.verbose,
            reference_variant=args.reference_variant,
            strong_threshold=args.strong_threshold,
            weak_threshold=args.weak_threshold,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    # Run structural risk adjustment if requested
    if args.structural_risk:
        try:
            from immunogenicity.structural_risk import run_structural_risk_assessment
            results = run_structural_risk_assessment(
                results,
                args.structure_dir,
                args.backmutation_dir,
                args.output,
                method=args.structural_method,
                verbose=args.verbose,
                numbering_json=args.numbering_json,
            )
            print("\n  Structural risk adjustment completed.")
        except ImportError as e:
            print(f"ERROR: Could not import structural_risk module: {e}",
                  file=sys.stderr)
            return 1
        except Exception as e:
            print(f"ERROR: Structural risk adjustment failed: {e}",
                  file=sys.stderr)
            return 1

    # Print text report
    title = f"Immunogenicity Analysis: {os.path.basename(results['input_fasta'])}"
    text_report = format_text_report(results, title)
    print(text_report)

    # Save reports
    stem = Path(results["input_fasta"]).stem
    if args.json or args.all_formats:
        json_path = os.path.join(args.output, f"{stem}_immunogenicity.json")
        save_json_report(results, json_path)
        print(f"\n  JSON report saved: {json_path}")

    if args.markdown or args.all_formats:
        md_path = os.path.join(args.output, f"{stem}_immunogenicity.md")
        save_markdown_report(results, md_path, title)
        print(f"  Markdown report saved: {md_path}")

    if args.word or args.all_formats:
        docx_path = os.path.join(args.output, f"{stem}_immunogenicity.docx")
        save_word_report(results, docx_path, title)
        print(f"  Word report saved: {docx_path}")

    if args.csv or args.all_formats:
        csv_path = os.path.join(args.output, f"{stem}_epitope_heatmap.csv")
        save_epitope_heatmap_csv(results, csv_path)
        print(f"  CSV heatmap saved: {csv_path}")

    errors = []
    total = 0
    for cdata in results.get("chains", {}).values():
        for vname, vdata in cdata.get("variants", {}).items():
            total += 1
            if vdata.get("error") or vdata.get("immunogenicity_score") is None:
                errors.append(f"{vname}: {vdata.get('error') or 'not evaluated'}")
    if errors:
        print(f"\n  WARNING: {len(errors)}/{total} variant(s) could not be evaluated:",
              file=sys.stderr)
        for msg in errors[:10]:
            print(f"    {msg}", file=sys.stderr)
        if len(errors) == total:
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
