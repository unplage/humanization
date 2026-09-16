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
    BackendResult, DEFAULT_ALLELES, STRONG_RANK_THRESHOLD,
    generate_peptides, auto_select_backend, check_backends, get_backend,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRONG_THRESHOLD = 0.5   # for heuristic/sigmoid scores
STRONG_RANK = 2.0        # %Rank < 2 = strong binder
WEAK_RANK = 10.0         # %Rank < 10 = weak binder


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

def group_by_chain(entries: List[Tuple[str, str, str]]) -> Dict[str, Dict[str, str]]:
    """Group sequences by chain type (H/L) and variant name.

    Returns: {chain_type: {variant_name: sequence}}
    """
    chains: Dict[str, Dict[str, str]] = defaultdict(dict)
    for hdr, seq, name in entries:
        # Determine chain type from header or name
        hdr_upper = hdr.upper()
        if "|H_" in hdr_upper or hdr_upper.startswith("H") or "VH" in hdr_upper:
            chain_type = "H"
        elif "|L_" in hdr_upper or hdr_upper.startswith("L") or "VL" in hdr_upper or "K" in hdr_upper:
            chain_type = "L"
        else:
            # Guess from length / content
            if len(seq) > 200 or seq[:4] in ("EVQL", "QVQL", "EVKL"):
                chain_type = "H"
            else:
                chain_type = "L"
        # Extract variant name from header if available
        parts = hdr.split("|")
        variant_name = parts[1].strip() if len(parts) > 1 else name
        chains[chain_type][variant_name] = seq
    return dict(chains)


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
    max_epitope: float = 0.0
    immunogenicity_score: float = 0.0   # 0-100
    per_peptide_hits: list = field(default_factory=list)
    positions_of_interest: List[Dict] = field(default_factory=list)


def compute_epitope_metrics(per_residue: Dict[int, float],
                            seq_len: int) -> Dict:
    """Compute aggregate metrics from per-residue epitope scores."""
    if not per_residue:
        return {"n_strong": 0, "n_weak": 0, "density": 0, "mean": 0, "max": 0}
    scores = list(per_residue.values())
    n_strong = sum(1 for s in scores if s >= STRONG_THRESHOLD)
    n_weak = sum(1 for s in scores if STRONG_RANK / 100 <= s < STRONG_THRESHOLD)
    # Epitope density: fraction of residues in at least one strong/weak window
    n_annotated = len(per_residue)
    density = n_annotated / max(1, seq_len)
    mean_load = sum(scores) / max(1, len(scores))
    max_score = max(scores)
    return {
        "n_strong": n_strong,
        "n_weak": n_weak,
        "density": round(density, 4),
        "mean": round(mean_load, 4),
        "max": round(max_score, 4),
    }


def compute_immunogenicity_score(metrics: Dict, donor_metrics: Optional[Dict] = None) -> float:
    """Compute 0-100 immunogenicity score.

    Lower is better (less immunogenic). Score = 100 * (variant_load / donor_load),
    capped at 100. If no donor, use absolute metrics.
    """
    variant_load = metrics["mean"]
    if donor_metrics and donor_metrics.get("mean", 0) > 0:
        ratio = variant_load / donor_metrics["mean"]
        return round(min(100.0, ratio * 100.0), 1)
    # No donor: use mean load directly (0-1 scale → 0-100)
    return round(min(100.0, variant_load * 100.0), 1)


def find_positions_of_interest(donor_seq: str, variant_seq: str,
                               donor_epi: Dict[int, float],
                               variant_epi: Dict[int, float]) -> List[Dict]:
    """Find positions where epitope score changes between donor and variant."""
    positions = []
    # Check positions that differ between donor and variant
    max_len = max(len(donor_seq), len(variant_seq))
    for i in range(max_len):
        d_aa = donor_seq[i] if i < len(donor_seq) else "-"
        v_aa = variant_seq[i] if i < len(variant_seq) else "-"
        if d_aa != v_aa:
            d_score = donor_epi.get(i, 0.0)
            v_score = variant_epi.get(i, 0.0)
            delta = v_score - d_score
            if abs(delta) > 0.05:  # meaningful change
                positions.append({
                    "position": i + 1,  # 1-based
                    "donor_aa": d_aa,
                    "variant_aa": v_aa,
                    "donor_epitope": round(d_score, 3),
                    "variant_epitope": round(v_score, 3),
                    "delta": round(delta, 3),
                    "direction": "reduced" if delta < 0 else "increased",
                })
    # Sort by absolute delta descending
    positions.sort(key=lambda x: abs(x["delta"]), reverse=True)
    return positions[:30]  # top 30


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
) -> Dict:
    """Run immunogenicity analysis on variants.

    Returns dict with per-chain, per-variant results and metrics.
    """
    if alleles is None:
        alleles = DEFAULT_ALLELES

    # Discover input FASTA
    fasta_path = input_path
    if not fasta_path and outdir:
        fasta_path = discover_variants(outdir)
    if not fasta_path or not os.path.isfile(fasta_path):
        raise FileNotFoundError(
            f"No variants.fasta found. Provide --input or --outdir containing variants.fasta")

    # Parse sequences
    entries = parse_fasta(fasta_path)
    if not entries:
        raise ValueError(f"No sequences found in {fasta_path}")

    # Group by chain
    chains = group_by_chain(entries)

    # Load donor if provided
    donor_chains = {}
    if donor_path and os.path.isfile(donor_path):
        donor_entries = parse_fasta(donor_path)
        donor_chains = group_by_chain(donor_entries)

    # Select MHC-II backend
    if backend_name == "auto":
        mhc_backend = auto_select_backend()
    else:
        mhc_backend = get_backend(backend_name)
        if not mhc_backend.is_available():
            raise RuntimeError(f"Backend '{backend_name}' is not available")

    if verbose:
        print(f"  Backend: {mhc_backend.describe()}")
        print(f"  Alleles: {', '.join(alleles)}")
        print(f"  Peptide length: {pep_len}")

    # Run analysis per chain
    all_results = {
        "input_fasta": fasta_path,
        "alleles": alleles,
        "pep_len": pep_len,
        "mhc_backend": mhc_backend.name,
        "chains": {},
    }

    for chain_type in sorted(chains.keys()):
        chain_variants = chains[chain_type]
        if verbose:
            print(f"\n  Chain {chain_type}: {len(chain_variants)} variants")

        # Run MHC-II backend on all variants in this chain
        mhc_results = mhc_backend.score(chain_variants, alleles, pep_len=pep_len)

        # Compute metrics per variant
        variant_analyses = {}
        donor_variant_name = None
        for vname, seq in chain_variants.items():
            # Identify donor (look for V0 or donor in name)
            is_donor = ("V0" in vname or "donor" in vname.lower()
                        or "graft" in vname.lower())
            if is_donor:
                donor_variant_name = vname

        # Use donor chain if available
        donor_seq = ""
        if chain_type in donor_chains:
            # Pick first (or V0) sequence
            for dn, ds in donor_chains[chain_type].items():
                donor_seq = ds
                break

        donor_epi = {}
        donor_metrics = None
        if donor_seq and donor_variant_name:
            donor_backend_res = mhc_results.get(donor_variant_name)
            if donor_backend_res and not donor_backend_res.error:
                donor_epi = donor_backend_res.per_residue
                donor_metrics = compute_epitope_metrics(donor_epi, len(donor_seq))
        elif donor_seq:
            # Score donor sequence separately
            donor_res = mhc_backend.score(
                {"donor": donor_seq}, alleles, pep_len=pep_len)
            donor_backend_res = donor_res.get("donor")
            if donor_backend_res and not donor_backend_res.error:
                donor_epi = donor_backend_res.per_residue
                donor_metrics = compute_epitope_metrics(donor_epi, len(donor_seq))

        # Score all variants (including V0 if not already done via donor)
        for vname, seq in chain_variants.items():
            backend_res = mhc_results.get(vname)
            if backend_res is None:
                variant_analyses[vname] = VariantAnalysis(
                    name=vname, chain_type=chain_type, sequence=seq)
                continue
            if backend_res.error:
                if verbose:
                    print(f"    {vname}: ERROR - {backend_res.error}")
                variant_analyses[vname] = VariantAnalysis(
                    name=vname, chain_type=chain_type, sequence=seq)
                continue

            epi = backend_res.per_residue
            metrics = compute_epitope_metrics(epi, len(seq))
            imm_score = compute_immunogenicity_score(metrics, donor_metrics)

            # Positions of interest
            poi = find_positions_of_interest(donor_seq, seq, donor_epi, epi) if donor_seq else []

            va = VariantAnalysis(
                name=vname,
                chain_type=chain_type,
                sequence=seq,
                epitope_scores=epi,
                n_strong=metrics["n_strong"],
                n_weak=metrics["n_weak"],
                n_total_windows=metrics["n_strong"] + metrics["n_weak"],
                epitope_density=metrics["density"],
                mean_epitope_load=metrics["mean"],
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
        }

    return all_results


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def format_text_report(results: Dict, title: str = "Immunogenicity Analysis") -> str:
    """Format results as text report (stdout)."""
    lines = []
    lines.append("=" * 80)
    lines.append(f"  {title}")
    lines.append("=" * 80)
    lines.append(f"  Input:    {results['input_fasta']}")
    lines.append(f"  Backend:  {results['mhc_backend']}")
    lines.append(f"  Alleles:  {', '.join(results['alleles'][:3])}..."
                 if len(results['alleles']) > 3 else f"  Alleles:  {', '.join(results['alleles'])}")
    lines.append(f"  Pep len:  {results['pep_len']}")
    lines.append("")

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        donor_m = cdata.get("donor_metrics")
        lines.append(f"  Chain {chain_type} ({len(variants)} variants)")
        lines.append("  " + "-" * 76)

        # Summary table
        lines.append(f"  {'Variant':<20} {'Score':>7} {'Strong':>7} {'Weak':>7} "
                     f"{'Density':>9} {'Mean EPI':>10}")
        lines.append("  " + "-" * 70)
        for vname in sorted(variants.keys()):
            v = variants[vname]
            lines.append(
                f"  {vname:<20} {v['immunogenicity_score']:>6.1f}% "
                f"{v['n_strong']:>7} {v['n_weak']:>7} "
                f"{v['epitope_density']:>8.3f} {v['mean_epitope_load']:>10.3f}"
            )
        lines.append("")

        # Donor comparison
        if donor_m:
            lines.append(f"  Donor reference: mean_epitope={donor_m['mean']:.3f}, "
                        f"strong={donor_m['n_strong']}, density={donor_m['density']:.3f}")
            lines.append("")

        # Top positions of interest
        for vname in sorted(variants.keys()):
            v = variants[vname]
            if v.get("positions_of_interest"):
                lines.append(f"  Top epitope-changing positions ({vname}):")
                lines.append(f"  {'Pos':>5} {'Donor':>6} {'Variant':>8} {'ΔEPI':>8} {'Dir':>10}")
                for p in v["positions_of_interest"][:10]:
                    lines.append(
                        f"  {p['position']:>5} {p['donor_aa']:>6} {p['variant_aa']:>8} "
                        f"{p['delta']:>+7.3f} {p['direction']:>10}"
                    )
                lines.append("")

        # Top strong-binding peptides per variant
        for vname in sorted(variants.keys()):
            v = variants[vname]
            peptides = v.get("per_peptide_hits", [])
            if not peptides:
                continue
            # Sort by score descending (higher = stronger binder)
            top_peptides = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top_peptides if p.get("score", 0) >= 0.5]
            if not strong:
                continue
            lines.append(f"  Top predicted epitopes ({vname}):")
            lines.append(f"  {'Peptide':<20} {'Score':>7} {'Rank':>7} {'Allele':<20}")
            lines.append("  " + "-" * 60)
            for p in strong:
                pep = p.get("peptide", "")
                score = p.get("score", 0)
                rank = p.get("rank", 0)
                allele = p.get("allele", "")
                lines.append(f"  {pep:<20} {score:>7.3f} {rank:>7.3f} {allele:<20}")
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
    lines.append(f"**Alleles:** {', '.join(results['alleles'][:5])}"
                 + (" ..." if len(results['alleles']) > 5 else ""))
    lines.append(f"**Peptide length:** {results['pep_len']}")
    lines.append("")

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        lines.append(f"## Chain {chain_type} ({len(variants)} variants)\n")

        # Summary table
        lines.append("| Variant | Score | Strong | Weak | Density | Mean EPI |")
        lines.append("|---------|-------|--------|------|---------|----------|")
        for vname in sorted(variants.keys()):
            v = variants[vname]
            lines.append(
                f"| {vname} | {v['immunogenicity_score']:.1f}% | "
                f"{v['n_strong']} | {v['n_weak']} | "
                f"{v['epitope_density']:.3f} | {v['mean_epitope_load']:.3f} |"
            )
        lines.append("")

        # Top positions of interest
        for vname in sorted(variants.keys()):
            v = variants[vname]
            if v.get("positions_of_interest"):
                lines.append(f"### Epitope-changing positions ({vname})\n")
                lines.append("| Pos | Donor | Variant | ΔEPI | Direction |")
                lines.append("|-----|-------|---------|------|-----------|")
                for p in v["positions_of_interest"][:15]:
                    lines.append(
                        f"| {p['position']} | {p['donor_aa']} | {p['variant_aa']} | "
                        f"{p['delta']:+.3f} | {p['direction']} |"
                    )
                lines.append("")

        # Top strong-binding peptides per variant
        for vname in sorted(variants.keys()):
            v = variants[vname]
            peptides = v.get("per_peptide_hits", [])
            if not peptides:
                continue
            top_peptides = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top_peptides if p.get("score", 0) >= 0.5]
            if not strong:
                continue
            lines.append(f"### Top predicted epitopes ({vname})\n")
            lines.append("| Peptide | Score | Rank | Allele |")
            lines.append("|---------|-------|------|--------|")
            for p in strong:
                pep = p.get("peptide", "")
                score = p.get("score", 0)
                rank = p.get("rank", 0)
                allele = p.get("allele", "")
                lines.append(f"| `{pep}` | {score:.3f} | {rank:.3f} | {allele} |")
            lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def save_word_report(results: Dict, output_path: str,
                     title: str = "Immunogenicity Analysis"):
    """Save results as Word document."""
    try:
        from docx import Document
        from docx.shared import Inches, Pt, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.enum.table import WD_TABLE_ALIGNMENT
    except ImportError:
        print("  WARNING: python-docx not installed. Skipping Word report.",
              file=sys.stderr)
        return

    doc = Document()
    heading = doc.add_heading(title, 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_paragraph(f"Input: {results['input_fasta']}", style="Intense Quote")
    doc.add_paragraph(f"Backend: {results['mhc_backend']}")
    doc.add_paragraph(f"Alleles: {', '.join(results['alleles'][:5])}"
                     + (" ..." if len(results['alleles']) > 5 else ""))

    for chain_type, cdata in results["chains"].items():
        variants = cdata["variants"]
        doc.add_heading(f"Chain {chain_type} ({len(variants)} variants)", level=1)

        # Summary table
        table = doc.add_table(rows=1 + len(variants), cols=6)
        table.style = "Table Grid"
        headers = ["Variant", "Score", "Strong", "Weak", "Density", "Mean EPI"]
        for i, h in enumerate(headers):
            cell = table.rows[0].cells[i]
            cell.text = h
            for p in cell.paragraphs:
                for run in p.runs:
                    run.bold = True
        for ri, vname in enumerate(sorted(variants.keys()), 1):
            v = variants[vname]
            row = table.rows[ri]
            vals = [vname, f"{v['immunogenicity_score']:.1f}%",
                    str(v['n_strong']), str(v['n_weak']),
                    f"{v['epitope_density']:.3f}", f"{v['mean_epitope_load']:.3f}"]
            for ci, val in enumerate(vals):
                row.cells[ci].text = val

        # Top positions of interest
        for vname in sorted(variants.keys()):
            v = variants[vname]
            if v.get("positions_of_interest"):
                doc.add_heading(f"Epitope-changing positions ({vname})", level=2)
                table2 = doc.add_table(rows=1 + min(15, len(v["positions_of_interest"])), cols=5)
                table2.style = "Table Grid"
                for i, h in enumerate(["Pos", "Donor", "Variant", "ΔEPI", "Direction"]):
                    cell = table2.rows[0].cells[i]
                    cell.text = h
                    for p in cell.paragraphs:
                        for run in p.runs:
                            run.bold = True
                for ri, p in enumerate(v["positions_of_interest"][:15], 1):
                    row = table2.rows[ri]
                    vals = [str(p['position']), p['donor_aa'], p['variant_aa'],
                            f"{p['delta']:+.3f}", p['direction']]
                    for ci, val in enumerate(vals):
                        row.cells[ci].text = val
                    delta_cell = row.cells[3]
                    for para in delta_cell.paragraphs:
                        for run in para.runs:
                            if p['delta'] < 0:
                                run.font.color.rgb = RGBColor(0, 128, 0)
                            elif p['delta'] > 0:
                                run.font.color.rgb = RGBColor(255, 0, 0)

        # Top predicted epitopes per variant
        for vname in sorted(variants.keys()):
            v = variants[vname]
            peptides = v.get("per_peptide_hits", [])
            if not peptides:
                continue
            top_peptides = sorted(peptides, key=lambda x: x.get("score", 0), reverse=True)[:10]
            strong = [p for p in top_peptides if p.get("score", 0) >= 0.5]
            if not strong:
                continue
            doc.add_heading(f"Top predicted epitopes ({vname})", level=2)
            table3 = doc.add_table(rows=1 + len(strong), cols=4)
            table3.style = "Table Grid"
            for i, h in enumerate(["Peptide", "Score", "Rank", "Allele"]):
                cell = table3.rows[0].cells[i]
                cell.text = h
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.bold = True
            for ri, p in enumerate(strong, 1):
                row = table3.rows[ri]
                vals = [p.get("peptide", ""), f"{p.get('score', 0):.3f}",
                        f"{p.get('rank', 0):.3f}", p.get("allele", "")]
                for ci, val in enumerate(vals):
                    row.cells[ci].text = val

    doc.save(output_path)


def save_epitope_heatmap_csv(results: Dict, output_path: str):
    """Save per-position epitope scores as CSV heatmap."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        # Collect all positions
        all_positions = set()
        for chain_type, cdata in results["chains"].items():
            for vname, vdata in cdata["variants"].items():
                all_positions.update(vdata.get("epitope_scores", {}).keys())
        if not all_positions:
            return
        sorted_pos = sorted(all_positions, key=lambda x: int(x) if isinstance(x, str) else x)
        # Header
        header = ["Chain", "Variant"] + [str(int(p) + 1) if isinstance(p, str) else str(p + 1) for p in sorted_pos]
        writer.writerow(header)
        # Data
        for chain_type, cdata in sorted(results["chains"].items()):
            for vname, vdata in sorted(cdata["variants"].items()):
                epi = vdata.get("epitope_scores", {})
                # Convert string keys to int for lookup
                epi_int = {int(k): v for k, v in epi.items()}
                row = [chain_type, vname] + [f"{epi_int.get(int(p) if isinstance(p, str) else p, 0.0):.3f}" for p in sorted_pos]
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
                        help="HLA-II alleles (default: 9 common DRB1)")
    parser.add_argument("--pep-len", type=int, default=15,
                        help="Peptide window length (default: 15)")
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
    parser.add_argument("--structural-method", default="relSASA",
                        choices=["relSASA", "binary"],
                        help="Structural risk adjustment method (default: relSASA)")

    args = parser.parse_args()

    # Check mode
    if args.check:
        print("Backend availability:")
        avail = check_backends()
        for name, ok in avail.items():
            status = "✓ available" if ok else "✗ NOT installed"
            print(f"  {name:<15} {status}")
        return 0

    # Validate input
    if not args.input and not args.outdir:
        parser.error("Either --input or --outdir must be specified")

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Run analysis
    try:
        results = analyze_immunogenicity(
            input_path=args.input,
            outdir=args.outdir,
            donor_path=args.donor,
            backend_name=args.backend,
            alleles=args.alleles,
            pep_len=args.pep_len,
            verbose=args.verbose,
        )
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    # Run structural risk adjustment if requested
    if args.structural_risk and args.structure_dir:
        try:
            from immunogenicity.structural_risk import run_structural_risk_assessment
            results = run_structural_risk_assessment(
                results,
                args.structure_dir,
                args.backmutation_dir,
                args.output,
                method=args.structural_method,
                verbose=args.verbose,
            )
            print("\n  Structural risk adjustment completed.")
        except ImportError as e:
            print(f"  WARNING: Could not import structural_risk module: {e}",
                  file=sys.stderr)
        except Exception as e:
            print(f"  WARNING: Structural risk adjustment failed: {e}",
                  file=sys.stderr)

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

    return 0


if __name__ == "__main__":
    sys.exit(main())
