#!/usr/bin/env python3
"""CDR Homology Analysis — query clinical antibody database.

Usage:
    # Single query (VH + VL)
    python3 tools/cdr_homology/query.py \
        --vh QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYWMHWVRQAPGQGLEWIG... \
        --vl DIQMTQSPSSLSASVGDRVTITCRASQGISSYLAWYQQKPGKAPKLLI... \
        --top 5 --output outputs/cdr_homology/

    # From FASTA (pairs of VH/VL by header suffix)
    python3 tools/cdr_homology/query.py \
        --input query.fasta --top 5 --output outputs/cdr_homology/
"""

import argparse
import csv
import json
import os
import sys
import textwrap
from collections import Counter, defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# ANARCI CDR extraction (self-contained, no pipeline imports)
# ---------------------------------------------------------------------------

_IMGT_CDR_RANGES = {
    "CDR1": (27, 38),
    "CDR2": (56, 65),
    "CDR3": (105, 117),
}


def _extract_cdrs_from_posmap(posmap: dict, chain_type: str) -> dict:
    """Extract CDR1/CDR2/CDR3 from an IMGT position map."""
    cdrs = {}
    for cdr_name, (start, end) in _IMGT_CDR_RANGES.items():
        seq = ""
        for pos in range(start, end + 1):
            aa = posmap.get(f"{chain_type}{pos}", None)
            if aa and aa != "-":
                seq += aa
            for letter in "ABCDEFGHIJK":
                ins_aa = posmap.get(f"{chain_type}{pos}{letter}", None)
                if ins_aa:
                    seq += ins_aa
        cdrs[cdr_name] = seq
    return cdrs


def _number_anarci(sequence: str, chain_type: str):
    """Number a sequence using ANARCI IMGT, return posmap."""
    try:
        from anarci import anarci
    except ImportError:
        return None

    seq = sequence.upper().strip()
    try:
        result = anarci([("q", seq)], scheme="imgt", output=False)
    except Exception:
        return None

    try:
        numbering_list = result[0][0][0][0]
    except (IndexError, TypeError):
        return None

    if not numbering_list:
        return None

    posmap = {}
    for (pos_num, ins_flag), aa in numbering_list:
        if aa == "-":
            continue
        label = f"{chain_type}{pos_num}{ins_flag}" if ins_flag and ins_flag != " " else f"{chain_type}{pos_num}"
        posmap[label] = aa.upper()
    return posmap


def extract_cdrs(sequence: str, chain_type: str):
    """Full CDR extraction: ANARCI numbering + posmap extraction.

    Returns dict {CDR1, CDR2, CDR3} or None on failure.
    """
    posmap = _number_anarci(sequence, chain_type)
    if not posmap:
        return None
    return _extract_cdrs_from_posmap(posmap, chain_type)


# ---------------------------------------------------------------------------
# FASTA parsing (self-contained)
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> list:
    """Parse FASTA file into [(header, sequence), ...]."""
    records = []
    name, lines = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(lines)))
                name = line[1:].strip()
                lines = []
            else:
                lines.append(line)
    if name is not None:
        records.append((name, "".join(lines)))
    return records


# ---------------------------------------------------------------------------
# Alignment engine (self-contained, uses Biopython)
# ---------------------------------------------------------------------------

_aligner = None


def _get_aligner(mode="global"):
    """Lazily build a Biopython PairwiseAligner with BLOSUM62."""
    global _aligner
    from Bio.Align import PairwiseAligner, substitution_matrices

    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.mode = mode
    return a


def align_sequences(query: str, subject: str, is_local: bool = False) -> dict:
    """Align two sequences, return identity / score / aligned strings.

    Uses BioPython PairwiseAligner with BLOSUM62.
    All CDRs use global alignment (Needleman-Wunsch) to ensure full-length comparison.
    Identity is calculated over the aligned length (including gaps).

    Note: BioPython aligner.align() takes (target, query) order.
    """
    if not query or not subject:
        return {"identity": 0.0, "score": 0.0, "matches": 0, "aligned_len": 0,
                "query_aligned": "", "subject_aligned": ""}

    # Always use global alignment for proper full-length comparison
    aligner = _get_aligner("global")
    try:
        # BioPython: aligner.align(target, query)
        alignments = aligner.align(subject, query)
    except Exception:
        return {"identity": 0.0, "score": 0.0, "matches": 0, "aligned_len": 0,
                "query_aligned": "", "subject_aligned": ""}

    if not alignments:
        return {"identity": 0.0, "score": 0.0, "matches": 0, "aligned_len": 0,
                "query_aligned": "", "subject_aligned": ""}

    best = alignments[0]
    coords = best.coordinates  # shape (2, num_segments+1)
    seqs = best.sequences
    # seqs[0] = target = subject, seqs[1] = query

    # Build aligned strings from segment coordinates
    q_parts = []
    s_parts = []
    for i in range(coords.shape[1] - 1):
        t_start = int(coords[0, i])
        t_end = int(coords[0, i + 1])
        q_start = int(coords[1, i])
        q_end = int(coords[1, i + 1])

        if t_start == t_end:
            # Gap in target (subject) — query has insertion
            q_parts.append(seqs[1][q_start:q_end])
            s_parts.append("-" * (q_end - q_start))
        elif q_start == q_end:
            # Gap in query — subject has deletion
            q_parts.append("-" * (t_end - t_start))
            s_parts.append(seqs[0][t_start:t_end])
        else:
            # Match/mismatch block
            q_parts.append(seqs[1][q_start:q_end])
            s_parts.append(seqs[0][t_start:t_end])

    aln_q = "".join(q_parts)  # query aligned
    aln_s = "".join(s_parts)  # subject aligned

    # Count matches (exclude gaps)
    matches = sum(1 for a, b in zip(aln_q, aln_s) if a == b and a not in ("-", "."))
    aligned_len = max(len(aln_q), 1)
    identity = matches / aligned_len

    return {
        "identity": round(identity, 4),
        "score": round(best.score, 1),
        "matches": matches,
        "aligned_len": aligned_len,
        "query_aligned": aln_q,
        "subject_aligned": aln_s,
    }


# ---------------------------------------------------------------------------
# Database loading
# ---------------------------------------------------------------------------

def load_database(db_path: str) -> list:
    """Load preprocessed CDR JSON database."""
    with open(db_path) as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# CDR homology search
# ---------------------------------------------------------------------------

def search_cdr(query_cdr: str, database: list, chain_key: str, cdr_name: str,
               is_cdr3: bool, top_n: int = 5) -> list:
    """Find top N matches for a single CDR from the database."""
    mode = "local" if is_cdr3 else "global"
    results = []

    for entry in database:
        for label_key in (chain_key,):
            chain_data = entry.get("chains", {}).get(label_key)
            if not chain_data:
                continue
            subject_cdr = chain_data.get(cdr_name, "")
            if not subject_cdr:
                continue

            aln = align_sequences(query_cdr, subject_cdr, is_local=is_cdr3)
            results.append({
                "therapeutic": entry["name"],
                "identity": aln["identity"],
                "score": aln["score"],
                "matches": aln["matches"],
                "aligned_len": aln["aligned_len"],
                "subject_cdr": subject_cdr,
                "query_aligned": aln["query_aligned"],
                "subject_aligned": aln["subject_aligned"],
                "format": entry.get("format", ""),
                "genetics": entry.get("genetics", ""),
                "target": entry.get("target", ""),
                "highest_trial": entry.get("highest_trial", ""),
            })

        # also check bispecific chains
        bispec = entry.get("bispecific", {})
        if chain_key in bispec:
            subject_cdr = bispec[chain_key].get(cdr_name, "")
            if not subject_cdr:
                continue
            aln = align_sequences(query_cdr, subject_cdr, is_local=is_cdr3)
            results.append({
                "therapeutic": f"{entry['name']} (bispec)",
                "identity": aln["identity"],
                "score": aln["score"],
                "matches": aln["matches"],
                "aligned_len": aln["aligned_len"],
                "subject_cdr": subject_cdr,
                "query_aligned": aln["query_aligned"],
                "subject_aligned": aln["subject_aligned"],
                "format": entry.get("format", ""),
                "genetics": entry.get("genetics", ""),
                "target": entry.get("target", ""),
                "highest_trial": entry.get("highest_trial", ""),
            })

    results.sort(key=lambda x: (-x["identity"], -x["score"]))
    return results[:top_n]


def run_query(queries: list, database: list, top_n: int = 5) -> list:
    """Run CDR homology search for all queries.

    queries: [(name, vh_seq, vl_seq), ...]
    Returns per-query result dicts.
    """
    all_results = []
    for name, vh_seq, vl_seq in queries:
        print(f"  Processing {name} ...")
        result = {"name": name, "vh_seq": vh_seq, "vl_seq": vl_seq, "cdrs": {}, "matches": {}}

        # Extract query CDRs
        for chain_key, seq in [("VH", vh_seq), ("VL", vl_seq)]:
            if not seq:
                continue
            cdrs = extract_cdrs(seq, chain_key[1])  # H or L
            if cdrs:
                result["cdrs"][chain_key] = cdrs
                for cdr_name, cdr_seq in cdrs.items():
                    full_key = f"{chain_key}_{cdr_name}"
                    is_cdr3 = "CDR3" in cdr_name
                    top_matches = search_cdr(cdr_seq, database, chain_key, cdr_name, is_cdr3, top_n)
                    result["matches"][full_key] = {
                        "query_cdr": cdr_seq,
                        "query_len": len(cdr_seq),
                        "top_match": top_matches[0] if top_matches else None,
                        "top_matches": top_matches,
                    }
        all_results.append(result)
    return all_results


# ---------------------------------------------------------------------------
# Report generation (Markdown + CSV)
# ---------------------------------------------------------------------------

def generate_report(all_results: list, database_stats: dict, output_dir: str, top_n: int):
    """Generate Markdown summary + per-query CSV detail files."""
    os.makedirs(output_dir, exist_ok=True)

    # ---- Summary Markdown ----
    lines = []
    lines.append("# CDR Homology Analysis Report\n")
    lines.append(f"**Queries analysed:** {len(all_results)}\n")
    lines.append(f"**Database size:** {database_stats.get('total_entries', '?')} clinical antibodies\n")
    lines.append("")

    # Summary table
    lines.append("## Summary Table\n")
    header = "| # | Antibody | VH_CDR3 best | VH_CDR3 id | VL_CDR3 best | VL_CDR3 id | Mean id |"
    sep    = "|---|----------|-------------|-----------|-------------|-----------|---------|"
    lines.append(header)
    lines.append(sep)
    for i, r in enumerate(all_results, 1):
        vh3 = r["matches"].get("VH_CDR3", {}).get("top_match")
        vl3 = r["matches"].get("VL_CDR3", {}).get("top_match")
        vh3_name = (vh3["therapeutic"][:25] + "..") if vh3 and len(vh3["therapeutic"]) > 25 else (vh3["therapeutic"] if vh3 else "N/A")
        vh3_id = f"{vh3['identity']:.1%}" if vh3 else "N/A"
        vl3_name = (vl3["therapeutic"][:25] + "..") if vl3 and len(vl3["therapeutic"]) > 25 else (vl3["therapeutic"] if vl3 else "N/A")
        vl3_id = f"{vl3['identity']:.1%}" if vl3 else "N/A"
        # overall mean
        ids = [v["top_match"]["identity"] for v in r["matches"].values() if v.get("top_match")]
        mean_id = f"{sum(ids)/len(ids):.1%}" if ids else "N/A"
        lines.append(f"| {i} | {r['name'][:30]} | {vh3_name} | {vh3_id} | {vl3_name} | {vl3_id} | {mean_id} |")
    lines.append("")

    # Per-query detail
    for i, r in enumerate(all_results, 1):
        lines.append(f"---\n## {i}. {r['name']}\n")

        # CDR table for this query
        lines.append("### CDR Sequences (IMGT numbering)\n")
        lines.append("| Region | Query CDR | Len | Top Match | Identity | Score |")
        lines.append("|--------|-----------|-----|-----------|----------|-------|")
        for full_key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(full_key)
            if not m:
                lines.append(f"| {full_key} | - | - | - | - | - |")
                continue
            q = m["query_cdr"]
            qlen = m["query_len"]
            tm = m["top_match"]
            if tm:
                tname = tm["therapeutic"][:30]
                tid = f"{tm['identity']:.1%}"
                tscore = f"{tm['score']:.0f}"
            else:
                tname = tid = tscore = "N/A"
            lines.append(f"| {full_key} | `{q}` | {qlen} | {tname} | {tid} | {tscore} |")
        lines.append("")

        # Top 5 per CDR
        for full_key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(full_key)
            if not m:
                continue
            lines.append(f"### {full_key} — Top {top_n} Matches\n")
            lines.append("| Rank | Therapeutic | Identity | Score | Subject CDR | Alignment |")
            lines.append("|------|-------------|----------|-------|-------------|-----------|")
            for j, hit in enumerate(m["top_matches"][:top_n], 1):
                aln_str = ""
                if hit["query_aligned"] and hit["subject_aligned"]:
                    aln_str = f"`{hit['query_aligned']}` / `{hit['subject_aligned']}`"
                lines.append(
                    f"| {j} | {hit['therapeutic'][:30]} | {hit['identity']:.1%} "
                    f"| {hit['score']:.0f} | `{hit['subject_cdr']}` | {aln_str} |"
                )
            lines.append("")

    # ---- Appendix: Database Statistics ----
    lines.append("---\n")
    lines.append("# Appendix: Clinical Antibody Database Summary (TheraSAbDab)\n")

    # Format counts
    lines.append("## A. Antibody Format Distribution\n")
    lines.append("| Format | Count |")
    lines.append("|--------|-------|")
    for fmt, cnt in sorted(database_stats.get("format_counts", {}).items(), key=lambda x: -x[1]):
        lines.append(f"| {fmt} | {cnt} |")
    lines.append("")

    # Genetics
    lines.append("## B. Genetics / Humanization Strategy\n")
    lines.append("| Strategy | Count |")
    lines.append("|----------|-------|")
    for g, cnt in sorted(database_stats.get("genetics_counts", {}).items(), key=lambda x: -x[1]):
        lines.append(f"| {g} | {cnt} |")
    lines.append("")

    # Light chain type
    lines.append("## C. Light Chain Type\n")
    lines.append("| Type | Count |")
    lines.append("|------|-------|")
    for v, cnt in sorted(database_stats.get("vd_lc_counts", {}).items(), key=lambda x: -x[1]):
        lines.append(f"| {v} | {cnt} |")
    lines.append("")

    # Clinical trial status
    lines.append("## D. Highest Clinical Trial Phase\n")
    lines.append("| Phase | Count |")
    lines.append("|-------|-------|")
    for t, cnt in sorted(database_stats.get("trial_counts", {}).items(), key=lambda x: -x[1]):
        lines.append(f"| {t} | {cnt} |")
    lines.append("")

    # CDR length distribution
    lines.append("## E. CDR Length Distribution (IMGT numbering)\n")
    lines.append("| CDR | Min | Max | Mean | Median | N |")
    lines.append("|-----|-----|-----|------|--------|---|")
    for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
        s = database_stats.get("cdr_length_summary", {}).get(key)
        if s:
            lines.append(f"| {key} | {s['min']} | {s['max']} | {s['mean']} | {s['median']} | {s['count']} |")
        else:
            lines.append(f"| {key} | - | - | - | - | 0 |")
    lines.append("")

    # Top targets
    lines.append("## F. Top 30 Therapeutic Targets\n")
    lines.append("| Target | Count |")
    lines.append("|--------|-------|")
    tc = database_stats.get("target_counts", {})
    for t, cnt in sorted(tc.items(), key=lambda x: -x[1])[:30]:
        lines.append(f"| {t} | {cnt} |")
    lines.append("")

    # Top 30 drugs
    lines.append("## G. Top 30 Drugs by Entry Count\n")
    drug_counts = Counter(e["name"] for e in (load_database_stats_from_results(all_results) or []))
    lines.append("| Drug | Entries |")
    lines.append("|------|---------|")
    for d, cnt in drug_counts.most_common(30):
        lines.append(f"| {d} | {cnt} |")
    lines.append("")

    # Write markdown
    md_path = os.path.join(output_dir, "cdr_homology_report.md")
    with open(md_path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Report written: {md_path}")

    # ---- CSV detail files ----
    csv_dir = os.path.join(output_dir, "detail")
    os.makedirs(csv_dir, exist_ok=True)
    for r in all_results:
        safe_name = r["name"].replace("/", "_").replace(" ", "_")[:60]
        csv_path = os.path.join(csv_dir, f"{safe_name}.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([
                "CDR", "Query CDR", "Query Len", "Rank",
                "Therapeutic", "Identity", "Score", "Subject CDR",
                "Format", "Genetics", "Target", "Highest Trial",
            ])
            for full_key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
                m = r["matches"].get(full_key)
                if not m:
                    continue
                for j, hit in enumerate(m["top_matches"], 1):
                    writer.writerow([
                        full_key, m["query_cdr"], m["query_len"], j,
                        hit["therapeutic"], f"{hit['identity']:.4f}", hit["score"],
                        hit["subject_cdr"], hit["format"], hit["genetics"],
                        hit["target"], hit["highest_trial"],
                    ])
        print(f"  CSV written: {csv_path}")

    return md_path


def generate_docx_report(all_results: list, database_stats: dict, output_dir: str, top_n: int):
    """Generate Word (.docx) report with summary + appendix."""
    from docx import Document
    from docx.shared import Inches, Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml.ns import qn

    def set_cell_shading(cell, color):
        shading = cell._element.get_or_add_tcPr().makeelement(qn("w:shd"), {
            qn("w:val"): "clear",
            qn("w:color"): "auto",
            qn("w:fill"): color,
        })
        cell._element.get_or_add_tcPr().append(shading)

    def add_table_row(table, values, bold=False, header=False):
        row = table.add_row()
        for i, v in enumerate(values):
            cell = row.cells[i]
            cell.text = str(v)
            for p in cell.paragraphs:
                p.style = doc.styles["Normal"]
                for run in p.runs:
                    run.font.size = Pt(8)
                    if bold:
                        run.bold = True
            if header:
                set_cell_shading(cell, "2F5496")
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.font.color.rgb = RGBColor(255, 255, 255)
        return row

    doc = Document()

    # ---- Title ----
    title = doc.add_heading("CDR Homology Analysis Report", level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER

    doc.add_paragraph(f"Queries analysed: {len(all_results)}")
    doc.add_paragraph(f"Database size: {database_stats.get('total_entries', '?')} clinical antibodies")
    doc.add_paragraph()

    # ---- Summary Table ----
    doc.add_heading("Summary Table", level=1)
    tbl = doc.add_table(rows=1, cols=7)
    tbl.style = "Table Grid"
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    hdr_cells = tbl.rows[0].cells
    for i, h in enumerate(["#", "Antibody", "VH_CDR3 Best", "VH_CDR3 Id", "VL_CDR3 Best", "VL_CDR3 Id", "Mean Id"]):
        hdr_cells[i].text = h
        set_cell_shading(hdr_cells[i], "2F5496")
        for p in hdr_cells[i].paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)

    for i, r in enumerate(all_results, 1):
        vh3 = r["matches"].get("VH_CDR3", {}).get("top_match")
        vl3 = r["matches"].get("VL_CDR3", {}).get("top_match")
        vh3_name = (vh3["therapeutic"][:25] + "..") if vh3 and len(vh3["therapeutic"]) > 25 else (vh3["therapeutic"] if vh3 else "N/A")
        vh3_id = f"{vh3['identity']:.1%}" if vh3 else "N/A"
        vl3_name = (vl3["therapeutic"][:25] + "..") if vl3 and len(vl3["therapeutic"]) > 25 else (vl3["therapeutic"] if vl3 else "N/A")
        vl3_id = f"{vl3['identity']:.1%}" if vl3 else "N/A"
        ids = [v["top_match"]["identity"] for v in r["matches"].values() if v.get("top_match")]
        mean_id = f"{sum(ids)/len(ids):.1%}" if ids else "N/A"
        add_table_row(tbl, [i, r["name"][:30], vh3_name, vh3_id, vl3_name, vl3_id, mean_id])
    doc.add_paragraph()

    # ---- Per-query detail ----
    for i, r in enumerate(all_results, 1):
        doc.add_heading(f"{i}. {r['name']}", level=1)

        # CDR summary table
        doc.add_heading("CDR Sequences (IMGT numbering)", level=2)
        tbl2 = doc.add_table(rows=1, cols=6)
        tbl2.style = "Table Grid"
        for j, h in enumerate(["Region", "Query CDR", "Len", "Top Match", "Identity", "Score"]):
            tbl2.rows[0].cells[j].text = h
            set_cell_shading(tbl2.rows[0].cells[j], "2F5496")
            for p in tbl2.rows[0].cells[j].paragraphs:
                for run in p.runs:
                    run.bold = True
                    run.font.color.rgb = RGBColor(255, 255, 255)
                    run.font.size = Pt(8)

        for full_key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(full_key)
            if not m:
                add_table_row(tbl2, [full_key, "-", "-", "-", "-", "-"])
                continue
            tm = m["top_match"]
            if tm:
                add_table_row(tbl2, [
                    full_key, m["query_cdr"], m["query_len"],
                    tm["therapeutic"][:30], f"{tm['identity']:.1%}", f"{tm['score']:.0f}"
                ])
            else:
                add_table_row(tbl2, [full_key, m["query_cdr"], m["query_len"], "N/A", "N/A", "N/A"])
        doc.add_paragraph()

        # Top N per CDR
        for full_key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(full_key)
            if not m:
                continue
            doc.add_heading(f"{full_key} — Top {top_n} Matches", level=2)
            tbl3 = doc.add_table(rows=1, cols=6)
            tbl3.style = "Table Grid"
            for j, h in enumerate(["Rank", "Therapeutic", "Identity", "Score", "Subject CDR", "Alignment"]):
                tbl3.rows[0].cells[j].text = h
                set_cell_shading(tbl3.rows[0].cells[j], "2F5496")
                for p in tbl3.rows[0].cells[j].paragraphs:
                    for run in p.runs:
                        run.bold = True
                        run.font.color.rgb = RGBColor(255, 255, 255)
                        run.font.size = Pt(8)

            for j, hit in enumerate(m["top_matches"][:top_n], 1):
                aln_str = ""
                if hit["query_aligned"] and hit["subject_aligned"]:
                    aln_str = f"{hit['query_aligned']} / {hit['subject_aligned']}"
                add_table_row(tbl3, [
                    j, hit["therapeutic"][:30], f"{hit['identity']:.1%}",
                    f"{hit['score']:.0f}", hit["subject_cdr"], aln_str
                ])
            doc.add_paragraph()

    # ---- Appendix ----
    doc.add_page_break()
    doc.add_heading("Appendix: Clinical Antibody Database Summary (TheraSAbDab)", level=1)

    # A. Format
    doc.add_heading("A. Antibody Format Distribution", level=2)
    tbl_a = doc.add_table(rows=1, cols=2)
    tbl_a.style = "Table Grid"
    tbl_a.rows[0].cells[0].text = "Format"
    tbl_a.rows[0].cells[1].text = "Count"
    for cell in tbl_a.rows[0].cells:
        set_cell_shading(cell, "2F5496")
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    for fmt, cnt in sorted(database_stats.get("format_counts", {}).items(), key=lambda x: -x[1]):
        add_table_row(tbl_a, [fmt, cnt])
    doc.add_paragraph()

    # B. Genetics
    doc.add_heading("B. Genetics / Humanization Strategy", level=2)
    tbl_b = doc.add_table(rows=1, cols=2)
    tbl_b.style = "Table Grid"
    tbl_b.rows[0].cells[0].text = "Strategy"
    tbl_b.rows[0].cells[1].text = "Count"
    for cell in tbl_b.rows[0].cells:
        set_cell_shading(cell, "2F5496")
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    for g, cnt in sorted(database_stats.get("genetics_counts", {}).items(), key=lambda x: -x[1]):
        add_table_row(tbl_b, [g, cnt])
    doc.add_paragraph()

    # C. Light chain
    doc.add_heading("C. Light Chain Type", level=2)
    tbl_c = doc.add_table(rows=1, cols=2)
    tbl_c.style = "Table Grid"
    tbl_c.rows[0].cells[0].text = "Type"
    tbl_c.rows[0].cells[1].text = "Count"
    for cell in tbl_c.rows[0].cells:
        set_cell_shading(cell, "2F5496")
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    for v, cnt in sorted(database_stats.get("vd_lc_counts", {}).items(), key=lambda x: -x[1]):
        add_table_row(tbl_c, [v, cnt])
    doc.add_paragraph()

    # D. Clinical trial
    doc.add_heading("D. Highest Clinical Trial Phase", level=2)
    tbl_d = doc.add_table(rows=1, cols=2)
    tbl_d.style = "Table Grid"
    tbl_d.rows[0].cells[0].text = "Phase"
    tbl_d.rows[0].cells[1].text = "Count"
    for cell in tbl_d.rows[0].cells:
        set_cell_shading(cell, "2F5496")
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    for t, cnt in sorted(database_stats.get("trial_counts", {}).items(), key=lambda x: -x[1]):
        add_table_row(tbl_d, [t, cnt])
    doc.add_paragraph()

    # E. CDR length
    doc.add_heading("E. CDR Length Distribution (IMGT numbering)", level=2)
    tbl_e = doc.add_table(rows=1, cols=6)
    tbl_e.style = "Table Grid"
    for j, h in enumerate(["CDR", "Min", "Max", "Mean", "Median", "N"]):
        tbl_e.rows[0].cells[j].text = h
        set_cell_shading(tbl_e.rows[0].cells[j], "2F5496")
        for p in tbl_e.rows[0].cells[j].paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
        s = database_stats.get("cdr_length_summary", {}).get(key)
        if s:
            add_table_row(tbl_e, [key, s["min"], s["max"], s["mean"], s["median"], s["count"]])
        else:
            add_table_row(tbl_e, [key, "-", "-", "-", "-", 0])
    doc.add_paragraph()

    # F. Top targets
    doc.add_heading("F. Top 30 Therapeutic Targets", level=2)
    tbl_f = doc.add_table(rows=1, cols=2)
    tbl_f.style = "Table Grid"
    tbl_f.rows[0].cells[0].text = "Target"
    tbl_f.rows[0].cells[1].text = "Count"
    for cell in tbl_f.rows[0].cells:
        set_cell_shading(cell, "2F5496")
        for p in cell.paragraphs:
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8)
    tc = database_stats.get("target_counts", {})
    for t, cnt in sorted(tc.items(), key=lambda x: -x[1])[:30]:
        add_table_row(tbl_f, [t, cnt])
    doc.add_paragraph()

    # Save
    docx_path = os.path.join(output_dir, "cdr_homology_report.docx")
    doc.save(docx_path)
    print(f"  Word report written: {docx_path}")
    return docx_path


def load_database_stats_from_results(results):
    """Dummy: not used in report, kept for structure."""
    return None


def load_stats(stats_path: str) -> dict:
    """Load precomputed statistics JSON."""
    if os.path.exists(stats_path):
        with open(stats_path) as f:
            return json.load(f)
    return {"total_entries": 0, "format_counts": {}, "genetics_counts": {},
            "vd_lc_counts": {}, "trial_counts": {}, "target_counts": {},
            "cdr_length_summary": {}}


# ---------------------------------------------------------------------------
# Input parsing helpers
# ---------------------------------------------------------------------------

def parse_query_inputs(args) -> list:
    """Return [(name, vh_seq, vl_seq), ...] from CLI args."""
    queries = []
    if args.input:
        records = parse_fasta(args.input)
        # Group by name: VH + VL pairs
        # Convention: header contains VH/VL/Light/Heavy or suffix _VH/_VL
        vh_map = {}
        vl_map = {}
        for hdr, seq in records:
            h = hdr.upper()
            if "VH" in h or "HEAVY" in h:
                vh_map[hdr] = seq
            elif "VL" in h or "LIGHT" in h:
                vl_map[hdr] = seq
            elif len(seq) > 130:
                vh_map[hdr] = seq
            else:
                vl_map[hdr] = seq

        # Pair by base name
        all_names = list(dict.fromkeys(list(vh_map.keys()) + list(vl_map.keys())))
        used = set()
        for hdr in all_names:
            if hdr in used:
                continue
            base = hdr.split("_")[0].split()[0]
            vh = vh_map.get(hdr)
            vl = vl_map.get(hdr)
            # try to find partner
            if vh and not vl:
                for vhdr in vl_map:
                    if vhdr.split("_")[0].split()[0] == base and vhdr not in used:
                        vl = vl_map[vhdr]
                        used.add(vhdr)
                        break
            elif vl and not vh:
                for vhdr in vh_map:
                    if vhdr.split("_")[0].split()[0] == base and vhdr not in used:
                        vh = vh_map[vhdr]
                        used.add(vhdr)
                        break
            queries.append((hdr, vh, vl))
            used.add(hdr)
    elif args.vh or args.vl:
        name = args.name or "query_1"
        queries.append((name, args.vh, args.vl))
    return queries


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="CDR Homology Analysis: compare antibody CDRs against clinical database"
    )
    parser.add_argument("--database", default="data/therasabdab_cdrs.json",
                        help="Preprocessed CDR database JSON")
    parser.add_argument("--stats", default="data/therasabdab_stats.json",
                        help="Precomputed statistics JSON")
    parser.add_argument("--input", "-i", help="FASTA file with query sequences")
    parser.add_argument("--vh", help="VH amino acid sequence")
    parser.add_argument("--vl", help="VL amino acid sequence")
    parser.add_argument("--name", help="Query name (for single --vh/--vl mode)")
    parser.add_argument("--top", type=int, default=5, help="Top N matches per CDR (default 5)")
    parser.add_argument("--output", "-o", default="outputs/cdr_homology",
                        help="Output directory")
    args = parser.parse_args()

    if not args.vh and not args.vl and not args.input:
        parser.error("Provide --input FASTA or --vh/--vl sequences")

    # Load database
    db_path = args.database
    if not os.path.exists(db_path):
        print(f"Error: database not found: {db_path}")
        print("Run preprocess.py first:  python3 tools/cdr_homology/preprocess.py")
        sys.exit(1)

    print(f"Loading database from {db_path} ...")
    database = load_database(db_path)
    print(f"  {len(database)} clinical antibodies loaded")

    # Parse queries
    queries = parse_query_inputs(args)
    if not queries:
        print("Error: no queries found")
        sys.exit(1)
    print(f"  {len(queries)} query(ies) to process")

    # Run search
    print("Running CDR homology search ...")
    all_results = run_query(queries, database, top_n=args.top)

    # Load stats for appendix
    stats = load_stats(args.stats)

    # Generate report
    print("Generating report ...")
    md_path = generate_report(all_results, stats, args.output, args.top)

    # Generate Word report
    print("Generating Word report ...")
    generate_docx_report(all_results, stats, args.output, args.top)

    # Also save raw JSON
    json_path = os.path.join(args.output, "cdr_homology_results.json")
    with open(json_path, "w") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)
    print(f"  Raw JSON written: {json_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
