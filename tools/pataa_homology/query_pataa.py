#!/usr/bin/env python3
"""CDR homology analysis against the PATAA patent antibody database.

Build the database first:
    python3 tools/pataa_homology/preprocess_pataa.py all

Then query:
    # single pair
    python3 tools/pataa_homology/query_pataa.py \
        --vh EVQLQQSGPELVKPGASVKMSCKAS... \
        --vl DIQMTQTTSSLSASLGDRVTISCRAS... \
        --top 5 --output outputs/pataa_homology/

    # from FASTA (VH/VL paired by header suffix)
    python3 tools/pataa_homology/query_pataa.py \
        --input outputs/.../variants.fasta --top 5 \
        --output outputs/pataa_homology/

Mirrors tools/cdr_homology/query.py (ANARCI IMGT CDR extraction + BLOSUM62
alignment) but reports patent / accession / germline metadata.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from collections import Counter


_IMGT_CDR_RANGES = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR3": (105, 117)}


# ---------------------------------------------------------------------------
# ANARCI CDR extraction
# ---------------------------------------------------------------------------

def _extract_cdrs_from_posmap(posmap: dict, chain_type: str) -> dict:
    cdrs = {}
    for cdr, (start, end) in _IMGT_CDR_RANGES.items():
        seq = ""
        for pos in range(start, end + 1):
            aa = posmap.get(f"{chain_type}{pos}")
            if aa and aa != "-":
                seq += aa
            for letter in "ABCDEFGHIJK":
                iaa = posmap.get(f"{chain_type}{pos}{letter}")
                if iaa:
                    seq += iaa
        cdrs[cdr] = seq
    return cdrs


def extract_cdrs(sequence: str, chain_type: str):
    """ANARCI (IMGT) numbering then CDR1/2/3 extraction; None on failure."""
    try:
        from anarci import anarci
    except ImportError:
        return None
    try:
        res = anarci([("q", sequence.upper().strip())], scheme="imgt", output=False)
        numbering = res[0][0][0][0]
    except (IndexError, TypeError, Exception):
        return None
    if not numbering:
        return None
    posmap = {}
    for (pos_num, ins_flag), aa in numbering:
        if aa == "-":
            continue
        label = (f"{chain_type}{pos_num}{ins_flag}"
                 if ins_flag and ins_flag != " " else f"{chain_type}{pos_num}")
        posmap[label] = aa.upper()
    return _extract_cdrs_from_posmap(posmap, chain_type)


# ---------------------------------------------------------------------------
# FASTA parsing
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> list:
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


# ---------------------------------------------------------------------------
# alignment (BLOSUM62, mirror cdr_homology: global for CDR1/2, local for CDR3)
# ---------------------------------------------------------------------------

_aligner = None


def _get_aligner(mode="global"):
    global _aligner
    from Bio.Align import PairwiseAligner, substitution_matrices
    a = PairwiseAligner()
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.mode = mode
    return a


def align_sequences(query: str, subject: str) -> dict:
    if not query or not subject:
        return {"identity": 0.0, "score": 0.0}
    aligner = _get_aligner("global")
    try:
        alns = aligner.align(subject, query)
    except Exception:
        return {"identity": 0.0, "score": 0.0}
    if not alns:
        return {"identity": 0.0, "score": 0.0}
    best = alns[0]
    coords = best.coordinates
    seqs = best.sequences
    q_parts, s_parts = [], []
    for i in range(coords.shape[1] - 1):
        ts, te = int(coords[0, i]), int(coords[0, i + 1])
        qs, qe = int(coords[1, i]), int(coords[1, i + 1])
        if ts == te:
            q_parts.append(seqs[1][qs:qe]); s_parts.append("-" * (qe - qs))
        elif qs == qe:
            q_parts.append("-" * (te - ts)); s_parts.append(seqs[0][ts:te])
        else:
            q_parts.append(seqs[1][qs:qe]); s_parts.append(seqs[0][ts:te])
    aln_q, aln_s = "".join(q_parts), "".join(s_parts)
    matches = sum(1 for a, b in zip(aln_q, aln_s) if a == b and a not in "-.")
    alen = max(len(aln_q), 1)
    return {"identity": round(matches / alen, 4), "score": round(best.score, 1),
            "matches": matches, "aligned_len": alen,
            "query_aligned": aln_q, "subject_aligned": aln_s}


# ---------------------------------------------------------------------------
# database search
# ---------------------------------------------------------------------------

def load_database(path: str) -> list:
    with open(path) as fh:
        return json.load(fh)


def build_index(database: list) -> dict:
    """[(entry_idx, cdr_str)] lists per (chain, cdr), for fast prefiltering."""
    index = {}
    for ck in ("VH", "VL"):
        for cdr in ("CDR1", "CDR2", "CDR3"):
            entries = []
            for i, e in enumerate(database):
                ch = e.get("chains", {}).get(ck)
                if ch and ch.get(cdr):
                    entries.append((i, ch[cdr]))
            index[(ck, cdr)] = entries
    return index


def _kmers(s: str, k: int) -> set:
    if len(s) < k:
        return {s}
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def search_cdr(query_cdr: str, database: list, index: dict, chain_key: str,
               cdr_name: str, top_n: int = 5, prefilter: int = 400) -> list:
    entries = index.get((chain_key, cdr_name), [])
    q = query_cdr
    L = len(q)
    if not q or not entries:
        return []

    # ---- cheap k-mer + length prefilter ----
    if L < 2:
        cand = [(1.0, i) for i, s in entries if s == q]
    else:
        k = 3 if L >= 6 else 2
        qk = _kmers(q, k)
        nq = len(qk)
        if cdr_name == "CDR3":
            lo, hi = L - 5, L + 5
        else:
            lo, hi = L - 2, L + 2
        overlap = []
        for i, s in entries:
            ls = len(s)
            if ls < lo or ls > hi:
                continue
            inter = len(qk & _kmers(s, k))
            if inter:
                overlap.append((inter / nq, i))
        overlap.sort(key=lambda x: -x[0])
        cand = overlap[:prefilter]

    results = []
    for frac, i in cand:
        e = database[i]
        s = e["chains"][chain_key][cdr_name]
        aln = align_sequences(q, s)
        results.append({
            "name": e.get("name", ""),
            "patents": e.get("patents", []),
            "species": e.get("species", ""),
            "v_germline": e.get("v_germline", ""),
            "v_identity": e.get("v_identity", 0.0),
            "v_bitscore": e.get("v_bitscore", 0.0),
            "identity": aln["identity"],
            "score": aln.get("score", 0.0),
            "subject_cdr": s,
            "query_aligned": aln.get("query_aligned", ""),
            "subject_aligned": aln.get("subject_aligned", ""),
        })
    results.sort(key=lambda x: (-x["identity"], -x["score"]))
    return results[:top_n]


def run_query(queries: list, database: list, top_n: int = 5) -> list:
    index = build_index(database)
    out = []
    for name, vh, vl in queries:
        print(f"  Processing {name} ...")
        res = {"name": name, "matches": {}}
        for chain_key, seq in (("VH", vh), ("VL", vl)):
            if not seq:
                continue
            cdrs = extract_cdrs(seq, chain_key[1])
            if not cdrs:
                continue
            for cdr_name, cdr_seq in cdrs.items():
                key = f"{chain_key}_{cdr_name}"
                hits = search_cdr(cdr_seq, database, index, chain_key, cdr_name, top_n)
                res["matches"][key] = {"query_cdr": cdr_seq, "top_matches": hits,
                                       "top_match": hits[0] if hits else None}
        out.append(res)
    return out


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------

def _pat(entry) -> str:
    pats = entry.get("patents") or []
    return "; ".join(pats[:3]) if pats else "-"


def generate_markdown(all_results, stats, output_dir, top_n):
    L = []
    L.append("# PATAA Patent Antibody CDR Homology Report\n")
    L.append(f"**Queries analysed:** {len(all_results)}\n")
    L.append(f"**Database:** {stats.get('total_after_dedup', '?')} unique patent "
             f"antibody sequences ({stats.get('total_positive_records', '?')} records "
             f"before dedup)\n")
    L.append(f"**Detection:** {stats.get('detection', '')}; "
             f"**Numbering:** {stats.get('numbering', '')}\n")

    for i, r in enumerate(all_results, 1):
        L.append(f"\n---\n## {i}. {r['name']}\n")
        L.append("### CDR table (IMGT)\n")
        L.append("| Region | Query CDR | Len | Top accession | Patents | CDR identity |")
        L.append("|--------|-----------|-----|---------------|---------|--------------|")
        for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(key)
            if not m:
                continue
            q = m["query_cdr"]
            tm = m["top_match"]
            if tm:
                L.append(f"| {key} | `{q}` | {len(q)} | {tm['name']} | {_pat(tm)} "
                         f"| {tm['identity']:.1%} |")
            else:
                L.append(f"| {key} | `{q}` | {len(q)} | - | - | - |")
        # top hits per CDR
        for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(key)
            if not m:
                continue
            L.append(f"\n#### {key} — top {top_n}\n")
            L.append("| # | Accession | Identity | V-germline | V-id | Patents | Subject CDR |")
            L.append("|---|-----------|----------|------------|------|---------|-------------|")
            for j, h in enumerate(m["top_matches"], 1):
                L.append(f"| {j} | {h['name']} | {h['identity']:.1%} | {h['v_germline']} "
                         f"| {h['v_identity']:.0f}% | {_pat(h)} | `{h['subject_cdr']}` |")
    # appendix
    L.append("\n---\n# Appendix: PATAA database summary\n")
    cd = stats.get("chain_distribution", {})
    if cd:
        L.append("## Chain distribution\n")
        L.append("| Chain | Count |\n|-------|-------|")
        for k, v in cd.items():
            L.append(f"| {k} | {v} |")
    L.append("\n## CDR length summary\n")
    L.append("| CDR | Count | Min | Max | Mean |")
    L.append("|-----|-------|-----|-----|------|")
    for k, s in stats.get("cdr_length_summary", {}).items():
        L.append(f"| {k} | {s['count']} | {s['min']} | {s['max']} | {s['mean']} |")
    L.append("\n## Top 20 V germlines in database\n")
    L.append("| Germline | Count |\n|----------|-------|")
    for g, c in list(stats.get("top_germlines", {}).items())[:20]:
        L.append(f"| {g} | {c} |")
    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "pataa_homology_report.md")
    with open(path, "w") as fh:
        fh.write("\n".join(L))
    print(f"  Markdown report: {path}")
    return path


def generate_csv(all_results, output_dir):
    csv_dir = os.path.join(output_dir, "detail")
    os.makedirs(csv_dir, exist_ok=True)
    for r in all_results:
        safe = r["name"].replace("/", "_").replace(" ", "_")[:60]
        path = os.path.join(csv_dir, f"{safe}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["CDR", "Query CDR", "Rank", "Accession", "Identity",
                        "Score", "V-germline", "V-identity", "Subject CDR", "Patents"])
            for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
                m = r["matches"].get(key)
                if not m:
                    continue
                for j, h in enumerate(m["top_matches"], 1):
                    w.writerow([key, m["query_cdr"], j, h["name"],
                                f"{h['identity']:.4f}", h["score"], h["v_germline"],
                                h["v_identity"], h["subject_cdr"], _pat(h)])
    print(f"  CSV detail dir: {csv_dir}")


def generate_docx_report(all_results, stats, output_dir, top_n):
    """Word report mirroring tools/cdr_homology/query.py styling."""
    from docx import Document
    from docx.shared import Pt, RGBColor
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml.ns import qn

    def shade(cell, color):
        el = cell._element.get_or_add_tcPr().makeelement(qn("w:shd"), {
            qn("w:val"): "clear", qn("w:color"): "auto", qn("w:fill"): color})
        cell._element.get_or_add_tcPr().append(el)

    def add_row(table, values, header=False):
        row = table.add_row()
        for i, v in enumerate(values):
            cell = row.cells[i]
            cell.text = str(v)
            for p in cell.paragraphs:
                p.style = doc.styles["Normal"]
                for run in p.runs:
                    run.font.size = Pt(8)
            if header:
                shade(cell, "2F5496")
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.bold = True
                        run.font.color.rgb = RGBColor(255, 255, 255)

    def head(table, headers):
        for i, h in enumerate(headers):
            table.rows[0].cells[i].text = h
        for cell in table.rows[0].cells:
            shade(cell, "2F5496")
            for p in cell.paragraphs:
                for run in p.runs:
                    run.bold = True
                    run.font.color.rgb = RGBColor(255, 255, 255)
                    run.font.size = Pt(8)

    doc = Document()
    title = doc.add_heading("PATAA Patent Antibody CDR Homology Report", level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.add_paragraph(f"Queries analysed: {len(all_results)}")
    doc.add_paragraph(
        f"Database: {stats.get('total_after_dedup', '?')} unique patent antibody "
        f"sequences ({stats.get('total_positive_records', '?')} records before dedup)")
    doc.add_paragraph(f"Detection: {stats.get('detection', '')}")
    doc.add_paragraph(f"Numbering: {stats.get('numbering', '')}")
    doc.add_paragraph()

    # ---- Summary table ----
    doc.add_heading("Summary Table", level=1)
    tbl = doc.add_table(rows=1, cols=7)
    tbl.style = "Table Grid"
    tbl.alignment = WD_TABLE_ALIGNMENT.CENTER
    head(tbl, ["#", "Antibody", "VH_CDR3 Best", "VH_CDR3 Id",
               "VL_CDR3 Best", "VL_CDR3 Id", "Mean Id"])
    for i, r in enumerate(all_results, 1):
        vh3 = r["matches"].get("VH_CDR3", {}).get("top_match")
        vl3 = r["matches"].get("VL_CDR3", {}).get("top_match")
        def nm(m):
            return (m["name"] if m else "N/A")
        def pid(m):
            return f"{m['identity']:.1%}" if m else "N/A"
        ids = [v["top_match"]["identity"] for v in r["matches"].values() if v.get("top_match")]
        mean_id = f"{sum(ids)/len(ids):.1%}" if ids else "N/A"
        add_row(tbl, [i, r["name"][:30], nm(vh3), pid(vh3), nm(vl3), pid(vl3), mean_id])
    doc.add_paragraph()

    # ---- Per-query detail ----
    for i, r in enumerate(all_results, 1):
        doc.add_heading(f"{i}. {r['name']}", level=1)

        doc.add_heading("CDR Sequences (IMGT numbering)", level=2)
        t2 = doc.add_table(rows=1, cols=6)
        t2.style = "Table Grid"
        head(t2, ["Region", "Query CDR", "Len", "Top Accession", "Patents", "Identity"])
        for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(key)
            if not m:
                continue
            tm = m["top_match"]
            if tm:
                add_row(t2, [key, m["query_cdr"], len(m["query_cdr"]),
                             tm["name"], _pat(tm), f"{tm['identity']:.1%}"])
            else:
                add_row(t2, [key, m["query_cdr"], len(m["query_cdr"]), "-", "-", "-"])
        doc.add_paragraph()

        for key in ["VH_CDR1", "VH_CDR2", "VH_CDR3", "VL_CDR1", "VL_CDR2", "VL_CDR3"]:
            m = r["matches"].get(key)
            if not m:
                continue
            doc.add_heading(f"{key} - Top {top_n} Matches", level=2)
            t3 = doc.add_table(rows=1, cols=7)
            t3.style = "Table Grid"
            head(t3, ["Rank", "Accession", "Identity", "V-germline",
                      "V-id", "Patents", "Subject CDR"])
            for j, h in enumerate(m["top_matches"][:top_n], 1):
                add_row(t3, [j, h["name"], f"{h['identity']:.1%}", h["v_germline"],
                             f"{h['v_identity']:.0f}%", _pat(h), h["subject_cdr"]])
            doc.add_paragraph()

    # ---- Appendix ----
    doc.add_page_break()
    doc.add_heading("Appendix: PATAA Database Summary", level=1)

    doc.add_heading("A. Chain Distribution", level=2)
    ta = doc.add_table(rows=1, cols=2)
    ta.style = "Table Grid"
    head(ta, ["Chain", "Count"])
    for k, v in stats.get("chain_distribution", {}).items():
        add_row(ta, [k, v])
    doc.add_paragraph()

    doc.add_heading("B. CDR Length Distribution (IMGT)", level=2)
    tb = doc.add_table(rows=1, cols=5)
    tb.style = "Table Grid"
    head(tb, ["CDR", "Count", "Min", "Max", "Mean"])
    for k, s in stats.get("cdr_length_summary", {}).items():
        add_row(tb, [k, s["count"], s["min"], s["max"], s["mean"]])
    doc.add_paragraph()

    doc.add_heading("C. Top 20 V Germlines in Database", level=2)
    tc = doc.add_table(rows=1, cols=2)
    tc.style = "Table Grid"
    head(tc, ["Germline", "Count"])
    for g, c in list(stats.get("top_germlines", {}).items())[:20]:
        add_row(tc, [g, c])

    path = os.path.join(output_dir, "pataa_homology_report.docx")
    doc.save(path)
    print(f"  Word report: {path}")
    return path


def parse_query_inputs(args) -> list:
    queries = []
    if args.input:
        recs = parse_fasta(args.input)
        vh_map, vl_map = {}, {}
        for hdr, seq in recs:
            h = hdr.upper()
            if "VH" in h or "HEAVY" in h or "_H" in h:
                vh_map[hdr] = seq
            elif "VL" in h or "LIGHT" in h or "_L" in h:
                vl_map[hdr] = seq
            elif len(seq) > 130:
                vh_map[hdr] = seq
            else:
                vl_map[hdr] = seq
        def base(n):
            for suf in ("_Heavy", "_Light", "_VH", "_VL"):
                if n.endswith(suf):
                    return n[:-len(suf)]
            return n
        all_names = list(dict.fromkeys(list(vh_map) + list(vl_map)))
        used = set()
        for hdr in all_names:
            if hdr in used:
                continue
            b = base(hdr)
            vh, vl = vh_map.get(hdr), vl_map.get(hdr)
            if vh and not vl:
                for cand, seq in vl_map.items():
                    if cand not in used and base(cand) == b:
                        vl = seq; used.add(cand); break
            elif vl and not vh:
                for cand, seq in vh_map.items():
                    if cand not in used and base(cand) == b:
                        vh = seq; used.add(cand); break
            queries.append((hdr, vh, vl))
            used.add(hdr)
    elif args.vh or args.vl:
        queries.append((args.name or "query_1", args.vh, args.vl))
    return queries


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--database", default="data/pataa_cdrs.json")
    ap.add_argument("--stats", default="data/pataa_stats.json")
    ap.add_argument("--input", "-i")
    ap.add_argument("--vh")
    ap.add_argument("--vl")
    ap.add_argument("--name")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--no-docx", action="store_true", help="skip Word report")
    ap.add_argument("--output", "-o", default="outputs/pataa_homology")
    args = ap.parse_args()

    if not args.vh and not args.vl and not args.input:
        ap.error("provide --input FASTA or --vh/--vl")
    if not os.path.exists(args.database):
        sys.exit(f"database not found: {args.database}\n"
                 f"run: python3 tools/pataa_homology/preprocess_pataa.py all")

    print(f"Loading PATAA database: {args.database}")
    db = load_database(args.database)
    print(f"  {len(db)} unique patent antibody sequences")
    stats = {}
    if os.path.exists(args.stats):
        with open(args.stats) as fh:
            stats = json.load(fh)

    queries = parse_query_inputs(args)
    print(f"  {len(queries)} query(ies)")
    results = run_query(queries, db, args.top)

    generate_markdown(results, stats, args.output, args.top)
    generate_csv(results, args.output)
    if not args.no_docx:
        try:
            generate_docx_report(results, stats, args.output, args.top)
        except ImportError:
            print("  [WARN] python-docx not installed; skipped Word report")
    json_path = os.path.join(args.output, "pataa_homology_results.json")
    os.makedirs(args.output, exist_ok=True)
    with open(json_path, "w") as fh:
        json.dump(results, fh, ensure_ascii=False, indent=2)
    print(f"  JSON: {json_path}")
    print("Done.")


if __name__ == "__main__":
    main()
