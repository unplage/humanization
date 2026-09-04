#!/usr/bin/env python3
"""Evaluate humanization degree of antibody sequences.

Compares input sequences against human germline V/J genes and reports:
  - Best matching germline (by FR identity)
  - FR identity, CDR identity, overall identity (Kabat and IMGT numbering)
  - Humanization level per USAN/WHO standards (IMGT-based)
  - Number of positions differing from germline (potential back-mutation candidates)

Usage:
    python3 -m scripts.humanize.evaluate_humanness --input seq.fasta
    python3 -m scripts.humanize.evaluate_humanness --vh EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS --vl DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
"""

import argparse
import sys
import os
from typing import Dict

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def evaluate_sequence(chain_type: str, sequence: str, db_dir: str = ""):
    """Evaluate a single chain against human germline."""
    from .numbering import number_heavy, number_light
    from .germline import load_germline_db, compare_to_germline

    # Number the sequence
    if chain_type == "H":
        numbered = number_heavy(sequence)
    else:
        numbered = number_light(sequence)

    # Load germline DB
    if not db_dir:
        db_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "germline")
    db = load_germline_db(db_dir)

    # Get human V genes for this chain type
    v_genes = db.human(chain_type)

    # Score each germline (Kabat-based)
    scored = []
    for gene in v_genes:
        scores = compare_to_germline(numbered, gene)
        scored.append((gene, scores))

    # Sort by FR identity (primary) and CDR identity (secondary)
    scored.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))

    return scored, numbered


def evaluate_sequence_imgt(chain_type: str, sequence: str, db_dir: str = ""):
    """Evaluate a single chain against human germline using IMGT numbering.

    Uses AbRSA with IMGT scheme to number the query sequence directly,
    then loads IMGT germline database for comparison.
    """
    from .imgt_numbering import (
        number_with_abrsa_imgt, load_imgt_germline, compare_imgt_posmaps_direct
    )

    # Number the sequence using AbRSA with IMGT scheme
    numbered = number_with_abrsa_imgt(sequence, chain_type)

    if numbered is not None:
        # AbRSA succeeded - use IMGT-numbered posmap directly
        query_imgt = numbered.posmap()
    else:
        # AbRSA failed - fallback to Kabat numbering + conversion
        from .numbering import number_heavy, number_light
        from .imgt_numbering import kabat_posmap_to_imgt_posmap

        if chain_type == "H":
            numbered_kabat = number_heavy(sequence)
        else:
            numbered_kabat = number_light(sequence)

        query_imgt = kabat_posmap_to_imgt_posmap(numbered_kabat.posmap(), chain_type)
        numbered = numbered_kabat  # For display purposes

    # Load IMGT germline DB (abnumber_human_imgt.json)
    imgt_genes = load_imgt_germline(db_dir)

    if not imgt_genes:
        # Fallback: use Kabat germline with conversion
        from .germline import load_germline_db
        from .imgt_numbering import kabat_posmap_to_imgt_posmap

        if not db_dir:
            db_dir = os.path.join(
                os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                "data", "germline"
            )
        db = load_germline_db(db_dir)
        v_genes = db.human(chain_type)

        scored = []
        for gene in v_genes:
            if gene.numbered is None:
                continue
            germline_imgt = kabat_posmap_to_imgt_posmap(gene.numbered.posmap(), chain_type)
            scores = compare_imgt_posmaps_direct(query_imgt, germline_imgt)
            scored.append((gene, scores))
    else:
        # Use IMGT germline directly
        scored = []
        for gene in imgt_genes:
            # Filter by chain type
            if chain_type == "H" and not gene.gene_id.startswith("IGHV"):
                continue
            if chain_type == "L" and not gene.gene_id.startswith(("IGKV", "IGLV")):
                continue

            # Compare using direct IMGT comparison
            scores = compare_imgt_posmaps_direct(query_imgt, gene._imgt_posmap)
            scored.append((gene, scores))

    # Sort by FR identity (primary) and CDR identity (secondary)
    scored.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))

    return scored, numbered


def classify_humanization_level(fr_identity: float) -> str:
    """Classify humanization level based on FR identity.
    
    NOTE: This is an ACADEMIC ESTIMATE based on sequence identity.
    The official USAN/INN naming is based on the technology used to create
    the antibody (transgenic mouse, CDR grafting, phage display, etc.),
    NOT on sequence identity thresholds.
    
    Academic conventions (commonly used in literature):
    - Human (-umab): FR identity >= 95% (typically from transgenic mice/phage display)
    - Humanized (-zumab): FR identity >= 85% (typically CDR grafting)
    - Chimeric (-ximab): FR identity >= 70% (typically mouse V + human C)
    - Murine (-omab): FR identity < 70%
    
    Official USAN/INN classification requires knowledge of the manufacturing process.
    """
    if fr_identity >= 0.95:
        return "Academic: likely Human-like"
    elif fr_identity >= 0.85:
        return "Academic: likely Humanized-like"
    elif fr_identity >= 0.70:
        return "Academic: likely Chimeric-like"
    else:
        return "Academic: likely Murine-like"


def format_sequence_comparison_table(numbered, germline_numbered, chain_type: str):
    """Generate a table showing sequences under different numbering schemes.

    Args:
        numbered: NumberedChain object for query sequence
        germline_numbered: NumberedChain object for germline sequence
        chain_type: "H" or "L"

    Returns:
        Formatted string with sequence comparison table
    """
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  Sequence Comparison Table (Query vs Germline)")
    lines.append(f"{'='*70}")

    # Get posmaps
    q_map = numbered.posmap()
    g_map = germline_numbered.posmap() if germline_numbered else {}

    def get_pos_sort_key(pos: str) -> tuple:
        """Extract sort key from position label like 'H31' or 'H100A'.
        
        Returns (numeric_part, letter_part) tuple for proper sorting.
        H100A -> (100, 'A')
        H100B -> (100, 'B')
        H31 -> (31, '')
        """
        num_str = ""
        letter_str = ""
        for c in pos[1:]:
            if c.isdigit():
                num_str += c
            else:
                letter_str += c
        return (int(num_str) if num_str else 0, letter_str)

    def get_region_seq(posmap: Dict, start_num: int, end_num: int, prefix: str) -> str:
        """Extract sequence for a region given start/end position numbers."""
        seq = ""
        for pos in sorted(posmap.keys(), key=lambda x: get_pos_sort_key(x)):
            if pos.startswith(prefix):
                num = get_pos_sort_key(pos)[0]
                if start_num <= num <= end_num:
                    seq += posmap[pos]
        return seq

    # Kabat numbering regions
    if chain_type == "H":
        kabat_regions = {
            "FR1": (1, 30),
            "CDR1": (31, 35),
            "FR2": (36, 49),
            "CDR2": (50, 65),
            "FR3": (66, 94),
            "CDR3": (95, 102),
            "FR4": (103, 113),
        }
        imgt_regions = {
            "FR1": (1, 26),
            "CDR1": (27, 38),
            "FR2": (39, 55),
            "CDR2": (56, 65),
            "FR3": (66, 104),
            "CDR3": (105, 117),
            "FR4": (118, 128),
        }
        prefix = "H"
    else:  # VL
        kabat_regions = {
            "FR1": (1, 23),
            "CDR1": (24, 34),
            "FR2": (35, 49),
            "CDR2": (50, 56),
            "FR3": (57, 88),
            "CDR3": (89, 97),
            "FR4": (98, 107),
        }
        imgt_regions = {
            "FR1": (1, 23),
            "CDR1": (24, 34),
            "FR2": (35, 55),
            "CDR2": (56, 65),
            "FR3": (66, 88),
            "CDR3": (89, 104),
            "FR4": (105, 118),
        }
        prefix = "L"

    # Kabat numbering table
    lines.append(f"\n  Kabat numbering:")
    lines.append(f"  {'Region':<8} {'Positions':<15} {'Query':<30} {'Germline':<30}")
    lines.append(f"  {'-'*8} {'-'*15} {'-'*30} {'-'*30}")

    for region, (start_num, end_num) in kabat_regions.items():
        q_seq = get_region_seq(q_map, start_num, end_num, prefix)
        g_seq = get_region_seq(g_map, start_num, end_num, prefix)

        # Truncate if too long
        q_display = q_seq[:28] + ".." if len(q_seq) > 30 else q_seq
        g_display = g_seq[:28] + ".." if len(g_seq) > 30 else g_seq

        lines.append(f"  {region:<8} {prefix}{start_num}-{prefix}{end_num:<12} {q_display:<30} {g_display:<30}")

    # IMGT numbering table
    lines.append(f"\n  IMGT numbering:")
    lines.append(f"  {'Region':<8} {'Positions':<15} {'Query':<30} {'Germline':<30}")
    lines.append(f"  {'-'*8} {'-'*15} {'-'*30} {'-'*30}")

    for region, (start_num, end_num) in imgt_regions.items():
        q_seq = get_region_seq(q_map, start_num, end_num, prefix)
        g_seq = get_region_seq(g_map, start_num, end_num, prefix)

        # Truncate if too long
        q_display = q_seq[:28] + ".." if len(q_seq) > 30 else q_seq
        g_display = g_seq[:28] + ".." if len(g_seq) > 30 else g_seq

        lines.append(f"  {region:<8} {prefix}{start_num}-{prefix}{end_num:<12} {q_display:<30} {g_display:<30}")

    # Sequence identity summary
    lines.append(f"\n  Sequence identity summary:")
    lines.append(f"  {'Region':<10} {'Identity':<12} {'Query Length':<15} {'Germline Length':<15}")
    lines.append(f"  {'-'*10} {'-'*12} {'-'*15} {'-'*15}")

    for region, (start_num, end_num) in kabat_regions.items():
        q_seq = get_region_seq(q_map, start_num, end_num, prefix)
        g_seq = get_region_seq(g_map, start_num, end_num, prefix)
        q_count = len(q_seq)
        g_count = len(g_seq)
        match_count = sum(1 for q_aa, g_aa in zip(q_seq, g_seq) if q_aa == g_aa and g_aa != "-")
        identity = match_count / q_count if q_count > 0 else 0.0
        lines.append(f"  {region:<10} {identity:.1%}     {q_count:<15} {g_count:<15}")

    return "\n".join(lines)


def format_results(chain_type: str, sequence: str, scored, numbered, top_n: int = 10):
    """Format Kabat-based results for display."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  {chain_type} chain ({len(sequence)} residues) - Kabat numbering")
    lines.append(f"{'='*70}")

    # Sequence with region annotation
    lines.append(f"\n  Sequence (first 50 aa): {sequence[:50]}...")

    # Top matches
    lines.append(f"\n  Top {min(top_n, len(scored))} human germline matches (Kabat):")
    lines.append(f"  {'Rank':<6} {'Gene':<20} {'FR id':<10} {'CDR id':<10} {'Overall':<10}")
    lines.append(f"  {'-'*6} {'-'*20} {'-'*10} {'-'*10} {'-'*10}")

    for i, (gene, scores) in enumerate(scored[:top_n]):
        lines.append(
            f"  {i+1:<6} {gene.gene_id:<20} "
            f"{scores['fr_identity']:<10.4f} "
            f"{scores['cdr_identity']:<10.4f} "
            f"{scores['all_identity']:<10.4f}"
        )

    # Best match details
    best_gene, best_scores = scored[0]
    lines.append(f"\n  Best match (Kabat): {best_gene.gene_id}")
    lines.append(f"    FR identity:  {best_scores['fr_identity']:.4f} ({best_scores['fr_identity']*100:.1f}%)")
    lines.append(f"    CDR identity: {best_scores['cdr_identity']:.4f} ({best_scores['cdr_identity']*100:.1f}%)")
    lines.append(f"    Overall:      {best_scores['all_identity']:.4f} ({best_scores['all_identity']*100:.1f}%)")

    # Detailed region breakdown
    lines.append(f"\n  Detailed region breakdown (Kabat):")
    lines.append(f"  {'Region':<10} {'Identity':<12} {'Match':<8} {'Total':<8}")
    lines.append(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*8}")
    lines.append(f"  {'FR1':<10} {best_scores.get('fr1_identity', 0):.1%}     {best_scores.get('n_fr1', 0):<8}")
    lines.append(f"  {'FR2':<10} {best_scores.get('fr2_identity', 0):.1%}     {best_scores.get('n_fr2', 0):<8}")
    lines.append(f"  {'FR3':<10} {best_scores.get('fr3_identity', 0):.1%}     {best_scores.get('n_fr3', 0):<8}")
    lines.append(f"  {'CDR1':<10} {best_scores.get('cdr1_identity', 0):.1%}     {best_scores.get('n_cdr1', 0):<8}")
    lines.append(f"  {'CDR2':<10} {best_scores.get('cdr2_identity', 0):.1%}     {best_scores.get('n_cdr2', 0):<8}")

    # Note: USAN classification is only shown in IMGT mode

    # Positions differing from best germline (potential back-mutations)
    q_map = numbered.posmap()
    g_map = best_gene.numbered.posmap() if best_gene.numbered else {}
    diffs = []
    for pos in q_map:
        if pos in g_map and q_map[pos] != g_map[pos]:
            region = numbered.region_of(pos) or "?"
            diffs.append((pos, q_map[pos], g_map[pos], region))

    if diffs:
        lines.append(f"\n  Positions differing from {best_gene.gene_id} ({len(diffs)} total):")
        lines.append(f"  {'Position':<12} {'Query':<8} {'Germline':<10} {'Region':<8}")
        lines.append(f"  {'-'*12} {'-'*8} {'-'*10} {'-'*8}")
        for pos, q_aa, g_aa, region in diffs[:30]:  # Show first 30
            lines.append(f"  {pos:<12} {q_aa:<8} {g_aa:<10} {region:<8}")
        if len(diffs) > 30:
            lines.append(f"  ... and {len(diffs)-30} more positions")

    return "\n".join(lines)


def format_results_imgt(chain_type: str, sequence: str, scored_imgt, numbered, top_n: int = 10):
    """Format IMGT-based results for display."""
    from .imgt_numbering import format_imgt_region_alignment

    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  {chain_type} chain ({len(sequence)} residues) - IMGT numbering")
    lines.append(f"  (Per WHO/INN/USAN standards for drug naming)")
    lines.append(f"{'='*70}")

    if not scored_imgt:
        lines.append(f"\n  [No IMGT-scored germlines available]")
        return "\n".join(lines)

    # Top matches
    lines.append(f"\n  Top {min(top_n, len(scored_imgt))} human germline matches (IMGT):")
    lines.append(f"  {'Rank':<6} {'Gene':<20} {'FR id':<10} {'CDR id':<10} {'Overall':<10}")
    lines.append(f"  {'-'*6} {'-'*20} {'-'*10} {'-'*10} {'-'*10}")

    for i, (gene, scores) in enumerate(scored_imgt[:top_n]):
        lines.append(
            f"  {i+1:<6} {gene.gene_id:<20} "
            f"{scores['fr_identity']:<10.4f} "
            f"{scores['cdr_identity']:<10.4f} "
            f"{scores['all_identity']:<10.4f}"
        )

    # Best match details
    best_gene, best_scores = scored_imgt[0]
    lines.append(f"\n  Best match (IMGT): {best_gene.gene_id}")
    lines.append(f"    FR identity:  {best_scores['fr_identity']:.4f} ({best_scores['fr_identity']*100:.1f}%)")
    lines.append(f"    CDR identity: {best_scores['cdr_identity']:.4f} ({best_scores['cdr_identity']*100:.1f}%)")
    lines.append(f"    Overall:      {best_scores['all_identity']:.4f} ({best_scores['all_identity']*100:.1f}%)")

    # IMGT region breakdown
    if "imgt_region_stats" in best_scores:
        lines.append(f"\n  IMGT region breakdown:")
        lines.append(f"  {'Region':<10} {'Identity':<12} {'Match':<8} {'Total':<8}")
        lines.append(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*8}")
        for region in ["FR1", "CDR1", "FR2", "CDR2", "FR3"]:
            stats = best_scores["imgt_region_stats"].get(region, {})
            identity = stats.get("identity", 0)
            count = stats.get("count", 0)
            lines.append(f"  {region:<10} {identity:.1%}     {count:<8}")

    # Humanization classification (Academic estimate based on sequence identity)
    # NOTE: Official USAN/INN naming is based on the technology used (transgenic mouse,
    # CDR grafting, phage display, etc.), NOT on sequence identity thresholds.
    # This is an academic convention commonly used in literature.
    fr_id = best_scores["fr_identity"]
    if fr_id >= 0.95:
        usan_class = "Academic: likely Human-like"
        usan_desc = "FR identity >= 95%"
    elif fr_id >= 0.85:
        usan_class = "Academic: likely Humanized-like"
        usan_desc = "FR identity >= 85%"
    elif fr_id >= 0.70:
        usan_class = "Academic: likely Chimeric-like"
        usan_desc = "FR identity >= 70%"
    else:
        usan_class = "Academic: likely Murine-like"
        usan_desc = "FR identity < 70%"
    lines.append(f"\n    Humanization estimate: {usan_class} ({usan_desc})")
    lines.append(f"    (Note: Official USAN/INN naming requires knowledge of manufacturing process)")

    # Region alignment
    if best_gene.numbered:
        lines.append(f"\n  Region alignment vs {best_gene.gene_id} (IMGT):")
        lines.append(f"  (X = differs from germline, . = matches)")
        alignment = format_imgt_region_alignment(
            numbered.posmap(),
            best_gene.numbered.posmap(),
            chain_type,
        )
        lines.append(alignment)

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Evaluate antibody humanization degree")
    parser.add_argument("--input", "-i", help="FASTA file with sequences")
    parser.add_argument("--vh", help="VH sequence (amino acids)")
    parser.add_argument("--vl", help="VL sequence (amino acids)")
    parser.add_argument("--db-dir", help="Germline database directory")
    parser.add_argument("--top", type=int, default=10, help="Number of top germlines to show")
    args = parser.parse_args()

    if not args.vh and not args.vl and not args.input:
        parser.error("Provide --input (FASTA) or --vh/--vl sequences")

    # Load sequences
    sequences = {}
    if args.input:
        from .germline import _parse_fasta_text
        with open(args.input) as f:
            for name, seq in _parse_fasta_text(f.read()):
                # Auto-detect chain type by name or length
                name_upper = name.upper()
                if "VH" in name_upper or "HEAVY" in name_upper:
                    sequences[name] = ("H", seq)
                elif "VL" in name_upper or "LIGHT" in name_upper or "KAPPA" in name_upper or "LAMBDA" in name_upper:
                    sequences[name] = ("L", seq)
                elif len(seq) > 130:
                    sequences[name] = ("H", seq)
                else:
                    sequences[name] = ("L", seq)
    else:
        if args.vh:
            sequences["query_VH"] = ("H", args.vh)
        if args.vl:
            sequences["query_VL"] = ("L", args.vl)

    # Evaluate each sequence with both numbering schemes
    for name, (chain_type, seq) in sequences.items():
        # Kabat-based evaluation
        scored, numbered = evaluate_sequence(chain_type, seq, args.db_dir)
        result = format_results(chain_type, seq, scored, numbered, args.top)
        print(result)

        # Sequence comparison table (Kabat)
        best_gene = scored[0][0]
        if best_gene.numbered:
            table = format_sequence_comparison_table(numbered, best_gene.numbered, chain_type)
            print(table)

        # IMGT-based evaluation
        scored_imgt, numbered_imgt = evaluate_sequence_imgt(chain_type, seq, args.db_dir)
        result_imgt = format_results_imgt(chain_type, seq, scored_imgt, numbered_imgt, args.top)
        print(result_imgt)

    # Combined summary for Fab
    if "H" in [c for c, _ in sequences.values()] and "L" in [c for c, _ in sequences.values()]:
        print(f"\n{'='*70}")
        print(f"  SUMMARY (Fab humanization)")
        print(f"{'='*70}")

        print(f"\n  {'Chain':<8} {'Best Kabat Gene':<20} {'FR id':<10} {'Best IMGT Gene':<20} {'FR id':<10}")
        print(f"  {'-'*8} {'-'*20} {'-'*10} {'-'*20} {'-'*10}")

        for name, (chain_type, seq) in sequences.items():
            scored, numbered = evaluate_sequence(chain_type, seq, args.db_dir)
            scored_imgt, _ = evaluate_sequence_imgt(chain_type, seq, args.db_dir)
            best_kabat, best_scores_kabat = scored[0]
            if scored_imgt:
                best_imgt, best_scores_imgt = scored_imgt[0]
                print(f"  {chain_type+' chain':<8} {best_kabat.gene_id:<20} {best_scores_kabat['fr_identity']:<10.4f} "
                      f"{best_imgt.gene_id:<20} {best_scores_imgt['fr_identity']:<10.4f}")
            else:
                print(f"  {chain_type+' chain':<8} {best_kabat.gene_id:<20} {best_scores_kabat['fr_identity']:<10.4f} "
                      f"{'N/A':<20} {'N/A':<10}")

        print(f"\n  Note: IMGT-based evaluation is used for WHO/INN/USAN drug naming.")
        print(f"  Kabat-based evaluation is for internal pipeline use.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
