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
    """Evaluate a single chain against human germline using IMGT numbering."""
    from .numbering import number_heavy, number_light
    from .germline import load_germline_db, compare_to_germline
    from .imgt_numbering import compare_to_germline_imgt

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

    # Score each germline using IMGT-based comparison
    scored = []
    for gene in v_genes:
        if gene.numbered is None:
            continue
        scores = compare_to_germline_imgt(
            numbered.posmap(),
            gene.numbered.posmap(),
            chain_type,
        )
        scored.append((gene, scores))

    # Sort by FR identity (primary) and CDR identity (secondary)
    scored.sort(key=lambda t: (-t[1]["fr_identity"], -t[1]["cdr_identity"]))

    return scored, numbered


def classify_humanization_level(fr_identity: float) -> str:
    """Classify humanization level based on FR identity."""
    if fr_identity >= 0.95:
        return "Very High (>95%) - Human (-umab)"
    elif fr_identity >= 0.90:
        return "High (90-95%) - Humanized (-zumab)"
    elif fr_identity >= 0.85:
        return "Moderate-High (85-90%) - Humanized (-zumab)"
    elif fr_identity >= 0.80:
        return "Moderate (80-85%) - Humanized (-zumab) [borderline]"
    elif fr_identity >= 0.70:
        return "Low (70-80%) - Chimeric (-xi-)"
    elif fr_identity >= 0.60:
        return "Very Low (60-70%) - Chimeric (-xi-)"
    else:
        return "Extremely Low (<60%) - Murine (-o-)"


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

    # Humanization classification
    fr_id = best_scores["fr_identity"]
    level = classify_humanization_level(fr_id)
    lines.append(f"    Humanization level: {level}")

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

    # Humanization classification (USAN standard)
    fr_id = best_scores["fr_identity"]
    if fr_id >= 0.85:
        usan_class = "Humanized (-zumab)"
    elif fr_id >= 0.70:
        usan_class = "Chimeric (-xi-)"
    else:
        usan_class = "Murine (-o-)"
    lines.append(f"    USAN naming class: {usan_class} (FR identity >= 85% = humanized)")

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

        # IMGT-based evaluation
        scored_imgt, numbered = evaluate_sequence_imgt(chain_type, seq, args.db_dir)
        result_imgt = format_results_imgt(chain_type, seq, scored_imgt, numbered, args.top)
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
