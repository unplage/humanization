#!/usr/bin/env python3
"""Evaluate humanization degree of antibody sequences.

Compares input sequences against human germline V genes and reports:
  - Best matching germline (by FR identity)
  - FR identity, CDR identity, overall identity (Kabat and IMGT numbering)
  - Humanization level per USAN/WHO standards (IMGT-based)
  - Number of positions differing from germline (potential back-mutation candidates)

Usage:
    python3 tools/humanness_score/evaluate_humanness.py --input seq.fasta
    python3 tools/humanness_score/evaluate_humanness.py --vh EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS --vl DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
"""

import argparse
import re
import sys
import os
from typing import Dict, List, Optional, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Import from scripts.humanize (absolute imports)
from scripts.humanize.numbering import number_heavy, number_light
from scripts.humanize.germline import load_germline_db, compare_to_germline
from scripts.humanize.anarci_adapter import number_with_anarci_imgt
from scripts.humanize.imgt_numbering import (
    load_imgt_germline, compare_imgt_posmaps_direct
)


# ---------------------------------------------------------------------------
# Tie-breaking for germline selection
# ---------------------------------------------------------------------------
# When several germlines give identical (FR, CDR) identity, the previous sort
# fell through to the germline database file order - effectively an arbitrary
# choice. Prefer well-established, druggable frameworks instead. Cascade:
#   1. curated druggable priority (lists below, in the listed order);
#   2. germline usage frequency from scripts/humanize/germline_frequency
#      (Adimab/OAS-style therapeutic usage data - a broader druggability
#      signal that also covers genes outside the curated lists);
#   3. gene id (deterministic; never order-dependent).
# Extend the lists as new clinical/developability evidence appears.
DRUGGABLE_VH_PRIORITY = ("IGHV1-69", "IGHV3-23", "IGHV1-46")
DRUGGABLE_VK_PRIORITY = ("IGKV1-39", "IGKV3-20", "IGKV4-1")


def _germline_tiebreak_key(chain_type: str, gene_id: str) -> tuple:
    """Sort key appended AFTER (FR identity, CDR identity); smaller wins.

    Only ever consulted for exact identity ties, so the primary ranking is
    unchanged.
    """
    from scripts.humanize.germline_frequency import get_frequency
    family = gene_id.split("*")[0]
    priority = DRUGGABLE_VH_PRIORITY if chain_type == "H" else DRUGGABLE_VK_PRIORITY
    try:
        rank = priority.index(family)
    except ValueError:
        rank = len(priority)
    freq = get_frequency(chain_type, gene_id) or 0.0
    return (rank, -freq, gene_id)


def _sort_scored(chain_type: str, scored: list) -> list:
    """Sort (gene, scores) pairs by FR identity, then CDR identity, then the
    druggability tie-break (see _germline_tiebreak_key)."""
    scored.sort(key=lambda t: (
        -t[1]["fr_identity"],
        -t[1]["cdr_identity"],
        *_germline_tiebreak_key(chain_type, t[0].gene_id),
    ))
    return scored


def evaluate_sequence(chain_type: str, sequence: str, db_dir: str = "",
                     numbered=None):
    """Evaluate a single chain against human germline.

    ``numbered`` may be supplied when the caller already numbered the
    sequence (e.g. from classify_sequences), avoiding a second numbering pass.
    """
    # Number the sequence unless the caller already did it
    if numbered is None:
        numbered = number_heavy(sequence) if chain_type == "H" else number_light(sequence)

    # Load germline DB
    if not db_dir:
        db_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "data", "germline")
    db = load_germline_db(db_dir)

    # Get human V genes for this chain type
    v_genes = db.human(chain_type)
    if not v_genes:
        raise ValueError(f"no human V genes found for chain type {chain_type!r}")

    # Score each germline (Kabat-based)
    scored = []
    for gene in v_genes:
        scores = compare_to_germline(numbered, gene)
        scored.append((gene, scores))

    # Sort by FR identity (primary), CDR identity (secondary), then the
    # druggability tie-break for exact ties (germline selection must be
    # deterministic and prefer developable frameworks)
    _sort_scored(chain_type, scored)

    return scored, numbered


def evaluate_sequence_imgt(chain_type: str, sequence: str, db_dir: str = "",
                          imgt_genes=None):
    """Evaluate a single chain against human germline using IMGT numbering.

    Uses ANARCI with IMGT scheme to number the query sequence directly,
    then loads the IMGT germline database for comparison. Like the main
    pipeline's ``--numbering imgt`` mode, this requires ANARCI and the
    bundled IMGT germline DB; when either is missing the function warns and
    returns ``(None, None)`` instead of emitting approximate converted
    (Kabat-to-IMGT) numbering.

    ``imgt_genes`` may be supplied by a batch caller to avoid reloading the
    IMGT germline JSON for every sequence.
    """
    # Number the sequence using ANARCI with IMGT scheme. The main pipeline
    # requires ANARCI for IMGT mode (there is no portable fallback), so the
    # tool does the same instead of emitting approximate converted numbers.
    numbered = number_with_anarci_imgt(sequence, chain_type)
    if numbered is None:
        print(
            "  [WARN] ANARCI IMGT numbering unavailable or failed; skipping "
            "IMGT scoring (install anarci for IMGT results)",
            file=sys.stderr,
        )
        return None, None

    # Use IMGT-numbered posmap directly
    query_imgt = numbered.posmap()

    # Load IMGT germline DB (abnumber_human_imgt.json). The pipeline's
    # --numbering imgt mode also requires this bundled IMGT database.
    if imgt_genes is None:
        imgt_genes = load_imgt_germline(db_dir)
        if not imgt_genes:
            print(
                "  [WARN] IMGT germline database unavailable; skipping IMGT "
                "scoring (expected data/germline/abnumber_human_imgt.json)",
                file=sys.stderr,
            )
            return None, None
    if not imgt_genes:
        return None, None

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

    # Sort by FR identity (primary), CDR identity (secondary), then the
    # druggability tie-break for exact ties (same cascade as the Kabat path)
    _sort_scored(chain_type, scored)

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


def _format_imgt_alignment(numbered, germline_posmap: Dict[str, str],
                           chain_type: str) -> str:
    """Format a true IMGT alignment using the query's numbered residue order.

    Unlike a numeric range walk, this preserves ANARCI insertion order
    (e.g. 111, 111A.., 112C..112A, 112) and therefore does not drop IMGT
    insertion residues from CDR1/CDR2/CDR3.
    """
    regions_order = ["FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4"]
    lines = []
    for region in regions_order:
        query_residues = [r for r in numbered.residues if r.region == region]
        if not query_residues:
            continue
        q_str = "".join(r.aa for r in query_residues)
        g_str = "".join(germline_posmap.get(r.pos, "-") for r in query_residues)
        diff = "".join(
            "X" if germline_posmap.get(r.pos, "-") not in ("-", r.aa) else "."
            for r in query_residues
        )
        lines.append(f"  {region:5s}: {q_str}  |  {g_str}  |  {diff}")
    return "\n".join(lines)


def format_sequence_comparison_table(numbered, germline_numbered, chain_type: str):
    """Generate a Kabat comparison table and a label-based diagnostic summary.

    Both NumberedChain inputs are Kabat-numbered (the caller passes the Kabat
    result of evaluate_sequence). True IMGT alignment is emitted separately in
    format_results_imgt() from real IMGT posmaps; this function no longer
    pretends to show IMGT by re-slicing Kabat labels.

    Args:
        numbered: NumberedChain object for query sequence (Kabat)
        germline_numbered: NumberedChain object for germline sequence (Kabat)
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

    # Canonical Kabat definitions, identical to KABAT_VH_REGIONS /
    # KABAT_VL_REGIONS in scripts/humanize/imgt_numbering.py.
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
        prefix = "L"

    # Kabat table: Kabat labels x Kabat boundaries (exact)
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

    # Raw per-region diagnostic over ALL 7 Kabat regions. This INCLUDES CDR3
    # and FR4, unlike the headline FR/CDR scores; it is compared by position
    # label (not by zipping concatenated strings) so insertion/deletion
    # offsets cannot shift the alignment. Identity uses the same common-
    # position denominator as compare_to_germline.
    lines.append(f"\n  Sequence identity summary (all regions, incl. CDR3+FR4):")
    lines.append(f"  {'Region':<10} {'Identity':<10} {'Match':<8} {'Compared':<10} {'QueryOnly':<11} {'GermOnly':<10}")
    lines.append(f"  {'-'*10} {'-'*10} {'-'*8} {'-'*10} {'-'*11} {'-'*10}")

    for region in kabat_regions:
        q_positions = [p for p in q_map if numbered.region_of(p) == region]
        g_positions = [p for p in g_map if germline_numbered and germline_numbered.region_of(p) == region]
        common = [p for p in q_positions if p in g_map]
        match_count = sum(1 for p in common if q_map[p] == g_map[p])
        identity = f"{match_count / len(common):.1%}" if common else "N/A"
        query_only = sum(1 for p in q_positions if p not in g_map)
        germline_only = sum(1 for p in g_positions if p not in q_map)
        lines.append(
            f"  {region:<10} {identity:<10} {match_count:<8} "
            f"{len(common):<10} {query_only:<11} {germline_only:<10}"
        )

    return "\n".join(lines)


def format_results(chain_type: str, sequence: str, scored, numbered, top_n: int = 10):
    """Format Kabat-based results for display."""
    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  {chain_type} chain ({len(sequence)} residues) - Kabat numbering")
    lines.append(f"{'='*70}")

    # Sequence with region annotation
    lines.append(f"\n  Sequence (first 50 aa): {sequence[:50]}...")

    if not scored:
        lines.append("\n  [No Kabat-scored germlines available]")
        return "\n".join(lines)

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
    lines.append(
        f"    Compared:     {best_scores.get('n_compared', 0)}/{best_scores.get('n_query', 0)} "
        f"query positions (coverage {best_scores.get('coverage', 0):.1%}); "
        f"query-only {best_scores.get('n_query_only', 0)}, "
        f"germline-only {best_scores.get('n_germline_only', 0)}"
    )

    # Detailed region breakdown
    lines.append(f"\n  Detailed region breakdown (Kabat):")
    lines.append(f"  {'Region':<10} {'Identity':<12} {'Match':<8} {'Total':<8}")
    lines.append(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*8}")
    for region, key in (("FR1", "fr1"), ("FR2", "fr2"), ("FR3", "fr3"),
                        ("CDR1", "cdr1"), ("CDR2", "cdr2")):
        lines.append(
            f"  {region:<10} {best_scores.get(key + '_identity', 0):.1%}     "
            f"{best_scores.get(key + '_match', 0):<8} {best_scores.get('n_' + key, 0):<8}"
        )

    lines.append(f"\n    Humanization estimate: {classify_humanization_level(best_scores['fr_identity'])}")
    lines.append("    (Academic estimate only; official USAN/INN naming requires manufacturing-process knowledge)")

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
    lines.append(
        f"    Compared:     {best_scores.get('n_compared', 0)}/{best_scores.get('n_query', 0)} "
        f"query positions (coverage {best_scores.get('coverage', 0):.1%}); "
        f"query-only {best_scores.get('n_query_only', 0)}, "
        f"germline-only {best_scores.get('n_germline_only', 0)}"
    )

    # IMGT region breakdown
    if "imgt_region_stats" in best_scores:
        lines.append(f"\n  IMGT region breakdown:")
        lines.append(f"  {'Region':<10} {'Identity':<12} {'Match':<8} {'Total':<8}")
        lines.append(f"  {'-'*10} {'-'*12} {'-'*8} {'-'*8}")
        # FR1+CDR1+FR2+CDR2+FR3 only (IMGT_EVAL_REGIONS); CDR3 and FR4 are
        # excluded here for consistency with the headline fr/cdr identities
        # above (CDR3 = donor-grafted, FR4 = J-derived).
        for region in ["FR1", "CDR1", "FR2", "CDR2", "FR3"]:
            stats = best_scores["imgt_region_stats"].get(region, {})
            identity = stats.get("identity", 0)
            match = stats.get("match", 0)
            count = stats.get("count", 0)
            lines.append(f"  {region:<10} {identity:.1%}     {match:<8} {count:<8}")

    lines.append(f"\n    Humanization estimate: {classify_humanization_level(best_scores['fr_identity'])}")
    lines.append("    (Academic estimate only; official USAN/INN naming requires manufacturing-process knowledge)")

    # Region alignment (true IMGT posmaps; preserves insertion order)
    germline_posmap = getattr(best_gene, "_imgt_posmap", None)
    if numbered is not None and germline_posmap:
        lines.append(f"\n  Region alignment vs {best_gene.gene_id} (IMGT):")
        lines.append(f"  (X = differs from germline, . = matches)")
        lines.append(_format_imgt_alignment(numbered, germline_posmap, chain_type))

    return "\n".join(lines)


def load_sequences(input_path: Optional[str] = None,
                   vh: Optional[str] = None,
                   vl: Optional[str] = None,
                   skip_invalid: bool = False):
    """Load and classify input sequences.

    Returns (sequences, errors), where each sequence is
    (name, chain_type, sequence, numbered_or_None). Duplicate FASTA names are
    made unique. With ``skip_invalid=True`` a bad record is reported in
    ``errors`` and skipped; otherwise the first bad record raises.
    """
    from scripts.humanize.sequences import parse_fasta, classify_sequences

    sequences: List[Tuple[str, str, str, object]] = []
    errors: List[Tuple[str, str]] = []
    seen: Dict[str, int] = {}

    def _add(name: str, chain_type: str, sequence: str, numbered) -> None:
        if name in seen:
            seen[name] += 1
            name = f"{name}_{seen[name]}"
        else:
            seen[name] = 0
        sequences.append((name, chain_type, sequence, numbered))

    def _clean(seq: str) -> str:
        return re.sub(r"[^A-Za-z]", "", seq or "").upper()

    if input_path:
        records = parse_fasta(input_path)
        if not records:
            raise ValueError(f"no sequences found in {input_path}")
        for raw_name, raw_seq in records:
            try:
                chain = classify_sequences([(raw_name, raw_seq)], numbering="kabat")[0]
            except ValueError as exc:
                message = str(exc) or exc.__class__.__name__
                errors.append((raw_name, message))
                if skip_invalid:
                    continue
                raise
            _add(chain.name, chain.chain_type, chain.sequence, chain.numbered)
    else:
        for name, raw_seq in (("query_VH", vh), ("query_VL", vl)):
            if not raw_seq:
                continue
            seq = _clean(raw_seq)
            if not seq:
                raise ValueError(f"{name} sequence is empty")
            _add(name, "H" if name == "query_VH" else "L", seq, None)

    return sequences, errors


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

    try:
        sequences, errors = load_sequences(
            input_path=args.input, vh=args.vh, vl=args.vl, skip_invalid=True)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    for name, message in errors:
        print(f"  [ERROR] {name}: {message}", file=sys.stderr)
    if not sequences:
        print("ERROR: no valid antibody sequences to evaluate", file=sys.stderr)
        return 1

    # Check ANARCI availability (required for IMGT)
    try:
        from scripts.humanize.anarci_adapter import number_with_anarci_imgt
        test_result = number_with_anarci_imgt("EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARMDYWGQGTLVTVSS", "H")
        if test_result is None:
            print("ERROR: ANARCI IMGT numbering failed. Please install ANARCI for IMGT evaluation:", file=sys.stderr)
            print("  pip install anarci", file=sys.stderr)
            return 1
    except ImportError:
        print("ERROR: ANARCI not installed. Required for IMGT numbering (WHO/INN/USAN standard).", file=sys.stderr)
        print("  Install with: pip install anarci", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"ERROR: ANARCI IMGT numbering failed: {e}", file=sys.stderr)
        return 1

    # Load IMGT germline database
    imgt_genes = load_imgt_germline(args.db_dir or "")
    if not imgt_genes:
        print(
            "ERROR: IMGT germline database unavailable. Required for IMGT evaluation.",
            file=sys.stderr,
        )
        print("  Expected: data/germline/abnumber_human_imgt.json", file=sys.stderr)
        return 1

    print("  [INFO] IMGT evaluation enabled (ANARCI + IMGT germline DB)")

    evaluated = []
    for name, chain_type, seq, numbered in sequences:
        try:
            scored, numbered = evaluate_sequence(
                chain_type, seq, args.db_dir, numbered=numbered)
        except ValueError as exc:
            print(f"  [ERROR] {name}: {exc}", file=sys.stderr)
            errors.append((name, str(exc)))
            continue

        print(format_results(chain_type, seq, scored, numbered, args.top))

        best_gene = scored[0][0]
        if best_gene.numbered:
            print(format_sequence_comparison_table(
                numbered, best_gene.numbered, chain_type))

        scored_imgt, numbered_imgt = evaluate_sequence_imgt(
            chain_type, seq, args.db_dir, imgt_genes=imgt_genes)
        if not scored_imgt or not numbered_imgt:
            print(f"  [ERROR] {name}: IMGT evaluation failed - ANARCI or IMGT germline DB unavailable", file=sys.stderr)
            errors.append((name, "IMGT evaluation failed"))
            continue

        print(format_results_imgt(
            chain_type, seq, scored_imgt, numbered_imgt, args.top))

        evaluated.append((name, chain_type, seq, scored, numbered, scored_imgt))

    # Combined summary for Fab
    if any(c == "H" for _, c, _, _, _, _ in evaluated) and \
            any(c == "L" for _, c, _, _, _, _ in evaluated):
        print(f"\n{'='*70}")
        print(f"  SUMMARY (Fab humanization)")
        print(f"{'='*70}")
        print(f"\n  {'Chain':<8} {'Best Kabat Gene':<20} {'FR id':<10} {'Best IMGT Gene':<20} {'FR id':<10}")
        print(f"  {'-'*8} {'-'*20} {'-'*10} {'-'*20} {'-'*10}")

        for name, chain_type, seq, scored, numbered, scored_imgt in evaluated:
            best_kabat, best_scores_kabat = scored[0]
            if scored_imgt:
                best_imgt, best_scores_imgt = scored_imgt[0]
                print(
                    f"  {chain_type+' chain':<8} {best_kabat.gene_id:<20} "
                    f"{best_scores_kabat['fr_identity']:<10.4f} "
                    f"{best_imgt.gene_id:<20} {best_scores_imgt['fr_identity']:<10.4f}"
                )
            else:
                print(
                    f"  {chain_type+' chain':<8} {best_kabat.gene_id:<20} "
                    f"{best_scores_kabat['fr_identity']:<10.4f} {'N/A':<20} {'N/A':<10}"
                )

        print(f"\n  Note: IMGT-based evaluation is used for WHO/INN/USAN drug naming.")
        print(f"  Kabat-based evaluation is for internal pipeline use.")

    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
