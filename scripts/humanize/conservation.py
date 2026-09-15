"""Functional conservation scoring for antibody positions.

P1 improvement: Calculate functional conservation based on:
1. Germline frequency (weighted by therapeutic antibody usage)
2. Evolutionary rate (from germline variation)
3. Structural constraint (from B-factor/plDDT, SASA)

Reference:
- Meaney et al., MAbs 2024: "Vernier, canonical, and interface residues are
  functionally constrained; mutations here are poorly tolerated"
- Anbarasu et al., J Struct Biol 2024: "Conservation analysis identifies
  structurally important residues"
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Dict, List, Optional, Tuple

if TYPE_CHECKING:
    from .numbering import NumberedChain
    from .backmut import BackMutationCandidate

from .germline_frequency import get_frequency


def calculate_functional_conservation(
    pos: str,
    chain_type: str,
    germline_db,
    donor: 'NumberedChain',
    structure_data=None,
    immunogenicity: Optional[Dict[str, float]] = None,
) -> Dict[str, float]:
    """Calculate functional conservation score for a position.
    
    Returns:
        Dictionary with:
        - germline_conservation: weighted germline conservation (0-1)
        - evolutionary_rate: estimated evolutionary rate (0-1, lower = more conserved)
        - structural_constraint: structural constraint score (0-1)
        - functional_conservation: combined score (0-1)
    """
    # 1. Germline conservation (weighted by therapeutic antibody usage)
    germline_cons = _calculate_germline_conservation(
        pos, chain_type, germline_db, donor
    )
    
    # 2. Evolutionary rate (from germline variation)
    evo_rate = _calculate_evolutionary_rate(
        pos, chain_type, germline_db, donor
    )
    
    # 3. Structural constraint (from structure data)
    struct_constraint = _calculate_structural_constraint(
        pos, chain_type, structure_data
    )
    
    # 4. Combined functional conservation
    # Weighted combination based on literature
    w_germline = 0.4   # Germline conservation
    w_evolution = 0.3  # Evolutionary rate
    w_struct = 0.3     # Structural constraint
    
    functional_cons = (
        w_germline * germline_cons +
        w_evolution * (1 - evo_rate) +  # Lower rate = higher conservation
        w_struct * struct_constraint
    )
    
    return {
        "germline_conservation": germline_cons,
        "evolutionary_rate": evo_rate,
        "structural_constraint": struct_constraint,
        "functional_conservation": functional_cons,
    }


def _calculate_germline_conservation(
    pos: str,
    chain_type: str,
    germline_db,
    donor: 'NumberedChain',
) -> float:
    """Calculate germline conservation weighted by therapeutic antibody usage.
    
    Returns:
        Conservation score (0-1): 1 = highly conserved, 0 = not conserved
    """
    # Get human germline genes for this chain type
    v_genes = germline_db.human(chain_type)
    
    if not v_genes:
        return 0.5  # Default medium conservation
    
    # Collect amino acid frequencies at this position
    aa_freq = {}
    total_freq = 0.0
    
    for gene in v_genes:
        # Get frequency for this gene (weighted by therapeutic usage)
        freq = get_frequency(chain_type, gene.gene_id) or 0.01  # Default low frequency
        
        # Get amino acid at this position
        if gene.numbered:
            posmap = gene.numbered.posmap()
            aa = posmap.get(pos, None)
            if aa and aa != '-':
                aa_freq[aa] = aa_freq.get(aa, 0.0) + freq
                total_freq += freq
    
    if total_freq == 0:
        return 0.5
    
    # Normalize frequencies
    for aa in aa_freq:
        aa_freq[aa] /= total_freq
    
    # Conservation = frequency of most common amino acid
    if aa_freq:
        max_freq = max(aa_freq.values())
        return max_freq
    
    return 0.5


def _calculate_evolutionary_rate(
    pos: str,
    chain_type: str,
    germline_db,
    donor: 'NumberedChain',
) -> float:
    """Calculate evolutionary rate from germline variation.
    
    Returns:
        Evolutionary rate (0-1): 0 = highly conserved, 1 = rapidly evolving
    """
    # Get human germline genes for this chain type
    v_genes = germline_db.human(chain_type)
    
    if not v_genes:
        return 0.5  # Default medium rate
    
    # Collect amino acid diversity at this position
    aa_count = {}
    total = 0
    
    for gene in v_genes:
        if gene.numbered:
            posmap = gene.numbered.posmap()
            aa = posmap.get(pos, None)
            if aa and aa != '-':
                aa_count[aa] = aa_count.get(aa, 0) + 1
                total += 1
    
    if total == 0:
        return 0.5
    
    # Shannon entropy as a measure of diversity
    import math
    entropy = 0.0
    for count in aa_count.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    
    # Normalize entropy (max entropy = log2(20) for 20 amino acids)
    max_entropy = math.log2(20)
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
    
    return normalized_entropy


def _calculate_structural_constraint(
    pos: str,
    chain_type: str,
    structure_data=None,
) -> float:
    """Calculate structural constraint from structure data.
    
    Returns:
        Constraint score (0-1): 1 = highly constrained, 0 = flexible
    """
    if structure_data is None:
        return 0.5  # Default medium constraint
    
    # Get structural metrics
    plddt = structure_data.plddt(pos)
    rel_sasa = structure_data.rel_sasa(pos)
    buried = structure_data.buried(chain_type, pos)
    
    # Calculate constraint from available data
    constraint = 0.5  # Default
    
    if plddt is not None:
        # High pLDDT = high confidence = high constraint
        # pLDDT ranges from 0-100, normalize to 0-1
        plddt_constraint = plddt / 100.0
        constraint = max(constraint, plddt_constraint)
    
    if rel_sasa is not None:
        # Low SASA = buried = high constraint
        # SASA ranges from 0-1, invert for constraint
        sasa_constraint = 1.0 - rel_sasa
        constraint = max(constraint, sasa_constraint)
    
    if buried is True:
        # Buried positions are highly constrained
        constraint = max(constraint, 0.8)
    elif buried is False:
        # Exposed positions are less constrained
        constraint = min(constraint, 0.6)
    
    return constraint


def annotate_positions_with_conservation(
    candidates: List['BackMutationCandidate'],
    chain_type: str,
    germline_db,
    donor: 'NumberedChain',
    structure_data=None,
) -> List['BackMutationCandidate']:
    """Annotate all candidates with functional conservation scores.
    
    Modifies candidates in place and returns them.
    """
    for candidate in candidates:
        scores = calculate_functional_conservation(
            candidate.position,
            chain_type,
            germline_db,
            donor,
            structure_data,
        )
        
        # Add conservation scores to rationale
        if scores["functional_conservation"] > 0.8:
            candidate.rationale.append(
                f"functional conservation: {scores['functional_conservation']:.2f} "
                f"(germline: {scores['germline_conservation']:.2f}, "
                f"evolution: {scores['evolutionary_rate']:.2f}, "
                f"struct: {scores['structural_constraint']:.2f})"
            )
        elif scores["functional_conservation"] < 0.3:
            candidate.rationale.append(
                f"low functional conservation ({scores['functional_conservation']:.2f})"
            )
    
    return candidates
