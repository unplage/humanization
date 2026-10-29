#!/usr/bin/env python3
"""Structural risk adjustment for immunogenicity prediction.

Two-step analysis workflow:
1. Sequence-based analysis (HLAIIPred/heuristic) → initial risk scores
2. Structure-based confirmation (AF3/IgFold) → adjust risk based on buriedness

Key insight: Buried residues are less accessible to MHC-II antigen processing,
so sequence-based predictions may overestimate risk for buried positions.

Reference:
- Tien et al. (2013) consensus 3D protein structure
- FreeSASA: Beer et al. (2016)
"""
from __future__ import annotations

import csv
import json
import os
import re
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Try to import freesasa
try:
    import freesasa
    HAS_FREESASA = True
except ImportError:
    HAS_FREESASA = False


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Reference SASA values: Miller et al. 1987 scale (A^2), often miscopied
# under a Tien et al. 2013 citation (Tien's recommended values are larger,
# e.g. ALA 129 / ASP 193; cf. scripts/humanize/structure.py). NOTE: this
# constant is currently unused (dead code).
REFERENCE_SASA = {
    'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0, 'CYS': 140.0,
    'GLU': 183.0, 'GLN': 179.0, 'GLY': 85.0, 'HIS': 194.0, 'ILE': 182.0,
    'LEU': 180.0, 'LYS': 211.0, 'MET': 203.0, 'PHE': 218.0, 'PRO': 143.0,
    'SER': 122.0, 'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0,
}

# Buriedness thresholds (single source of truth)
BURIED_THRESHOLD = 0.15    # relSASA < 0.15 = buried
EXPOSED_THRESHOLD = 0.25   # relSASA > 0.25 = exposed
# 0.15 <= relSASA <= 0.25 = uncertain zone

# Risk adjustment factors for the binary method
RISK_ADJUSTMENT = {
    'buried': 0.10,
    'exposed': 1.0,
    'uncertain': 0.5,
    'no_structure': 1.0,   # neutral: never silently halve risk
}


def _kabat_num(pos: str) -> Optional[int]:
    """Numeric part of a Kabat label, tolerant of insertion letters.

    "H5" -> 5, "H82A" -> 82, "L100B" -> 100, "H103" -> 103. Returns None when
    no digits are present. Kabat insertion letters (82A/82B/82C, 100A..K,
    27A..) are common in FR3/CDR loops; a plain int() would raise on them.
    """
    m = re.search(r"(\d+)", pos or "")
    return int(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PositionBuriedness:
    """Buriedness data for a single residue position."""
    position: str           # e.g. "H5"
    chain: str              # "H" or "L"
    seq_index: int          # 0-based index in sequence
    resname: str            # 3-letter residue name
    abs_sasa: float         # absolute SASA (A^2)
    rel_sasa: float         # relative SASA (0-1)
    buried: Optional[bool]  # True=buried, False=exposed, None=uncertain
    is_backmutation: bool   # True if this is a backmutation position
    donor_aa: Optional[str] = None
    human_aa: Optional[str] = None
    tier: Optional[str] = None  # T1/T2/T3
    cdr_contact: Optional[bool] = None
    antigen_contact: Optional[bool] = None

    @property
    def risk_class(self) -> str:
        """Risk classification based on buriedness."""
        if self.buried is True:
            return 'buried'
        elif self.buried is False:
            return 'exposed'
        else:
            return 'uncertain'

    @property
    def risk_factor(self) -> float:
        """Risk adjustment factor (0-1)."""
        return RISK_ADJUSTMENT.get(self.risk_class, 0.5)


@dataclass
class PeptideRiskAssessment:
    """Structural risk assessment for a single peptide."""
    peptide: str
    start_index: int        # 0-based start position in sequence
    end_index: int          # 0-based end position (exclusive)
    raw_score: float        # HLAIIPred/heuristic score
    positions: List[PositionBuriedness] = field(default_factory=list)
    adjustment_factor: float = 1.0
    adjusted_score: float = 0.0
    n_buried: int = 0
    n_exposed: int = 0
    n_uncertain: int = 0
    n_backmutation: int = 0
    rationale: List[str] = field(default_factory=list)

    def __post_init__(self):
        self.adjusted_score = self.raw_score * self.adjustment_factor


@dataclass
class VariantStructuralRisk:
    """Structural risk assessment for an entire variant."""
    variant_name: str
    chain_type: str
    sequence: str
    peptides: List[PeptideRiskAssessment] = field(default_factory=list)
    overall_adjustment: float = 1.0
    adjusted_mean: Optional[float] = None
    adjusted_immunogenicity_score: Optional[float] = None
    structure_ok: bool = True
    error: Optional[str] = None
    structure_source: str = ""  # e.g. "V2.pdb", "donor.pdb"


# ---------------------------------------------------------------------------
# SASA computation
# ---------------------------------------------------------------------------

_AA3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'SEC': 'U', 'PYL': 'O',
}


def _resseq_num(value) -> Optional[int]:
    m = re.search(r"(-?\d+)", str(value or ""))
    return int(m.group(1)) if m else None


def _ordered_chain(chain_areas) -> List[Tuple[int, object]]:
    ordered = []
    for resseq_str, area_obj in chain_areas.items():
        n = _resseq_num(resseq_str)
        if n is not None:
            ordered.append((n, area_obj))
    ordered.sort(key=lambda item: item[0])
    return ordered


def _chain_sequence(chain_areas) -> str:
    return "".join(
        _AA3.get(str(getattr(obj, "residueType", "")).upper(), "X")
        for _, obj in _ordered_chain(chain_areas)
    )


def _identity(a: str, b: str) -> float:
    if not a or not b:
        return 0.0
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / n


def _pick_chain(residue_areas, chain_id: str,
                expected_sequence: Optional[str]) -> str:
    if chain_id in residue_areas:
        return chain_id
    if not residue_areas:
        raise ValueError("structure contains no protein chains")
    if expected_sequence:
        scored = [
            (_identity(_chain_sequence(areas), expected_sequence.upper()), name)
            for name, areas in residue_areas.items()
        ]
        scored.sort(reverse=True)
        if scored and scored[0][0] >= 0.5:
            return scored[0][1]
    if len(residue_areas) == 1:
        return next(iter(residue_areas))
    raise ValueError(
        f"chain '{chain_id}' not found and no chain matches the variant sequence")


def compute_rel_sasa(pdb_path: str, chain_id: str = 'H',
                     expected_sequence: Optional[str] = None
                     ) -> Dict[str, PositionBuriedness]:
    """Compute relative SASA for all residues of one chain.

    If ``chain_id`` is absent (common for AF3/Boltz A/B chains), the chain
    whose sequence best matches ``expected_sequence`` is selected. A chain is
    never silently treated as all-missing.
    """
    if not HAS_FREESASA:
        raise ImportError("FreeSASA is required. Install with: pip install freesasa")

    import freesasa as _freesasa
    structure = _freesasa.Structure(pdb_path)
    result = _freesasa.calc(structure)

    residue_areas = result.residueAreas()
    chain_id = _pick_chain(residue_areas, chain_id, expected_sequence)
    chain_areas = residue_areas.get(chain_id, {})
    observed = _chain_sequence(chain_areas)
    if expected_sequence and observed:
        ident = _identity(observed, expected_sequence.upper())
        if ident < 0.80:
            raise ValueError(
                f"chain {chain_id} sequence identity {ident:.2f} is too low "
                f"for this variant; refusing to use a mismatched structure")

    rel_sasa = {}
    for seq_index, (resseq, area_obj) in enumerate(_ordered_chain(chain_areas)):
        total_area = area_obj.total
        rel_total = area_obj.relativeTotal
        resname = area_obj.residueType

        if rel_total < BURIED_THRESHOLD:
            buried = True
        elif rel_total > EXPOSED_THRESHOLD:
            buried = False
        else:
            buried = None

        position = f'{chain_id}{resseq}'
        rel_sasa[position] = PositionBuriedness(
            position=position,
            chain=chain_id,
            seq_index=seq_index,
            resname=resname,
            abs_sasa=total_area,
            rel_sasa=rel_total,
            buried=buried,
            is_backmutation=False,
        )

    return rel_sasa


def compute_all_chains_sasa(pdb_path: str) -> Dict[str, Dict[str, PositionBuriedness]]:
    """Compute relSASA for all chains in a PDB file.

    Returns:
        Dict mapping chain_id -> {position -> PositionBuriedness}
    """
    if not HAS_FREESASA:
        raise ImportError("FreeSASA is required")

    import freesasa as _freesasa
    structure = _freesasa.Structure(pdb_path)
    result = _freesasa.calc(structure)

    residue_areas = result.residueAreas()
    all_chains = {}

    for chain_id, chain_areas in residue_areas.items():
        chain_data = {}
        for seq_index, (resseq, area_obj) in enumerate(_ordered_chain(chain_areas)):
            total_area = area_obj.total
            rel_total = area_obj.relativeTotal
            resname = area_obj.residueType

            if rel_total < BURIED_THRESHOLD:
                buried = True
            elif rel_total > EXPOSED_THRESHOLD:
                buried = False
            else:
                buried = None

            position = f'{chain_id}{resseq}'
            chain_data[position] = PositionBuriedness(
                position=position,
                chain=chain_id,
                seq_index=seq_index,
                resname=resname,
                abs_sasa=total_area,
                rel_sasa=rel_total,
                buried=buried,
                is_backmutation=False,
            )
        all_chains[chain_id] = chain_data

    return all_chains


# ---------------------------------------------------------------------------
# Position mapping
# ---------------------------------------------------------------------------

def map_peptide_positions(
    peptide_start: int,
    peptide_end: int,
    sequence: str,
    chain_id: str,
    structure_data: Dict[str, PositionBuriedness],
    backmutations: Optional[List[Dict]] = None,
    position_labels: Optional[List[str]] = None,
) -> List[PositionBuriedness]:
    """Map a peptide range to its constituent positions with buriedness data.

    Args:
        peptide_start: 0-based start index of peptide in sequence
        peptide_end: 0-based end index (exclusive)
        sequence: full chain sequence
        chain_id: chain identifier ("H" or "L")
        structure_data: dict from compute_rel_sasa
        backmutations: list of backmutation dicts from pipeline output
        position_labels: optional per-residue labels aligned to ``sequence``
            (index 0 -> label of the first residue), in the run's coordinate
            space (Kabat or IMGT). When given, backmutations are matched by
            their exact label, which is the only way positions that carry
            insertion letters (Kabat FR3 82A-C / CDR 100A-K / 27A-, IMGT
            111A-112F) or FR4 can be matched. Without it the tool falls back
            to a sequential index, which is approximate whenever the position
            number != sequence index.

    Returns:
        List of PositionBuriedness for positions within the peptide
    """
    # Build backmutation lookups: exact label, plus numeric fallback.
    backmut_by_label: Dict[str, Dict] = {}
    backmut_by_num: Dict[int, Dict] = {}
    if backmutations:
        for bm in backmutations:
            pos_str = bm.get('position', '') or ''
            if not pos_str.startswith(chain_id):
                continue
            backmut_by_label[pos_str] = bm
            n = _kabat_num(pos_str)
            if n is not None:
                backmut_by_num[n] = bm

    def _label(seq_idx: int) -> str:
        if position_labels is not None and seq_idx < len(position_labels):
            return position_labels[seq_idx]
        return f'{chain_id}{seq_idx + 1}'

    positions = []
    for seq_idx in range(peptide_start, min(peptide_end, len(sequence))):
        position = _label(seq_idx)
        pos_num = _kabat_num(position)

        # Get structure data: by label first (re-keyed structures), then by
        # the raw chain+index key produced by compute_rel_sasa.
        struct_pos = structure_data.get(position)
        if struct_pos is None:
            struct_pos = structure_data.get(f'{chain_id}{seq_idx + 1}')

        # Check if this is a backmutation (exact label, then numeric fallback)
        bm_info = backmut_by_label.get(position)
        if bm_info is None and pos_num is not None:
            bm_info = backmut_by_num.get(pos_num)

        if struct_pos:
            # Update with backmutation info
            struct_pos.is_backmutation = bm_info is not None
            if bm_info:
                struct_pos.donor_aa = bm_info.get('donor_aa')
                struct_pos.human_aa = bm_info.get('human_aa')
                struct_pos.tier = bm_info.get('tier')
                struct_pos.cdr_contact = bm_info.get('cdr_contact')
                struct_pos.antigen_contact = bm_info.get('antigen_contact')
            positions.append(struct_pos)
        else:
            # No structure data for this position
            positions.append(PositionBuriedness(
                position=position,
                chain=chain_id,
                seq_index=seq_idx,
                resname='',
                abs_sasa=0.0,
                rel_sasa=0.0,
                buried=None,
                is_backmutation=bm_info is not None,
                donor_aa=bm_info.get('donor_aa') if bm_info else None,
                human_aa=bm_info.get('human_aa') if bm_info else None,
                tier=bm_info.get('tier') if bm_info else None,
                cdr_contact=bm_info.get('cdr_contact') if bm_info else None,
                antigen_contact=bm_info.get('antigen_contact') if bm_info else None,
            ))

    return positions


# ---------------------------------------------------------------------------
# Risk adjustment calculation
# ---------------------------------------------------------------------------

def _rel_sasa_weight(rel_sasa: float) -> float:
    if rel_sasa < 0.05:
        return 0.05
    if rel_sasa < 0.15:
        return 0.05 + 0.25 * (rel_sasa - 0.05) / 0.10
    if rel_sasa < 0.25:
        return 0.30 + 0.40 * (rel_sasa - 0.15) / 0.10
    if rel_sasa < 0.40:
        return 0.70 + 0.20 * (rel_sasa - 0.25) / 0.15
    return 1.0


def calculate_adjustment_factor(
    positions: List[PositionBuriedness],
    method: str = "relSASA",
    core_start: Optional[int] = None,
    core_end: Optional[int] = None,
) -> Tuple[float, List[str]]:
    """Calculate a peptide-level structural adjustment factor.

    The factor is a weighted average over all residues of the peptide, with
    core residues weighted 1.0 and flanking residues 0.5. Back-mutation tier is
    never used as a structural weight: tiers describe reversion importance, not
    the reliability of a buriedness measurement.
    """
    rationale = []
    if not positions:
        rationale.append("No positions available for structural assessment")
        return 1.0, rationale

    weighted_sum = 0.0
    weight_sum = 0.0
    n_missing = 0

    for pos in positions:
        weight = 1.0
        if core_start is not None and core_end is not None:
            weight = 1.0 if core_start <= pos.seq_index < core_end else 0.5

        if pos.rel_sasa > 0:
            struct_weight = _rel_sasa_weight(pos.rel_sasa)
        elif pos.buried is True:
            struct_weight = RISK_ADJUSTMENT['buried']
        elif pos.buried is False:
            struct_weight = RISK_ADJUSTMENT['exposed']
        else:
            struct_weight = RISK_ADJUSTMENT['no_structure']
            n_missing += 1

        if method == "binary":
            if pos.buried is True:
                struct_weight = RISK_ADJUSTMENT['buried']
            elif pos.buried is False:
                struct_weight = RISK_ADJUSTMENT['exposed']
            else:
                struct_weight = RISK_ADJUSTMENT['no_structure']

        if pos.cdr_contact and pos.rel_sasa > 0:
            struct_weight = max(struct_weight, 0.6)
            rationale.append(f"{pos.position}: CDR contact (structurally relevant)")

        weighted_sum += weight * struct_weight
        weight_sum += weight
        rationale.append(
            f"{pos.position}: relSASA={pos.rel_sasa:.3f}, "
            f"weight={struct_weight:.3f}, core={weight:.1f}"
        )

    if method not in ("relSASA", "binary"):
        raise ValueError(f"Unknown method: {method}")
    if weight_sum == 0:
        return 1.0, rationale

    factor = weighted_sum / weight_sum
    if n_missing == len(positions):
        rationale.append("No structure data for any peptide position; no adjustment applied")
        return 1.0, rationale
    return factor, rationale


# ---------------------------------------------------------------------------
# Main structural risk assessment
# ---------------------------------------------------------------------------

def _rekey_structure_by_label(
    structure_data: Dict[str, PositionBuriedness],
    chain_id: str,
    position_labels: List[str],
) -> Dict[str, PositionBuriedness]:
    """Re-key a chain+resseq structure map onto position labels (Kabat/IMGT).

    Residues are matched in resseq order to ``position_labels`` (index 0 =
    first residue), so this works whether the PDB is sequentially numbered
    (IgFold/AF3) or already carries another numbering. Missing residues shift
    the alignment only after the first gap.
    """
    def _resseq(key: str) -> int:
        m = re.search(r"(\d+)$", key or "")
        return int(m.group(1)) if m else 0

    ordered = sorted(structure_data.items(), key=lambda kv: _resseq(kv[0]))
    if len(ordered) != len(position_labels):
        print(
            f"  WARNING: structure has {len(ordered)} residues but "
            f"{len(position_labels)} labels; label alignment may be shifted",
            file=sys.stderr,
        )
    rekeyed: Dict[str, PositionBuriedness] = {}
    for (key, val), label in zip(ordered, position_labels):
        val.position = label
        rekeyed[label] = val
    return rekeyed


def filter_backmutations_for_variant(
    backmutations: Optional[List[Dict]],
    sequence: str,
    position_labels: Optional[List[str]],
    chain_type: str,
) -> List[Dict]:
    """Keep only candidate positions actually reverted in this variant."""
    if not backmutations:
        return []
    label_to_idx = {label: i for i, label in enumerate(position_labels or [])}
    out: List[Dict] = []
    for bm in backmutations:
        pos = str(bm.get('position', ''))
        if not pos.startswith(chain_type):
            continue
        idx = label_to_idx.get(pos)
        if idx is None and not label_to_idx:
            num = _kabat_num(pos)
            idx = num - 1 if num is not None else None
        if idx is None or not (0 <= idx < len(sequence)):
            continue
        donor = str(bm.get('donor_aa', '')).upper()
        human = str(bm.get('human_aa', '')).upper()
        if donor and donor != human and sequence[idx].upper() == donor:
            out.append(bm)
    return out


def assess_variant_structural_risk(
    variant_name: str,
    chain_type: str,
    sequence: str,
    per_peptide_hits: List[Dict],
    structure_path: str,
    backmutations: Optional[List[Dict]] = None,
    method: str = "relSASA",
    position_labels: Optional[List[str]] = None,
) -> VariantStructuralRisk:
    """Assess one variant against its own structure."""
    structure_data = compute_rel_sasa(
        structure_path, chain_type, expected_sequence=sequence)
    if position_labels:
        structure_data = _rekey_structure_by_label(
            structure_data, chain_type, position_labels)

    variant_backmuts = filter_backmutations_for_variant(
        backmutations, sequence, position_labels, chain_type)

    peptides = []
    search_cursor = 0
    for hit in per_peptide_hits:
        peptide = hit.get('peptide', '')
        score = float(hit.get('score', 0.0) or 0.0)
        if 'start_index' in hit and hit.get('start_index') is not None:
            start_idx = int(hit.get('start_index'))
        else:
            found = sequence.find(peptide, search_cursor)
            if found < 0:
                found = sequence.find(peptide)
            start_idx = found if found >= 0 else int(hit.get('core_pos', 0) or 0)
            search_cursor = max(search_cursor, start_idx + 1)

        core_pos = int(hit.get('core_pos', 0) or 0)
        core_start = start_idx + core_pos
        core_end = min(len(sequence), core_start + 9)
        end_idx = min(len(sequence), start_idx + len(peptide))

        positions = map_peptide_positions(
            start_idx, end_idx, sequence, chain_type,
            structure_data, variant_backmuts, position_labels=position_labels
        )
        factor, rationale = calculate_adjustment_factor(
            positions, method, core_start=core_start, core_end=core_end)
        adjusted_score = score * factor

        n_buried = sum(1 for p in positions if p.buried is True)
        n_exposed = sum(1 for p in positions if p.buried is False)
        n_uncertain = sum(1 for p in positions if p.buried is None)
        n_backmutation = sum(1 for p in positions if p.is_backmutation)

        peptides.append(PeptideRiskAssessment(
            peptide=peptide,
            start_index=start_idx,
            end_index=end_idx,
            raw_score=score,
            positions=positions,
            adjustment_factor=factor,
            adjusted_score=adjusted_score,
            n_buried=n_buried,
            n_exposed=n_exposed,
            n_uncertain=n_uncertain,
            n_backmutation=n_backmutation,
            rationale=rationale,
        ))

    total_raw = sum(p.raw_score for p in peptides)
    total_adjusted = sum(p.adjusted_score for p in peptides)
    overall_adjustment = total_adjusted / total_raw if total_raw > 0 else 1.0
    adjusted_mean = (total_adjusted / len(peptides)) if peptides else None

    return VariantStructuralRisk(
        variant_name=variant_name,
        chain_type=chain_type,
        sequence=sequence,
        peptides=peptides,
        overall_adjustment=overall_adjustment,
        adjusted_mean=adjusted_mean,
        structure_source=os.path.basename(structure_path),
    )


def load_numbering_json(path: str) -> Dict[str, List[str]]:
    """Load per-variant position labels: ``{variant_name: [label, ...]}``.

    Written by the pipeline as ``variants_numbering.json`` (next to
    ``variants.fasta``). Each list is aligned to that variant chain's sequence,
    so index 0 is the first residue. Enables exact FR4/backmutation matching.

    Two shapes are accepted:
      * legacy flat ``{variant: [label, ...]}`` (Kabat runs);
      * ``{"numbering_scheme": ..., "variants": {variant: [label, ...]}}``.
    Labels share the coordinate space of the run (Kabat or IMGT); matching is
    exact-by-label, so back-mutation CSV positions and these labels always
    agree within one run.
    """
    if not path or not os.path.isfile(path):
        return {}
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    if isinstance(data.get("variants"), dict):
        data = data["variants"]
    return {k: v for k, v in data.items() if isinstance(v, list)}


def _labels_for_variant(
    numbering: Dict[str, List[str]], vname: str,
    chain_type: Optional[str] = None,
) -> Optional[List[str]]:
    """Resolve a variant's label list.

    Tolerates plain keys (``combo_V0_V1``), chain-aware keys
    (``combo_V0_V1|H``), and ``name|variant`` style keys. When ``chain_type``
    is given, chain-qualified keys are preferred so the H and L chains of a
    paired variant resolve to their own numbering.
    """
    if not numbering:
        return None
    exact = []
    if chain_type:
        exact.append(f"{vname}|{chain_type}")
    exact.append(vname)
    for k in exact:
        if k in numbering:
            return numbering[k]
    for key, labels in numbering.items():
        head, _, tail = key.partition("|")
        if chain_type and tail and tail != chain_type:
            continue
        if key.split("|")[-1].strip() == vname or head.strip() == vname:
            return labels
    for key, labels in numbering.items():
        tail = key.split("|")[-1].strip()
        if tail.endswith(vname) or vname.endswith(tail):
            return labels
    return None


def load_backmutations_from_csv(csv_path: str) -> List[Dict]:
    """Load backmutation data from pipeline CSV output."""
    backmutations = []
    if not os.path.isfile(csv_path):
        return backmutations

    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                composite = float(row.get('composite') or 0)
            except (TypeError, ValueError):
                composite = 0.0
            backmutations.append({
                'position': row.get('position', ''),
                'donor_aa': row.get('donor_aa', ''),
                'human_aa': row.get('human_aa', ''),
                'tier': row.get('tier', ''),
                'buried': row.get('buried', '') == 'True',
                'cdr_contact': row.get('cdr_contact', '') == 'True',
                'antigen_contact': row.get('antigen_contact', '') == 'True',
                'composite': composite,
            })
    return backmutations


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_structural_risk_report(
    assessments: List[VariantStructuralRisk],
    output_path: str,
    title: str = "Structural Risk Assessment Report",
) -> None:
    """Generate detailed Markdown report of structural risk assessment."""
    lines = [
        f"# {title}",
        "",
        "## Overview",
        "",
        "This report applies structural risk adjustment to immunogenicity predictions.",
        "",
        "**Two-step workflow:**",
        "1. **Sequence-based analysis**: HLAIIPred/heuristic predicts MHC-II binding from sequence",
        "2. **Structure-based confirmation**: FreeSASA computes relSASA from actual variant structure",
        "",
        "**Key principle**: Buried residues (relSASA < 0.15) are less accessible to MHC-II",
        "antigen processing, so sequence-based predictions may overestimate risk.",
        "",
        "The adjusted mean is used to rewrite the final immunogenicity score in the "
        "main report (relative to the adjusted reference mean).",
        "",
        "---",
        "",
    ]

    for assessment in assessments:
        lines.append(f"## {assessment.variant_name} ({assessment.chain_type} chain)")
        lines.append("")
        lines.append(f"- **Structure**: {assessment.structure_source}")
        lines.append(f"- **Overall adjustment**: {assessment.overall_adjustment:.3f}")
        if assessment.adjusted_mean is not None:
            lines.append(f"- **Adjusted mean epitope load**: {assessment.adjusted_mean:.3f}")
        if assessment.adjusted_immunogenicity_score is not None:
            lines.append(
                f"- **Final structural immunogenicity score**: "
                f"{assessment.adjusted_immunogenicity_score:.1f}%")
        if assessment.error:
            lines.append(f"- **Error**: {assessment.error}")
        lines.append("")

        # Peptide summary table
        lines.append("### Peptide Risk Summary")
        lines.append("")
        lines.append("| Peptide | Raw Score | Adj Factor | Adj Score | Buried | Exposed | Risk Level |")
        lines.append("|---------|-----------|------------|-----------|--------|---------|------------|")

        for p in sorted(assessment.peptides, key=lambda x: x.adjusted_score, reverse=True):
            if p.raw_score < 0.3:  # skip low-confidence peptides
                continue
            risk_level = "HIGH" if p.adjusted_score > 0.5 else ("MEDIUM" if p.adjusted_score > 0.2 else "LOW")
            lines.append(
                f"| `{p.peptide[:20]}...` | {p.raw_score:.3f} | "
                f"{p.adjustment_factor:.3f} | {p.adjusted_score:.3f} | "
                f"{p.n_buried} | {p.n_exposed} | {risk_level} |"
            )
        lines.append("")

        # Detailed position analysis for high-risk peptides
        high_risk = [p for p in assessment.peptides if p.adjusted_score > 0.3]
        if high_risk:
            lines.append("### Detailed Position Analysis")
            lines.append("")
            for p in high_risk[:5]:  # top 5
                lines.append(f"#### Peptide: `{p.peptide}`")
                lines.append(f"- Raw score: {p.raw_score:.3f}")
                lines.append(f"- Adjusted score: {p.adjusted_score:.3f}")
                lines.append("")
                lines.append("| Position | relSASA | Buried | CDR Contact | Tier | Risk Factor |")
                lines.append("|----------|---------|--------|-------------|------|-------------|")
                for pos in p.positions:
                    if pos.is_backmutation:
                        buried_str = "Yes" if pos.buried else ("No" if pos.buried is False else "Uncertain")
                        cdr_str = "Yes" if pos.cdr_contact else "No"
                        lines.append(
                            f"| {pos.position} | {pos.rel_sasa:.3f} | {buried_str} | "
                            f"{cdr_str} | {pos.tier} | {pos.risk_factor:.3f} |"
                        )
                lines.append("")

        # Rationale
        all_rationale = []
        for p in assessment.peptides:
            all_rationale.extend(p.rationale)
        if all_rationale:
            lines.append("### Rationale")
            lines.append("")
            for r in set(all_rationale):  # unique rationale
                lines.append(f"- {r}")
            lines.append("")

        lines.append("---")
        lines.append("")

    # Summary comparison
    if len(assessments) > 1:
        lines.append("## Cross-Variant Comparison")
        lines.append("")
        lines.append("| Variant | Structure | Overall Adj | Adj Mean | Final Score | High-Risk Peptides |")
        lines.append("|---------|-----------|-------------|----------|-------------|-------------------|")
        for a in assessments:
            n_high = sum(1 for p in a.peptides if p.adjusted_score > 0.5)
            adj_mean = (f"{a.adjusted_mean:.3f}"
                        if a.adjusted_mean is not None else "N/A")
            score = (f"{a.adjusted_immunogenicity_score:.1f}%"
                     if a.adjusted_immunogenicity_score is not None else "N/A")
            lines.append(
                f"| {a.variant_name} | {a.structure_source} | "
                f"{a.overall_adjustment:.3f} | {adj_mean} | {score} | {n_high} |"
            )
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))


def generate_structural_risk_word(
    assessments: List[VariantStructuralRisk],
    output_path: str,
    title: str = "Structural Risk Assessment Report",
) -> None:
    """Generate Word document report of structural risk assessment."""
    try:
        from docx import Document
        from docx.shared import Inches, Pt, RGBColor
        from docx.enum.text import WD_ALIGN_PARAGRAPH
    except ImportError:
        print("  WARNING: python-docx not installed. Skipping Word report.",
              file=sys.stderr)
        return

    doc = Document()
    heading = doc.add_heading(title, 0)
    heading.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # Overview
    doc.add_heading("Overview", level=1)
    doc.add_paragraph(
        "This report applies structural risk adjustment to immunogenicity predictions."
    )
    doc.add_paragraph(
        "Two-step workflow: (1) Sequence-based analysis (HLAIIPred/heuristic) "
        "predicts MHC-II binding; (2) Structure-based confirmation (FreeSASA) "
        "adjusts risk based on actual buriedness from variant structures."
    )

    for assessment in assessments:
        doc.add_heading(f"{assessment.variant_name} ({assessment.chain_type} chain)", level=1)
        doc.add_paragraph(f"Structure: {assessment.structure_source}")
        doc.add_paragraph(f"Overall adjustment: {assessment.overall_adjustment:.3f}")
        if assessment.adjusted_mean is not None:
            doc.add_paragraph(f"Adjusted mean epitope load: {assessment.adjusted_mean:.3f}")
        if assessment.adjusted_immunogenicity_score is not None:
            doc.add_paragraph(
                "Final structural immunogenicity score: "
                f"{assessment.adjusted_immunogenicity_score:.1f}%")
        if assessment.error:
            doc.add_paragraph(f"Error: {assessment.error}")

        # Peptide summary table
        doc.add_heading("Peptide Risk Summary", level=2)
        high_risk_peptides = [p for p in assessment.peptides if p.raw_score >= 0.3]
        if high_risk_peptides:
            table = doc.add_table(rows=1 + len(high_risk_peptides), cols=7)
            table.style = "Table Grid"
            headers = ["Peptide", "Raw Score", "Adj Factor", "Adj Score",
                       "Buried", "Exposed", "Risk Level"]
            for i, h in enumerate(headers):
                cell = table.rows[0].cells[i]
                cell.text = h
                for p in cell.paragraphs:
                    for run in p.runs:
                        run.bold = True

            for ri, pep in enumerate(sorted(high_risk_peptides,
                                            key=lambda x: x.adjusted_score,
                                            reverse=True), 1):
                row = table.rows[ri]
                risk_level = "HIGH" if pep.adjusted_score > 0.5 else (
                    "MEDIUM" if pep.adjusted_score > 0.2 else "LOW")
                vals = [
                    pep.peptide[:25] + ("..." if len(pep.peptide) > 25 else ""),
                    f"{pep.raw_score:.3f}",
                    f"{pep.adjustment_factor:.3f}",
                    f"{pep.adjusted_score:.3f}",
                    str(pep.n_buried),
                    str(pep.n_exposed),
                    risk_level,
                ]
                for ci, val in enumerate(vals):
                    row.cells[ci].text = val

        # Detailed position analysis
        detailed_peptides = [p for p in assessment.peptides if p.adjusted_score > 0.3]
        if detailed_peptides:
            doc.add_heading("Detailed Position Analysis", level=2)
            for pep in detailed_peptides[:3]:
                doc.add_heading(f"Peptide: {pep.peptide}", level=3)
                doc.add_paragraph(f"Raw score: {pep.raw_score:.3f}")
                doc.add_paragraph(f"Adjusted score: {pep.adjusted_score:.3f}")

                backmut_positions = [p for p in pep.positions if p.is_backmutation]
                if backmut_positions:
                    table2 = doc.add_table(rows=1 + len(backmut_positions), cols=5)
                    table2.style = "Table Grid"
                    for i, h in enumerate(["Position", "relSASA", "Buried",
                                           "CDR Contact", "Tier"]):
                        cell = table2.rows[0].cells[i]
                        cell.text = h
                        for p in cell.paragraphs:
                            for run in p.runs:
                                run.bold = True

                    for ri, pos in enumerate(backmut_positions, 1):
                        row = table2.rows[ri]
                        buried_str = "Yes" if pos.buried else (
                            "No" if pos.buried is False else "Uncertain")
                        vals = [
                            pos.position,
                            f"{pos.rel_sasa:.3f}",
                            buried_str,
                            "Yes" if pos.cdr_contact else "No",
                            pos.tier or "",
                        ]
                        for ci, val in enumerate(vals):
                            row.cells[ci].text = val

                        # Color-code based on buriedness
                        if pos.buried is True:
                            row.cells[1].paragraphs[0].runs[0].font.color.rgb = (
                                RGBColor(0, 128, 0))  # green
                        elif pos.buried is False:
                            row.cells[1].paragraphs[0].runs[0].font.color.rgb = (
                                RGBColor(255, 0, 0))  # red

    # Cross-variant comparison
    if len(assessments) > 1:
        doc.add_heading("Cross-Variant Comparison", level=1)
        table3 = doc.add_table(rows=1 + len(assessments), cols=6)
        table3.style = "Table Grid"
        for i, h in enumerate(["Variant", "Structure", "Overall Adj",
                               "Adj Mean", "Final Score", "High-Risk Peptides"]):
            cell = table3.rows[0].cells[i]
            cell.text = h
            for p in cell.paragraphs:
                for run in p.runs:
                    run.bold = True

        for ri, a in enumerate(assessments, 1):
            row = table3.rows[ri]
            n_high = sum(1 for p in a.peptides if p.adjusted_score > 0.5)
            vals = [
                a.variant_name,
                a.structure_source,
                f"{a.overall_adjustment:.3f}",
                f"{a.adjusted_mean:.3f}" if a.adjusted_mean is not None else "N/A",
                (f"{a.adjusted_immunogenicity_score:.1f}%"
                 if a.adjusted_immunogenicity_score is not None else "N/A"),
                str(n_high),
            ]
            for ci, val in enumerate(vals):
                row.cells[ci].text = val

    doc.save(output_path)


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def _structure_names(vname: str) -> List[str]:
    variant_part = vname.split('_', 1)[1] if '_' in vname else vname
    bases = []
    for base in (vname, variant_part, vname.replace('_', ''),
                 variant_part.replace('_', '')):
        if base and base not in bases:
            bases.append(base)
    names: List[str] = []
    for base in bases:
        for prefix in ('', 'igfold_', 'af3_', 'boltz_', 'protenix_'):
            for suffix in ('.pdb', '.cif'):
                for value in (base, base.lower()):
                    names.append(f"{prefix}{value}{suffix}")
    seen = set()
    return [n for n in names if not (n in seen or seen.add(n))]


def _find_structure(structure_dir: str, vname: str) -> Optional[str]:
    names = _structure_names(vname)
    for name in names:
        candidate = os.path.join(structure_dir, name)
        if os.path.isfile(candidate):
            return candidate
    for root, _dirs, _files in os.walk(structure_dir):
        for name in names:
            candidate = os.path.join(root, name)
            if os.path.isfile(candidate):
                return candidate
    return None


def _load_chain_backmutations(csv_dir: Optional[str], chain_type: str) -> List[Dict]:
    if not csv_dir:
        return []
    import glob
    patterns = [
        f'backmutations_*.{chain_type}.csv',
        f'backmutations_{chain_type}.csv',
    ]
    for pattern in patterns:
        matches = sorted(glob.glob(os.path.join(csv_dir, pattern)))
        if matches:
            return load_backmutations_from_csv(matches[0])
    return []


def _recompute_counts(vdata: Dict, strong_threshold: float,
                      weak_threshold: float) -> None:
    scores = [float(p.get('adjusted_score', 0.0))
              for p in vdata.get('adjusted_peptides', [])]
    if not scores:
        return
    n_strong = sum(1 for s in scores if s >= strong_threshold)
    n_weak = sum(1 for s in scores if weak_threshold <= s < strong_threshold)
    vdata['n_strong'] = n_strong
    vdata['n_weak'] = n_weak
    vdata['n_total_windows'] = len(scores)
    vdata['epitope_density'] = round((n_strong + n_weak) / len(scores), 4)
    vdata['mean_epitope_load'] = round(sum(scores) / len(scores), 4)
    vdata['max_epitope'] = round(max(scores), 4)


def run_structural_risk_assessment(
    immunogenicity_results: Dict,
    structure_dir: str,
    backmutation_csv_dir: Optional[str] = None,
    output_dir: Optional[str] = None,
    method: str = "relSASA",
    verbose: bool = False,
    numbering_json: Optional[str] = None,
) -> Dict:
    """Run variant-specific structural adjustment and rewrite final scores."""
    if output_dir is None:
        output_dir = os.path.dirname(immunogenicity_results.get('input_fasta', '.'))
    os.makedirs(output_dir, exist_ok=True)

    if numbering_json is None and backmutation_csv_dir:
        cand = os.path.join(backmutation_csv_dir, 'variants_numbering.json')
        if os.path.isfile(cand):
            numbering_json = cand
    numbering = load_numbering_json(numbering_json) if numbering_json else {}
    if verbose and numbering:
        print(f"  Loaded position labels for {len(numbering)} variant(s) from "
              f"{numbering_json}")
    elif verbose:
        print("  WARNING: no variants_numbering.json found; backmutation/FR4 "
              "matching falls back to sequence index (approximate)")

    strong_threshold = float(immunogenicity_results.get('strong_threshold', 0.5))
    weak_threshold = float(immunogenicity_results.get('weak_threshold', 0.3))
    assessments: List[VariantStructuralRisk] = []

    for chain_type, cdata in immunogenicity_results.get('chains', {}).items():
        variants = cdata.get('variants', {})
        backmutations = _load_chain_backmutations(backmutation_csv_dir, chain_type)
        if verbose and backmutations:
            print(f"  Loaded {len(backmutations)} backmutation candidates for chain {chain_type}")

        assessed_names = []
        for vname, vdata in variants.items():
            vdata['structural_score_applied'] = False
            vdata['structural_adjustment'] = None
            structure_path = _find_structure(structure_dir, vname)
            if structure_path is None:
                vdata['structural_error'] = "no structure file found"
                if verbose:
                    print(f"  WARNING: No structure found for {vname}")
                continue

            labels = _labels_for_variant(numbering, vname, chain_type)
            variant_backmuts = filter_backmutations_for_variant(
                backmutations, vdata.get('sequence', ''), labels, chain_type)
            if verbose:
                print(f"  Assessing {vname} with {os.path.basename(structure_path)} "
                      f"({len(variant_backmuts)} reverted positions)")

            try:
                assessment = assess_variant_structural_risk(
                    variant_name=vname,
                    chain_type=chain_type,
                    sequence=vdata.get('sequence', ''),
                    per_peptide_hits=vdata.get('per_peptide_hits', []),
                    structure_path=structure_path,
                    backmutations=variant_backmuts,
                    method=method,
                    position_labels=labels,
                )
            except Exception as exc:
                vdata['structural_error'] = str(exc)
                vdata['structural_adjustment'] = None
                if verbose:
                    print(f"    structural assessment failed: {exc}")
                continue

            assessments.append(assessment)
            assessed_names.append(vname)
            vdata.pop('structural_error', None)
            vdata['structural_adjustment'] = round(assessment.overall_adjustment, 4)
            vdata['structural_source'] = assessment.structure_source
            vdata['structural_adjusted_mean'] = (
                round(assessment.adjusted_mean, 4)
                if assessment.adjusted_mean is not None else None)
            vdata['raw_mean_epitope_load'] = vdata.get('raw_mean_epitope_load',
                                                       vdata.get('mean_epitope_load'))
            vdata['adjusted_peptides'] = [
                {
                    'peptide': p.peptide,
                    'start_index': p.start_index,
                    'raw_score': p.raw_score,
                    'adjusted_score': p.adjusted_score,
                    'adjustment_factor': p.adjustment_factor,
                    'n_buried': p.n_buried,
                    'n_exposed': p.n_exposed,
                    'n_backmutation': p.n_backmutation,
                }
                for p in assessment.peptides
            ]
            _recompute_counts(vdata, strong_threshold, weak_threshold)

        ref_name = cdata.get('reference_name')
        ref_adj_mean = None
        if ref_name and ref_name in variants:
            ref_adj_mean = variants[ref_name].get('structural_adjusted_mean')
        donor_m = cdata.get('donor_metrics') or {}
        if ref_adj_mean is None:
            ref_adj_mean = donor_m.get('mean')
        cdata['structural_reference_name'] = ref_name
        cdata['structural_reference_mean'] = ref_adj_mean

        for vname in assessed_names:
            vdata = variants[vname]
            adj_mean = vdata.get('structural_adjusted_mean')
            if adj_mean is None:
                continue
            if ref_adj_mean and ref_adj_mean > 0:
                score = adj_mean / ref_adj_mean * 100.0
            else:
                score = adj_mean * 100.0
            vdata['immunogenicity_score'] = round(score, 1)
            vdata['structural_immunogenicity_score'] = round(score, 1)
            vdata['structural_score_applied'] = True
            assessment = next(a for a in assessments
                              if a.variant_name == vname and a.chain_type == chain_type)
            assessment.adjusted_immunogenicity_score = round(score, 1)
            if verbose:
                print(f"    {vname}: adjusted mean={adj_mean:.3f}, "
                      f"score={score:.1f}%")

    immunogenicity_results['structural_risk'] = {
        'method': method,
        'applied_to_final_score': True,
        'score_semantics': ('percent of structurally adjusted reference mean; '
                            '>100 means higher than reference'),
    }
    immunogenicity_results['score_semantics'] = (
        'structural: percent of structurally adjusted reference mean; '
        '>100 means higher than reference')

    if assessments:
        md_path = os.path.join(output_dir, 'structural_risk_assessment.md')
        generate_structural_risk_report(assessments, md_path)
        if verbose:
            print(f"  Generated Markdown report: {md_path}")
        docx_path = os.path.join(output_dir, 'structural_risk_assessment.docx')
        generate_structural_risk_word(assessments, docx_path)
        if verbose:
            print(f"  Generated Word report: {docx_path}")

    return immunogenicity_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Structural risk adjustment for immunogenicity prediction"
    )
    parser.add_argument("--immunogenicity-json", required=True,
                        help="Path to immunogenicity results JSON")
    parser.add_argument("--structure-dir", required=True,
                        help="Directory containing variant PDB files")
    parser.add_argument("--backmutation-dir",
                        help="Directory containing backmutation CSVs")
    parser.add_argument("--numbering-json",
                        help="variants_numbering.json (per-variant position labels); "
                             "auto-discovered in --backmutation-dir when omitted")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory for reports")
    parser.add_argument("--method", default="relSASA",
                        choices=["relSASA", "binary"],
                        help="Adjustment method")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    # Load immunogenicity results
    with open(args.immunogenicity_json) as f:
        results = json.load(f)

    # Run assessment
    updated = run_structural_risk_assessment(
        results,
        args.structure_dir,
        args.backmutation_dir,
        args.output_dir,
        method=args.method,
        verbose=args.verbose,
        numbering_json=args.numbering_json,
    )

    # Save updated results
    output_json = os.path.join(args.output_dir, 'immunogenicity_with_structural_risk.json')
    with open(output_json, 'w') as f:
        json.dump(updated, f, indent=2)
    print(f"Saved updated results to {output_json}")
