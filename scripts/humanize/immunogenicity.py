"""Optional per-position immunogenicity scoring (T-cell epitope propensity).

The back-mutation `benefit` term is by default a surface-exposure x germline
rarity proxy. This module lets a server replace/augment that proxy with a real
MHC-II binding prediction, so positions whose donor residue sits in a strong
predicted T-cell epitope get a higher humanization benefit (reverting them to
the human germline is more likely to reduce ADA risk).

Two input paths (both optional, both degrade gracefully):

  * ``load_scores_json``      - a precomputed {kabat_position: score} JSON map;
  * ``netmhciipan_scores``    - run NetMHCIIpan locally and parse its ``-xls``
                                output into per-position scores.

Scores are in [0, 1] (higher = stronger predicted epitope). No third-party
Python dependency; NetMHCIIpan is an external binary.
"""

from __future__ import annotations

import json
import os
import subprocess
from typing import Dict, Optional

from .numbering import NumberedChain


def load_scores_json(path: str) -> Dict[str, float]:
    """Load ``{kabat_position: score}`` (e.g. ``{"H71": 0.8, "L45": 0.3}``)."""
    with open(path) as fh:
        data = json.load(fh)
    if isinstance(data, dict) and "scores" in data:
        data = data["scores"]
    return {str(k): float(v) for k, v in data.items()}


def _norm_rank(pct_rank: float) -> float:
    """Map NetMHCIIpan %Rank to an epitope score.

    Strong binders (%Rank < 2) -> ~0.8-1.0; weak binders (%Rank < 10) -> >0;
    %Rank >= 10 -> 0.
    """
    return max(0.0, min(1.0, 1.0 - pct_rank / 10.0))


def parse_netmhciipan_xls(path: str, numbered: NumberedChain) -> Dict[str, float]:
    """Parse a NetMHCIIpan ``-xls`` file into {kabat_position: max score}.

    The parser is column-name driven (``Pos`` / ``Peptide`` / ``%Rank``), so it
    tolerates version-specific extra columns. A peptide at 1-based input
    position ``p`` of length ``L`` covers input residues ``p..p+L-1``.
    """
    idx_to_pos = {r.index: r.pos for r in numbered.residues}
    scores: Dict[str, float] = {}
    col: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.rstrip("\n").split()
            if not col:
                low = [p.lower() for p in parts]
                if "pos" in low and any("%rank" in p for p in low):
                    col = {name: i for i, name in enumerate(low)}
                continue
            pos_i = col.get("pos")
            rank_i = next((i for name, i in col.items() if "%rank" in name), None)
            pep_i = col.get("peptide")
            if pos_i is None or rank_i is None or pos_i >= len(parts):
                continue
            try:
                start = int(parts[pos_i])
                rank = float(parts[rank_i])
            except (ValueError, IndexError):
                continue
            length = len(parts[pep_i]) if (pep_i is not None and pep_i < len(parts)) else 15
            score = _norm_rank(rank)
            if score <= 0:
                continue
            for offset in range(max(1, length)):
                kabat = idx_to_pos.get(start - 1 + offset)
                if kabat and score > scores.get(kabat, 0.0):
                    scores[kabat] = score
    return scores


def _write_peptide_fasta(sequence: str, path: str, length: int = 15) -> int:
    n = 0
    with open(path, "w") as fh:
        for i in range(max(0, len(sequence) - length + 1)):
            fh.write(f">p{i + 1}\n{sequence[i:i + length]}\n")
            n += 1
    return n


def netmhciipan_scores(
    numbered: NumberedChain,
    binary: str,
    alleles: str = "HLA-DRB1*01:01,HLA-DRB1*04:01,HLA-DRB1*07:01,"
                    "HLA-DRB1*15:01,HLA-DRB3*01:01,HLA-DRB3*02:02,"
                    "HLA-DRB4*01:01,HLA-DRB5*01:01",
    env: Optional[str] = None,
    workdir: Optional[str] = None,
) -> Optional[Dict[str, float]]:
    """Run NetMHCIIpan and return {kabat_position: epitope score}.

    Returns None (with no exception) when the tool is unavailable or fails, so
    the pipeline silently falls back to the default immunogenicity proxy.
    """
    if not binary or not os.path.exists(binary):
        import warnings
        warnings.warn(f"[immunogenicity] NetMHCIIpan not found: {binary!r}; "
                      "using the exposure x rarity proxy")
        return None
    import tempfile
    workdir = workdir or tempfile.mkdtemp(prefix="netmhciipan_")
    os.makedirs(workdir, exist_ok=True)
    pep_fa = os.path.join(workdir, "peptides.fasta")
    out_xls = os.path.join(workdir, "netmhciipan.xls")
    _write_peptide_fasta(numbered.sequence, pep_fa, length=15)

    cmd = [binary, "-f", pep_fa, "-inptype", "1", "-a", alleles,
           "-xls", "-xlsfile", out_xls]
    if env:
        cmd = ["conda", "run", "-n", env] + cmd
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except (OSError, ValueError) as e:
        import warnings
        warnings.warn(f"[immunogenicity] NetMHCIIpan failed to start: {e}")
        return None
    if proc.returncode != 0 or not os.path.exists(out_xls):
        import warnings
        warnings.warn(
            f"[immunogenicity] NetMHCIIpan failed (rc={proc.returncode}): "
            f"{proc.stderr[-300:]}")
        return None
    return parse_netmhciipan_xls(out_xls, numbered)
