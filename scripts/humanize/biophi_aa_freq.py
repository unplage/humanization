"""BioPhi/OAS per-position amino-acid frequency lookup (frozen table).

Reads the frozen ``data/germline/biophi_germline_frequency.json`` produced by
``scripts/build_biophi_frequency.py``. Never requires BioPhi to be installed;
missing file or unknown position/gene simply yields ``None`` (no exception),
so report rendering never fails when BioPhi is absent.
"""
from __future__ import annotations

import json
import os
import re
from typing import Optional

_DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                          "..", "data", "germline",
                          "biophi_germline_frequency.json")

_DATA: Optional[dict] = None


def _load() -> Optional[dict]:
    global _DATA
    if _DATA is None:
        try:
            with open(_DATA_PATH) as fh:
                _DATA = json.load(fh)
        except (OSError, ValueError):
            _DATA = {}
    return _DATA


def _chain_key(chain_type: str) -> str:
    if chain_type == "H":
        return "H"
    if chain_type in ("K", "Kappa", "KAPPA"):
        return "K"
    return "L"


def _family_key(gene_id: str) -> str:
    """BioPhi germline-family table key, e.g. IGHV3-23*01 -> IGHV3."""
    base = gene_id.split("*")[0]
    match = re.match(r"^(IG[HKL]V\d+D?)", base)
    return match.group(1) if match else base


def get_aa_freq(chain_type: str, gene_id: Optional[str],
                imgt_pos: Optional[str], aa: Optional[str]) -> Optional[float]:
    """Return OAS frequency of residue ``aa`` at IMGT position ``imgt_pos``.

    Looks up the germline-family table first (best match for the chosen V
    gene), falling back to the chain-wide table. Returns ``None`` when the
    table is missing, the position is absent, or the residue is unobserved.
    """
    if aa is None or aa == "-" or not imgt_pos:
        return None
    data = _load()
    if not data:
        return None
    chain_table = data.get("chain", {}).get(_chain_key(chain_type or "H"), {})
    chain_pos = chain_table.get(imgt_pos)
    if gene_id:
        family = _family_key(gene_id)
        fam_pos = data.get("family", {}).get(family, {}).get(imgt_pos)
        if fam_pos and aa.upper() in fam_pos:
            return float(fam_pos[aa.upper()])
    if chain_pos and aa.upper() in chain_pos:
        return float(chain_pos[aa.upper()])
    return None


def get_top_aa_freqs(chain_type: str, gene_id: Optional[str],
                     imgt_pos: Optional[str], top_n: int = 5) -> list:
    """Return top-N (aa, freq) pairs at IMGT position ``imgt_pos``.

    Family table first, chain-wide table as fallback. Returns an empty list
    when the table/position is missing — never raises, so reports render fine
    even without BioPhi installed.
    """
    if not imgt_pos:
        return []
    data = _load()
    if not data:
        return []
    pos_table = None
    if gene_id:
        family = _family_key(gene_id)
        fam_pos = data.get("family", {}).get(family, {}).get(imgt_pos)
        if fam_pos:
            pos_table = fam_pos
    if pos_table is None:
        chain_table = data.get("chain", {}).get(_chain_key(chain_type or "H"), {})
        pos_table = chain_table.get(imgt_pos)
    if not pos_table:
        return []
    ranked = sorted(pos_table.items(), key=lambda kv: -kv[1])
    out = []
    for aa, freq in ranked[:top_n]:
        try:
            out.append((str(aa).upper(), float(freq)))
        except (TypeError, ValueError):
            continue
    return out
