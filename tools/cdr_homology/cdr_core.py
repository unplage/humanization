#!/usr/bin/env python3
"""Shared core for the CDR-homology tools.

Self-contained (no pipeline imports); stdlib + Biopython only.

Provides:
  * robust IMGT numbering via ANARCI (version-agnostic return-shape probing)
  * IMGT CDR extraction (27-38 / 56-65 / 105-117) in ANARCI sequence order
  * alignment with global *and* local modes plus coverage-aware metrics
  * VHH hallmark detection (IMGT 42/49/50/52)
  * filesystem-safe, collision-free names
  * provenance helpers (timestamps, SHA-256, tool/library versions)
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import sys
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

TOOL_NAME = "cdr_homology"
TOOL_VERSION = "2.0.0"

# IMGT CDR boundaries (Lefranc 2003).  Kept in sync with
# ``scripts/humanize/numbering.py::NATIVE_IMGT_CDRS`` and asserted by the tests.
IMGT_CDR_RANGES: Dict[str, Tuple[int, int]] = {
    "CDR1": (27, 38),
    "CDR2": (56, 65),
    "CDR3": (105, 117),
}

# Alignment mode per CDR.  CDR1/CDR2 are short and near-invariant, so a full
# length comparison is meaningful.  CDR3 varies from 3 to >30 residues, so a
# full length (global) alignment divides the identity by the length mismatch
# and systematically under-reports homology; the homologous *core* is scored
# with a local (Smith-Waterman) alignment instead.
DEFAULT_CDR_MODE: Dict[str, str] = {
    "CDR1": "global",
    "CDR2": "global",
    "CDR3": "local",
}

# Camelid VHH FR2 hallmark, IMGT numbering (Kabat 37/44/45/47).
VHH_HALLMARK_IMGT: Dict[int, Tuple[str, ...]] = {
    42: ("F", "Y"),
    49: ("E", "Q"),
    50: ("R",),
    52: ("G", "F", "S"),
}
VHH_MIN_MATCHES = 3


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

def parse_fasta(path: str) -> List[Tuple[str, str]]:
    """Return ``[(header, sequence), ...]``; the header keeps its description."""
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    seq: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line)
    if name is not None:
        records.append((name, "".join(seq)))
    return records


# ---------------------------------------------------------------------------
# ANARCI IMGT numbering
# ---------------------------------------------------------------------------

def _is_numbering_entry(obj) -> bool:
    """True for a single ANARCI numbering entry ``((pos, ins), aa)``."""
    return (
        isinstance(obj, tuple)
        and len(obj) == 2
        and isinstance(obj[0], tuple)
        and len(obj[0]) == 2
        and isinstance(obj[0][0], int)
        and isinstance(obj[1], str)
        and len(obj[1]) == 1
    )


def _find_numbering(obj, depth: int = 0):
    """Depth-limited search for the numbering list in an arbitrary return value."""
    if depth > 8:
        return None
    if isinstance(obj, list) and obj and _is_numbering_entry(obj[0]):
        return obj
    if isinstance(obj, (list, tuple)):
        for item in obj:
            if isinstance(item, (list, tuple)):
                found = _find_numbering(item, depth + 1)
                if found is not None:
                    return found
    return None


def extract_numbering(result):
    """Pull the IMGT numbering list out of an ANARCI ``output=False`` result.

    ANARCI has used two incompatible return shapes:

      * ``anarci >= 2024``  -> ``(results, hits, species)`` where
        ``results[0][0]`` is ``[(numbering, start, end)]``
      * ``anarci 0.1.x``    -> ``[(numbering, start, end), ...]``

    The old code indexed ``result[0][0][0][0]`` unconditionally, which yields
    an ``int``/``tuple`` on the legacy shape and then raises an *uncaught*
    ``TypeError`` while iterating.  The shape is therefore probed explicitly
    first and only then via a generic search.
    """
    candidates = []
    try:
        candidates.append(result[0][0][0][0])   # anarci >= 2024
    except (TypeError, IndexError, KeyError):
        pass
    try:
        candidates.append(result[0][0])          # anarci 0.1.x
    except (TypeError, IndexError, KeyError):
        pass
    for cand in candidates:
        if isinstance(cand, list) and cand and _is_numbering_entry(cand[0]):
            return cand
    return _find_numbering(result)


def anarci_available() -> bool:
    """True when the ANARCI package can be imported."""
    import importlib.util

    return importlib.util.find_spec("anarci") is not None


# One ANARCI run per unique sequence; the result is reused for numbering,
# chain-type detection and VHH detection.  Bounded to keep memory in check.
_NUMBER_CACHE: Dict[str, Tuple[Optional[list], Optional[str]]] = {}
_NUMBER_CACHE_MAX = 4096


def _anarci_raw(sequence: str) -> Tuple[Optional[list], Optional[str]]:
    """Run ANARCI once; return ``(numbering_list, chain_type)``.

    ``chain_type`` is ``'H'`` for heavy and ``'L'`` for kappa/lambda, taken
    from ANARCI's own HMM hit table -- this is what lets a headerless FASTA
    record be classified correctly instead of guessing from sequence length.
    """
    seq = sequence.upper().strip()
    if not seq:
        return None, None
    if seq in _NUMBER_CACHE:
        return _NUMBER_CACHE[seq]
    try:
        from anarci import anarci
    except ImportError:
        return None, None

    try:
        result = anarci([("query", seq)], scheme="imgt", output=False)
    except Exception:
        return None, None

    numbering = extract_numbering(result)
    chain_type: Optional[str] = None
    try:
        hits = result[1][0]
        if hits:
            called = str(hits[0].get("chain_type", "")).upper()
            chain_type = {"H": "H", "K": "L", "L": "L"}.get(called)
    except (IndexError, KeyError, TypeError, AttributeError):
        pass

    if len(_NUMBER_CACHE) >= _NUMBER_CACHE_MAX:
        _NUMBER_CACHE.clear()
    _NUMBER_CACHE[seq] = (numbering, chain_type)
    return numbering, chain_type


def number_sequence_imgt(sequence: str, chain_type: str) -> Optional[Dict[str, str]]:
    """Number ``sequence`` with ANARCI (IMGT) and return ``{label: aa}``.

    ``chain_type`` only selects the label prefix (``H``/``L``); the actual
    chain call comes from ANARCI.  Returns ``None`` when ANARCI is
    unavailable, raises, or cannot number the sequence.
    """
    numbering, _ = _anarci_raw(sequence)
    if not numbering:
        return None

    posmap: Dict[str, str] = {}
    for entry in numbering:
        if not _is_numbering_entry(entry):
            continue
        (pos_num, ins_flag), aa = entry
        if not aa or aa == "-":
            continue
        label = f"{chain_type}{pos_num}"
        if ins_flag and ins_flag != " ":
            label = f"{label}{ins_flag}"
        posmap[label] = aa.upper()
    return posmap or None


def detect_chain_type(sequence: str) -> Optional[str]:
    """Return ``'H'`` or ``'L'`` from ANARCI's chain call, or ``None``."""
    _, chain_type = _anarci_raw(sequence)
    return chain_type


def _position_number(label: str) -> Optional[int]:
    """Numeric IMGT position of a label like ``H105`` or ``H111A``."""
    digits = "".join(c for c in label[1:] if c.isdigit())
    return int(digits) if digits else None


def extract_cdrs_from_posmap(posmap: Dict[str, str], chain_type: str) -> Dict[str, str]:
    """Extract CDR1/CDR2/CDR3 from an IMGT position map.

    The posmap is walked in ANARCI sequence order.  IMGT places CDR3
    insertions as 111, 111A, 111B, ... then 112C, 112B, 112A and finally the
    112 base, so a numeric walk would transpose the loop.
    """
    cdrs = {name: "" for name in IMGT_CDR_RANGES}
    for label, aa in posmap.items():
        if not aa or aa == "-" or not label.startswith(chain_type):
            continue
        pos_num = _position_number(label)
        if pos_num is None:
            continue
        for cdr_name, (start, end) in IMGT_CDR_RANGES.items():
            if start <= pos_num <= end:
                cdrs[cdr_name] += aa
                break
    return cdrs


def extract_cdrs(sequence: str, chain_type: str) -> Optional[Dict[str, str]]:
    """Full CDR extraction.  Returns ``None`` when the chain cannot be numbered."""
    posmap = number_sequence_imgt(sequence, chain_type)
    if not posmap:
        return None
    return extract_cdrs_from_posmap(posmap, chain_type)


# ---------------------------------------------------------------------------
# VHH detection
# ---------------------------------------------------------------------------

def detect_vhh(posmap: Dict[str, str], chain_type: str = "H"):
    """Detect the camelid FR2 hallmark in an IMGT position map.

    Returns ``(is_vhh, n_matched, matched_labels)``.
    """
    if chain_type != "H" or not posmap:
        return False, 0, []
    matched = []
    for pos, allowed in VHH_HALLMARK_IMGT.items():
        aa = posmap.get(f"H{pos}")
        if aa and aa.upper() in allowed:
            matched.append(f"H{pos}{aa}")
    return len(matched) >= VHH_MIN_MATCHES, len(matched), matched


def chain_kind(vh_seq: str, vl_seq: str) -> Tuple[str, Dict[str, object]]:
    """Classify a query as ``fab`` / ``vhh`` / ``vhh_suspect`` / ``heavy_only``.

    A single heavy chain carrying the camelid hallmark is a VHH; a single
    heavy chain without it is flagged as a suspect rather than silently
    reported as an ordinary VH.
    """
    info: Dict[str, object] = {}
    if vl_seq:
        info["kind"] = "fab"
        return "fab", info
    if not vh_seq:
        info["kind"] = "unknown"
        return "unknown", info
    posmap = number_sequence_imgt(vh_seq, "H")
    if posmap:
        is_vhh, n_matched, matched = detect_vhh(posmap, "H")
        info["vhh_hallmark"] = matched
        info["vhh_hallmark_count"] = n_matched
        if is_vhh:
            info["kind"] = "vhh"
            return "vhh", info
    info["kind"] = "vhh_suspect"
    return "vhh_suspect", info


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

_ALIGNERS: Dict[str, Any] = {}
_SELF_SCORES: Dict[Tuple[str, str], float] = {}

BLOSUM62 = "BLOSUM62"
OPEN_GAP_SCORE = -10.0
EXTEND_GAP_SCORE = -0.5


def _get_aligner(mode: str = "global"):
    """Return a cached Biopython aligner for ``mode`` (``global``/``local``)."""
    if mode in _ALIGNERS:
        return _ALIGNERS[mode]
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load(BLOSUM62)
    aligner.open_gap_score = OPEN_GAP_SCORE
    aligner.extend_gap_score = EXTEND_GAP_SCORE
    _ALIGNERS[mode] = aligner
    return aligner


def self_score(sequence: str, mode: str = "global") -> float:
    """Score of ``sequence`` against itself -- the normaliser for raw scores."""
    key = (mode, sequence)
    if key not in _SELF_SCORES:
        aligner = _get_aligner(mode)
        try:
            _SELF_SCORES[key] = float(aligner.score(sequence, sequence))
        except Exception:
            _SELF_SCORES[key] = 0.0
    return _SELF_SCORES[key]


def _aligned_strings(best) -> Tuple[str, str]:
    """Render a pairwise alignment as ``(query_aligned, subject_aligned)``.

    Biopython calls the first argument the *target*; ``query.py`` passes the
    database CDR as target and the query CDR as query.
    """
    coords = best.coordinates
    seqs = best.sequences
    q_parts: List[str] = []
    s_parts: List[str] = []
    for i in range(coords.shape[1] - 1):
        t_start, t_end = int(coords[0, i]), int(coords[0, i + 1])
        q_start, q_end = int(coords[1, i]), int(coords[1, i + 1])
        if t_start == t_end:
            q_parts.append(seqs[1][q_start:q_end])
            s_parts.append("-" * (q_end - q_start))
        elif q_start == q_end:
            q_parts.append("-" * (t_end - t_start))
            s_parts.append(seqs[0][t_start:t_end])
        else:
            q_parts.append(seqs[1][q_start:q_end])
            s_parts.append(seqs[0][t_start:t_end])
    return "".join(q_parts), "".join(s_parts)


def align_sequences(query: str, subject: str, mode: str = "global") -> Dict[str, Any]:
    """Align two CDR sequences and return identity / coverage / score metrics.

    Metrics
    -------
    ``identity``
        Headline identity.  ``global`` mode: identical residues divided by the
        full alignment length (gaps included) -- the classic Needleman-Wunsch
        value.  ``local`` mode: identical residues divided by the length of
        the homologous core, i.e. the local alignment span.
    ``coverage``
        Fraction of the *query* CDR that participates in the alignment.  For a
        global alignment this is 1.0; for a local alignment it says how much
        of the query CDR3 is actually matched.
    ``identity_over_query``
        Identical residues divided by the query length.  For a short query
        CDR3 matched inside a long subject CDR3 this is the metric that
        answers "how much of my loop is shared".
    ``score_norm``
        Raw BLOSUM62 score divided by the self-alignment score of the query,
        so scores are comparable across queries of different lengths.
    """
    empty = {
        "identity": 0.0,
        "identity_global": 0.0,
        "identity_core": 0.0,
        "identity_over_query": 0.0,
        "coverage": 0.0,
        "score": 0.0,
        "score_norm": 0.0,
        "matches": 0,
        "aligned_len": 0,
        "query_aligned": "",
        "subject_aligned": "",
        "mode": mode,
    }
    if not query or not subject:
        return empty

    aligner = _get_aligner(mode)
    try:
        alignments = aligner.align(subject, query)
        best = alignments[0]
    except Exception:
        return empty

    aln_q, aln_s = _aligned_strings(best)
    if not aln_q and not aln_s:
        return empty

    matches = sum(1 for a, b in zip(aln_q, aln_s) if a == b and a not in ("-", "."))
    aligned_len = max(len(aln_q), 1)
    q_residues = sum(1 for c in aln_q if c not in ("-", "."))

    identity_core = matches / aligned_len
    identity_global = identity_core
    if mode == "local":
        # A local alignment may still place gaps inside the core; the global
        # equivalent is obtained by aligning the two full sequences.
        glob = align_sequences(query, subject, mode="global")
        identity_global = float(glob["identity_global"])

    coverage = q_residues / max(len(query), 1)
    identity_over_query = matches / max(len(query), 1)

    raw_score = float(getattr(best, "score", 0.0))
    norm = self_score(query, mode)
    score_norm = (raw_score / norm) if norm > 0 else 0.0

    return {
        "identity": round(identity_core, 4),
        "identity_global": round(identity_global, 4),
        "identity_core": round(identity_core, 4),
        "identity_over_query": round(identity_over_query, 4),
        "coverage": round(coverage, 4),
        "score": round(raw_score, 1),
        "score_norm": round(score_norm, 4),
        "matches": matches,
        "aligned_len": aligned_len,
        "query_aligned": aln_q,
        "subject_aligned": aln_s,
        "mode": mode,
    }


def cdr_mode(cdr_name: str, override: Optional[str] = None) -> str:
    """Alignment mode for a CDR (``auto`` -> per-CDR default)."""
    if override and override != "auto":
        return override
    return DEFAULT_CDR_MODE.get(cdr_name, "global")


# ---------------------------------------------------------------------------
# Names and paths
# ---------------------------------------------------------------------------

_UNSAFE = re.compile(r"[^A-Za-z0-9_.-]+")
_CHAIN_SUFFIXES = ("_VHH", "_VH", "_VL", "_HEAVY", "_LIGHT", "_H", "_L")


def sanitize(text: str) -> str:
    """Make a text token safe to use inside a FASTA record name."""
    return _UNSAFE.sub("_", str(text)).strip("_") or "query"


def base_name(header: str) -> str:
    """Strip the description and a chain suffix from a FASTA header.

    ``'4D5_VH mouse anti-HER2 heavy chain'`` -> ``'4D5'``;
    ``'AMG110_V0_VH'`` -> ``'AMG110_V0'``.
    """
    token = str(header).split()[0] if str(header).split() else str(header)
    token = token.strip()
    for suffix in _CHAIN_SUFFIXES:
        if token.upper().endswith(suffix):
            token = token[: -len(suffix)]
            break
    return sanitize(token)


def safe_filename(name: str, used: Optional[set] = None, max_len: int = 60) -> str:
    """Filesystem-safe file stem, collision-free against ``used``.

    Collisions are resolved by appending ``_2``, ``_3``, ... and by mixing a
    short digest into names that had to be truncated.
    """
    base = _UNSAFE.sub("_", str(name)).strip("._-") or "query"
    if len(base) > max_len:
        digest = hashlib.sha1(base.encode("utf-8")).hexdigest()[:8]
        base = f"{base[:max_len - 9]}_{digest}"
    if used is None:
        return base
    candidate = base
    counter = 2
    while candidate in used:
        candidate = f"{base}_{counter}"
        counter += 1
    used.add(candidate)
    return candidate


def markdown_cell(value) -> str:
    """Escape a value so it cannot break a Markdown table cell."""
    text = "" if value is None else str(value)
    text = text.replace("\\", "\\\\").replace("|", "\\|")
    text = text.replace("\r\n", " ").replace("\n", " ").replace("\r", " ")
    return text


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def file_sha256(path: str, chunk: int = 1 << 20) -> Optional[str]:
    """SHA-256 of a file, or ``None`` when it cannot be read."""
    try:
        digest = hashlib.sha256()
        with open(path, "rb") as fh:
            for block in iter(lambda: fh.read(chunk), b""):
                digest.update(block)
        return digest.hexdigest()
    except OSError:
        return None


def library_versions() -> Dict[str, str]:
    """Versions of the optional libraries the analysis depends on."""
    versions: Dict[str, str] = {}
    for module in ("anarci", "Bio", "pandas"):
        try:
            mod = __import__(module)
            versions[module] = str(getattr(mod, "__version__", "unknown"))
        except ImportError:
            versions[module] = "not installed"
    return versions


def timestamp() -> str:
    """UTC timestamp in ISO-8601 form."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_json(path: str, payload) -> str:
    """Write JSON with a stable encoding (UTF-8) and create parent dirs."""
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    return path


def read_json(path: str):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def python_version() -> str:
    return sys.version.split()[0]
