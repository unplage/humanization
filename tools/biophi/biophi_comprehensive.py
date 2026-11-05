#!/usr/bin/env python3
"""Full BioPhi capability reproduction via the Python API.

Reproduces **all** BioPhi analysis capabilities with maximal defaults:

- **Sapiens**: mean score, per-position scores, full 20-aa probability
  matrix, best-score sequence, masked sequence, residue + sequence
  embeddings, Sapiens-humanized FASTA, mutation proposals (Designer scores).
- **OASis humanness**: identity, percentile, germline content, V/J germline
  matches, peptide-level table, positional humanness, OASis curve, top
  germline frequencies, chain-type / germline-family residue frequencies,
  paired-antibody (VH+VL) humanness, full XLSX workbook.
- **Humanization**: Sapiens, CDR grafting, Manual/Designer, and coordinated
  full-antibody (VH+VL) modes, each with alignments, mutation lists and
  top mutation proposals; optional OASis before/after comparison.
- **Reports**: JSON / Markdown / CSV / heatmap / embeddings / humanized
  FASTA, plus cross-validation against pipeline results.

Runs inside the dedicated ``biophi`` conda environment: if BioPhi is not
importable the script re-execs itself with the environment interpreter.

Examples
--------
python3 tools/biophi/biophi_comprehensive.py check

python3 tools/biophi/biophi_comprehensive.py evaluate-all \\
    --input outputs/step2/variants.fasta \\
    --metadata outputs/step2/variants_numbering.json \\
    --pipeline-result outputs/step2/humanization_result.json \\
    --outdir outputs/biophi_full

python3 tools/biophi/biophi_comprehensive.py sapiens \\
    --input data/examples/mouse_4d5_fab.fasta --outdir outputs/biophi_full
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

TOOL_VERSION = "1.0.0"
DEFAULT_ENV = "biophi"
DEFAULT_DB = "data/biophi/OASis_9mers_v1.db"
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

# Pairing suffixes, longest first so "_VH" wins over "_H".
_PAIR_SUFFIXES: Tuple[Tuple[str, str], ...] = (
    ("_VH", "H"), ("_VL", "L"), ("_HC", "H"), ("_LC", "L"),
    (".VH", "H"), (".VL", "L"), (".HC", "H"), (".LC", "L"),
    ("_H", "H"), ("_L", "L"),
    (".H", "H"), (".L", "L"),
)

# ---------------------------------------------------------------------------
# Bootstrap: offline HF defaults + auto re-exec into the biophi conda env
# ---------------------------------------------------------------------------


def _argv_option(name: str, default: Optional[str] = None) -> Optional[str]:
    argv = sys.argv[1:]
    prefix = name + "="
    for i, arg in enumerate(argv):
        if arg == name and i + 1 < len(argv):
            return argv[i + 1]
        if arg.startswith(prefix):
            return arg[len(prefix):]
    return default


def _apply_hf_env() -> None:
    endpoint = _argv_option("--hf-endpoint")
    if endpoint:
        os.environ["HF_ENDPOINT"] = endpoint
        os.environ.pop("HF_HUB_OFFLINE", None)
        os.environ.pop("TRANSFORMERS_OFFLINE", None)
    else:
        # Server has no direct huggingface.co access; local cache is populated.
        os.environ.setdefault("HF_HUB_OFFLINE", "1")
        os.environ.setdefault("TRANSFORMERS_OFFLINE", "1")
        os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")


def _find_env_python(env: str) -> Optional[str]:
    candidates: List[Path] = []
    conda = shutil.which("conda")
    if conda:
        conda_root = Path(conda).resolve().parent.parent
        candidates.append(conda_root / "envs" / env / "bin" / "python")
        try:
            out = subprocess.run(
                [conda, "info", "--envs"], capture_output=True, text=True, timeout=60
            )
            for line in out.stdout.splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[0] == env:
                    candidates.append(Path(parts[1]) / "bin" / "python")
        except Exception:
            pass
    prefix = os.environ.get("CONDA_PREFIX", "")
    if prefix and Path(prefix).name == env:
        candidates.append(Path(prefix) / "bin" / "python")
    for cand in candidates:
        if cand.exists() and os.access(str(cand), os.X_OK):
            return str(cand)
    return None


def _bootstrap() -> None:
    """Re-exec into the biophi conda env when BioPhi is not importable."""
    _apply_hf_env()
    try:
        import sapiens  # noqa: F401
        import biophi  # noqa: F401
        import abnumber  # noqa: F401
        return
    except ImportError:
        pass

    env = _argv_option("--env", DEFAULT_ENV) or DEFAULT_ENV
    script = str(Path(__file__).resolve())
    env_python = _find_env_python(env)
    if env_python:
        os.execv(env_python, [env_python, script, *sys.argv[1:]])
    conda = shutil.which("conda")
    if conda:
        os.execv(conda, [conda, "run", "-n", env, "python", script, *sys.argv[1:]])
    sys.stderr.write(
        f"ERROR: cannot import BioPhi and conda env {env!r} was not found.\n"
        "Install the env as described in tools/biophi/README.md.\n"
    )
    raise SystemExit(2)


_bootstrap()

# Heavy imports only after bootstrap (guaranteed to run inside the env).
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from abnumber import Chain  # noqa: E402
from sapiens.predict import predict_scores  # noqa: E402

from biophi.humanization.methods.humanness import (  # noqa: E402
    ChainHumanness,
    OASisParams,
    get_antibody_humanness,
    get_chain_humanness,
)
from biophi.humanization.methods.humanization import (  # noqa: E402
    CDRGraftingHumanizationParams,
    HumanizationParams,
    ManualHumanizationParams,
    SapiensHumanizationParams,
    humanize_antibody,
    humanize_chain,
    cdr_grafting_humanize_chain,
    sapiens_humanize_chain,
)

# Shared stdlib helpers from the sibling CLI wrapper.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from biophi_analyze import (  # noqa: E402
    ToolError,
    _db_info,
    _env_info,
    _hf_cache_info,
    _md5_file,
    _now,
    _sha256_file,
    _write_csv,
    _write_json,
    read_fasta,
)


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return [_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    return str(obj)


def _round(value: Any, digits: int = 6) -> Any:
    try:
        if value is None:
            return None
        return round(float(value), digits)
    except (TypeError, ValueError):
        return value


def _pos_str(pos: Any) -> str:
    try:
        return pos.format()
    except Exception:
        return str(pos)


def _chain_type_for_sapiens(chain: Optional[Chain], record_id: str) -> str:
    """Sapiens accepts H / K / L; K and L share the light-chain model."""
    if chain is not None:
        return chain.chain_type if chain.chain_type in ("H", "K", "L") else "L"
    rid_upper = record_id.upper()
    for suffix, token in _PAIR_SUFFIXES:
        if rid_upper.endswith(suffix):
            return token
    if "|H_" in record_id or rid_upper.startswith("H"):
        return "H"
    return "L"


def _split_record_id(record_id: str) -> Dict[str, str]:
    """Parse pipeline FASTA ids: ``base|H_V2|description``."""
    parts = record_id.split("|")
    base = parts[0]
    variant = parts[1] if len(parts) > 1 else ""
    description = "|".join(parts[2:]) if len(parts) > 2 else ""
    chain_token = ""
    if variant.startswith(("H_", "L_")):
        chain_token = variant[0]
    return {
        "base": base,
        "variant": variant,
        "description": description,
        "chain_token": chain_token,
        "pair_base": base[:-3] if base.endswith(("_VH", "_VL", "_HC", "_LC")) else base,
        "pair_variant": variant[2:] if variant.startswith(("H_", "L_")) else variant,
    }


# Pairing suffixes live near the module constants (_PAIR_SUFFIXES).


def _pair_key(record_id: str, chain: Optional[Chain]) -> Optional[str]:
    """Key used to pair VH/VL records of the same antibody/variant.

    Recognises ``base|H_V2`` pipeline ids as well as the common FASTA
    suffix styles: ``name_VH/_VL``, ``name.VH/.VL``, ``name.H/.L``,
    ``name_H/_L``, ``name_HC/_LC``.
    """
    parts = _split_record_id(record_id)
    if "|" in record_id and parts["variant"]:
        return f"{parts['pair_base']}|{parts['pair_variant']}"
    for suffix, _token in _PAIR_SUFFIXES:
        if record_id.endswith(suffix) and len(record_id) > len(suffix):
            return record_id[: -len(suffix)]
    if chain is not None and chain.name:
        name = chain.name
        for suffix, _token in _PAIR_SUFFIXES:
            if name.endswith(suffix) and len(name) > len(suffix):
                return name[: -len(suffix)]
    # Single unpaired chain: analyse on its own.
    return None


def _number_records(
    records: Sequence[Tuple[str, str]], scheme: str, cdr_definition: str
) -> List[Dict[str, Any]]:
    numbered: List[Dict[str, Any]] = []
    for record_id, seq in records:
        entry: Dict[str, Any] = {
            "id": record_id,
            "sequence": seq,
            "length": len(seq),
            "chain": None,
            "numbering_error": None,
        }
        entry.update(_split_record_id(record_id))
        try:
            entry["chain"] = Chain(
                seq, scheme=scheme, cdr_definition=cdr_definition, name=record_id
            )
        except Exception as exc:  # numbering failure is non-fatal
            entry["numbering_error"] = str(exc)
        entry["chain_type"] = (
            entry["chain"].chain_type if entry["chain"] is not None else ""
        )
        entry["sapiens_chain_type"] = _chain_type_for_sapiens(
            entry["chain"], record_id
        )
        entry["pair_key"] = _pair_key(record_id, entry["chain"])
        numbered.append(entry)
    return numbered


def _group_pairs(
    numbered: Sequence[Dict[str, Any]]
) -> Tuple[List[Tuple[Optional[Dict], Optional[Dict]]], List[Dict[str, Any]]]:
    """Split numbered records into (VH, VL) pairs and leftover singles."""
    by_key: Dict[str, Dict[str, Dict[str, Any]]] = {}
    singles: List[Dict[str, Any]] = []
    for entry in numbered:
        key = entry.get("pair_key")
        if not key:
            singles.append(entry)
            continue
        slot = by_key.setdefault(key, {})
        ctype = entry.get("chain_type") or entry.get("chain_token") or ""
        side = "vh" if ctype in ("H", "") and entry.get("chain_token") != "L" else "vl"
        if side in slot:
            # Ambiguous: keep the first, treat the second as a single.
            singles.append(entry)
        else:
            slot[side] = entry
    pairs: List[Tuple[Optional[Dict], Optional[Dict]]] = []
    for key in sorted(by_key):
        slot = by_key[key]
        pairs.append((slot.get("vh"), slot.get("vl")))
    return pairs, singles


def _write_fasta(path: Path, records: Sequence[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for record_id, seq in records:
        lines.append(f">{record_id}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _spearman(xs: Sequence[Optional[float]], ys: Sequence[Optional[float]]) -> Optional[float]:
    pairs = [
        (float(x), float(y))
        for x, y in zip(xs, ys)
        if x is not None and y is not None
    ]
    n = len(pairs)
    if n < 3:
        return None

    def _rank(values: List[float]) -> List[float]:
        order = sorted(range(n), key=lambda i: values[i])
        ranks = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and values[order[j + 1]] == values[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                ranks[order[k]] = avg
            i = j + 1
        return ranks

    rx = _rank([p[0] for p in pairs])
    ry = _rank([p[1] for p in pairs])
    d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
    denom = n * (n * n - 1)
    if denom == 0:
        return None
    return _round(1.0 - 6.0 * d2 / denom, 4)


# ---------------------------------------------------------------------------
# Sapiens: full prediction suite
# ---------------------------------------------------------------------------


def run_sapiens_record(
    record_id: str,
    seq: str,
    sapiens_chain_type: str,
    with_embeddings: bool = True,
    metadata_positions: Optional[List[str]] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """One model pass -> mean / per-position / full matrix / best / masked / embeddings."""
    if verbose:
        print(f"  [sapiens] {record_id}", file=sys.stderr, flush=True)

    if with_embeddings:
        probs, embeddings = predict_scores(
            seq, chain_type=sapiens_chain_type, return_embeddings=True
        )
    else:
        probs = predict_scores(seq, chain_type=sapiens_chain_type)
        embeddings = None

    per_position: List[Dict[str, Any]] = []
    matrix_rows: List[Dict[str, Any]] = []
    observed: List[float] = []
    best_chars: List[str] = []
    for idx, residue in enumerate(seq):
        row = probs.iloc[idx]
        score = float(row[residue]) if residue in probs.columns else None
        if score is not None:
            observed.append(score)
        best_aa = str(row.idxmax())
        best_chars.append(best_aa)
        label = (
            metadata_positions[idx]
            if metadata_positions and idx < len(metadata_positions)
            else str(idx + 1)
        )
        entry = {
            "position": idx + 1,
            "numbering": label,
            "residue": residue,
            "score": _round(score),
        }
        per_position.append(entry)
        matrix_row: Dict[str, Any] = {
            "id": record_id,
            "chain": sapiens_chain_type,
            "position": idx + 1,
            "numbering": label,
            "input_aa": residue,
            "sapiens_score": _round(score),
        }
        for aa in AA_ORDER:
            matrix_row[aa] = _round(float(row[aa])) if aa in probs.columns else None
        matrix_rows.append(matrix_row)

    best_sequence = "".join(best_chars)
    # Same rule as sapiens.predict.predict_masked.
    masked_sequence = "".join(
        b if a in "*X" else a for a, b in zip(seq, best_sequence)
    )

    result: Dict[str, Any] = {
        "id": record_id,
        "chain_type": sapiens_chain_type,
        "length": len(seq),
        "sapiens_mean": _round(sum(observed) / len(observed)) if observed else None,
        "sapiens_min": _round(min(observed)) if observed else None,
        "sapiens_max": _round(max(observed)) if observed else None,
        "per_position": per_position,
        "matrix_rows": matrix_rows,
        "best_score_sequence": best_sequence,
        "masked_sequence": masked_sequence,
        "n_below_0_5": sum(1 for s in observed if s < 0.5),
    }
    if observed:
        hotspots = sorted(
            (p for p in per_position if p["score"] is not None),
            key=lambda p: p["score"],
        )[:10]
        result["hotspots"] = hotspots

    if embeddings is not None:
        emb = np.asarray(embeddings)
        if emb.ndim == 2:  # single layer -> (L, D)
            emb = emb[None, ...]
        residue_layers = [
            [[_round(float(v), 5) for v in vec] for vec in layer] for layer in emb
        ]
        seq_emb = emb.mean(axis=-2)
        result["embeddings"] = {
            "n_layers": int(emb.shape[0]),
            "residue_embedding_layers": residue_layers,
            "sequence_embedding_layers": [
                [_round(float(v), 6) for v in layer] for layer in seq_emb
            ],
        }
    return result


def run_sapiens_full(
    records: Sequence[Tuple[str, str]],
    numbered: Sequence[Dict[str, Any]],
    with_embeddings: bool = True,
    metadata: Optional[Dict[str, Any]] = None,
    verbose: bool = False,
) -> Dict[str, Dict[str, Any]]:
    numbering_map: Dict[str, List[str]] = {}
    if metadata:
        raw = metadata.get("variants", metadata)
        if isinstance(raw, dict):
            numbering_map = {
                k: v for k, v in raw.items() if isinstance(v, list)
            }
    results: Dict[str, Dict[str, Any]] = {}
    by_id = {e["id"]: e for e in numbered}
    for record_id, seq in records:
        entry = by_id.get(record_id, {})
        labels = numbering_map.get(entry.get("variant") or record_id)
        results[record_id] = run_sapiens_record(
            record_id,
            seq,
            entry.get("sapiens_chain_type") or _chain_type_for_sapiens(None, record_id),
            with_embeddings=with_embeddings,
            metadata_positions=labels,
            verbose=verbose,
        )
    return results


def sapiens_humanized_fasta(
    numbered: Sequence[Dict[str, Any]],
    scheme: str,
    cdr_definition: str,
    iterations: int = 1,
    humanize_cdrs: bool = False,
    backmutate_vernier: bool = True,
    model_version: str = "latest",
    verbose: bool = False,
) -> List[Tuple[str, str]]:
    if humanize_cdrs:
        backmutate_vernier = False  # BioPhi forbids combining both
    out: List[Tuple[str, str]] = []
    for entry in numbered:
        chain = entry.get("chain")
        if chain is None:
            continue
        if verbose:
            print(f"  [humanize:sapiens] {entry['id']}", file=sys.stderr, flush=True)
        params = SapiensHumanizationParams(
            scheme=scheme,
            cdr_definition=cdr_definition,
            model_version=model_version,
            humanize_cdrs=humanize_cdrs,
            backmutate_vernier=backmutate_vernier,
            iterations=iterations,
        )
        result = sapiens_humanize_chain(chain, params)
        out.append((entry["id"], result.humanized_chain.seq))
    return out


# ---------------------------------------------------------------------------
# OASis humanness: full analysis
# ---------------------------------------------------------------------------


def _peptide_to_dict(peptide: Any) -> Dict[str, Any]:
    return {
        "seq": peptide.seq,
        "num_oas_subjects": peptide.num_oas_subjects,
        "fraction_oas_subjects": _round(peptide.fraction_oas_subjects),
        "num_oas_occurrences": peptide.num_oas_occurrences,
    }


def run_oasis_chain(chain: Chain, oasis_params: OASisParams) -> Dict[str, Any]:
    mfs = oasis_params.min_fraction_subjects
    hum = get_chain_humanness(chain, params=oasis_params)

    positional = []
    for pos, aa, nonhuman in hum.get_positional_humanness(mfs):
        positional.append(
            {
                "position": _pos_str(pos),
                "aa": aa,
                "nonhuman_peptides": [_peptide_to_dict(p) for p in nonhuman],
            }
        )

    top_freqs = []
    for rank_row in hum.get_top_freqs(3):
        top_freqs.append(
            [
                {"position": _pos_str(pos), "aa": aa, "freq": _round(freq, 4)}
                for pos, aa, freq in rank_row
            ]
        )

    chain_freq = {
        _pos_str(pos): {aa: _round(f, 4) for aa, f in freqs.items()}
        for pos, freqs in hum.chain_type_residue_frequency.items()
        if freqs
    }
    family_freq = {
        _pos_str(pos): {aa: _round(f, 4) for aa, f in freqs.items()}
        for pos, freqs in hum.germline_family_residue_frequency.items()
        if freqs
    }

    summary = {
        "id": chain.name,
        "chain_type": chain.chain_type,
        "oasis_identity": _round(hum.get_oasis_identity(mfs)),
        "oasis_percentile": _round(hum.get_oasis_percentile(mfs)),
        "germline_content": _round(hum.get_germline_content()),
        "num_peptides": hum.get_num_peptides(),
        "num_human_peptides": hum.get_num_human_peptides(mfs),
        "num_nonhuman_peptides": hum.get_num_nonhuman_peptides(mfs),
        "peptide_length": hum.get_peptide_length(),
        "v_germline_names": list(hum.v_germline_names),
        "j_germline_names": list(hum.j_germline_names),
        "v_germline_family": hum.v_germline_family,
        "v_germline_suffix": hum.v_germline_suffix,
        "oasis_curve": {str(k): _round(v) for k, v in hum.get_oasis_curve().items()},
        "top_freqs": top_freqs,
        "positional_humanness": positional,
        "chain_type_residue_frequency": chain_freq,
        "germline_family_residue_frequency": family_freq,
    }
    return {"summary": summary, "humanness": hum, "peptide_df": hum.to_peptide_dataframe()}


def run_oasis_antibody(
    vh: Optional[Chain], vl: Optional[Chain], oasis_params: OASisParams
) -> Optional[Dict[str, Any]]:
    if vh is None and vl is None:
        return None
    mfs = oasis_params.min_fraction_subjects
    ab = get_antibody_humanness(vh=vh, vl=vl, params=oasis_params)
    return {
        "oasis_identity": _round(ab.get_oasis_identity(mfs)),
        "oasis_percentile": _round(ab.get_oasis_percentile(mfs)),
        "germline_content": _round(ab.get_germline_content()),
        "num_peptides": ab.get_num_peptides(),
        "num_human_peptides": ab.get_num_human_peptides(mfs),
        "num_nonhuman_peptides": ab.get_num_nonhuman_peptides(mfs),
        "oasis_curve": {str(k): _round(v) for k, v in ab.get_oasis_curve().items()},
    }


def run_humanness_full(
    numbered: Sequence[Dict[str, Any]],
    oasis_db: Optional[str],
    min_fraction_subjects: float = 0.1,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Per-chain + paired-antibody OASis with all BioPhi detail sheets."""
    oasis_params: Optional[OASisParams] = None
    if oasis_db:
        db_path = Path(oasis_db)
        if db_path.suffix == ".gz":
            raise ToolError(f"OASis database must be unzipped first: gunzip {db_path}")
        if not db_path.exists():
            raise ToolError(f"OASis database not found: {db_path}")
        oasis_params = OASisParams(
            oasis_db_path=str(db_path.resolve()),
            min_fraction_subjects=min_fraction_subjects,
        )

    per_chain: Dict[str, Dict[str, Any]] = {}
    peptide_frames: List[pd.DataFrame] = []
    if oasis_params is None:
        return {
            "available": False,
            "reason": "OASis database not provided or missing",
            "per_chain": {},
            "per_antibody": {},
            "peptide_frames": {},
        }

    for entry in numbered:
        chain = entry.get("chain")
        if chain is None:
            continue
        if verbose:
            print(f"  [oasis] {entry['id']}", file=sys.stderr, flush=True)
        try:
            data = run_oasis_chain(chain, oasis_params)
        except Exception as exc:
            print(f"  OASis error for {entry['id']}: {exc}", file=sys.stderr)
            continue
        summary = data["summary"]
        summary["id"] = entry["id"]
        summary["variant"] = entry.get("variant", "")
        summary["base"] = entry.get("base", "")
        per_chain[entry["id"]] = summary
        df = data["peptide_df"].copy()
        df.insert(0, "record_id", entry["id"])
        peptide_frames.append(df)

    pairs, singles = _group_pairs(numbered)
    per_antibody: Dict[str, Dict[str, Any]] = {}
    for vh_entry, vl_entry in pairs:
        vh = vh_entry.get("chain") if vh_entry else None
        vl = vl_entry.get("chain") if vl_entry else None
        if vh is None and vl is None:
            continue
        anchor = vh_entry or vl_entry
        if anchor is None:
            continue
        key = anchor.get("pair_key") or anchor["id"]
        if verbose:
            print(f"  [oasis:antibody] {key}", file=sys.stderr, flush=True)
        try:
            ab_summary = run_oasis_antibody(vh, vl, oasis_params)
        except Exception as exc:
            print(f"  OASis antibody error for {key}: {exc}", file=sys.stderr)
            continue
        if ab_summary:
            ab_summary["pair_key"] = key
            ab_summary["vh_id"] = vh_entry["id"] if vh_entry else None
            ab_summary["vl_id"] = vl_entry["id"] if vl_entry else None
            per_antibody[key] = ab_summary
    for entry in singles:
        chain = entry.get("chain")
        if chain is None:
            continue
        try:
            ab_summary = run_oasis_antibody(
                chain if chain.is_heavy_chain() else None,
                chain if chain.is_light_chain() else None,
                oasis_params,
            )
        except Exception as exc:
            print(f"  OASis antibody error for {entry['id']}: {exc}", file=sys.stderr)
            continue
        if ab_summary:
            ab_summary["pair_key"] = None
            ab_summary["vh_id"] = entry["id"] if chain.is_heavy_chain() else None
            ab_summary["vl_id"] = entry["id"] if chain.is_light_chain() else None
            per_antibody[entry["id"]] = ab_summary

    return {
        "available": True,
        "min_fraction_subjects": min_fraction_subjects,
        "db": _db_info(Path(oasis_db), with_md5=False) if oasis_db else None,
        "per_chain": per_chain,
        "per_antibody": per_antibody,
        "peptide_frames": peptide_frames,
    }


# ---------------------------------------------------------------------------
# Humanization: Sapiens / CDR grafting / Manual / full antibody
# ---------------------------------------------------------------------------


def _mutations_between(
    parental: Chain, humanized: Chain
) -> List[Dict[str, Any]]:
    mutations: List[Dict[str, Any]] = []
    for pos, aa_p in parental.positions.items():
        aa_h = humanized.positions.get(pos)
        if aa_h is None:
            continue
        if aa_h != aa_p:
            mutations.append(
                {
                    "position": _pos_str(pos),
                    "parental": aa_p,
                    "humanized": aa_h,
                    "region": pos.get_region() if hasattr(pos, "get_region") else "",
                    "in_cdr": bool(pos.is_in_cdr()) if hasattr(pos, "is_in_cdr") else None,
                }
            )
    return mutations


def _top_rank_scores(result: Any, rank: int = 0, limit: int = 0,
                     next_scores: bool = False) -> List[Dict[str, Any]]:
    try:
        rows = result.get_top_scores(n=rank + 1, next=next_scores)[rank]
    except Exception:
        return []
    out = []
    for pos, aa, score in rows:
        out.append(
            {"position": _pos_str(pos), "aa": aa, "score": _round(score, 4)}
        )
    if limit:
        out = out[:limit]
    return out


def _humanization_entry(
    entry: Dict[str, Any],
    result: Any,
    mode: str,
    oasis_params: Optional[OASisParams] = None,
    germline_info: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    parental = result.parental_chain
    humanized = result.humanized_chain
    data: Dict[str, Any] = {
        "mode": mode,
        "id": entry["id"],
        "chain_type": parental.chain_type,
        "parental_sequence": parental.seq,
        "humanized_sequence": humanized.seq,
        "num_mutations": int(result.num_mutations()),
        "mutations": _mutations_between(parental, humanized),
        "alignment": result.get_alignment_string(),
        "top_scores": _top_rank_scores(result, rank=0),
        "second_scores": _top_rank_scores(result, rank=1),
        "third_scores": _top_rank_scores(result, rank=2),
        # Proposals computed on the *humanized* sequence (iterate further).
        "next_top_scores": _top_rank_scores(result, rank=0, next_scores=True),
        "next_second_scores": _top_rank_scores(result, rank=1, next_scores=True),
        "changed": humanized.seq != parental.seq,
    }
    if germline_info:
        data["germline"] = germline_info
    if oasis_params is not None and data["changed"]:
        try:
            hum_parent = run_oasis_chain(parental, oasis_params)["summary"]
            hum_after = run_oasis_chain(humanized, oasis_params)["summary"]
            data["oasis_before"] = {
                k: hum_parent[k]
                for k in ("oasis_identity", "oasis_percentile", "germline_content",
                          "num_nonhuman_peptides")
            }
            data["oasis_after"] = {
                k: hum_after[k]
                for k in ("oasis_identity", "oasis_percentile", "germline_content",
                          "num_nonhuman_peptides")
            }
        except Exception as exc:
            data["oasis_error"] = str(exc)
    return data


def _germline_candidates(chain: Chain, limit: int = 5) -> Dict[str, Any]:
    try:
        v_chains, j_chains = chain.find_human_germlines(limit=limit)
        return {
            "v_candidates": [c.name for c in v_chains],
            "j_candidates": [c.name for c in j_chains],
            "v_top": v_chains[0].name if v_chains else None,
            "j_top": j_chains[0].name if j_chains else None,
        }
    except Exception as exc:
        return {"error": str(exc)}


def run_humanization_modes(
    numbered: Sequence[Dict[str, Any]],
    modes: Sequence[str],
    scheme: str,
    cdr_definition: str,
    iterations: int = 1,
    humanize_cdrs: bool = False,
    backmutate_vernier: bool = True,
    model_version: str = "latest",
    heavy_v_germline: str = "auto",
    light_v_germline: str = "auto",
    sapiens_iterations: int = 0,
    oasis_params: Optional[OASisParams] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    modes = list(modes)
    per_record: Dict[str, Dict[str, Any]] = {}
    antibody_results: List[Dict[str, Any]] = []

    for entry in numbered:
        chain = entry.get("chain")
        if chain is None:
            continue
        record_modes: Dict[str, Any] = {}

        if "sapiens" in modes:
            if verbose:
                print(f"  [humanize:sapiens] {entry['id']}", file=sys.stderr, flush=True)
            params = SapiensHumanizationParams(
                scheme=scheme,
                cdr_definition=cdr_definition,
                model_version=model_version,
                humanize_cdrs=humanize_cdrs,
                backmutate_vernier=(backmutate_vernier and not humanize_cdrs),
                iterations=iterations,
            )
            result = sapiens_humanize_chain(chain, params)
            record_modes["sapiens"] = _humanization_entry(
                entry, result, "sapiens", oasis_params
            )

        if "cdr_grafting" in modes:
            if verbose:
                print(f"  [humanize:cdr_grafting] {entry['id']}", file=sys.stderr, flush=True)
            v_gene = heavy_v_germline if chain.is_heavy_chain() else light_v_germline
            germline_info = _germline_candidates(chain)
            if v_gene != "auto":
                germline_info["requested_v_germline"] = v_gene
            params = CDRGraftingHumanizationParams(
                scheme=scheme,
                cdr_definition=cdr_definition,
                heavy_v_germline=heavy_v_germline,
                light_v_germline=light_v_germline,
                backmutate_vernier=backmutate_vernier,
                sapiens_iterations=sapiens_iterations,
            )
            result = cdr_grafting_humanize_chain(chain, params)
            record_modes["cdr_grafting"] = _humanization_entry(
                entry, result, "cdr_grafting", oasis_params, germline_info
            )

        if "manual" in modes:
            if verbose:
                print(f"  [humanize:manual] {entry['id']}", file=sys.stderr, flush=True)
            params = ManualHumanizationParams(
                scheme=scheme, cdr_definition=cdr_definition
            )
            result = humanize_chain(chain, params)
            record_modes["manual"] = _humanization_entry(
                entry, result, "manual", oasis_params
            )
            record_modes["manual"]["note"] = (
                "Designer mode: sequence unchanged; top_scores are ranked "
                "humanizing proposals for manual selection"
            )

        per_record[entry["id"]] = record_modes

    if "full_antibody" in modes:
        pairs, _singles = _group_pairs(numbered)
        for vh_entry, vl_entry in pairs:
            if vh_entry is None or vl_entry is None:
                continue
            if verbose:
                print(
                    f"  [humanize:antibody] {vh_entry['id']} + {vl_entry['id']}",
                    file=sys.stderr,
                    flush=True,
                )
            try:
                params = SapiensHumanizationParams(
                    scheme=scheme,
                    cdr_definition=cdr_definition,
                    model_version=model_version,
                    humanize_cdrs=humanize_cdrs,
                    backmutate_vernier=(backmutate_vernier and not humanize_cdrs),
                    iterations=iterations,
                )
                ab = humanize_antibody(
                    vh=vh_entry["chain"], vl=vl_entry["chain"], params=params
                )
                item: Dict[str, Any] = {
                    "pair_key": vh_entry.get("pair_key"),
                    "vh_id": vh_entry["id"],
                    "vl_id": vl_entry["id"],
                    "alignment": ab.get_alignment_string(),
                    "vh": _humanization_entry(vh_entry, ab.vh, "full_antibody", oasis_params),
                    "vl": _humanization_entry(vl_entry, ab.vl, "full_antibody", oasis_params),
                }
                antibody_results.append(item)
            except Exception as exc:
                print(
                    f"  full-antibody error for {vh_entry['id']}: {exc}",
                    file=sys.stderr,
                )

    return {"per_record": per_record, "full_antibody": antibody_results}


def humanize_fasta_from_results(
    results: Dict[str, Any], mode: str
) -> List[Tuple[str, str]]:
    out: List[Tuple[str, str]] = []
    for record_id, modes in results.get("per_record", {}).items():
        data = modes.get(mode)
        if data:
            out.append((record_id, data["humanized_sequence"]))
    if mode == "full_antibody":
        for item in results.get("full_antibody", []):
            for side in ("vh", "vl"):
                data = item.get(side)
                if data:
                    out.append((data["id"], data["humanized_sequence"]))
    return out


# ---------------------------------------------------------------------------
# OASis XLSX workbook (replicates `biophi oasis` report sheets)
# ---------------------------------------------------------------------------


def build_oasis_workbook(
    numbered: Sequence[Dict[str, Any]],
    humanness: Dict[str, Any],
    out_path: Path,
    min_fraction_subjects: float,
    humanized_seqs: Optional[Dict[str, str]] = None,
) -> Optional[Path]:
    if not humanness.get("available"):
        return None
    per_chain = humanness.get("per_chain", {})
    if not per_chain:
        return None

    threshold = f"{min_fraction_subjects:.0%}"
    overview_rows = []
    for key, ab in humanness.get("per_antibody", {}).items():
        row: Dict[str, Any] = {
            "Antibody": key or ab.get("vh_id") or ab.get("vl_id"),
            "Threshold": threshold,
            "OASis Percentile": ab.get("oasis_percentile"),
            "OASis Identity": ab.get("oasis_identity"),
            "Germline Content": ab.get("germline_content"),
        }
        vh_summary = per_chain.get(ab.get("vh_id") or "")
        vl_summary = per_chain.get(ab.get("vl_id") or "")
        if vh_summary:
            row.update(
                {
                    "Heavy V Germline": (vh_summary.get("v_germline_names") or [""])[0],
                    "Heavy J Germline": (vh_summary.get("j_germline_names") or [""])[0],
                    "Heavy OASis Percentile": vh_summary.get("oasis_percentile"),
                    "Heavy OASis Identity": vh_summary.get("oasis_identity"),
                    "Heavy Non-human peptides": vh_summary.get("num_nonhuman_peptides"),
                    "Heavy Germline Content": vh_summary.get("germline_content"),
                }
            )
        if vl_summary:
            row.update(
                {
                    "Light V Germline": (vl_summary.get("v_germline_names") or [""])[0],
                    "Light J Germline": (vl_summary.get("j_germline_names") or [""])[0],
                    "Light OASis Percentile": vl_summary.get("oasis_percentile"),
                    "Light OASis Identity": vl_summary.get("oasis_identity"),
                    "Light Non-human peptides": vl_summary.get("num_nonhuman_peptides"),
                    "Light Germline Content": vl_summary.get("germline_content"),
                }
            )
        overview_rows.append(row)

    overview = pd.DataFrame(overview_rows).set_index("Antibody")

    # Chain sequence sheets + OASis curves + per-antibody peptide sheets.
    vh_hums, vl_hums = [], []
    curve_index, curve_rows = [], []
    antibody_sheets: Dict[str, pd.DataFrame] = {}
    for entry in numbered:
        if entry["id"] not in per_chain:
            continue
        chain = entry.get("chain")
        if chain is None:
            continue
        try:
            hum = get_chain_humanness(
                chain, params=_oasis_params_from(humanness)
            )
        except Exception:
            continue
        if chain.is_heavy_chain():
            vh_hums.append(hum)
        elif chain.is_light_chain():
            vl_hums.append(hum)
        curve_index.append(entry["id"])
        curve_rows.append(hum.get_oasis_curve())
        antibody_sheets[entry["id"]] = hum.to_peptide_dataframe()

    sheets: Dict[str, pd.DataFrame] = {"Overview": overview}
    if vh_hums:
        sheets["VH"] = ChainHumanness.to_sequence_dataframe(vh_hums, species=False)
    if vl_hums:
        sheets["VL"] = ChainHumanness.to_sequence_dataframe(vl_hums, species=False)
    if curve_rows:
        curves = pd.DataFrame(curve_rows, index=curve_index)
        curves.index.name = "Antibody"
        sheets["OASis Curves"] = curves
    for name, df in antibody_sheets.items():
        sheets[name] = df

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name[:31])
    return out_path


def _oasis_params_from(humanness: Dict[str, Any]) -> OASisParams:
    db = humanness.get("db") or {}
    return OASisParams(
        oasis_db_path=str(db.get("path", DEFAULT_DB)),
        min_fraction_subjects=float(
            humanness.get("min_fraction_subjects", 0.1)
        ),
    )


# ---------------------------------------------------------------------------
# Cross-validation: BioPhi vs pipeline
# ---------------------------------------------------------------------------


def _pipeline_index(pipeline_result: Optional[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    index: Dict[str, Dict[str, Any]] = {}
    if not pipeline_result:
        return index
    for chain in pipeline_result.get("chains", []):
        backmutations = {
            bm.get("position"): bm
            for bm in chain.get("backmutations", [])
            if isinstance(bm, dict)
        }
        for variant in chain.get("variants", []):
            name = variant.get("name")
            if not name:
                continue
            bm_positions = [
                p for p in (variant.get("backmutations") or []) if isinstance(p, str)
            ]
            composites: List[float] = []
            for p in bm_positions:
                bm = backmutations.get(p)
                comp = bm.get("composite") if bm else None
                if comp is not None:
                    composites.append(float(comp))
            tiers: Dict[str, int] = {}
            for p in bm_positions:
                tier = (backmutations.get(p) or {}).get("tier") or "unknown"
                tiers[tier] = tiers.get(tier, 0) + 1
            index[name] = {
                "chain_name": chain.get("name"),
                "chain_type": chain.get("chain_type"),
                "human_likeness": chain.get("human_likeness", {}),
                "cvi_homology": chain.get("cvi_homology"),
                "germline_v": (chain.get("germline") or {}).get("v_gene"),
                "variant_description": variant.get("description", ""),
                "n_backmutations": len(bm_positions),
                "backmutations": bm_positions,
                "backmutation_composite_mean": (
                    _round(sum(composites) / len(composites)) if composites else None
                ),
                "tier_distribution": tiers,
            }
    return index


def build_cross_validation(
    sapiens_results: Dict[str, Dict[str, Any]],
    humanness: Dict[str, Any],
    pipeline_result: Optional[Dict[str, Any]],
    numbered: Sequence[Dict[str, Any]],
) -> Dict[str, Any]:
    per_variant: List[Dict[str, Any]] = []
    pipeline_index = _pipeline_index(pipeline_result)
    per_chain = humanness.get("per_chain", {}) if humanness else {}

    for entry in numbered:
        rid = entry["id"]
        sapiens = sapiens_results.get(rid, {})
        oasis = per_chain.get(rid, {})
        variant_name = entry.get("variant") or ""
        pipe = pipeline_index.get(variant_name) or pipeline_index.get(rid)
        row: Dict[str, Any] = {
            "id": rid,
            "variant": variant_name,
            "chain_type": entry.get("chain_type") or entry.get("chain_token"),
            "sapiens_mean": sapiens.get("sapiens_mean"),
            "oasis_identity": oasis.get("oasis_identity"),
            "oasis_percentile": oasis.get("oasis_percentile"),
            "oasis_germline_content": oasis.get("germline_content"),
            "pipeline": None,
        }
        if pipe:
            row["pipeline"] = {
                "chain_name": pipe["chain_name"],
                "human_likeness": pipe["human_likeness"],
                "cvi_homology": pipe["cvi_homology"],
                "germline_v": pipe["germline_v"],
                "n_backmutations": pipe["n_backmutations"],
                "backmutation_composite_mean": pipe["backmutation_composite_mean"],
                "tier_distribution": pipe["tier_distribution"],
            }
        per_variant.append(row)

    sapiens_vals = [r["sapiens_mean"] for r in per_variant]
    oasis_vals = [r["oasis_identity"] for r in per_variant]
    bm_vals = [
        (r["pipeline"] or {}).get("n_backmutations") for r in per_variant
    ]
    cvi_vals = [(r["pipeline"] or {}).get("cvi_homology") for r in per_variant]

    notes = [
        "Sapiens mean and OASis identity are BioPhi humanness metrics; "
        "pipeline human_likeness / CVI are germline-identity metrics - "
        "they are different estimators and are reported side by side.",
        "Spearman rho measures rank agreement across variants (1 = same order).",
    ]
    spearman = {
        "sapiens_vs_oasis": _spearman(sapiens_vals, oasis_vals),
        "sapiens_vs_pipeline_n_backmutations": _spearman(sapiens_vals, bm_vals),
        "oasis_vs_pipeline_n_backmutations": _spearman(oasis_vals, bm_vals),
        "sapiens_vs_pipeline_cvi": _spearman(sapiens_vals, cvi_vals),
    }

    # Flag rank disagreements between Sapiens and OASis (both BioPhi metrics).
    ranked = [
        (r["id"], r["sapiens_mean"], r["oasis_identity"])
        for r in per_variant
        if r["sapiens_mean"] is not None and r["oasis_identity"] is not None
    ]
    disagreements = []
    if len(ranked) >= 3:
        by_sapiens = [r[0] for r in sorted(ranked, key=lambda x: -x[1])]
        by_oasis = [r[0] for r in sorted(ranked, key=lambda x: -x[2])]
        for i, rid in enumerate(by_sapiens):
            j = by_oasis.index(rid)
            if abs(i - j) >= max(2, len(ranked) // 3):
                disagreements.append(
                    {
                        "id": rid,
                        "sapiens_rank": i + 1,
                        "oasis_rank": j + 1,
                    }
                )

    return {
        "per_variant": per_variant,
        "spearman": spearman,
        "rank_disagreements": disagreements,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Report writers
# ---------------------------------------------------------------------------


def _variant_summary_rows(
    numbered: Sequence[Dict[str, Any]],
    sapiens_results: Dict[str, Dict[str, Any]],
    humanness: Dict[str, Any],
    pipeline_result: Optional[Dict[str, Any]],
    humanization: Optional[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    pipeline_index = _pipeline_index(pipeline_result)
    per_chain = humanness.get("per_chain", {}) if humanness else {}
    rows: List[Dict[str, Any]] = []
    for entry in numbered:
        rid = entry["id"]
        sapiens = sapiens_results.get(rid, {})
        oasis = per_chain.get(rid, {})
        variant_name = entry.get("variant") or ""
        pipe = pipeline_index.get(variant_name)
        modes = (humanization or {}).get("per_record", {}).get(rid, {})
        row = {
            "id": rid,
            "base": entry.get("base", ""),
            "variant": variant_name,
            "chain_type": entry.get("chain_type") or entry.get("chain_token"),
            "length": entry.get("length"),
            "sapiens_mean": sapiens.get("sapiens_mean"),
            "sapiens_min": sapiens.get("sapiens_min"),
            "oasis_identity": oasis.get("oasis_identity"),
            "oasis_percentile": oasis.get("oasis_percentile"),
            "germline_content": oasis.get("germline_content"),
            "v_germline": (oasis.get("v_germline_names") or [None])[0],
            "pipeline_n_backmutations": (
                pipe.get("n_backmutations") if pipe else None
            ),
            "pipeline_cvi_homology": pipe.get("cvi_homology") if pipe else None,
            "sapiens_humanize_mutations": (
                modes.get("sapiens", {}).get("num_mutations") if modes else None
            ),
            "cdr_grafting_mutations": (
                modes.get("cdr_grafting", {}).get("num_mutations") if modes else None
            ),
        }
        rows.append(row)
    return rows


def _fmt(value: Any) -> str:
    if value is None or value == "":
        return "N/A"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def _md_cell(value: Any) -> str:
    """Escape pipes so FASTA ids like ``4D5_VH|H_V0|pure`` stay one cell."""
    if value is None or value == "":
        return "N/A"
    if isinstance(value, float):
        text = f"{value:.4f}"
    else:
        text = str(value)
    return text.replace("|", "\\|")


def write_markdown_report(
    path: Path,
    stem: str,
    input_path: Path,
    numbered: Sequence[Dict[str, Any]],
    sapiens_results: Dict[str, Dict[str, Any]],
    humanness: Dict[str, Any],
    humanization: Optional[Dict[str, Any]],
    cross_validation: Optional[Dict[str, Any]],
    parameters: Dict[str, Any],
    environment: Dict[str, Any],
    outputs: Dict[str, str],
) -> None:
    summary_rows = _variant_summary_rows(
        numbered,
        sapiens_results,
        humanness,
        parameters.get("pipeline_result_obj"),
        humanization,
    )
    lines: List[str] = [
        f"# BioPhi Comprehensive Analysis: {input_path.name}",
        "",
        f"- Generated: {_now()}",
        f"- Tool: tools/biophi/biophi_comprehensive.py v{TOOL_VERSION}",
        f"- BioPhi: {environment.get('biophi', 'unknown')} | "
        f"Sapiens: {environment.get('sapiens', 'unknown')} | "
        f"abnumber: {environment.get('abnumber', 'unknown')}",
        f"- Numbering: {parameters.get('scheme')} / CDR definition: "
        f"{parameters.get('cdr_definition')}",
        f"- Records: {len(numbered)}",
        "",
        "## Variant Summary",
        "",
        "| ID | Variant | Chain | Sapiens mean | OASis identity | OASis percentile "
        "| Germline content | Pipeline back-mutations | Sapiens humanize #mut | "
        "CDR graft #mut |",
        "|----|---------|-------|--------------|----------------|-------------------"
        "|------------------|-------------------------|----------------------|"
        "----------------|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {_md_cell(row['id'])} | {_md_cell(row['variant'])} | {row['chain_type']} | "
            f"{_fmt(row['sapiens_mean'])} | {_fmt(row['oasis_identity'])} | "
            f"{_fmt(row['oasis_percentile'])} | {_fmt(row['germline_content'])} | "
            f"{_fmt(row['pipeline_n_backmutations'])} | "
            f"{_fmt(row['sapiens_humanize_mutations'])} | "
            f"{_fmt(row['cdr_grafting_mutations'])} |"
        )

    lines += ["", "## Per-Variant Details", ""]
    per_chain = humanness.get("per_chain", {}) if humanness else {}
    for entry in numbered:
        rid = entry["id"]
        sapiens = sapiens_results.get(rid, {})
        oasis = per_chain.get(rid, {})
        lines.append(f"### {rid}")
        if entry.get("description"):
            lines.append(f"- Description: {entry['description']}")
        lines.append(f"- Length: {entry.get('length')} | Chain type: {entry.get('chain_type')}")
        lines.append(
            f"- Sapiens: mean={_fmt(sapiens.get('sapiens_mean'))} "
            f"min={_fmt(sapiens.get('sapiens_min'))} max={_fmt(sapiens.get('sapiens_max'))} "
            f"positions<0.5={sapiens.get('n_below_0_5')}"
        )
        if sapiens.get("hotspots"):
            hot = ", ".join(
                f"{h['numbering']}({h['residue']})={h['score']}"
                for h in sapiens["hotspots"][:5]
            )
            lines.append(f"- Sapiens hotspots (lowest): {hot}")
        if oasis:
            v_names = ", ".join((oasis.get("v_germline_names") or [])[:3])
            j_names = ", ".join((oasis.get("j_germline_names") or [])[:3])
            lines.append(
                f"- OASis: identity={_fmt(oasis.get('oasis_identity'))} "
                f"percentile={_fmt(oasis.get('oasis_percentile'))} "
                f"germline_content={_fmt(oasis.get('germline_content'))} "
                f"non-human peptides={oasis.get('num_nonhuman_peptides')}/"
                f"{oasis.get('num_peptides')}"
            )
            lines.append(f"- V germline: {v_names} | J germline: {j_names}")
        modes = (humanization or {}).get("per_record", {}).get(rid, {})
        for mode_name, data in modes.items():
            if not isinstance(data, dict) or "num_mutations" not in data:
                continue
            lines.append(
                f"- Humanization [{mode_name}]: {data['num_mutations']} mutations"
            )
            if data.get("germline", {}).get("v_top"):
                lines.append(
                    f"  - germline used (top): {data['germline'].get('v_top')} / "
                    f"{data['germline'].get('j_top')}"
                )
            if data.get("oasis_before") and data.get("oasis_after"):
                b, a = data["oasis_before"], data["oasis_after"]
                lines.append(
                    f"  - OASis before/after: "
                    f"{_fmt(b.get('oasis_identity'))} -> {_fmt(a.get('oasis_identity'))} "
                    f"(percentile {_fmt(b.get('oasis_percentile'))} -> "
                    f"{_fmt(a.get('oasis_percentile'))})"
                )
        lines.append("")

    ab_results = (humanization or {}).get("full_antibody", [])
    if ab_results:
        lines += ["## Full-Antibody Humanization (VH+VL)", ""]
        for item in ab_results:
            lines.append(f"### {item.get('pair_key') or item.get('vh_id')}")
            lines.append(f"- VH: {item['vh_id']} ({item['vh']['num_mutations']} mut)")
            lines.append(f"- VL: {item['vl_id']} ({item['vl']['num_mutations']} mut)")
            lines.append("")

    if cross_validation and cross_validation.get("per_variant"):
        lines += [
            "## Cross-Validation: BioPhi vs Pipeline",
            "",
            "| ID | Sapiens | OASis identity | Pipeline human-likeness (kabat) "
            "| Pipeline CVI | Pipeline back-mutations | Composite mean |",
            "|----|---------|----------------|-------------------------------"
            "|--------------|-------------------------|----------------|",
        ]
        for row in cross_validation["per_variant"]:
            pipe = row.get("pipeline") or {}
            hl = (pipe.get("human_likeness") or {}).get("kabat")
            lines.append(
                f"| {_md_cell(row['id'])} | {_fmt(row['sapiens_mean'])} | "
                f"{_fmt(row['oasis_identity'])} | {_fmt(hl)} | "
                f"{_fmt(pipe.get('cvi_homology'))} | "
                f"{_fmt(pipe.get('n_backmutations'))} | "
                f"{_fmt(pipe.get('backmutation_composite_mean'))} |"
            )
        lines += ["", "### Spearman rank correlations", ""]
        for key, value in (cross_validation.get("spearman") or {}).items():
            lines.append(f"- {key}: {_fmt(value)}")
        if cross_validation.get("rank_disagreements"):
            lines += ["", "### Sapiens vs OASis rank disagreements", ""]
            for d in cross_validation["rank_disagreements"]:
                lines.append(
                    f"- {d['id']}: Sapiens rank {d['sapiens_rank']} vs "
                    f"OASis rank {d['oasis_rank']}"
                )
        lines += ["", "### Notes", ""]
        for note in cross_validation.get("notes", []):
            lines.append(f"- {note}")
        lines.append("")

    lines += ["## Output Files", ""]
    for kind, target in sorted(outputs.items()):
        lines.append(f"- {kind}: `{target}`")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_heatmap(
    sapiens_results: Dict[str, Dict[str, Any]]
) -> Dict[str, Any]:
    heatmap: Dict[str, Any] = {}
    for rid, data in sapiens_results.items():
        heatmap[rid] = [
            {
                "position": p["position"],
                "numbering": p.get("numbering"),
                "residue": p["residue"],
                "score": p["score"],
            }
            for p in data.get("per_position", [])
            if p.get("score") is not None
        ]
    return {"heatmap": heatmap, "generated": _now(), "tool_version": TOOL_VERSION}


def build_embeddings_payload(
    sapiens_results: Dict[str, Dict[str, Any]]
) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for rid, data in sapiens_results.items():
        if data.get("embeddings"):
            out[rid] = data["embeddings"]
    return {"embeddings": out, "generated": _now(), "tool_version": TOOL_VERSION}


# ---------------------------------------------------------------------------
# Shared command helpers
# ---------------------------------------------------------------------------


def _load_inputs(args: argparse.Namespace) -> Tuple[Path, List[Tuple[str, str]]]:
    input_path = Path(args.input)
    records = read_fasta(input_path)
    limit = getattr(args, "limit", None)
    if limit:
        records = records[:limit]
    return input_path, records


def _load_json(path: Optional[str]) -> Optional[Dict[str, Any]]:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        print(f"WARNING: JSON not found: {p}", file=sys.stderr)
        return None
    with p.open() as fh:
        return json.load(fh)


def _oasis_params_from_args(args: argparse.Namespace) -> Optional[OASisParams]:
    db = getattr(args, "oasis_db", None)
    run = getattr(args, "run_oasis", True)
    if not run or not db:
        return None
    p = Path(db)
    if p.suffix == ".gz":
        raise ToolError(f"OASis database must be unzipped first: gunzip {p}")
    if not p.exists():
        raise ToolError(f"OASis database not found: {p}")
    return OASisParams(
        oasis_db_path=str(p.resolve()),
        min_fraction_subjects=float(getattr(args, "min_fraction_subjects", 0.1)),
    )


def _write_scores_csv(path: Path, sapiens_results: Dict[str, Dict[str, Any]]) -> None:
    """Long-format full 20-aa matrix (BioPhi ``--scores-only`` equivalent)."""
    fieldnames = ["id", "chain", "position", "numbering", "input_aa",
                  "sapiens_score", *AA_ORDER]
    rows: List[Dict[str, Any]] = []
    for data in sapiens_results.values():
        # Rebuild matrix rows from per_position + we stored matrix only inside
        # run_sapiens_record - regenerate from per_position is impossible, so
        # keep a parallel list when available.
        rows.extend(data.get("matrix_rows", []))
    if not rows:
        return
    _write_csv(path, rows, fieldnames)


def _write_mean_csv(path: Path, sapiens_results: Dict[str, Dict[str, Any]]) -> None:
    rows = [
        {
            "id": rid,
            "chain": data.get("chain_type"),
            "sapiens_score": data.get("sapiens_mean"),
        }
        for rid, data in sapiens_results.items()
    ]
    _write_csv(path, rows, ["id", "chain", "sapiens_score"])


def _write_peptides_csv(path: Path, humanness: Dict[str, Any]) -> None:
    frames = humanness.get("peptide_frames") or []
    if not frames:
        return
    df = pd.concat(frames, ignore_index=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def _environment_report(args: argparse.Namespace) -> Dict[str, Any]:
    return _env_info(getattr(args, "env", DEFAULT_ENV))


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------


def cmd_check(args: argparse.Namespace) -> int:
    problems: List[str] = []
    warnings: List[str] = []
    report: Dict[str, Any] = {
        "tool": "tools/biophi/biophi_comprehensive.py",
        "tool_version": TOOL_VERSION,
        "timestamp": _now(),
        "python": sys.executable,
        "capabilities": {
            "sapiens": True,
            "biophi": True,
            "abnumber": True,
            "pandas": True,
            "openpyxl": True,
        },
    }
    try:
        import biophi
        import sapiens
        import abnumber
        report["versions"] = {
            "biophi": getattr(biophi, "__version__", "?"),
            "sapiens": getattr(sapiens, "__version__", "?"),
            "abnumber": getattr(abnumber, "__version__", "?"),
            "pandas": pd.__version__,
            "numpy": np.__version__,
        }
    except ImportError as exc:
        problems.append(f"BioPhi stack not importable: {exc}")

    hf = _hf_cache_info()
    report["sapiens_cache"] = hf
    if not hf["complete"]:
        warnings.append(
            "Sapiens VH/VL/tokenizer not fully cached; see README (hf-mirror)"
        )

    missing_hmmer = [
        name for name in ("hmmsearch", "hmmbuild", "hmmscan") if shutil.which(name) is None
    ]
    report["hmmer"] = {"missing": missing_hmmer}
    if missing_hmmer:
        warnings.append(
            "HMMER missing (needed for CDR-grafting germline search): "
            + ", ".join(missing_hmmer)
        )

    db = Path(getattr(args, "oasis_db", DEFAULT_DB))
    report["oasis_db"] = _db_info(db, with_md5=getattr(args, "md5", False))
    if not db.exists():
        warnings.append(f"OASis database not found: {db} (OASis steps will fail)")
    elif db.suffix == ".gz":
        problems.append(f"OASis database is still gzipped: {db} (run gunzip)")

    report["problems"] = problems
    report["warnings"] = warnings
    print(json.dumps(_jsonable(report), indent=2, ensure_ascii=False))
    return 1 if problems else 0


def cmd_sapiens(args: argparse.Namespace) -> int:
    input_path, records = _load_inputs(args)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    metadata = _load_json(getattr(args, "metadata", None))

    print("[1/3] Numbering sequences...", file=sys.stderr, flush=True)
    numbered = _number_records(records, args.scheme, args.cdr_definition)

    print("[2/3] Running full Sapiens predictions...", file=sys.stderr, flush=True)
    sapiens_results = run_sapiens_full(
        records,
        numbered,
        with_embeddings=not args.no_embeddings,
        metadata=metadata,
        verbose=args.verbose,
    )

    outputs: Dict[str, str] = {}
    print("[3/3] Writing outputs...", file=sys.stderr, flush=True)

    payload = {
        "tool": "tools/biophi/biophi_comprehensive.py",
        "tool_version": TOOL_VERSION,
        "mode": "sapiens",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
            "n_records": len(records),
        },
        "parameters": {
            "scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
            "embeddings": not args.no_embeddings,
        },
        "environment": _environment_report(args),
        "results": [_jsonable(v) for v in sapiens_results.values()],
    }
    json_path = outdir / f"{stem}_biophi_full_sapiens.json"
    _write_json(json_path, payload)
    outputs["sapiens json"] = str(json_path)

    mean_path = outdir / f"{stem}_biophi_sapiens.csv"
    _write_mean_csv(mean_path, sapiens_results)
    outputs["mean csv"] = str(mean_path)

    scores_path = outdir / f"{stem}_biophi_scores.csv"
    _write_scores_csv(scores_path, sapiens_results)
    if scores_path.exists():
        outputs["full score matrix csv"] = str(scores_path)

    if not args.no_embeddings:
        emb_path = outdir / f"{stem}_biophi_embeddings.json"
        _write_json(emb_path, build_embeddings_payload(sapiens_results))
        outputs["embeddings json"] = str(emb_path)

    if not args.no_humanize:
        print("  Generating Sapiens-humanized FASTA...", file=sys.stderr, flush=True)
        humanized = sapiens_humanized_fasta(
            numbered,
            args.scheme,
            args.cdr_definition,
            iterations=args.iterations,
            humanize_cdrs=args.humanize_cdrs,
            backmutate_vernier=not args.no_backmutate_vernier,
            model_version=args.model_version,
            verbose=args.verbose,
        )
        fa_path = outdir / f"{stem}_biophi_humanized.fa"
        _write_fasta(fa_path, humanized)
        outputs["humanized fasta"] = str(fa_path)

    heatmap_path = outdir / f"{stem}_biophi_heatmap.json"
    _write_json(heatmap_path, build_heatmap(sapiens_results))
    outputs["heatmap json"] = str(heatmap_path)

    md_path = outdir / f"{stem}_biophi_full_sapiens.md"
    lines = [
        f"# BioPhi Sapiens (full) report: {input_path.name}",
        "",
        f"- Generated: {payload['timestamp']}",
        f"- BioPhi: {payload['environment'].get('biophi', 'unknown')}",
        f"- Numbering: {args.scheme} / CDR: {args.cdr_definition}",
        "",
        "| ID | Chain | Length | Mean | Min | Max | #<0.5 |",
        "|----|-------|--------|------|-----|-----|-------|",
    ]
    for rid, data in sapiens_results.items():
        lines.append(
            f"| {_md_cell(rid)} | {data.get('chain_type')} | {data.get('length')} | "
            f"{_fmt(data.get('sapiens_mean'))} | {_fmt(data.get('sapiens_min'))} | "
            f"{_fmt(data.get('sapiens_max'))} | {data.get('n_below_0_5')} |"
        )
    lines += ["", "## Output files", ""]
    for kind, target in sorted(outputs.items()):
        lines.append(f"- {kind}: `{target}`")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")
    outputs["markdown"] = str(md_path)

    for target in sorted(set(outputs.values())):
        print(f"Wrote {target}", file=sys.stderr, flush=True)
    return 0


def cmd_humanness(args: argparse.Namespace) -> int:
    input_path, records = _load_inputs(args)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem

    print("[1/3] Numbering sequences...", file=sys.stderr, flush=True)
    numbered = _number_records(records, args.scheme, args.cdr_definition)

    print("[2/3] Running full OASis humanness analysis...", file=sys.stderr, flush=True)
    humanness = run_humanness_full(
        numbered,
        oasis_db=args.oasis_db if args.run_oasis else None,
        min_fraction_subjects=args.min_fraction_subjects,
        verbose=args.verbose,
    )
    if not humanness.get("available"):
        raise ToolError(humanness.get("reason", "OASis unavailable"))

    print("[3/3] Writing outputs...", file=sys.stderr, flush=True)
    outputs: Dict[str, str] = {}

    payload = {
        "tool": "tools/biophi/biophi_comprehensive.py",
        "tool_version": TOOL_VERSION,
        "mode": "humanness",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
        },
        "parameters": {
            "scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
            "min_fraction_subjects": args.min_fraction_subjects,
        },
        "oasis_db": humanness.get("db"),
        "environment": _environment_report(args),
        "per_chain": _jsonable(humanness.get("per_chain")),
        "per_antibody": _jsonable(humanness.get("per_antibody")),
    }
    json_path = outdir / f"{stem}_biophi_humanness.json"
    _write_json(json_path, payload)
    outputs["humanness json"] = str(json_path)

    peptides_path = outdir / f"{stem}_biophi_peptides.csv"
    _write_peptides_csv(peptides_path, humanness)
    if peptides_path.exists():
        outputs["peptide-level csv"] = str(peptides_path)

    xlsx_path = outdir / f"{stem}_biophi_oasis.xlsx"
    written = build_oasis_workbook(
        numbered, humanness, xlsx_path, args.min_fraction_subjects
    )
    if written:
        outputs["oasis xlsx workbook"] = str(written)

    # CSV summary of per-chain + per-antibody metrics
    summary_path = outdir / f"{stem}_biophi_humanness.csv"
    rows = []
    for rid, summary in humanness.get("per_chain", {}).items():
        rows.append(
            {
                "scope": "chain",
                "id": rid,
                "chain_type": summary.get("chain_type"),
                "oasis_identity": summary.get("oasis_identity"),
                "oasis_percentile": summary.get("oasis_percentile"),
                "germline_content": summary.get("germline_content"),
                "num_peptides": summary.get("num_peptides"),
                "num_nonhuman_peptides": summary.get("num_nonhuman_peptides"),
                "v_germline": (summary.get("v_germline_names") or [None])[0],
                "j_germline": (summary.get("j_germline_names") or [None])[0],
            }
        )
    for key, summary in humanness.get("per_antibody", {}).items():
        rows.append(
            {
                "scope": "antibody",
                "id": key,
                "chain_type": "paired",
                "oasis_identity": summary.get("oasis_identity"),
                "oasis_percentile": summary.get("oasis_percentile"),
                "germline_content": summary.get("germline_content"),
                "num_peptides": summary.get("num_peptides"),
                "num_nonhuman_peptides": summary.get("num_nonhuman_peptides"),
                "v_germline": None,
                "j_germline": None,
            }
        )
    _write_csv(
        summary_path,
        rows,
        [
            "scope", "id", "chain_type", "oasis_identity", "oasis_percentile",
            "germline_content", "num_peptides", "num_nonhuman_peptides",
            "v_germline", "j_germline",
        ],
    )
    outputs["summary csv"] = str(summary_path)

    md_path = outdir / f"{stem}_biophi_humanness.md"
    lines = [
        f"# BioPhi OASis humanness (full): {input_path.name}",
        "",
        f"- Generated: {payload['timestamp']}",
        f"- Threshold: >={args.min_fraction_subjects:.0%} of OAS subjects",
        f"- DB: `{(payload.get('oasis_db') or {}).get('path')}`",
        "",
        "| ID | Chain | OASis identity | Percentile | Germline content "
        "| Non-human peptides | Top V germline |",
        "|----|-------|----------------|------------|------------------"
        "|-------------------|---------------|",
    ]
    for rid, summary in humanness.get("per_chain", {}).items():
        lines.append(
            f"| {_md_cell(rid)} | {summary.get('chain_type')} | "
            f"{_fmt(summary.get('oasis_identity'))} | "
            f"{_fmt(summary.get('oasis_percentile'))} | "
            f"{_fmt(summary.get('germline_content'))} | "
            f"{summary.get('num_nonhuman_peptides')}/{summary.get('num_peptides')} | "
            f"{(summary.get('v_germline_names') or ['N/A'])[0]} |"
        )
    lines += ["", "## Paired antibody", "",
              "| Pair | OASis identity | Percentile | Germline content |",
              "|------|----------------|------------|------------------|"]
    for key, summary in humanness.get("per_antibody", {}).items():
        lines.append(
            f"| {_md_cell(key)} | {_fmt(summary.get('oasis_identity'))} | "
            f"{_fmt(summary.get('oasis_percentile'))} | "
            f"{_fmt(summary.get('germline_content'))} |"
        )
    lines += ["", "## Output files", ""]
    for kind, target in sorted(outputs.items()):
        lines.append(f"- {kind}: `{target}`")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")
    outputs["markdown"] = str(md_path)

    for target in sorted(set(outputs.values())):
        print(f"Wrote {target}", file=sys.stderr, flush=True)
    return 0


def cmd_humanize(args: argparse.Namespace) -> int:
    input_path, records = _load_inputs(args)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem

    modes = _resolve_modes(args.modes)
    print("[1/3] Numbering sequences...", file=sys.stderr, flush=True)
    numbered = _number_records(records, args.scheme, args.cdr_definition)

    oasis_params = None
    if args.run_oasis:
        try:
            oasis_params = _oasis_params_from_args(args)
        except ToolError as exc:
            print(f"WARNING: {exc} (skipping OASis before/after)", file=sys.stderr)
            oasis_params = None

    print(f"[2/3] Running humanization modes: {', '.join(modes)}...",
          file=sys.stderr, flush=True)
    humanization = run_humanization_modes(
        numbered,
        modes=modes,
        scheme=args.scheme,
        cdr_definition=args.cdr_definition,
        iterations=args.iterations,
        humanize_cdrs=args.humanize_cdrs,
        backmutate_vernier=not args.no_backmutate_vernier,
        model_version=args.model_version,
        heavy_v_germline=args.heavy_v_germline,
        light_v_germline=args.light_v_germline,
        sapiens_iterations=args.sapiens_iterations,
        oasis_params=oasis_params,
        verbose=args.verbose,
    )

    print("[3/3] Writing outputs...", file=sys.stderr, flush=True)
    outputs: Dict[str, str] = {}
    payload = {
        "tool": "tools/biophi/biophi_comprehensive.py",
        "tool_version": TOOL_VERSION,
        "mode": "humanize",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
        },
        "parameters": {
            "scheme": args.scheme,
            "cdr_definition": args.cdr_definition,
            "modes": modes,
            "iterations": args.iterations,
            "humanize_cdrs": args.humanize_cdrs,
            "backmutate_vernier": not args.no_backmutate_vernier,
            "model_version": args.model_version,
            "heavy_v_germline": args.heavy_v_germline,
            "light_v_germline": args.light_v_germline,
            "sapiens_iterations": args.sapiens_iterations,
        },
        "environment": _environment_report(args),
        "humanization": _jsonable(humanization),
    }
    json_path = outdir / f"{stem}_biophi_humanize.json"
    _write_json(json_path, payload)
    outputs["humanize json"] = str(json_path)

    for mode in modes:
        records_out = humanize_fasta_from_results(humanization, mode)
        if not records_out:
            continue
        fa_path = outdir / f"{stem}_biophi_humanized_{mode}.fa"
        _write_fasta(fa_path, records_out)
        outputs[f"humanized fasta ({mode})"] = str(fa_path)

    # Mutation summary CSV
    csv_rows = []
    for rid, mode_dict in humanization.get("per_record", {}).items():
        for mode, data in mode_dict.items():
            if not isinstance(data, dict) or "num_mutations" not in data:
                continue
            csv_rows.append(
                {
                    "id": rid,
                    "mode": mode,
                    "num_mutations": data["num_mutations"],
                    "changed": data.get("changed"),
                    "oasis_before": (data.get("oasis_before") or {}).get(
                        "oasis_identity"
                    ),
                    "oasis_after": (data.get("oasis_after") or {}).get(
                        "oasis_identity"
                    ),
                    "v_germline": (data.get("germline") or {}).get("v_top"),
                }
            )
    for item in humanization.get("full_antibody", []):
        for side in ("vh", "vl"):
            data = item[side]
            csv_rows.append(
                {
                    "id": data["id"],
                    "mode": f"full_antibody:{side}",
                    "num_mutations": data["num_mutations"],
                    "changed": data.get("changed"),
                    "oasis_before": (data.get("oasis_before") or {}).get(
                        "oasis_identity"
                    ),
                    "oasis_after": (data.get("oasis_after") or {}).get(
                        "oasis_identity"
                    ),
                    "v_germline": None,
                }
            )
    csv_path = outdir / f"{stem}_biophi_humanize.csv"
    _write_csv(
        csv_path,
        csv_rows,
        ["id", "mode", "num_mutations", "changed", "oasis_before",
         "oasis_after", "v_germline"],
    )
    outputs["mutation summary csv"] = str(csv_path)

    md_path = outdir / f"{stem}_biophi_humanize.md"
    lines = [
        f"# BioPhi humanization: {input_path.name}",
        "",
        f"- Generated: {payload['timestamp']}",
        f"- Modes: {', '.join(modes)}",
        f"- Numbering: {args.scheme} / CDR: {args.cdr_definition}",
        "",
        "| ID | Mode | Mutations | OASis before | OASis after | Germline |",
        "|----|------|-----------|--------------|-------------|----------|",
    ]
    for row in csv_rows:
        lines.append(
            f"| {_md_cell(row['id'])} | {row['mode']} | {row['num_mutations']} | "
            f"{_fmt(row['oasis_before'])} | {_fmt(row['oasis_after'])} | "
            f"{row['v_germline'] or 'N/A'} |"
        )
    lines += ["", "## Output files", ""]
    for kind, target in sorted(outputs.items()):
        lines.append(f"- {kind}: `{target}`")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")
    outputs["markdown"] = str(md_path)

    for target in sorted(set(outputs.values())):
        print(f"Wrote {target}", file=sys.stderr, flush=True)
    return 0


def cmd_cdr_grafting(args: argparse.Namespace) -> int:
    args.modes = "cdr_grafting"
    return cmd_humanize(args)


def cmd_full_antibody(args: argparse.Namespace) -> int:
    args.modes = "full_antibody"
    return cmd_humanize(args)


def _resolve_modes(modes_arg: str) -> List[str]:
    all_modes = ["sapiens", "cdr_grafting", "manual", "full_antibody"]
    if modes_arg in ("all", ""):
        return all_modes
    modes = [m.strip() for m in modes_arg.split(",") if m.strip()]
    unknown = [m for m in modes if m not in all_modes]
    if unknown:
        raise ToolError(
            f"Unknown humanization mode(s): {', '.join(unknown)}. "
            f"Choose from: {', '.join(all_modes)} or 'all'"
        )
    return modes


def cmd_evaluate_all(args: argparse.Namespace) -> int:
    """Everything at once: Sapiens + OASis + all humanization modes + reports."""
    input_path, records = _load_inputs(args)
    metadata = _load_json(getattr(args, "metadata", None))
    pipeline_result = _load_json(getattr(args, "pipeline_result", None))
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = input_path.stem
    modes = _resolve_modes(args.modes)
    outputs: Dict[str, str] = {}

    print(f"[1/7] Numbering {len(records)} sequences...", file=sys.stderr, flush=True)
    numbered = _number_records(records, args.scheme, args.cdr_definition)
    n_err = sum(1 for e in numbered if e.get("numbering_error"))
    if n_err:
        print(f"  WARNING: numbering failed for {n_err} record(s)",
              file=sys.stderr, flush=True)

    print("[2/7] Running full Sapiens predictions...", file=sys.stderr, flush=True)
    sapiens_results = run_sapiens_full(
        records,
        numbered,
        with_embeddings=not args.no_embeddings,
        metadata=metadata,
        verbose=args.verbose,
    )

    print("[3/7] Running full OASis humanness analysis...", file=sys.stderr, flush=True)
    humanness: Dict[str, Any] = {
        "available": False,
        "reason": "disabled via --no-oasis",
        "per_chain": {},
        "per_antibody": {},
        "peptide_frames": [],
    }
    if args.run_oasis:
        try:
            humanness = run_humanness_full(
                numbered,
                oasis_db=args.oasis_db,
                min_fraction_subjects=args.min_fraction_subjects,
                verbose=args.verbose,
            )
        except ToolError as exc:
            print(f"WARNING: {exc} (continuing without OASis)",
                  file=sys.stderr, flush=True)
            humanness = {
                "available": False,
                "reason": str(exc),
                "per_chain": {},
                "per_antibody": {},
                "peptide_frames": [],
            }

    humanization: Optional[Dict[str, Any]] = None
    if not args.no_humanization:
        print(f"[4/7] Running humanization modes: {', '.join(modes)}...",
              file=sys.stderr, flush=True)
        oasis_params = None
        if humanness.get("available"):
            try:
                oasis_params = _oasis_params_from_args(args)
            except ToolError as exc:
                print(f"WARNING: {exc}", file=sys.stderr)
        humanization = run_humanization_modes(
            numbered,
            modes=modes,
            scheme=args.scheme,
            cdr_definition=args.cdr_definition,
            iterations=args.iterations,
            humanize_cdrs=args.humanize_cdrs,
            backmutate_vernier=not args.no_backmutate_vernier,
            model_version=args.model_version,
            heavy_v_germline=args.heavy_v_germline,
            light_v_germline=args.light_v_germline,
            sapiens_iterations=args.sapiens_iterations,
            oasis_params=oasis_params,
            verbose=args.verbose,
        )
    else:
        print("[4/7] Skipping humanization (--no-humanization)",
              file=sys.stderr, flush=True)

    print("[5/7] Building cross-validation...", file=sys.stderr, flush=True)
    cross_validation = build_cross_validation(
        sapiens_results, humanness, pipeline_result, numbered
    )

    print("[6/7] Writing data files...", file=sys.stderr, flush=True)
    mean_path = outdir / f"{stem}_biophi_sapiens.csv"
    _write_mean_csv(mean_path, sapiens_results)
    outputs["sapiens mean csv"] = str(mean_path)

    scores_path = outdir / f"{stem}_biophi_scores.csv"
    _write_scores_csv(scores_path, sapiens_results)
    if scores_path.exists():
        outputs["full score matrix csv"] = str(scores_path)

    heatmap_path = outdir / f"{stem}_biophi_heatmap.json"
    _write_json(heatmap_path, build_heatmap(sapiens_results))
    outputs["heatmap json"] = str(heatmap_path)

    if not args.no_embeddings:
        emb_path = outdir / f"{stem}_biophi_embeddings.json"
        _write_json(emb_path, build_embeddings_payload(sapiens_results))
        outputs["embeddings json"] = str(emb_path)

    if humanness.get("available"):
        peptides_path = outdir / f"{stem}_biophi_peptides.csv"
        _write_peptides_csv(peptides_path, humanness)
        if peptides_path.exists():
            outputs["peptide-level csv"] = str(peptides_path)
        xlsx_path = outdir / f"{stem}_biophi_oasis.xlsx"
        written = build_oasis_workbook(
            numbered, humanness, xlsx_path, args.min_fraction_subjects
        )
        if written:
            outputs["oasis xlsx workbook"] = str(written)

    if humanization:
        for mode in modes:
            recs_out = humanize_fasta_from_results(humanization, mode)
            if recs_out:
                fa_path = outdir / f"{stem}_biophi_humanized_{mode}.fa"
                _write_fasta(fa_path, recs_out)
                outputs[f"humanized fasta ({mode})"] = str(fa_path)

    print("[7/7] Writing unified report...", file=sys.stderr, flush=True)
    environment = _environment_report(args)
    parameters = {
        "scheme": args.scheme,
        "cdr_definition": args.cdr_definition,
        "modes": modes,
        "run_oasis": args.run_oasis,
        "oasis_db": args.oasis_db,
        "min_fraction_subjects": args.min_fraction_subjects,
        "iterations": args.iterations,
        "humanize_cdrs": args.humanize_cdrs,
        "backmutate_vernier": not args.no_backmutate_vernier,
        "model_version": args.model_version,
        "heavy_v_germline": args.heavy_v_germline,
        "light_v_germline": args.light_v_germline,
        "sapiens_iterations": args.sapiens_iterations,
        "embeddings": not args.no_embeddings,
        "metadata": str(Path(args.metadata).resolve()) if args.metadata else None,
        "pipeline_result": (
            str(Path(args.pipeline_result).resolve()) if args.pipeline_result else None
        ),
        "pipeline_result_obj": pipeline_result,
    }

    unified: Dict[str, Any] = {
        "tool": "tools/biophi/biophi_comprehensive.py",
        "tool_version": TOOL_VERSION,
        "mode": "evaluate-all",
        "timestamp": _now(),
        "input": {
            "path": str(input_path.resolve()),
            "sha256": _sha256_file(input_path),
            "n_records": len(records),
            "metadata": parameters["metadata"],
            "pipeline_result": parameters["pipeline_result"],
        },
        "parameters": {k: v for k, v in parameters.items() if k != "pipeline_result_obj"},
        "environment": environment,
        "variants": [
            {
                "id": e["id"],
                "base": e.get("base"),
                "variant": e.get("variant"),
                "description": e.get("description"),
                "chain_type": e.get("chain_type"),
                "length": e.get("length"),
                "sequence": e.get("sequence"),
                "numbering_error": e.get("numbering_error"),
                "sapiens": _jsonable(sapiens_results.get(e["id"])),
                "oasis": _jsonable((humanness.get("per_chain") or {}).get(e["id"])),
                "humanization": _jsonable(
                    (humanization or {}).get("per_record", {}).get(e["id"])
                ),
            }
            for e in numbered
        ],
        "antibody_humanness": _jsonable(humanness.get("per_antibody")),
        "full_antibody_humanization": _jsonable(
            (humanization or {}).get("full_antibody")
        ),
        "cross_validation": _jsonable(cross_validation),
        "outputs": outputs,
        "summary": {
            "n_variants": len(numbered),
            "numbering_errors": n_err,
            "oasis_available": bool(humanness.get("available")),
            "humanization_modes": modes if humanization else [],
            "sapiens_mean_range": _sapiens_range(sapiens_results),
        },
    }
    json_path = outdir / f"{stem}_biophi_comprehensive.json"
    _write_json(json_path, unified)
    outputs["comprehensive json"] = str(json_path)

    csv_rows = _variant_summary_rows(
        numbered, sapiens_results, humanness, pipeline_result, humanization
    )
    csv_path = outdir / f"{stem}_biophi_comprehensive.csv"
    _write_csv(
        csv_path,
        csv_rows,
        [
            "id", "base", "variant", "chain_type", "length",
            "sapiens_mean", "sapiens_min", "oasis_identity", "oasis_percentile",
            "germline_content", "v_germline", "pipeline_n_backmutations",
            "pipeline_cvi_homology", "sapiens_humanize_mutations",
            "cdr_grafting_mutations",
        ],
    )
    outputs["summary csv"] = str(csv_path)

    md_path = outdir / f"{stem}_biophi_comprehensive.md"
    write_markdown_report(
        md_path,
        stem,
        input_path,
        numbered,
        sapiens_results,
        humanness,
        humanization,
        cross_validation,
        parameters,
        environment,
        outputs,
    )
    outputs["markdown report"] = str(md_path)

    # Refresh JSON with the complete output index.
    unified["outputs"] = outputs
    _write_json(json_path, unified)

    for target in sorted(set(outputs.values())):
        print(f"Wrote {target}", file=sys.stderr, flush=True)
    return 0


def _sapiens_range(sapiens_results: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, float]]:
    values = [
        v["sapiens_mean"]
        for v in sapiens_results.values()
        if v.get("sapiens_mean") is not None
    ]
    if not values:
        return None
    return {"min": _round(min(values)), "max": _round(max(values))}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Full BioPhi capability reproduction (Sapiens + OASis + all "
            "humanization modes). Runs inside the biophi conda env; all "
            "analyses default to the fullest output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tool-version", action="version", version=TOOL_VERSION)
    sub = parser.add_subparsers(dest="mode", required=True)

    def add_common(p: argparse.ArgumentParser) -> None:
        p.add_argument("--env", default=DEFAULT_ENV,
                       help=f"conda env containing biophi (default: {DEFAULT_ENV})")
        p.add_argument("--hf-endpoint", default=None,
                       help="HuggingFace endpoint/mirror; omit for offline cache mode")
        p.add_argument("--verbose", action="store_true",
                       help="Stream progress for each record")
        p.add_argument("--md5", action="store_true",
                       help="Compute MD5 of the OASis DB (slow for ~22 GB)")

    def add_numbering(p: argparse.ArgumentParser) -> None:
        p.add_argument("--scheme", default="kabat",
                       choices=["kabat", "chothia", "imgt", "aho"])
        p.add_argument("--cdr-definition", default="kabat",
                       choices=["kabat", "chothia", "imgt", "north"])

    def add_sapiens_opts(p: argparse.ArgumentParser) -> None:
        p.add_argument("--no-embeddings", dest="no_embeddings",
                       action="store_true",
                       help="Skip residue/sequence embeddings (default: on)")
        p.add_argument("--iterations", type=int, default=1,
                       help="Sapiens humanization iterations (default: 1)")
        p.add_argument("--humanize-cdrs", action="store_true",
                       help="Allow humanizing mutations inside CDRs")
        p.add_argument("--no-backmutate-vernier",
                       dest="no_backmutate_vernier", action="store_true",
                       help="Disable Vernier back-mutations")
        p.add_argument("--model-version", default="latest",
                       help="Sapiens model version (default: latest)")

    def add_oasis_opts(p: argparse.ArgumentParser) -> None:
        p.add_argument("--oasis-db", default=DEFAULT_DB,
                       help=f"OASis 9-mer DB (default: {DEFAULT_DB})")
        p.add_argument("--min-fraction-subjects", type=float, default=0.1,
                       help="Minimum fraction of OAS subjects (default: 0.1 = 10%%)")
        p.add_argument("--no-oasis", dest="run_oasis", action="store_false",
                       default=True, help="Skip OASis analysis")

    def add_humanize_opts(p: argparse.ArgumentParser) -> None:
        p.add_argument("--modes", default="all",
                       help="Comma list of sapiens,cdr_grafting,manual,"
                            "full_antibody or 'all' (default: all)")
        p.add_argument("--heavy-v-germline", default="auto",
                       help="CDR-grafting heavy V germline (default: auto)")
        p.add_argument("--light-v-germline", default="auto",
                       help="CDR-grafting light V germline (default: auto)")
        p.add_argument("--sapiens-iterations", type=int, default=0,
                       help="Extra Sapiens pass after CDR grafting (default: 0)")

    check = sub.add_parser("check", help="Check env, cache, HMMER and OASis DB")
    check.add_argument("--oasis-db", default=DEFAULT_DB)
    add_common(check)
    check.set_defaults(func=cmd_check)

    sapiens = sub.add_parser(
        "sapiens",
        help="Full Sapiens analysis (mean/matrix/best/masked/embeddings/humanize)",
    )
    sapiens.add_argument("--input", required=True)
    sapiens.add_argument("--outdir", default="outputs/biophi_full")
    sapiens.add_argument("--metadata", default=None,
                         help="variants_numbering.json for Kabat-mapped outputs")
    sapiens.add_argument("--no-humanize", dest="no_humanize", action="store_true",
                         help="Skip Sapiens-humanized FASTA (default: on)")
    sapiens.add_argument("--limit", type=int, default=None,
                         help="Process only the first N records")
    add_numbering(sapiens)
    add_sapiens_opts(sapiens)
    add_common(sapiens)
    sapiens.set_defaults(func=cmd_sapiens)

    humanness = sub.add_parser(
        "humanness",
        help="Full OASis humanness (peptides, curves, frequencies, XLSX)",
    )
    humanness.add_argument("--input", required=True)
    humanness.add_argument("--outdir", default="outputs/biophi_full")
    humanness.add_argument("--limit", type=int, default=None)
    add_numbering(humanness)
    add_oasis_opts(humanness)
    add_common(humanness)
    humanness.set_defaults(func=cmd_humanness)

    humanize = sub.add_parser(
        "humanize",
        help="All humanization modes (sapiens, cdr_grafting, manual, full_antibody)",
    )
    humanize.add_argument("--input", required=True)
    humanize.add_argument("--outdir", default="outputs/biophi_full")
    humanize.add_argument("--limit", type=int, default=None)
    add_numbering(humanize)
    add_sapiens_opts(humanize)
    add_oasis_opts(humanize)
    add_humanize_opts(humanize)
    add_common(humanize)
    humanize.set_defaults(func=cmd_humanize)

    cdr = sub.add_parser("cdr-grafting", help="CDR grafting with full parameters")
    cdr.add_argument("--input", required=True)
    cdr.add_argument("--outdir", default="outputs/biophi_full")
    cdr.add_argument("--limit", type=int, default=None)
    add_numbering(cdr)
    add_sapiens_opts(cdr)
    add_oasis_opts(cdr)
    add_humanize_opts(cdr)
    add_common(cdr)
    cdr.set_defaults(func=cmd_cdr_grafting)

    antibody = sub.add_parser(
        "full-antibody", help="Coordinated VH+VL (full antibody) humanization"
    )
    antibody.add_argument("--input", required=True)
    antibody.add_argument("--outdir", default="outputs/biophi_full")
    antibody.add_argument("--limit", type=int, default=None)
    add_numbering(antibody)
    add_sapiens_opts(antibody)
    add_oasis_opts(antibody)
    add_humanize_opts(antibody)
    add_common(antibody)
    antibody.set_defaults(func=cmd_full_antibody)

    evaluate = sub.add_parser(
        "evaluate-all",
        help="Complete analysis: Sapiens + OASis + all humanization + reports",
    )
    evaluate.add_argument("--input", required=True,
                          help="FASTA (variants.fasta or any antibody FASTA)")
    evaluate.add_argument("--metadata", default=None,
                          help="variants_numbering.json (Kabat mapping)")
    evaluate.add_argument("--pipeline-result", default=None,
                          help="humanization_result.json for cross-validation")
    evaluate.add_argument("--outdir", default="outputs/biophi_full")
    evaluate.add_argument("--no-humanization", dest="no_humanization",
                          action="store_true",
                          help="Skip humanization modes (default: run all)")
    evaluate.add_argument("--limit", type=int, default=None)
    add_numbering(evaluate)
    add_sapiens_opts(evaluate)
    add_oasis_opts(evaluate)
    add_humanize_opts(evaluate)
    add_common(evaluate)
    evaluate.set_defaults(func=cmd_evaluate_all)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    mfs = getattr(args, "min_fraction_subjects", None)
    if mfs is not None and not (0 < mfs <= 0.9):
        parser.error("--min-fraction-subjects must be in (0, 0.9]")
    try:
        return int(args.func(args))
    except ToolError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
