#!/usr/bin/env python3
"""Immunogenicity prediction backends.

Pluggable MHC-II T-cell epitope backends.
Each backend: name, is_available(), score().

Backend priority (auto): HLAIIPred > heuristic.

All backends produce a BackendResult with:
  - per_peptide: List[{peptide, allele, rank_or_score, core}]
  - per_residue: Dict[int, float]  (raw sequence index → epitope score 0-1)
  - meta: Dict  (tool version, allele count, etc.)
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class PeptideHit:
    peptide: str
    allele: str
    rank: float           # compatibility field: 1 - presentation probability
    score: float = 0.0    # presentation probability (0-1)
    core: str = ""
    core_pos: int = 0     # 0-based position of core start in peptide
    start_index: int = 0  # 0-based position of peptide start in the chain


@dataclass
class BackendResult:
    backend: str
    per_peptide: List[PeptideHit] = field(default_factory=list)
    per_residue: Dict[int, float] = field(default_factory=dict)
    per_residue_n_epitopes: Dict[int, int] = field(default_factory=dict)
    meta: Dict = field(default_factory=dict)
    error: Optional[str] = None


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Default 9 common DRB1 alleles (HLAIIPred reference panel)
DEFAULT_ALLELES = [
    "DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*07:01",
    "DRB1*08:01", "DRB1*09:01", "DRB1*11:01", "DRB1*13:01",
    "DRB1*15:01",
]

# Broader panel: DRB1 + DRB3/4/5 + DP/DQ heterodimers. Heterodimers use the
# "alpha-beta" spelling accepted by HLAIIPred (see its predict_peptide_example).
EXTENDED_ALLELES = [
    "DRB1*01:01", "DRB1*03:01", "DRB1*04:01", "DRB1*04:05",
    "DRB1*07:01", "DRB1*08:01", "DRB1*09:01", "DRB1*11:01",
    "DRB1*12:01", "DRB1*13:01", "DRB1*15:01", "DRB1*16:01",
    "DRB3*01:01", "DRB3*02:02", "DRB4*01:01", "DRB5*01:01",
    "DPA1*01:03-DPB1*02:01", "DPA1*02:01-DPB1*01:01",
    "DPA1*01:03-DPB1*04:01", "DPA1*02:01-DPB1*05:01",
    "DQA1*01:02-DQB1*06:02", "DQA1*01:01-DQB1*05:01",
    "DQA1*03:01-DQB1*03:02", "DQA1*05:01-DQB1*02:01",
]

ALLELE_PANELS = {
    "9common": DEFAULT_ALLELES,
    "extended": EXTENDED_ALLELES,
}

# HLAIIPred accepts at most 14 alleles per call.
MAX_ALLELES_PER_CALL = 14

# HLAIIPred returns a presentation probability (sigmoid). The 0.5 cutoff is the
# model's natural binary decision boundary; 0.3 is a softer "weak" cutoff.
DEFAULT_STRONG_THRESHOLD = 0.5
DEFAULT_WEAK_THRESHOLD = 0.3

# Peptide lengths
DEFAULT_PEP_LEN = 15
DEFAULT_PEP_STEP = 1
MIN_PEP_LEN = 9
MAX_PEP_LEN = 30


# ---------------------------------------------------------------------------
# Peptide generation
# ---------------------------------------------------------------------------

def generate_peptides(seq: str, length: int = DEFAULT_PEP_LEN,
                      step: int = DEFAULT_PEP_STEP) -> List[Tuple[str, int]]:
    """Generate overlapping peptides: list of (peptide, start_index_0based)."""
    out = []
    for i in range(0, len(seq) - length + 1, step):
        out.append((seq[i:i + length], i))
    return out


# ---------------------------------------------------------------------------
# Backend base
# ---------------------------------------------------------------------------

class Backend(ABC):
    name: str = "base"
    requires_structure: bool = False
    is_neural: bool = False
    license_type: str = "unknown"

    @abstractmethod
    def is_available(self) -> bool:
        """Return True if this backend can run (env, binaries, weights)."""
        ...

    @abstractmethod
    def score(self, sequences: Dict[str, str], alleles: List[str],
              pep_len: int = DEFAULT_PEP_LEN, **kwargs) -> Dict[str, BackendResult]:
        """Score peptides from sequences.

        Args:
            sequences: {chain_id: sequence}
            alleles: list of HLA-II alleles
            pep_len: peptide length for windowing
        Returns:
            {chain_id: BackendResult}
        """
        ...

    def describe(self) -> str:
        return f"{self.name} ({self.license_type})"


# ---------------------------------------------------------------------------
# Heuristic backend (pure Python, zero deps)
# ---------------------------------------------------------------------------

class HeuristicBackend(Backend):
    """Allele-agnostic MHC-II anchor-propensity heuristic.

    This is an explicitly non-ML, allele-agnostic fallback. It scores each
    9-mer register inside a peptide by the hydrophobicity/aromaticity of the
    P1/P4/P6/P9 anchors and a small non-self (human proteome rarity) term.
    It is useful for relative ranking only; it is NOT a substitute for a
    trained MHC-II presentation model.
    """
    name = "heuristic"
    requires_structure = False
    is_neural = False
    license_type = "MIT"

    _ANCHOR_PREF = {
        "A": 0.55, "C": 0.45, "D": 0.15, "E": 0.15, "F": 0.85,
        "G": 0.35, "H": 0.30, "I": 0.80, "K": 0.20, "L": 0.85,
        "M": 0.75, "N": 0.25, "P": 0.30, "Q": 0.25, "R": 0.25,
        "S": 0.30, "T": 0.35, "V": 0.80, "W": 0.85, "Y": 0.85,
    }
    _AA_FREQ = {
        "A": 0.0825, "C": 0.0137, "D": 0.0545, "E": 0.0675,
        "F": 0.0386, "G": 0.0707, "H": 0.0227, "I": 0.0596,
        "K": 0.0584, "L": 0.0966, "M": 0.0242, "N": 0.0406,
        "P": 0.0470, "Q": 0.0393, "R": 0.0553, "S": 0.0656,
        "T": 0.0534, "V": 0.0687, "W": 0.0110, "Y": 0.0292,
    }
    _ANCHOR_WEIGHTS = (0.35, 0.25, 0.20, 0.20)

    def is_available(self) -> bool:
        return True

    def _core_score(self, core: str) -> float:
        anchor = sum(
            w * self._ANCHOR_PREF.get(aa, 0.35)
            for w, aa in zip(self._ANCHOR_WEIGHTS, (core[0], core[3], core[5], core[8]))
        )
        max_freq = max(self._AA_FREQ.values())
        rarity = sum(1.0 - self._AA_FREQ.get(aa, 0.05) / max_freq for aa in core) / len(core)
        return min(1.0, 0.75 * anchor + 0.25 * rarity)

    def score(self, sequences: Dict[str, str], alleles: List[str],
              pep_len: int = DEFAULT_PEP_LEN, **kwargs) -> Dict[str, BackendResult]:
        chain_types = kwargs.get("chain_types") or {}
        results = {}
        for chain_id, seq in sequences.items():
            chain_type = chain_types.get(chain_id)
            if chain_type not in ("H", "L"):
                chain_type = "H" if chain_id.startswith(("H", "VH")) else "L"
            peptides = generate_peptides(seq, pep_len)
            per_residue: Dict[int, float] = {}
            per_peptide: List[PeptideHit] = []
            for pep, start in peptides:
                best_score = 0.0
                best_core_pos = 0
                best_core = pep[:9]
                for core_pos in range(0, len(pep) - 8):
                    core = pep[core_pos:core_pos + 9]
                    sc = self._core_score(core)
                    if sc > best_score:
                        best_score, best_core_pos, best_core = sc, core_pos, core
                core_start = start + best_core_pos
                for j in range(9):
                    gidx = core_start + j
                    if gidx not in per_residue or best_score > per_residue[gidx]:
                        per_residue[gidx] = round(best_score, 4)
                per_peptide.append(PeptideHit(
                    peptide=pep,
                    allele="pan-allele",
                    rank=round(1.0 - best_score, 4),
                    score=round(best_score, 4),
                    core=best_core,
                    core_pos=best_core_pos,
                    start_index=start,
                ))

            results[chain_id] = BackendResult(
                backend=self.name,
                per_peptide=per_peptide,
                per_residue=per_residue,
                meta={"pep_len": pep_len, "peptides": len(peptides),
                      "note": "allele-agnostic anchor heuristic (not ML)"},
            )
        return results


# ---------------------------------------------------------------------------
# HLAIIPred backend
# ---------------------------------------------------------------------------

class HLAIIPredBackend(Backend):
    """HLAIIPred (Pfizer, Apache-2.0): Transformer pan-allele MHC-II predictor.

    Requires: conda env 'hlapred' with hlaiipred package + models dir.
    """
    name = "hlaiipred"
    requires_structure = False
    is_neural = True
    license_type = "apache-2.0"

    def __init__(self, env: str = "hlapred", models_dir: str = "models",
                 runner_path: Optional[str] = None):
        self.env = env
        self.models_dir = self._resolve_models_dir(models_dir)
        self.runner_path = runner_path or os.path.join(
            os.path.dirname(__file__), "runners", "hlaiipred_runner.py")

    @staticmethod
    def _resolve_models_dir(models_dir: str) -> str:
        candidates = [
            models_dir,
            os.path.join(os.path.dirname(__file__), models_dir),
            os.path.join(os.path.dirname(__file__), "models"),
            os.path.expanduser("~/.hlaiipred/models"),
        ]
        for candidate in candidates:
            if candidate and os.path.isdir(candidate):
                return os.path.abspath(candidate)
        return models_dir

    def _models_ready(self) -> bool:
        required = ("epT_0.pt", "epT_1.pt", "config.yaml")
        return (
            os.path.isdir(self.models_dir)
            and all(os.path.isfile(os.path.join(self.models_dir, f)) for f in required)
            and os.path.isdir(os.path.join(self.models_dir, "mhcII"))
            and os.path.isfile(self.runner_path)
        )

    def is_available(self) -> bool:
        try:
            if not self._models_ready():
                return False
            proc = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True, text=True, timeout=15
            )
            if proc.returncode != 0:
                return False
            envs = json.loads(proc.stdout).get("envs", [])
            env_names = [Path(p).name for p in envs]
            if self.env not in env_names:
                return False
            # Verify hlapred is importable in that env
            py = self._find_env_python()
            if not py:
                return False
            proc2 = subprocess.run(
                [py, "-c", "import hlapred"],
                capture_output=True, timeout=10
            )
            return proc2.returncode == 0
        except Exception:
            return False

    def _find_env_python(self) -> str:
        """Find the python binary inside the hlapred conda env."""
        try:
            proc = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True, text=True, timeout=15
            )
            if proc.returncode == 0:
                envs = json.loads(proc.stdout).get("envs", [])
                for env_path in envs:
                    if env_path.endswith(self.env) or env_path.endswith(f"/{self.env}"):
                        py = os.path.join(env_path, "bin", "python")
                        if os.path.isfile(py):
                            return py
        except Exception:
            pass
        # Fallback: use conda run to find python
        return sys.executable

    def score(self, sequences: Dict[str, str], alleles: List[str],
              pep_len: int = 15, **kwargs) -> Dict[str, BackendResult]:
        results = {}
        if not alleles:
            for chain_id in sequences:
                results[chain_id] = BackendResult(
                    backend=self.name, error="no HLA-II alleles supplied")
            return results
        if not (MIN_PEP_LEN <= pep_len <= MAX_PEP_LEN):
            for chain_id in sequences:
                results[chain_id] = BackendResult(
                    backend=self.name,
                    error=f"peptide length {pep_len} outside {MIN_PEP_LEN}-{MAX_PEP_LEN}")
            return results
        input_data = {
            "sequences": sequences,
            "alleles": alleles,
            "pep_len": pep_len,
            "models_dir": self.models_dir,
        }
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(input_data, f)
            input_path = f.name
        try:
            hlapred_python = self._find_env_python()
            timeout = max(600, 120 * max(1, len(sequences)))
            proc = subprocess.run(
                [hlapred_python, self.runner_path, input_path],
                capture_output=True, text=True, timeout=timeout
            )
            if proc.returncode != 0:
                err = proc.stderr[-800:] if proc.stderr else "unknown error"
                for chain_id in sequences:
                    results[chain_id] = BackendResult(
                        backend=self.name, error=f"HLAIIPred failed: {err}")
                return results
            output = json.loads(proc.stdout)
            for chain_id in sequences:
                res_dict = output.get(chain_id)
                if not res_dict:
                    results[chain_id] = BackendResult(
                        backend=self.name, error="HLAIIPred returned no result for chain")
                    continue
                meta = res_dict.get("meta", {}) or {}
                err = res_dict.get("error") or meta.get("error")
                if err:
                    results[chain_id] = BackendResult(
                        backend=self.name, error=str(err))
                    continue
                per_peptide = [
                    PeptideHit(
                        peptide=h.get("peptide", ""),
                        allele=h.get("allele", ""),
                        rank=float(h.get("rank", 0.0)),
                        score=float(h.get("score", 0.0)),
                        core=h.get("core", ""),
                        core_pos=int(h.get("core_pos", 0)),
                        start_index=int(h.get("start_index", 0)),
                    )
                    for h in res_dict.get("per_peptide", [])
                ]
                raw_per_residue = res_dict.get("per_residue", {})
                per_residue = {int(k): v for k, v in raw_per_residue.items()}
                results[chain_id] = BackendResult(
                    backend=self.name,
                    per_peptide=per_peptide,
                    per_residue=per_residue,
                    meta=meta,
                )
        except subprocess.TimeoutExpired:
            for chain_id in sequences:
                results[chain_id] = BackendResult(
                    backend=self.name, error="HLAIIPred timed out")
        except Exception as e:
            for chain_id in sequences:
                results[chain_id] = BackendResult(
                    backend=self.name, error=f"HLAIIPred error: {e}")
        finally:
            try:
                os.unlink(input_path)
            except OSError:
                pass
        return results


# ---------------------------------------------------------------------------
# Backend registry
# ---------------------------------------------------------------------------

ALL_BACKENDS = {
    "hlaiipred": HLAIIPredBackend,
    "heuristic": HeuristicBackend,
}

# Auto priority for MHC-II epitope prediction
AUTO_PRIORITY = ["hlaiipred", "heuristic"]


def get_backend(name: str, **kwargs) -> Backend:
    """Get a backend instance by name."""
    cls = ALL_BACKENDS.get(name)
    if cls is None:
        raise ValueError(f"Unknown backend: {name}. Available: {list(ALL_BACKENDS.keys())}")
    return cls(**kwargs)


def auto_select_backend(preferred: Optional[str] = None,
                        prefer_online: bool = False) -> Backend:
    """Auto-select the best available MHC-II backend.

    preferred: force a specific backend
    prefer_online: if True, skip heuristic even if others unavailable
    """
    if preferred:
        b = get_backend(preferred)
        if b.is_available():
            return b
        raise RuntimeError(f"Preferred backend '{preferred}' is not available")

    for name in AUTO_PRIORITY:
        if name == "heuristic" and prefer_online:
            continue
        try:
            b = get_backend(name)
            if b.is_available():
                return b
        except Exception:
            continue
    return HeuristicBackend()


def check_backends() -> Dict[str, bool]:
    """Check availability of all backends."""
    out = {}
    for name, cls in ALL_BACKENDS.items():
        try:
            b = cls()
            out[name] = b.is_available()
        except Exception:
            out[name] = False
    return out
