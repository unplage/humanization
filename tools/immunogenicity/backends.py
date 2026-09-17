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
import re
import shutil
import subprocess
import sys
import tempfile
import zlib
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
    rank: float           # sigmoid score (HLAIIPred) or %Rank
    score: float = 0.0    # raw score
    core: str = ""
    core_pos: int = 0     # 0-based position of core start in peptide


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

# Strong binder threshold
STRONG_RANK_THRESHOLD = 2.0   # %Rank < 2
WEAK_RANK_THRESHOLD = 10.0    # %Rank < 10

# Peptide lengths
DEFAULT_PEP_LEN = 15
DEFAULT_PEP_STEP = 1


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
    """Heuristic: germline rarity × surface exposure proxy.

    Always available. No ML, no binaries. Use as fallback.
    Approximates epitope propensity as (1 - human_germline_identity).
    """
    name = "heuristic"
    requires_structure = False
    is_neural = False
    license_type = "MIT"

    def is_available(self) -> bool:
        return True

    def _germline_rarity(self, aa: str, chain_type: str, pos_idx: int) -> float:
        """Approximate germline rarity from common human germline amino acids."""
        common_human = {
            "H": {0: "Q", 1: "V", 2: "Q", 3: "L", 4: "V", 5: "E", 6: "S", 7: "G",
                  8: "G", 9: "A", 10: "E", 11: "V", 12: "K", 13: "K", 14: "P"},
            "L": {0: "D", 1: "I", 2: "V", 3: "M", 4: "T", 5: "Q", 6: "S", 7: "P"},
        }
        chain_common = common_human.get(chain_type, {})
        if aa == chain_common.get(pos_idx % len(chain_common), ""):
            return 0.1
        return 0.5 + 0.3 * (zlib.crc32(f"{chain_type}{pos_idx}{aa}".encode()) % 10) / 10

    def _surface_exposure_proxy(self, pos_idx: int, seq_len: int) -> float:
        """Simple exposure proxy based on sequence position."""
        frac = pos_idx / max(1, seq_len)
        if frac < 0.27 or 0.33 < frac < 0.58 or frac > 0.84:
            return 0.8
        return 0.3

    def score(self, sequences: Dict[str, str], alleles: List[str],
              pep_len: int = DEFAULT_PEP_LEN, **kwargs) -> Dict[str, BackendResult]:
        results = {}
        for chain_id, seq in sequences.items():
            chain_type = "H" if chain_id.startswith("H") or chain_id.startswith("VH") else "L"
            peptides = generate_peptides(seq, pep_len)
            per_residue: Dict[int, float] = {}
            for pep, start in peptides:
                pep_score = 0.0
                for j, aa in enumerate(pep):
                    global_idx = start + j
                    rarity = self._germline_rarity(aa, chain_type, global_idx)
                    exposure = self._surface_exposure_proxy(global_idx, len(seq))
                    pep_score += rarity * exposure
                pep_score /= max(1, len(pep))
                for j in range(len(pep)):
                    global_idx = start + j
                    if global_idx not in per_residue or pep_score > per_residue[global_idx]:
                        per_residue[global_idx] = min(1.0, pep_score)

            results[chain_id] = BackendResult(
                backend=self.name,
                per_residue=per_residue,
                meta={"pep_len": pep_len, "note": "heuristic proxy (not ML)"},
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
        self.models_dir = models_dir
        self.runner_path = runner_path or os.path.join(
            os.path.dirname(__file__), "runners", "hlaiipred_runner.py")

    def is_available(self) -> bool:
        try:
            proc = subprocess.run(
                ["conda", "env", "list", "--json"],
                capture_output=True, text=True, timeout=15
            )
            if proc.returncode == 0:
                envs = json.loads(proc.stdout).get("envs", [])
                env_names = [Path(p).name for p in envs]
                if self.env not in env_names:
                    return False
            else:
                return False
            if not os.path.isdir(self.models_dir):
                for candidate in [
                    self.models_dir,
                    os.path.join(os.path.dirname(__file__), self.models_dir),
                    os.path.expanduser("~/.hlaiipred/models"),
                ]:
                    if os.path.isdir(candidate):
                        self.models_dir = candidate
                        return True
                return False
            return True
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
            # Find hlapred env's python
            hlapred_python = self._find_env_python()
            proc = subprocess.run(
                [hlapred_python, self.runner_path, input_path],
                capture_output=True, text=True, timeout=600
            )
            if proc.returncode != 0:
                err = proc.stderr[-500:] if proc.stderr else "unknown error"
                for chain_id in sequences:
                    results[chain_id] = BackendResult(
                        backend=self.name, error=f"HLAIIPred failed: {err}")
                return results
            output = json.loads(proc.stdout)
            for chain_id, res_dict in output.items():
                per_peptide = [PeptideHit(**h) for h in res_dict.get("per_peptide", [])]
                # Normalize per_residue keys to int (JSON deserializes dict keys as strings)
                raw_per_residue = res_dict.get("per_residue", {})
                per_residue = {int(k): v for k, v in raw_per_residue.items()}
                results[chain_id] = BackendResult(
                    backend=self.name,
                    per_peptide=per_peptide,
                    per_residue=per_residue,
                    meta=res_dict.get("meta", {}),
                )
        except subprocess.TimeoutExpired:
            for chain_id in sequences:
                results[chain_id] = BackendResult(
                    backend=self.name, error="HLAIIPred timed out (600s)")
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
