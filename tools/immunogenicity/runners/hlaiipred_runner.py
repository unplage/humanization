#!/usr/bin/env python3
"""HLAIIPred runner script.

Runs inside the 'hlapred' conda environment. Reads JSON input, predicts
MHC-II peptide presentation, writes JSON to stdout.

The HLAIIPred API returns:
  * ``y_pred``: best presentation score over all alleles and 9-mer registers
    for each input peptide;
  * ``scores``: per-peptide, per-allele, per-9-mer scores.

This runner uses ``scores`` to recover the best allele and the core register,
instead of blindly labelling every peptide with the first allele.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import torch
from hlapred.predict import HLAIIPredict

MAX_ALLELES_PER_CALL = 14
MIN_PEP_LEN = 9
MAX_PEP_LEN = 30


def _to_numpy(value):
    if value is None:
        return None
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def _predictors(models_dir, mhc2_root, device):
    predictors = []
    for fold_idx in range(2):
        try:
            predictors.append(HLAIIPredict(models_dir, fold_idx, device, mhc2_root))
        except Exception as exc:
            print(f"Warning: Failed to load fold {fold_idx}: {exc}", file=sys.stderr)
    if not predictors:
        raise RuntimeError("No HLAIIPred models loaded")
    return predictors


def run_hlaiipred(sequences, alleles, pep_len=15, models_dir="models"):
    if not (MIN_PEP_LEN <= int(pep_len) <= MAX_PEP_LEN):
        raise ValueError(
            f"peptide length {pep_len} outside supported range "
            f"{MIN_PEP_LEN}-{MAX_PEP_LEN}")

    device = torch.device("cpu" if not torch.cuda.is_available() else "cuda:0")

    mhc2_root = None
    for candidate in [
        os.path.join(models_dir, "mhcII"),
        os.path.join(os.path.dirname(models_dir), "mhcII"),
        os.path.expanduser("~/.hlaiipred/mhcII"),
    ]:
        if os.path.isdir(candidate):
            mhc2_root = candidate
            break
    if mhc2_root is None:
        raise RuntimeError(f"mhcII pseudosequence directory not found under {models_dir}")

    predictors = _predictors(models_dir, mhc2_root, device)
    alleles = [str(a) for a in alleles if a not in (None, "", 0)]
    if not alleles:
        raise ValueError("no HLA-II alleles supplied")
    chunks = [alleles[i:i + MAX_ALLELES_PER_CALL]
              for i in range(0, len(alleles), MAX_ALLELES_PER_CALL)]

    results = {}
    for chain_id, seq in sequences.items():
        seq = str(seq).strip().upper()
        if len(seq) < pep_len:
            results[chain_id] = {
                "per_peptide": [],
                "per_residue": {},
                "error": (f"sequence length {len(seq)} is shorter than "
                          f"peptide length {pep_len}"),
                "meta": {"pep_len": pep_len, "sequence_length": len(seq)},
            }
            continue

        peptides = [seq[i:i + pep_len] for i in range(0, len(seq) - pep_len + 1)]
        starts = list(range(0, len(seq) - pep_len + 1))
        n_pep = len(peptides)
        n_registers = pep_len - 8

        best_scores = np.full(n_pep, -1.0, dtype=float)
        best_alleles = [""] * n_pep
        best_core_pos = [0] * n_pep

        for chunk in chunks:
            padded = list(chunk) + [0] * (MAX_ALLELES_PER_CALL - len(chunk))
            allele_batch = [padded] * n_pep
            fold_y = []
            fold_scores = []
            for predictor in predictors:
                inputs = predictor.prepare_input(peptides, allele_batch)
                y_pred, scores = predictor.predict(
                    inputs, batch_size=300, sigmoid=True)
                y_arr = _to_numpy(y_pred).reshape(-1)
                fold_y.append(y_arr)
                s_arr = _to_numpy(scores)
                if s_arr is None:
                    continue
                if s_arr.ndim == 3:
                    s_arr = s_arr[:, :len(chunk), :]
                else:
                    s_arr = s_arr.reshape(n_pep, MAX_ALLELES_PER_CALL, -1)
                    s_arr = s_arr[:, :len(chunk), :]
                fold_scores.append(s_arr)

            avg_y = np.mean(np.stack(fold_y, axis=0), axis=0)
            avg_scores = None
            if fold_scores:
                avg_scores = np.mean(np.stack(fold_scores, axis=0), axis=0)

            for i in range(n_pep):
                score = float(avg_y[i])
                if score <= best_scores[i]:
                    continue
                best_scores[i] = score
                best_alleles[i] = chunk[0]
                best_core_pos[i] = 0
                if avg_scores is not None and avg_scores.shape[1] > 0:
                    sub = avg_scores[i]
                    a_idx, c_idx = np.unravel_index(int(np.argmax(sub)), sub.shape)
                    best_alleles[i] = chunk[int(a_idx)] if int(a_idx) < len(chunk) else chunk[0]
                    best_core_pos[i] = int(c_idx)

        per_peptide = []
        per_residue = {}
        for i, (pep, start) in enumerate(zip(peptides, starts)):
            score = max(0.0, min(1.0, float(best_scores[i])))
            core_pos = max(0, min(best_core_pos[i], n_registers - 1))
            core = pep[core_pos:core_pos + 9]
            per_peptide.append({
                "peptide": pep,
                "allele": best_alleles[i],
                "rank": round(1.0 - score, 6),
                "score": round(score, 6),
                "core": core,
                "core_pos": core_pos,
                "start_index": start,
            })
            for j in range(len(core)):
                gidx = start + core_pos + j
                if gidx not in per_residue or score > per_residue[gidx]:
                    per_residue[gidx] = round(score, 4)

        results[chain_id] = {
            "per_peptide": per_peptide,
            "per_residue": per_residue,
            "meta": {
                "alleles": alleles,
                "allele_chunks": len(chunks),
                "pep_len": pep_len,
                "n_folds": len(predictors),
                "n_peptides": n_pep,
                "model": "HLAIIPred",
                "score_type": "presentation_probability",
            },
        }

    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: hlaiipred_runner.py <input.json>", file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    with open(input_path) as fh:
        data = json.load(fh)

    sequences = data["sequences"]
    alleles = data.get("alleles", ["DRB1*01:01"])
    pep_len = data.get("pep_len", 15)
    models_dir = data.get("models_dir", "models")

    try:
        results = run_hlaiipred(sequences, alleles, pep_len, models_dir)
        json.dump(results, sys.stdout, indent=2)
    except Exception as exc:
        import traceback
        traceback.print_exc(file=sys.stderr)
        error_results = {}
        for chain_id in sequences:
            error_results[chain_id] = {
                "per_peptide": [],
                "per_residue": {},
                "error": str(exc),
                "meta": {},
            }
        json.dump(error_results, sys.stdout, indent=2)
        sys.exit(1)


if __name__ == "__main__":
    main()
