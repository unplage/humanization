#!/usr/bin/env python3
"""HLAIIPred runner script.

Runs inside the 'hlapred' conda environment.
Reads JSON input file, predicts MHC-II peptide presentation, outputs JSON.

Usage:
    conda run -n hlapred python tools/immunogenicity/runners/hlaiipred_runner.py input.json > output.json
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import torch
from hlapred.predict import HLAIIPredict


def run_hlaiipred(sequences, alleles, pep_len=15, models_dir="models"):
    """Run HLAIIPred on sequences.

    Args:
        sequences: {chain_id: amino_acid_sequence}
        alleles: list of HLA-II allele strings
        pep_len: peptide window length
        models_dir: path to HLAIIPred models directory
    Returns:
        {chain_id: {per_peptide: [...], per_residue: {...}, meta: {...}}}
    """
    device = torch.device("cpu" if not torch.cuda.is_available() else "cuda:0")

    # Find MHC-II data directory (allele definitions)
    mhc2_root = None
    for candidate in [
        os.path.join(models_dir, "mhcII"),
        os.path.join(os.path.dirname(models_dir), "mhcII"),
        os.path.expanduser("~/.hlaiipred/mhcII"),
    ]:
        if os.path.isdir(candidate):
            mhc2_root = candidate
            break

    # Load 2-fold models
    predictors = []
    for fold_idx in range(2):
        model_path = os.path.join(models_dir)
        try:
            p = HLAIIPredict(model_path, fold_idx, device, mhc2_root)
            predictors.append(p)
        except Exception as e:
            print(f"Warning: Failed to load fold {fold_idx}: {e}", file=sys.stderr)

    if not predictors:
        raise RuntimeError("No HLAIIPred models loaded")

    # Alleles should be in format like DRB1*01:01 (matching FASTA IDs)
    # No normalization needed - HLAIIPred uses the original format
    norm_alleles = list(alleles)

    results = {}
    for chain_id, seq in sequences.items():
        # Generate peptides
        peptides = []
        pep_starts = []
        for i in range(0, len(seq) - pep_len + 1):
            peptides.append(seq[i:i + pep_len])
            pep_starts.append(i)

        if not peptides:
            results[chain_id] = {
                "per_peptide": [],
                "per_residue": {},
                "meta": {"error": "sequence too short"},
            }
            continue

        # Prepare allele input (pad to 14 alleles as required by HLAIIPred)
        allele_input = []
        for a in norm_alleles:
            allele_input.append(a)
        allele_input += [0] * (14 - len(allele_input))
        allele_input_batch = [allele_input] * len(peptides)

        # Run predictions (average over 2 folds)
        all_scores = None
        for p in predictors:
            try:
                inputs = p.prepare_input(peptides, allele_input_batch)
                y_pred, scores = p.predict(inputs, batch_size=300, sigmoid=True)
                if all_scores is None:
                    all_scores = y_pred
                else:
                    all_scores = all_scores + y_pred
            except Exception as e:
                import traceback
                tb = traceback.format_exc()
                print(f"Warning: Prediction failed for fold: {e}", file=sys.stderr)
                print(tb, file=sys.stderr)
                raise  # Re-raise to be caught by outer handler

        if all_scores is None:
            results[chain_id] = {
                "per_peptide": [],
                "per_residue": {},
                "meta": {"error": "prediction failed"},
            }
            continue

        all_scores = all_scores / len(predictors)  # average folds

        # Build per-peptide results
        per_peptide = []
        per_residue = {}
        for i, (pep, start) in enumerate(zip(peptides, pep_starts)):
            score = float(all_scores[i].item()) if i < len(all_scores) else 0.0
            # Find best allele
            best_allele = norm_alleles[0] if norm_alleles else ""

            per_peptide.append({
                "peptide": pep,
                "allele": best_allele,
                "rank": 1.0 - score,  # convert sigmoid to rank-like (lower = stronger)
                "score": score,
                "core": "",
                "core_pos": 0,
            })

            # Map to per-residue
            for j in range(len(pep)):
                gidx = start + j
                if gidx not in per_residue or score > per_residue[gidx]:
                    per_residue[gidx] = round(score, 4)

        results[chain_id] = {
            "per_peptide": per_peptide,
            "per_residue": per_residue,
            "meta": {
                "alleles": norm_alleles,
                "pep_len": pep_len,
                "n_folds": len(predictors),
                "model": "HLAIIPred",
            },
        }

    return results


def main():
    if len(sys.argv) < 2:
        print("Usage: hlaiipred_runner.py <input.json>", file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    with open(input_path) as f:
        data = json.load(f)

    sequences = data["sequences"]
    alleles = data.get("alleles", ["DRB1*01:01"])
    pep_len = data.get("pep_len", 15)
    models_dir = data.get("models_dir", "models")

    try:
        results = run_hlaiipred(sequences, alleles, pep_len, models_dir)
        json.dump(results, sys.stdout, indent=2)
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"ERROR in run_hlaiipred: {e}", file=sys.stderr)
        print(tb, file=sys.stderr)
        # Output error as JSON
        error_results = {}
        for chain_id in sequences:
            error_results[chain_id] = {
                "per_peptide": [],
                "per_residue": {},
                "meta": {"error": str(e)},
            }
        json.dump(error_results, sys.stdout, indent=2)
        sys.exit(1)


if __name__ == "__main__":
    main()
