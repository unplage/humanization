#!/usr/bin/env python3
"""Run ThermoMPNN inference for an explicit list of point mutations.

This helper is executed with the ThermoMPNN Python environment. It monkey-
patches ``SSM.get_ssm_mutations`` so that ThermoMPNN's own
``analysis/custom_inference.py`` evaluates only the requested
``(wildtype, seq_pos, mutant)`` entries instead of the full 20 x L
site-saturation library.

It intentionally imports nothing outside the Python standard library and the
ThermoMPNN repository, so it can run in the isolated ThermoMPNN environment.
"""

from __future__ import annotations

import argparse
import json
import os
import runpy
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--thermompnn-dir", required=True)
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--model-path", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--targets", required=True,
                        help="JSON list of [wildtype, seq_pos, mutant]")
    args = parser.parse_args()

    tm_dir = os.path.abspath(os.path.expanduser(args.thermompnn_dir))
    analysis_dir = os.path.join(tm_dir, "analysis")
    for path in (analysis_dir, tm_dir):
        if path not in sys.path:
            sys.path.insert(0, path)

    with open(args.targets, "r", encoding="utf-8") as fh:
        raw_targets = json.load(fh)

    targets = []
    for item in raw_targets:
        wt, pos, mut = item
        targets.append((str(wt).strip().upper(), int(pos), str(mut).strip().upper()))

    import SSM  # type: ignore

    def _patched_get_ssm_mutations(_pdb_dict):
        return [f"{wt}{pos}{mut}" for wt, pos, mut in targets]

    SSM.get_ssm_mutations = _patched_get_ssm_mutations

    script = os.path.join(analysis_dir, "custom_inference.py")
    if not os.path.exists(script):
        print(f"custom_inference.py not found: {script}", file=sys.stderr)
        return 2

    sys.argv = [
        "custom_inference.py",
        "--pdb", args.pdb,
        "--chain", args.chain,
        "--model_path", args.model_path,
        "--out_dir", args.out_dir,
    ]
    runpy.run_path(script, run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
