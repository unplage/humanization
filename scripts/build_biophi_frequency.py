#!/usr/bin/env python3
"""Extract BioPhi's OAS-derived germline residue frequencies.

BioPhi ships ``GERMLINE_FREQUENCY`` in
``biophi/humanization/methods/stats.py``: per-IMGT-position amino-acid
frequencies observed in human antibody repertoires (Observed Antibody Space),
both chain-wide (H/K/L) and per germline family (IGHV1, IGKV1D, ...).

The humanization pipeline is stdlib-only, so this script freezes that table
into ``data/germline/biophi_germline_frequency.json``. Re-run it when BioPhi
is upgraded:

    ~/miniconda3/envs/biophi/bin/python scripts/build_biophi_frequency.py

BioPhi is MIT-licensed (Copyright (c) 2021 Merck Sharp & Dohme Corp.).
"""
from __future__ import annotations

import json
import os

from biophi.__version__ import __version__ as biophi_version
from biophi.humanization.methods.stats import GERMLINE_FREQUENCY

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DST = os.path.join(ROOT, "data", "germline", "biophi_germline_frequency.json")

CHAIN_KEYS = ("H", "K", "L")
VERSION = 1


def main() -> None:
    chain = {key: value for key, value in GERMLINE_FREQUENCY.items()
             if key in CHAIN_KEYS}
    family = {key: value for key, value in GERMLINE_FREQUENCY.items()
              if key not in CHAIN_KEYS}
    missing = [key for key in CHAIN_KEYS if key not in chain]
    if missing:
        raise SystemExit(f"BioPhi table is missing chain keys: {missing}")
    if not family:
        raise SystemExit("BioPhi table has no germline-family keys")

    blob = {
        "version": VERSION,
        "source": "BioPhi GERMLINE_FREQUENCY "
                  "(biophi/humanization/methods/stats.py)",
        "biophi_version": biophi_version,
        "license": "MIT",
        "citation": "Prihoda et al. (2022) BioPhi, mAbs 14:1",
        "note": "IMGT-numbered per-position residue frequencies in OAS human "
                "repertoires; chain-wide H/K/L plus germline-family tables.",
        "chain": chain,
        "family": family,
    }
    os.makedirs(os.path.dirname(DST), exist_ok=True)
    tmp = DST + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(blob, fh, sort_keys=True, separators=(",", ":"))
        fh.write("\n")
    os.replace(tmp, DST)
    print(f"wrote {DST} (BioPhi {biophi_version}, "
          f"{len(chain)} chain tables, {len(family)} family tables)")


if __name__ == "__main__":
    main()
