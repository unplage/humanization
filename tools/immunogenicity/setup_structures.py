#!/usr/bin/env python3
"""Create alias symlinks for structural-risk structure lookup.

The immunogenicity structural-risk tool locates a variant's structure by
filename. This utility links each predicted Fv PDB under the alternate names the
tool probes (see ``tools/immunogenicity/structural_risk.py``).

Supports both naming conventions:
  * ``combo_<H>_<L>.pdb``      (scripts/humanize/combinatorial_igfold.py)
  * ``H_<H>_L_<L>.pdb``        (older runs)

Usage:
    python3 tools/immunogenicity/setup_structures.py \
        --input-dir af3/amg110_step3_combinatorial \
        --output-dir outputs/.../immunogenicity_input/structures
"""

from __future__ import annotations

import argparse
import glob
import os
import re

H_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]
L_VARIANTS = ["V0", "V1", "V2a", "V2b", "V2", "V3", "Vmin"]

# repo root = parent of tools/immunogenicity/
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def find_source(input_dir: str, h: str, l: str) -> str:
    """Locate the predicted PDB for one H/L pair, tolerating both namings.

    ``combinatorial_igfold.py`` names same-index pairs ``combo_<X>.pdb`` and
    cross pairs ``combo_<H>_<L>.pdb``; older runs use ``H_<H>_L_<L>.pdb``.
    """
    candidates = [f"combo_{h}_{l}.pdb", f"H_{h}_L_{l}.pdb"]
    if h == l:
        candidates.insert(0, f"combo_{h}.pdb")
    for name in candidates:
        path = os.path.join(input_dir, name)
        if os.path.isfile(path):
            return path
    return ""


def setup_structures(input_dir: str, output_dir: str) -> int:
    """Create alias symlinks; return the number created."""
    if not os.path.isdir(input_dir):
        raise SystemExit(f"input directory not found: {input_dir}")
    os.makedirs(output_dir, exist_ok=True)

    created = 0
    for h in H_VARIANTS:
        for l in L_VARIANTS:
            src = find_source(input_dir, h, l)
            if not src:
                print(f"WARNING: no PDB for H_{h} + L_{l} in {input_dir}")
                continue
            for alias in (f"H{h}_L{l}_VH.pdb", f"H{h}_L{l}_VL.pdb", f"combo_{h}_{l}.pdb"):
                link = os.path.join(output_dir, alias)
                if not os.path.exists(link):
                    os.symlink(os.path.abspath(src), link)
                    created += 1

    print(f"Structure aliases created: {created}")
    print(f"Output directory: {output_dir}")
    return created


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create alias symlinks for structural-risk structure lookup")
    parser.add_argument("--input-dir", "-i",
                        default=os.path.join(REPO_ROOT, "af3", "amg110_step3_combinatorial"),
                        help="directory with predicted Fv PDBs")
    parser.add_argument("--output-dir", "-o", required=True,
                        help="directory to create alias symlinks in")
    args = parser.parse_args()
    setup_structures(args.input_dir, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
