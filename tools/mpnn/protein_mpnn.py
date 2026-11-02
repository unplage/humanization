#!/usr/bin/env python3
"""Deprecated alias for the ProteinMPNN adapter.

Historically this file was a standalone wrapper that looked for a
``protein_mpnn.py`` script and, failing that, wrote a placeholder
``consensus.fasta``. That interface no longer exists upstream (the entry point
is now ``protein_mpnn_run.py``), so the old wrapper could never produce real
designs.

It is kept only for backward compatibility: it now forwards its arguments
verbatim to ``proteinmpnn_adapter.py`` in the same directory. New code should
call the adapter directly.

    --mpnn-script tools/mpnn/proteinmpnn_adapter.py
"""

import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ADAPTER = os.path.join(HERE, "proteinmpnn_adapter.py")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    return subprocess.call([sys.executable, ADAPTER, *argv])


if __name__ == "__main__":
    sys.exit(main())
