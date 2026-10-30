#!/usr/bin/env python3
"""Compatibility entry point for the ThermoMPNN stability analyzer.

All logic now lives in :mod:`thermompnn_analyzer`. This module exists so the
documented command ``python3 tools/stability_predictor/stability_analyzer.py``
and the pipeline CLI (``humanize stability``) keep working, and so the script
can be imported either as part of the package or run directly.
"""

from __future__ import annotations

import sys

try:  # imported as tools.stability_predictor.stability_analyzer
    from .thermompnn_analyzer import (  # noqa: F401
        BackMutationCandidate,
        MutationResult,
        StabilityReport,
        analyze_batch,
        analyze_structure,
        main,
    )
except ImportError:  # executed directly as a script
    from thermompnn_analyzer import (  # type: ignore  # noqa: F401
        BackMutationCandidate,
        MutationResult,
        StabilityReport,
        analyze_batch,
        analyze_structure,
        main,
    )


if __name__ == "__main__":
    sys.exit(main())
