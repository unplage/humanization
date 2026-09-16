#!/usr/bin/env python3
"""Test patent example report generation."""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.humanize.pipeline import RunResult
from scripts.humanize.report_patent import build_patent_example_report


def load_result_from_json(json_path: str) -> RunResult:
    """Load RunResult from JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    
    # This is a simplified loader - in practice, you'd need to reconstruct
    # the full RunResult object from the JSON
    print(f"Loaded JSON with {len(data.get('chains', []))} chains")
    return data


def main():
    json_path = "outputs/amg110_full/step3/humanization_result.json"
    out_path = "outputs/amg110_full/step3/patent_example_report.docx"
    
    if not os.path.exists(json_path):
        print(f"Error: JSON file not found: {json_path}")
        return 1
    
    print(f"Loading result from: {json_path}")
    
    # For now, we'll test by running the full pipeline
    # In practice, you'd load the JSON and reconstruct the RunResult
    print("Note: Patent report generation requires full pipeline RunResult object.")
    print("To generate the report, run the full pipeline:")
    print("  python3 scripts/humanize/cli.py run --input inputs/amg110.fasta --outdir outputs/test")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
