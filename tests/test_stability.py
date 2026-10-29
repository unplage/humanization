#!/usr/bin/env python3
"""Tests for stability_predictor ThermoMPNN analyzer.

Tests CSV parsing, position filtering, region annotation, and output formatting
without requiring ThermoMPNN installation.
"""

import os
import sys
import tempfile
import json
import csv

# Add repo root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tools.stability_predictor.thermompnn_analyzer import (
    parse_thermompnn_csv,
    MutationResult,
    StabilityReport,
    format_text_report,
    save_json_report,
    save_markdown_report,
    parse_positions,
)
from tools.stability_predictor.stability_analyzer import (
    StabilityResult as StabilityAnalyzerResult,
    BatchStabilityResult,
    format_result_json,
    format_batch_json,
)


def test_parse_positions():
    """Test Kabat position parsing."""
    # Without chain prefix
    assert parse_positions("37,44,45") == ["H37", "H44", "H45"]
    # With chain prefix
    assert parse_positions("H37,H44,H45") == ["H37", "H44", "H45"]
    # Mixed
    assert parse_positions("37,H44,45") == ["H37", "H44", "H45"]
    # Light chain
    assert parse_positions("L24,L27") == ["L24", "L27"]
    # Whitespace
    assert parse_positions(" 37 , 44 , 45 ") == ["H37", "H44", "H45"]
    # Empty
    assert parse_positions("") == []


def test_parse_thermompnn_csv():
    """Test parsing ThermoMPNN output CSV."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.DictWriter(f, fieldnames=[
            'Model', 'Dataset', 'ddG_pred', 'position', 'wildtype', 'mutation', 'pdb', 'chain'
        ])
        writer.writeheader()
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '-1.23',
            'position': '37', 'wildtype': 'V', 'mutation': 'F', 'pdb': 'test', 'chain': 'H'
        })
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '0.45',
            'position': '44', 'wildtype': 'G', 'mutation': 'E', 'pdb': 'test', 'chain': 'H'
        })
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '2.10',
            'position': '45', 'wildtype': 'L', 'mutation': 'R', 'pdb': 'test', 'chain': 'H'
        })
        csv_path = f.name
    
    try:
        # Parse all
        mutations = parse_thermompnn_csv(csv_path, "test.pdb", "H")
        assert len(mutations) == 3
        assert mutations[0].position_kabat == "H37"
        assert mutations[0].ddG == -1.23
        assert mutations[0].is_stabilizing == True
        assert mutations[0].is_significant == True
        
        assert mutations[1].position_kabat == "H44"
        assert mutations[1].ddG == 0.45
        assert mutations[1].is_stabilizing == False
        assert mutations[1].is_significant == False  # |0.45| < 1.0
        
        assert mutations[2].position_kabat == "H45"
        assert mutations[2].ddG == 2.10
        assert mutations[2].is_stabilizing == False
        assert mutations[2].is_significant == True
        
        # Parse with position filter
        mutations_filtered = parse_thermompnn_csv(csv_path, "test.pdb", "H", positions_kabat=["H37", "H45"])
        assert len(mutations_filtered) == 2
        assert mutations_filtered[0].position_kabat == "H37"
        assert mutations_filtered[1].position_kabat == "H45"
    finally:
        os.unlink(csv_path)


def test_mutation_result_properties():
    """Test MutationResult computed properties."""
    # Stabilizing
    mut = MutationResult(
        pdb_file="test.pdb", chain="H", position_pdb=37, position_kabat="H37",
        wildtype="V", mutation="F", ddG=-1.5,
        region_kabat="FR2", region_imgt="FR2", is_stabilizing=True, is_significant=True
    )
    assert mut.is_stabilizing == True
    assert mut.is_significant == True
    
    # Destabilizing but not significant
    mut2 = MutationResult(
        pdb_file="test.pdb", chain="H", position_pdb=44, position_kabat="H44",
        wildtype="G", mutation="E", ddG=0.5,
        region_kabat="FR2", region_imgt="FR2", is_stabilizing=False, is_significant=False
    )
    assert mut2.is_stabilizing == False
    assert mut2.is_significant == False
    
    # Destabilizing and significant
    mut3 = MutationResult(
        pdb_file="test.pdb", chain="H", position_pdb=45, position_kabat="H45",
        wildtype="L", mutation="R", ddG=2.0,
        region_kabat="FR2", region_imgt="FR2", is_stabilizing=False, is_significant=True
    )
    assert mut3.is_stabilizing == False
    assert mut3.is_significant == True


def test_format_text_report():
    """Test text report formatting."""
    report = StabilityReport(
        pdb_file="test.pdb",
        chain="H",
        mutations=[
            MutationResult(
                pdb_file="test.pdb", chain="H", position_pdb=37, position_kabat="H37",
                wildtype="V", mutation="F", ddG=-1.23,
                region_kabat="FR2", region_imgt="FR2", is_stabilizing=True, is_significant=True
            ),
            MutationResult(
                pdb_file="test.pdb", chain="H", position_pdb=44, position_kabat="H44",
                wildtype="G", mutation="E", ddG=0.45,
                region_kabat="FR2", region_imgt="FR2", is_stabilizing=False, is_significant=False
            ),
        ],
        summary={
            "total": 2, "stabilizing": 1, "destabilizing": 1,
            "significant": 1, "mean_ddg": -0.39, "min_ddg": -1.23, "max_ddg": 0.45
        }
    )
    
    text = format_text_report(report)
    assert "ThermoMPNN Stability Analysis" in text
    assert "H37" in text
    assert "V" in text and "F" in text  # WT and Mut in table
    assert "STABILIZING" in text
    assert "DESTABILIZING" in text
    assert "SUMMARY" in text
    assert "Total mutations analyzed:  2" in text


def test_format_json_report():
    """Test JSON report formatting using stability_analyzer format."""
    # Use stability_analyzer's StabilityResult for JSON formatting
    result = StabilityAnalyzerResult(
        pdb_file="test.pdb",
        chain_id="H",
        ddg=-1.23,
        model="ThermoMPNN",
        error=None,
    )
    
    json_data = format_result_json(result)
    assert json_data["pdb_file"] == "test.pdb"
    assert json_data["chain_id"] == "H"
    assert json_data["ddg"] == -1.23
    assert json_data["is_stabilizing"] == True
    # is_significant not in this class
    
    # Test batch JSON
    batch = BatchStabilityResult()
    batch.add_result(result)
    
    batch_json = format_batch_json(batch)
    assert batch_json["summary"]["n_mutations"] == 1
    assert len(batch_json["results"]) == 1
    assert batch_json["results"][0]["ddg"] == -1.23
    
    # Test save
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json_path = f.name
    
    try:
        with open(json_path, 'w') as f:
            json.dump(batch_json, f, indent=2)
        with open(json_path) as f:
            loaded = json.load(f)
        assert loaded["summary"]["n_mutations"] == 1
        assert loaded["results"][0]["ddg"] == -1.23
    finally:
        os.unlink(json_path)


def test_save_markdown_report():
    """Test Markdown report saving."""
    report = StabilityReport(
        pdb_file="test.pdb",
        chain="H",
        mutations=[
            MutationResult(
                pdb_file="test.pdb", chain="H", position_pdb=37, position_kabat="H37",
                wildtype="V", mutation="F", ddG=-1.23,
                region_kabat="FR2", region_imgt="FR2", is_stabilizing=True, is_significant=True
            ),
        ],
        summary={
            "total": 1, "stabilizing": 1, "destabilizing": 0,
            "significant": 1, "mean_ddg": -1.23, "min_ddg": -1.23, "max_ddg": -1.23
        }
    )
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False) as f:
        md_path = f.name
    
    try:
        save_markdown_report(report, md_path)
        with open(md_path) as f:
            content = f.read()
        assert "# ThermoMPNN Stability Analysis" in content
        assert "H37" in content
        assert "-1.23" in content
        assert "STABILIZING" in content
    finally:
        os.unlink(md_path)


def test_csv_with_errors():
    """Test CSV parsing with malformed rows."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        writer = csv.DictWriter(f, fieldnames=[
            'Model', 'Dataset', 'ddG_pred', 'position', 'wildtype', 'mutation', 'pdb', 'chain'
        ])
        writer.writeheader()
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '-1.23',
            'position': '37', 'wildtype': 'V', 'mutation': 'F', 'pdb': 'test', 'chain': 'H'
        })
        # Malformed row (missing ddG)
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '',
            'position': '44', 'wildtype': 'G', 'mutation': 'E', 'pdb': 'test', 'chain': 'H'
        })
        # Malformed row (non-numeric position)
        writer.writerow({
            'Model': 'ThermoMPNN', 'Dataset': 'test', 'ddG_pred': '0.5',
            'position': 'abc', 'wildtype': 'L', 'mutation': 'R', 'pdb': 'test', 'chain': 'H'
        })
        csv_path = f.name
    
    try:
        mutations = parse_thermompnn_csv(csv_path, "test.pdb", "H")
        # Should only parse the valid row
        assert len(mutations) == 1
        assert mutations[0].position_kabat == "H37"
    finally:
        os.unlink(csv_path)


def test_batch_stability_result():
    """Test BatchStabilityResult aggregation using stability_analyzer classes."""
    batch = BatchStabilityResult()
    
    # Add stabilizing mutation
    batch.add_result(StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=-1.23, model="ThermoMPNN"
    ))
    
    # Add destabilizing mutation
    batch.add_result(StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=0.45, model="ThermoMPNN"
    ))
    
    # Add significant destabilizing
    batch.add_result(StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=2.0, model="ThermoMPNN"
    ))
    
    assert batch.summary["n_mutations"] == 3
    assert batch.summary["n_successful"] == 3
    assert batch.summary["n_stabilizing"] == 1
    assert batch.summary["n_destabilizing"] == 2
    assert batch.summary["n_significant"] == 2
    assert abs(batch.summary["ddg_mean"] - 0.407) < 0.01
    assert batch.summary["ddg_min"] == -1.23
    assert batch.summary["ddg_max"] == 2.0


def test_stability_analyzer_result_properties():
    """Test StabilityAnalyzerResult properties."""
    # Stabilizing
    result = StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=-1.5, model="ThermoMPNN"
    )
    assert result.is_stabilizing == True
    
    # Destabilizing
    result2 = StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=0.5, model="ThermoMPNN"
    )
    assert result2.is_stabilizing == False
    
    # None ddg defaults to not stabilizing (new behavior)
    result3 = StabilityAnalyzerResult(
        pdb_file="test.pdb", chain_id="H", ddg=None, model="ThermoMPNN"
    )
    assert result3.is_stabilizing == False  # New behavior: False when ddg is None


def run_all_tests():
    """Run all tests."""
    tests = [
        test_parse_positions,
        test_parse_thermompnn_csv,
        test_mutation_result_properties,
        test_format_text_report,
        test_format_json_report,
        test_save_markdown_report,
        test_csv_with_errors,
        test_batch_stability_result,
        test_stability_analyzer_result_properties,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            print(f"PASS: {test.__name__}")
            passed += 1
        except Exception as e:
            print(f"FAIL: {test.__name__}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print(f"\n{passed} passed, {failed} failed")
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)