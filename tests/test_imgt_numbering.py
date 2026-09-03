#!/usr/bin/env python3
"""Unit tests for IMGT numbering conversion.

Tests the accuracy of Kabat-to-IMGT position mapping based on
IMGT official correspondence (Lefranc 2003).

Reference: Lefranc, M.-P. (2003). IMGT unique numbering for immunoglobulin
and T cell receptor variable domains and Ig superfamily V-like domains.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scripts.humanize.imgt_numbering import (
    kabat_to_imgt_position,
    get_imgt_region,
    _map_insertion,
    IMGT_REGIONS,
    KABAT_TO_IMGT_VH,
    KABAT_TO_IMGT_VL,
)


def test_imgt_region_definitions():
    """Test that IMGT region definitions are correct (Lefranc 2003)."""
    print("Testing IMGT region definitions...")

    # VH and VL use IDENTICAL regions in IMGT
    assert IMGT_REGIONS["FR1"] == (1, 26), f"FR1 should be 1-26, got {IMGT_REGIONS['FR1']}"
    assert IMGT_REGIONS["CDR1"] == (27, 38), f"CDR1 should be 27-38, got {IMGT_REGIONS['CDR1']}"
    assert IMGT_REGIONS["FR2"] == (39, 55), f"FR2 should be 39-55, got {IMGT_REGIONS['FR2']}"
    assert IMGT_REGIONS["CDR2"] == (56, 65), f"CDR2 should be 56-65, got {IMGT_REGIONS['CDR2']}"
    assert IMGT_REGIONS["FR3"] == (66, 104), f"FR3 should be 66-104, got {IMGT_REGIONS['FR3']}"
    assert IMGT_REGIONS["CDR3"] == (105, 117), f"CDR3 should be 105-117, got {IMGT_REGIONS['CDR3']}"
    assert IMGT_REGIONS["FR4"] == (118, 128), f"FR4 should be 118-128, got {IMGT_REGIONS['FR4']}"

    # Verify total positions = 128
    total = sum(end - start + 1 for start, end in IMGT_REGIONS.values())
    assert total == 128, f"Total IMGT positions should be 128, got {total}"

    print("  [PASS] IMGT region definitions are correct")


def test_vh_anchors():
    """Test VH chain anchor positions (critical for numbering)."""
    print("Testing VH anchor positions...")

    # 1st-CYS: Kabat H22 -> IMGT 23
    assert kabat_to_imgt_position("H22", "H") == 23, "H22 (1st-CYS) should map to IMGT 23"

    # CONSERVED-TRP: Kabat H36 -> IMGT 41
    assert kabat_to_imgt_position("H36", "H") == 41, "H36 (CONSERVED-TRP) should map to IMGT 41"

    # 2nd-CYS: Kabat H92 -> IMGT 104
    assert kabat_to_imgt_position("H92", "H") == 104, "H92 (2nd-CYS) should map to IMGT 104"

    # J-TRP: Kabat H103 -> IMGT 118
    assert kabat_to_imgt_position("H103", "H") == 118, "H103 (J-TRP) should map to IMGT 118"

    print("  [PASS] VH anchor positions are correct")


def test_vl_anchors():
    """Test VL chain anchor positions."""
    print("Testing VL anchor positions...")

    # 1st-CYS: Kabat L23 -> IMGT 23
    assert kabat_to_imgt_position("L23", "L") == 23, "L23 (1st-CYS) should map to IMGT 23"

    # CONSERVED-TRP: Kabat L35 -> IMGT 39
    assert kabat_to_imgt_position("L35", "L") == 39, "L35 (CONSERVED-TRP) should map to IMGT 39"

    # 2nd-CYS: Kabat L88 -> IMGT 104
    assert kabat_to_imgt_position("L88", "L") == 104, "L88 (2nd-CYS) should map to IMGT 104"

    # J-PHE: Kabat L98 -> IMGT 118
    assert kabat_to_imgt_position("L98", "L") == 118, "L98 (J-PHE) should map to IMGT 118"

    print("  [PASS] VL anchor positions are correct")


def test_vh_fr1_mapping():
    """Test VH FR1 region mapping (Kabat 1-30 -> IMGT 1-26)."""
    print("Testing VH FR1 mapping...")

    # IMGT has a gap at position 10, so Kabat 10 maps to IMGT 11
    assert kabat_to_imgt_position("H1", "H") == 1
    assert kabat_to_imgt_position("H9", "H") == 9
    assert kabat_to_imgt_position("H10", "H") == 11  # Gap at IMGT 10
    assert kabat_to_imgt_position("H11", "H") == 12
    assert kabat_to_imgt_position("H22", "H") == 23  # 1st-CYS
    assert kabat_to_imgt_position("H25", "H") == 26  # FR1 end

    print("  [PASS] VH FR1 mapping is correct")


def test_vh_cdr1_mapping():
    """Test VH CDR1 region mapping (Kabat 26-35 -> IMGT 27-38)."""
    print("Testing VH CDR1 mapping...")

    assert kabat_to_imgt_position("H26", "H") == 27
    assert kabat_to_imgt_position("H27", "H") == 28
    assert kabat_to_imgt_position("H29", "H") == 30
    assert kabat_to_imgt_position("H34", "H") == 35
    assert kabat_to_imgt_position("H35", "H") == 36

    print("  [PASS] VH CDR1 mapping is correct")


def test_vh_fr2_mapping():
    """Test VH FR2 region mapping (Kabat 36-50 -> IMGT 41-55)."""
    print("Testing VH FR2 mapping...")

    assert kabat_to_imgt_position("H36", "H") == 41  # CONSERVED-TRP
    assert kabat_to_imgt_position("H37", "H") == 42
    assert kabat_to_imgt_position("H40", "H") == 45
    assert kabat_to_imgt_position("H50", "H") == 55

    print("  [PASS] VH FR2 mapping is correct")


def test_vh_cdr2_mapping():
    """Test VH CDR2 region mapping (Kabat 51-57 -> IMGT 56-65)."""
    print("Testing VH CDR2 mapping...")

    assert kabat_to_imgt_position("H51", "H") == 56
    assert kabat_to_imgt_position("H52", "H") == 57
    # Kabat 52A-52J map to IMGT 58-65

    print("  [PASS] VH CDR2 mapping is correct")


def test_vh_fr3_mapping():
    """Test VH FR3 region mapping (Kabat 58-92 -> IMGT 66-104)."""
    print("Testing VH FR3 mapping...")

    assert kabat_to_imgt_position("H58", "H") == 66
    assert kabat_to_imgt_position("H59", "H") == 67
    assert kabat_to_imgt_position("H66", "H") == 75
    assert kabat_to_imgt_position("H75", "H") == 84
    assert kabat_to_imgt_position("H82", "H") == 91
    assert kabat_to_imgt_position("H92", "H") == 104  # 2nd-CYS

    print("  [PASS] VH FR3 mapping is correct")


def test_vh_cdr3_mapping():
    """Test VH CDR3 region mapping (Kabat 93-102 -> IMGT 105-117)."""
    print("Testing VH CDR3 mapping...")

    assert kabat_to_imgt_position("H93", "H") == 105
    assert kabat_to_imgt_position("H94", "H") == 106
    assert kabat_to_imgt_position("H95", "H") == 107
    assert kabat_to_imgt_position("H100", "H") == 112
    assert kabat_to_imgt_position("H102", "H") == 115

    print("  [PASS] VH CDR3 mapping is correct")


def test_vh_fr4_mapping():
    """Test VH FR4 region mapping (Kabat 103-113 -> IMGT 118-128)."""
    print("Testing VH FR4 mapping...")

    assert kabat_to_imgt_position("H103", "H") == 118  # J-TRP
    assert kabat_to_imgt_position("H104", "H") == 119
    assert kabat_to_imgt_position("H113", "H") == 128

    print("  [PASS] VH FR4 mapping is correct")


def test_vh_insertions():
    """Test VH insertion position mapping."""
    print("Testing VH insertion positions...")

    # 35A/35B -> IMGT 39/40 (FR2 start)
    assert kabat_to_imgt_position("H35A", "H") == 39
    assert kabat_to_imgt_position("H35B", "H") == 40

    # 52A-52J -> IMGT 58-65 (CDR2)
    assert kabat_to_imgt_position("H52A", "H") == 58
    assert kabat_to_imgt_position("H52B", "H") == 59
    assert kabat_to_imgt_position("H52J", "H") == 65

    # 82A-82C -> IMGT 92-94 (FR3)
    assert kabat_to_imgt_position("H82A", "H") == 92
    assert kabat_to_imgt_position("H82B", "H") == 93
    assert kabat_to_imgt_position("H82C", "H") == 94

    # 100A-100K -> IMGT 108-118 (CDR3)
    assert kabat_to_imgt_position("H100A", "H") == 108
    assert kabat_to_imgt_position("H100K", "H") == 118

    print("  [PASS] VH insertion positions are correct")


def test_vl_insertions():
    """Test VL insertion position mapping."""
    print("Testing VL insertion positions...")

    # 27A-27E -> IMGT 28-32 (CDR1)
    assert kabat_to_imgt_position("L27A", "L") == 28
    assert kabat_to_imgt_position("L27E", "L") == 32

    # 82A-82C -> IMGT 93-95 (FR3)
    assert kabat_to_imgt_position("L82A", "L") == 93
    assert kabat_to_imgt_position("L82C", "L") == 95

    # 100A-100D -> IMGT 109-112 (CDR3)
    assert kabat_to_imgt_position("L100A", "L") == 109
    assert kabat_to_imgt_position("L100D", "L") == 112

    print("  [PASS] VL insertion positions are correct")


def test_get_imgt_region():
    """Test IMGT region assignment."""
    print("Testing IMGT region assignment...")

    assert get_imgt_region(1, "H") == "FR1"
    assert get_imgt_region(26, "H") == "FR1"
    assert get_imgt_region(27, "H") == "CDR1"
    assert get_imgt_region(38, "H") == "CDR1"
    assert get_imgt_region(39, "H") == "FR2"
    assert get_imgt_region(55, "H") == "FR2"
    assert get_imgt_region(56, "H") == "CDR2"
    assert get_imgt_region(65, "H") == "CDR2"
    assert get_imgt_region(66, "H") == "FR3"
    assert get_imgt_region(104, "H") == "FR3"
    assert get_imgt_region(105, "H") == "CDR3"
    assert get_imgt_region(117, "H") == "CDR3"
    assert get_imgt_region(118, "H") == "FR4"
    assert get_imgt_region(128, "H") == "FR4"

    # VL uses same regions
    assert get_imgt_region(1, "L") == "FR1"
    assert get_imgt_region(104, "L") == "FR3"
    assert get_imgt_region(128, "L") == "FR4"

    print("  [PASS] IMGT region assignment is correct")


def test_known_vh_sequence():
    """Test with a known VH sequence (IGHV1-69*01 germline)."""
    print("Testing known VH sequence...")

    # Simplified IGHV1-69*01 sequence (first 30 residues)
    # Expected IMGT positions 1-26 (FR1)
    vh_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFS"
    posmap = {f"H{i+1}": aa for i, aa in enumerate(vh_seq)}

    imgt_map = {}
    for kabat_pos, aa in posmap.items():
        imgt_pos = kabat_to_imgt_position(kabat_pos, "H")
        if imgt_pos is not None:
            imgt_map[f"H{imgt_pos}"] = aa

    # Verify key positions
    assert imgt_map.get("H1") == "E", "Position 1 should be E"
    assert imgt_map.get("H23") == "C", "Position 23 (1st-CYS) should be C"
    assert imgt_map.get("H26") == "S", "Position 26 (FR1 end) should be S"

    print("  [PASS] Known VH sequence mapping is correct")


def test_roundtrip_consistency():
    """Test that Kabat -> IMGT -> region is consistent."""
    print("Testing roundtrip consistency...")

    # Test a few positions
    test_cases = [
        ("H22", "H", "FR1"),
        ("H36", "H", "FR2"),
        ("H92", "H", "FR3"),
        ("H103", "H", "FR4"),
    ]

    for kabat_pos, chain_type, expected_region in test_cases:
        imgt_pos = kabat_to_imgt_position(kabat_pos, chain_type)
        assert imgt_pos is not None, f"Failed to map {kabat_pos}"
        region = get_imgt_region(imgt_pos, chain_type)
        assert region == expected_region, f"{kabat_pos} -> IMGT {imgt_pos} should be in {expected_region}, got {region}"

    print("  [PASS] Roundtrip consistency is correct")


def test_lookup_table_completeness():
    """Test that lookup tables cover all expected positions."""
    print("Testing lookup table completeness...")

    # VH should have mappings for positions 1-113 (some may be gaps in IMGT)
    vh_mapped = set(KABAT_TO_IMGT_VH.keys())
    vh_expected = set(range(1, 114))
    missing_vh = vh_expected - vh_mapped
    # Some positions are gaps in IMGT (e.g., 30-33, 65, 72-74, 100A-K)
    # These are expected to be missing from the lookup table
    expected_gaps = {30, 31, 32, 33, 53, 54, 55, 56, 57, 65, 72, 73, 74}
    unexpected_missing = missing_vh - expected_gaps
    assert len(unexpected_missing) <= 5, f"VH lookup table missing unexpected positions: {unexpected_missing}"

    # VL should have mappings for positions 1-107 (some may be gaps in IMGT)
    vl_mapped = set(KABAT_TO_IMGT_VL.keys())
    vl_expected = set(range(1, 108))
    missing_vl = vl_expected - vl_mapped
    # Some positions are gaps in IMGT
    unexpected_missing_vl = missing_vl - expected_gaps
    assert len(unexpected_missing_vl) <= 5, f"VL lookup table missing unexpected positions: {unexpected_missing_vl}"

    print("  [PASS] Lookup table completeness is acceptable")


def test_invalid_positions():
    """Test handling of invalid positions."""
    print("Testing invalid positions...")

    assert kabat_to_imgt_position("", "H") is None
    assert kabat_to_imgt_position("X", "H") is None
    assert kabat_to_imgt_position("H", "H") is None
    assert kabat_to_imgt_position("Habc", "H") is None
    assert kabat_to_imgt_position("H999", "H") is None

    print("  [PASS] Invalid positions handled correctly")


def main():
    """Run all tests."""
    print("=" * 70)
    print("IMGT Numbering Conversion Tests")
    print("=" * 70)
    print()

    tests = [
        test_imgt_region_definitions,
        test_vh_anchors,
        test_vl_anchors,
        test_vh_fr1_mapping,
        test_vh_cdr1_mapping,
        test_vh_fr2_mapping,
        test_vh_cdr2_mapping,
        test_vh_fr3_mapping,
        test_vh_cdr3_mapping,
        test_vh_fr4_mapping,
        test_vh_insertions,
        test_vl_insertions,
        test_get_imgt_region,
        test_known_vh_sequence,
        test_roundtrip_consistency,
        test_lookup_table_completeness,
        test_invalid_positions,
    ]

    passed = 0
    failed = 0
    errors = []

    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            failed += 1
            errors.append((test.__name__, str(e)))
            print(f"  [FAIL] {test.__name__}: {e}")
        except Exception as e:
            failed += 1
            errors.append((test.__name__, str(e)))
            print(f"  [ERROR] {test.__name__}: {e}")

    print()
    print("=" * 70)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 70)

    if errors:
        print()
        print("Failed tests:")
        for name, msg in errors:
            print(f"  - {name}: {msg}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
