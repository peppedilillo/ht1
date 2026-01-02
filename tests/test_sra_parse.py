"""
Test script for sra_parse function.
"""

from pathlib import Path

import pytest

from fm1trig.search import sra_parse, _NQUADRANTS, _NBANDS, InvalidSRA


def test_sra_parse_all_files():
    """Test sra_parse on all SRA files in test data, checking abts and data are not empty."""
    test_data_dir = Path(__file__).parent / "data"
    sra_files = list(test_data_dir.glob("**/SRA/*.raw"))

    assert len(sra_files) > 0, "No SRA files found in test data directory"

    print(f"\nTesting sra_parse on {len(sra_files)} SRA files:")

    for sra_file in sorted(sra_files):
        data, abts = sra_parse(sra_file)

        # Check that abts is not empty
        assert len(abts) > 0, f"ABTs list is empty for {sra_file}"

        # Check that data for each quadrant and each band is not empty
        for band in range(_NBANDS):
            for quad in range(_NQUADRANTS):
                assert len(data[band][quad]) > 0, f"Data is empty for band {band}, quadrant {quad} in {sra_file}"

        print(f"{sra_file.parent.parent.name}/{sra_file.parent.name}/{sra_file.name}: "
              f"{len(abts)} ABTs, {len(data[0][0])} data points per band/quad")

    print(f"PASS: All {len(sra_files)} files parsed successfully with valid data")


def test_sra_parse_empty_file():
    """Empty file should raise ValueError about file size smaller than header."""
    empty_file = Path(__file__).parent / "data" / "invalid" / "empty"
    with pytest.raises(InvalidSRA) as exc_info:
        sra_parse(empty_file)
    assert "File size smaller than header size" in str(exc_info.value)


def test_sra_parse_invalid_header():
    """File with invalid header (ASCII art) should raise ValueError about invalid header fsize."""
    invalid_file = Path(__file__).parent / "data" / "invalid" / "invalid.txt"
    with pytest.raises(InvalidSRA) as exc_info:
        sra_parse(invalid_file)
    assert "Invalid header fsize" in str(exc_info.value)


if __name__ == "__main__":
    test_sra_parse_all_files()
    test_sra_parse_empty_file()
    test_sra_parse_invalid_header()
    print("\nAll tests passed!")
