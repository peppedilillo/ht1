"""
Test script for sra_parse function.
"""

from pathlib import Path
from fm1trig.search import sra_parse, NQUADRANTS, NBANDS


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
        for band in range(NBANDS):
            for quad in range(NQUADRANTS):
                assert len(data[band][quad]) > 0, f"Data is empty for band {band}, quadrant {quad} in {sra_file}"

        print(f"{sra_file.parent.parent.name}/{sra_file.parent.name}/{sra_file.name}: "
              f"{len(abts)} ABTs, {len(data[0][0])} data points per band/quad")

    print(f"PASS: All {len(sra_files)} files parsed successfully with valid data")


if __name__ == "__main__":
    test_sra_parse_all_files()
    print("\nAll tests passed!")
