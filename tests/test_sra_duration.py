"""
Test script for sra_duration function.
"""

from pathlib import Path
from fm1trig.extras import sra_duration

def test_sra_duration_all_files():
    """Test sra_duration on all SRA files in test data, and checks their duration to be an integer."""
    test_data_dir = Path(__file__).parent / "data"
    sra_files = list(test_data_dir.glob("**/SRA/*.raw"))

    assert len(sra_files) > 0, "No SRA files found in test data directory"

    print(f"\nTesting sra_duration on {len(sra_files)} SRA files:")

    for sra_file in sorted(sra_files):
        duration = sra_duration(sra_file)
        assert duration.is_integer(), "SRA files duration is not integer"
        print(f"{sra_file.parent.parent.name}/{sra_file.parent.name}/{sra_file.name}: {duration:.2f} seconds")

    print(f"PASS: All {len(sra_files)} files returned valid durations")


if __name__ == "__main__":
    test_sra_duration_all_files()
    print("\nAll tests passed!")
