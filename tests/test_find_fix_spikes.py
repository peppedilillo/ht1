#!/usr/bin/env python3
"""
Test script for find_fix_spikes function on real SRA data.
"""

from pathlib import Path
from fm1trig.search import sra_parse, find_fix_spikes


def test_find_fix_spikes_all_files():
    """Test find_fix_spikes on all SRA files in test data without raising errors."""
    test_data_dir = Path(__file__).parent / "data"
    sra_files = list(test_data_dir.glob("**/SRA/*.raw"))

    assert len(sra_files) > 0, "No SRA files found in test data directory"

    print(f"\nTesting find_fix_spikes on {len(sra_files)} SRA files:")

    for sra_file in sorted(sra_files):
        data, abts = sra_parse(sra_file)

        # Run spike detection and fixing - should not raise any errors
        find_fix_spikes(data, thr=10.0)

        print(f"{sra_file.parent.parent.name}/{sra_file.parent.name}/{sra_file.name}: OK")

    print(f"PASS: All {len(sra_files)} files processed successfully")


if __name__ == "__main__":
    test_find_fix_spikes_all_files()
    print("\nAll tests passed!")
