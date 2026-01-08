"""
Tests for search_qbdata, search_data, and search_filepath functions.
"""

from pathlib import Path

import pytest

from ht1.ht1 import EnBand
from ht1.ht1 import InvalidSRA
from ht1.ht1 import Quadrant
from ht1.ht1 import search_data
from ht1.ht1 import search_filepath
from ht1.ht1 import search_qbdata

# --- search_qbdata tests ---


def test_search_qbdata_empty_data():
    """Empty data returns empty hits."""
    xs = []
    hits = search_qbdata(xs, size=5, maxtest=4, threshold=5.0)
    assert hits == []


def test_search_qbdata_all_zeros():
    """All-zero data returns empty hits (invalid range)."""
    xs = [0] * 100
    hits = search_qbdata(xs, size=5, maxtest=4, threshold=5.0)
    assert hits == []


def test_search_qbdata_data_too_short():
    """Data too short for moving average window returns empty hits."""
    xs = [0, 10, 10, 0]  # only 2 non-zero, size=5 needs more
    hits = search_qbdata(xs, size=5, maxtest=4, threshold=5.0)
    assert hits == []


def test_search_qbdata_quiet_data():
    """Constant quiet data returns no hits."""
    xs = [10] * 500
    hits = search_qbdata(xs, size=20, maxtest=8, threshold=5.0)
    assert hits == []


def test_search_qbdata_with_transient():
    """Data with obvious transient returns hits."""
    xs = [10] * 500
    xs[300] = 500  # big spike
    hits = search_qbdata(xs, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0
    # Hit should be near the spike
    assert any(abs(t - 300) <= 8 for t, _ in hits)


# --- search_data tests ---


def _make_test_data(npoints, transient_idx=None, transient_value=500):
    """Helper to create test Data structure with optional transient."""
    # Data structure: data[band][quadrant] = list of counts
    # 3 bands, 3 quadrants
    data = tuple(tuple([10] * npoints for _ in range(3)) for _ in range(3))
    if transient_idx is not None:
        # Insert transient in all bands and quadrants
        for band in range(3):
            for quad in range(3):
                data[band][quad][transient_idx] = transient_value
    return data


def test_search_data_different_transient_times():
    """Transients at different times in different checks don't coincide."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))
    # Transient at different times in different quadrants
    data[EnBand.MID][Quadrant.B][300] = 500
    data[EnBand.MID][Quadrant.C][350] = 500  # different time
    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    # No coincident hits
    assert hits == []


def test_search_data_only_4_channels_insufficient():
    """Very obvious transient in only 4 quadrant/band combinations is missed."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))

    # Insert very strong transient (10x background) in only 4 channels
    # MID band: all 3 quadrants (3 channels)
    data[EnBand.MID][Quadrant.B][300] = 200
    data[EnBand.MID][Quadrant.C][300] = 200
    data[EnBand.MID][Quadrant.D][300] = 200

    # HIGH band: only 1 quadrant (1 channel)
    data[EnBand.HIGH][Quadrant.B][300] = 200
    # HIGH/C and HIGH/D remain at background level

    # Total: 4 channels triggered, need 5
    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    assert hits == [], f"Expected no hits with only 4 channels, got {len(hits)} hits"


def test_search_data_5_channels_detected():
    """Less obvious transient in 5 quadrant/band combinations is detected."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))

    # Insert moderate transient (enough to trigger at 5 sigma) in 5 channels
    # Use a smaller excess than previous test, but over 5 channels
    transient_idx = 300
    transient_value = 100  # 10x background, should trigger at 5 sigma

    # MID band: all 3 quadrants
    data[EnBand.MID][Quadrant.B][transient_idx] = transient_value
    data[EnBand.MID][Quadrant.C][transient_idx] = transient_value
    data[EnBand.MID][Quadrant.D][transient_idx] = transient_value

    # HIGH band: 2 quadrants (total: 5 channels)
    data[EnBand.HIGH][Quadrant.B][transient_idx] = transient_value
    data[EnBand.HIGH][Quadrant.C][transient_idx] = transient_value
    # HIGH/D remains at background

    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0, "Expected hits with 5 channels triggering"

    # Verify the hit is at or near the transient time
    hit_times = [t for t, _ in hits]
    assert any(abs(t - transient_idx) <= 8 for t in hit_times), \
        f"Expected hit near index {transient_idx}, got hits at {hit_times}"


def test_search_data_all_6_channels_detected():
    """Transient in all 6 quadrant/band combinations is detected."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))

    transient_idx = 250
    transient_value = 120

    # All 6 channels trigger
    data[EnBand.MID][Quadrant.B][transient_idx] = transient_value
    data[EnBand.MID][Quadrant.C][transient_idx] = transient_value
    data[EnBand.MID][Quadrant.D][transient_idx] = transient_value
    data[EnBand.HIGH][Quadrant.B][transient_idx] = transient_value
    data[EnBand.HIGH][Quadrant.C][transient_idx] = transient_value
    data[EnBand.HIGH][Quadrant.D][transient_idx] = transient_value

    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0, "Expected hits with all 6 channels triggering"

    hit_times = [t for t, _ in hits]
    assert any(abs(t - transient_idx) <= 8 for t in hit_times), \
        f"Expected hit near index {transient_idx}, got hits at {hit_times}"


def test_search_data_quiet_in_all_channels():
    """Quiet data in all channels returns no hits."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))
    # No transients, just constant background
    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    assert hits == []


def test_search_data_5_channels_different_window_sizes():
    """Transient detected with 5 channels, verify window size reported."""
    data = tuple(tuple([10] * 500 for _ in range(3)) for _ in range(3))

    # Create a multi-bin transient to test window sizing
    transient_start = 300
    transient_duration = 4  # 4 bins
    transient_value = 80

    for i in range(transient_start, transient_start + transient_duration):
        data[EnBand.MID][Quadrant.B][i] = transient_value
        data[EnBand.MID][Quadrant.C][i] = transient_value
        data[EnBand.MID][Quadrant.D][i] = transient_value
        data[EnBand.HIGH][Quadrant.B][i] = transient_value
        data[EnBand.HIGH][Quadrant.C][i] = transient_value

    hits = search_data(data, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0, "Expected hits with extended transient"

    # Should have hits with window sizes that match transient duration
    window_sizes = [h for _, h in hits]
    assert any(w >= transient_duration for w in window_sizes), \
        f"Expected window size >= {transient_duration}, got {window_sizes}"


# --- search_filepath tests ---

TEST_DATA_DIR = Path(__file__).parent / "data"


def test_search_filepath_file_not_found():
    """Non-existent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        search_filepath(
            Path("/nonexistent/path/file.raw"), size=20, maxtest=8, threshold=5.0
        )


def test_search_filepath_invalid_sra(tmp_path):
    """Invalid SRA file raises InvalidSRA."""
    # Create a file that's too small
    bad_file = tmp_path / "bad.raw"
    bad_file.write_bytes(b"short")

    with pytest.raises(InvalidSRA):
        search_filepath(bad_file, size=20, maxtest=8, threshold=5.0)


def test_search_filepath_with_real_data():
    """Search on real SRA file completes without error."""
    # Use first available test file
    test_file = TEST_DATA_DIR / "20251014_160753" / "SRA" / "561937C6.raw"
    if not test_file.exists():
        pytest.skip("Test data file not available")

    # Should complete without raising
    hits = search_filepath(test_file, size=210, maxtest=8, threshold=5.0)
    assert isinstance(hits, list)


def test_search_filepath_accepts_string_path():
    """search_filepath accepts string path as well as Path."""
    test_file = TEST_DATA_DIR / "20251014_160753" / "SRA" / "561937C6.raw"
    if not test_file.exists():
        pytest.skip("Test data file not available")

    hits = search_filepath(str(test_file), size=210, maxtest=8, threshold=5.0)
    assert isinstance(hits, list)


if __name__ == "__main__":
    test_search_qbdata_empty_data()
    test_search_qbdata_all_zeros()
    test_search_qbdata_data_too_short()
    test_search_qbdata_quiet_data()
    test_search_qbdata_with_transient()
    test_search_data_different_transient_times()
    test_search_data_only_4_channels_insufficient()
    test_search_data_5_channels_detected()
    test_search_data_all_6_channels_detected()
    test_search_data_quiet_in_all_channels()
    test_search_data_5_channels_different_window_sizes()

    # Skip filepath tests in direct run (need pytest fixtures)
    print("All basic tests passed!")
