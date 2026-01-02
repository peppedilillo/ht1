"""
Tests for search_qbdata, search_data, and search_filepath functions.
"""

from pathlib import Path
import pytest

from ht1.ht1 import (
    search_qbdata,
    search_data,
    search_filepath,
    InvalidSRA,
    Quadrant,
    EnBand,
)


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
    data = tuple(
        tuple([10] * npoints for _ in range(3))
        for _ in range(3)
    )
    if transient_idx is not None:
        # Insert transient in all bands and quadrants
        for band in range(3):
            for quad in range(3):
                data[band][quad][transient_idx] = transient_value
    return data


def test_search_data_empty_checks():
    """Empty checks list returns empty hits."""
    data = _make_test_data(500, transient_idx=300)
    hits = search_data(data, checks=[], size=20, maxtest=8, threshold=5.0)
    assert hits == []


def test_search_data_single_check():
    """Single check returns hits from that band/quadrant."""
    data = _make_test_data(500, transient_idx=300)
    checks = [(EnBand.MID, Quadrant.B)]
    hits = search_data(data, checks, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0


def test_search_data_coincidence_all_match():
    """Hits present in all checks are returned."""
    data = _make_test_data(500, transient_idx=300)
    checks = [
        (EnBand.MID, Quadrant.B),
        (EnBand.MID, Quadrant.C),
        (EnBand.HIGH, Quadrant.B),
    ]
    hits = search_data(data, checks, size=20, maxtest=8, threshold=5.0)
    assert len(hits) > 0


def test_search_data_coincidence_partial_match():
    """Hits only in some checks are filtered out."""
    # Create data with transient only in one quadrant
    data = tuple(
        tuple([10] * 500 for _ in range(3))
        for _ in range(3)
    )
    # Insert transient only in quadrant B, band MID
    data[EnBand.MID][Quadrant.B][300] = 500

    checks = [
        (EnBand.MID, Quadrant.B),  # has transient
        (EnBand.MID, Quadrant.C),  # no transient
    ]
    hits = search_data(data, checks, size=20, maxtest=8, threshold=5.0)
    # No hits because transient not in all checks
    assert hits == []


def test_search_data_quiet_all_checks():
    """Quiet data in all checks returns no hits."""
    data = _make_test_data(500, transient_idx=None)
    checks = [
        (EnBand.MID, Quadrant.B),
        (EnBand.MID, Quadrant.C),
        (EnBand.MID, Quadrant.D),
    ]
    hits = search_data(data, checks, size=20, maxtest=8, threshold=5.0)
    assert hits == []


def test_search_data_different_transient_times():
    """Transients at different times in different checks don't coincide."""
    data = tuple(
        tuple([10] * 500 for _ in range(3))
        for _ in range(3)
    )
    # Transient at different times in different quadrants
    data[EnBand.MID][Quadrant.B][300] = 500
    data[EnBand.MID][Quadrant.C][350] = 500  # different time

    checks = [
        (EnBand.MID, Quadrant.B),
        (EnBand.MID, Quadrant.C),
    ]
    hits = search_data(data, checks, size=20, maxtest=8, threshold=5.0)
    # No coincident hits
    assert hits == []


# --- search_filepath tests ---

TEST_DATA_DIR = Path(__file__).parent / "data"


def test_search_filepath_file_not_found():
    """Non-existent file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        search_filepath(
            Path("/nonexistent/path/file.raw"),
            checks=[(EnBand.MID, Quadrant.B)],
            size=20,
            maxtest=8,
            threshold=5.0
        )


def test_search_filepath_invalid_sra(tmp_path):
    """Invalid SRA file raises InvalidSRA."""
    # Create a file that's too small
    bad_file = tmp_path / "bad.raw"
    bad_file.write_bytes(b"short")

    with pytest.raises(InvalidSRA):
        search_filepath(
            bad_file,
            checks=[(EnBand.MID, Quadrant.B)],
            size=20,
            maxtest=8,
            threshold=5.0
        )


def test_search_filepath_with_real_data():
    """Search on real SRA file completes without error."""
    # Use first available test file
    test_file = TEST_DATA_DIR / "20251014_160753" / "SRA" / "561937C6.raw"
    if not test_file.exists():
        pytest.skip("Test data file not available")

    checks = [
        (EnBand.MID, Quadrant.B),
        (EnBand.MID, Quadrant.C),
        (EnBand.MID, Quadrant.D),
    ]
    # Should complete without raising
    hits = search_filepath(test_file, checks, size=210, maxtest=8, threshold=5.0)
    assert isinstance(hits, list)


def test_search_filepath_accepts_string_path():
    """search_filepath accepts string path as well as Path."""
    test_file = TEST_DATA_DIR / "20251014_160753" / "SRA" / "561937C6.raw"
    if not test_file.exists():
        pytest.skip("Test data file not available")

    checks = [(EnBand.MID, Quadrant.B)]
    # Pass as string
    hits = search_filepath(str(test_file), checks, size=210, maxtest=8, threshold=5.0)
    assert isinstance(hits, list)


if __name__ == "__main__":
    test_search_qbdata_empty_data()
    test_search_qbdata_all_zeros()
    test_search_qbdata_data_too_short()
    test_search_qbdata_quiet_data()
    test_search_qbdata_with_transient()

    test_search_data_empty_checks()
    test_search_data_single_check()
    test_search_data_coincidence_all_match()
    test_search_data_coincidence_partial_match()
    test_search_data_quiet_all_checks()
    test_search_data_different_transient_times()

    # Skip filepath tests in direct run (need pytest fixtures)
    print("All basic tests passed!")
