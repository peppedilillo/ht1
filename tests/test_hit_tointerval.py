"""
Tests for hit_tointerval function.
"""

from fm1trig.search import hit_tointerval


def test_single_bin_hit():
    """Window length 1 returns interval containing just that bin."""
    hit = (10, 1)  # time_idx=10, window_len=1
    start, end = hit_tointerval(hit)
    assert (start, end) == (10, 11)


def test_multi_bin_hit():
    """Window length 4 returns interval of 4 bins ending at time_idx."""
    hit = (15, 4)  # time_idx=15, window_len=4
    start, end = hit_tointerval(hit)
    assert (start, end) == (12, 16)


def test_slice_captures_transient():
    """Slicing data with interval captures the transient bins."""
    # Data with transient at indices 10, 11, 12, 13
    xs = [5] * 20
    xs[10:14] = [100, 100, 100, 100]

    # Hit at end of transient with window covering it
    hit = (13, 4)
    start, end = hit_tointerval(hit)

    sliced = xs[start:end]
    assert sliced == [100, 100, 100, 100]


def test_slice_captures_spike():
    """Slicing with single-bin window captures the spike."""
    xs = [10] * 50
    xs[25] = 999

    hit = (25, 1)
    start, end = hit_tointerval(hit)

    assert xs[start:end] == [999]


if __name__ == "__main__":
    test_single_bin_hit()
    test_multi_bin_hit()
    test_slice_captures_transient()
    test_slice_captures_spike()
    print("All tests passed!")
