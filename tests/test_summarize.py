"""
Tests for summarize function.
"""

from ht1.ht1 import summarize, _INVALID_INTERVAL


def test_empty_hits():
    """Empty hits returns zero count and invalid interval."""
    nhits, interval = summarize([])
    assert nhits == 0
    assert interval == _INVALID_INTERVAL


def test_single_hit():
    """Single hit returns count 1 and its interval."""
    hits = [(10, 2)]  # time_idx=10, window_len=2 -> interval (9, 11)
    nhits, (start, end) = summarize(hits)
    assert nhits == 1
    assert (start, end) == (9, 11)


def test_multiple_hits_same_location():
    """Multiple hits at same time with different windows."""
    hits = [(10, 1), (10, 2), (10, 4)]
    nhits, (start, end) = summarize(hits)
    assert nhits == 3
    # Intervals: (10,11), (9,11), (7,11) -> union (7, 11)
    assert start == 7
    assert end == 11


def test_multiple_hits_span():
    """Multiple hits spanning a range."""
    hits = [(20, 2), (30, 4)]
    nhits, (start, end) = summarize(hits)
    assert nhits == 2
    # Intervals: (19,21), (27,31) -> span (19, 31)
    assert start == 19
    assert end == 31


if __name__ == "__main__":
    test_empty_hits()
    test_single_hit()
    test_multiple_hits_same_location()
    test_multiple_hits_span()
    print("All tests passed!")
