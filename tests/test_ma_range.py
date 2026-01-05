"""
Tests for data_valid_range function.
"""

from ht1.ht1 import _INVALID_INTERVAL
from ht1.ht1 import moving_average_range


def test_empty_data():
    """Empty data should return invalid range."""
    data = []
    size = 1
    expected = _INVALID_INTERVAL
    assert moving_average_range(data, size) == expected


def test_all_zeros():
    """Data with all zeros should return invalid range."""
    data = [0, 0, 0, 0, 0]
    size = 1
    expected = _INVALID_INTERVAL
    assert moving_average_range(data, size) == expected


def test_single_nonzero_too_short():
    """Single non-zero element is too short for any size > 0."""
    data = [0, 0, 1, 0, 0]
    size = 1
    expected = _INVALID_INTERVAL
    assert moving_average_range(data, size) == expected


def test_data_too_short_for_size():
    """Non-zero span too short for the given window size."""
    data = [0, 1, 2, 0]
    size = 1
    # first non-zero: i=1, last non-zero: j=2
    # after padding: i=2, j=2, so j-i=0 < 1
    expected = _INVALID_INTERVAL
    assert moving_average_range(data, size) == expected


def test_minimum_valid_length():
    """Minimum valid case: exactly enough data for one output point."""
    data = [0, 1, 2, 3, 0]
    size = 1
    # first non-zero: i=1, last non-zero: j=3
    # after padding: i=2, j=3, so j-i=1 >= 1
    expected = (2, 3)
    assert moving_average_range(data, size) == expected


def test_leading_zeros():
    """Data with leading zeros only."""
    data = [0, 0, 0, 1, 2, 3, 4, 5]
    size = 1
    # first non-zero: i=3, last non-zero: j=7
    # after padding: i=4, j=7
    expected = (4, 7)
    assert moving_average_range(data, size) == expected


def test_trailing_zeros():
    """Data with trailing zeros only."""
    data = [1, 2, 3, 4, 5, 0, 0, 0]
    size = 1
    # first non-zero: i=0, last non-zero: j=4
    # after padding: i=1, j=4
    expected = (1, 4)
    assert moving_average_range(data, size) == expected


def test_leading_and_trailing_zeros():
    """Data with both leading and trailing zeros."""
    data = [0, 0, 1, 2, 3, 4, 5, 0, 0]
    size = 1
    # first non-zero: i=2, last non-zero: j=6
    # after padding: i=3, j=6
    expected = (3, 6)
    assert moving_average_range(data, size) == expected


def test_no_zeros():
    """Data with no zeros at all."""
    data = [1, 2, 3, 4, 5]
    size = 1
    # first non-zero: i=0, last non-zero: j=4
    # after padding: i=1, j=4
    expected = (1, 4)
    assert moving_average_range(data, size) == expected


def test_larger_size():
    """Test with larger window size."""
    data = [0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 0]
    size = 3
    # first non-zero: i=2, last non-zero: j=11
    # after padding: i=5, j=9
    expected = (5, 9)
    assert moving_average_range(data, size) == expected


def test_larger_size_too_short():
    """Larger window size makes previously valid data invalid."""
    data = [0, 1, 2, 3, 4, 5, 0]
    size = 2
    # first non-zero: i=1, last non-zero: j=5
    # after padding: i=3, j=4, so j-i=1 >= 1
    expected = (3, 4)
    assert moving_average_range(data, size) == expected

    size = 3
    # after padding: i=4, j=3, so j-i=-1 < 1
    expected = _INVALID_INTERVAL
    assert moving_average_range(data, size) == expected


def test_size_zero():
    """Edge case: size=0 means no padding needed."""
    data = [0, 1, 0]
    size = 0
    # first non-zero: i=1, last non-zero: j=1
    # after padding: i=1, j=2, so j-i=1 >= 1
    expected = (1, 2)
    assert moving_average_range(data, size) == expected


def test_zeros_in_middle():
    """Zeros in middle of data don't affect range calculation."""
    data = [0, 1, 0, 0, 0, 2, 0]
    size = 1
    # first non-zero: i=1, last non-zero: j=5
    # after padding: i=2, j=5
    expected = (2, 5)
    assert moving_average_range(data, size) == expected


def test_range_iteration():
    """Verify returned range can be used correctly with range()."""
    data = [0, 0, 5, 6, 7, 8, 9, 0, 0]
    size = 1
    # first non-zero: i=2, last non-zero: j=6
    # after padding: i=3, j=6
    result = moving_average_range(data, size)
    assert result == (3, 6)
    # the range should yield 3 indices: 3, 4, 5
    assert list(range(*result)) == [3, 4, 5]


if __name__ == "__main__":
    test_empty_data()
    test_all_zeros()
    test_single_nonzero_too_short()
    test_data_too_short_for_size()
    test_minimum_valid_length()
    test_leading_zeros()
    test_trailing_zeros()
    test_leading_and_trailing_zeros()
    test_no_zeros()
    test_larger_size()
    test_larger_size_too_short()
    test_size_zero()
    test_zeros_in_middle()
    test_range_iteration()
    print("All tests passed!")
