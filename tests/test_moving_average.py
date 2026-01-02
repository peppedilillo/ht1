"""
Tests for moving_average function.
"""

from ht1.ht1 import moving_average, _INVALID_INTERVAL
import pytest


def test_empty_data():
    """Empty data should return empty sma and invalid range."""
    data = []
    size = 1
    sma, vrange = moving_average(data, size)
    assert sma == []
    assert vrange == _INVALID_INTERVAL


def test_all_zeros():
    """Data with all zeros should return empty sma and invalid range."""
    data = [0, 0, 0, 0, 0]
    size = 1
    sma, vrange = moving_average(data, size)
    assert sma == []
    assert vrange == _INVALID_INTERVAL


def test_data_too_short():
    """Data too short for window size returns empty sma and invalid range."""
    data = [0, 1, 2, 0]
    size = 1
    sma, vrange = moving_average(data, size)
    assert sma == []
    assert vrange == _INVALID_INTERVAL


def test_simple_sequence():
    """Simple ascending sequence with size=1 (window=3)."""
    data = [1, 2, 3, 4, 5]
    size = 1
    # vrange: first=0, last=4, padded: i=1, j=4
    # window_len = 3
    # i=1: avg(1,2,3) = 2.0
    # i=2: avg(2,3,4) = 3.0
    # i=3: avg(3,4,5) = 4.0
    expected_sma = [2.0, 3.0, 4.0]
    expected_vrange = (1, 4)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


def test_constant_values():
    """Constant values should yield constant moving average."""
    data = [5, 5, 5, 5, 5]
    size = 1
    # All averages should be 5.0
    expected_sma = [5.0, 5.0, 5.0]
    expected_vrange = (1, 4)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


def test_with_leading_trailing_zeros():
    """Data with leading and trailing zeros."""
    data = [0, 0, 1, 2, 3, 4, 5, 0, 0]
    size = 1
    # vrange: first=2, last=6, padded: i=3, j=6
    # i=3: avg(1,2,3) = 2.0
    # i=4: avg(2,3,4) = 3.0
    # i=5: avg(3,4,5) = 4.0
    expected_sma = [2.0, 3.0, 4.0]
    expected_vrange = (3, 6)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


def test_size_zero():
    """Size=0 means window of 1, so average equals the value itself."""
    data = [0, 3, 7, 2, 0]
    size = 0
    # vrange: first=1, last=3, padded: i=1, j=4
    # window_len = 1
    # i=1: avg(3) = 3.0
    # i=2: avg(7) = 7.0
    # i=3: avg(2) = 2.0
    expected_sma = [3.0, 7.0, 2.0]
    expected_vrange = (1, 4)
    with pytest.raises(ValueError) as exc_info:
        moving_average(data, size)
    assert "Parameter `size` should be a positive integer." in str(exc_info.value)


def test_larger_window():
    """Test with larger window size=2 (window=5)."""
    data = [1, 2, 3, 4, 5, 6, 7]
    size = 2
    # vrange: first=0, last=6, padded: i=2, j=5
    # window_len = 5
    # i=2: avg(1,2,3,4,5) = 15/5 = 3.0
    # i=3: avg(2,3,4,5,6) = 20/5 = 4.0
    # i=4: avg(3,4,5,6,7) = 25/5 = 5.0
    expected_sma = [3.0, 4.0, 5.0]
    expected_vrange = (2, 5)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


def test_sma_length_matches_range():
    """Verify sma length equals the number of indices in vrange."""
    data = [0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0]
    size = 2
    sma, vrange = moving_average(data, size)
    expected_length = vrange[1] - vrange[0]
    assert len(sma) == expected_length


def test_non_integer_averages():
    """Test case where averages are not whole numbers."""
    data = [1, 2, 4, 5, 7]
    size = 1
    # window_len = 3
    # i=1: avg(1,2,4) = 7/3 = 2.333...
    # i=2: avg(2,4,5) = 11/3 = 3.666...
    # i=3: avg(4,5,7) = 16/3 = 5.333...
    expected_vrange = (1, 4)
    sma, vrange = moving_average(data, size)
    assert vrange == expected_vrange
    assert len(sma) == 3
    assert abs(sma[0] - 7/3) < 1e-10
    assert abs(sma[1] - 11/3) < 1e-10
    assert abs(sma[2] - 16/3) < 1e-10


def test_zeros_in_middle():
    """Zeros in middle of data are included in averages."""
    data = [0, 3, 0, 0, 6, 0]
    size = 1
    # vrange: first=1, last=4, padded: i=2, j=4
    # i=2: avg(3,0,0) = 1.0
    # i=3: avg(0,0,6) = 2.0
    expected_sma = [1.0, 2.0]
    expected_vrange = (2, 4)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


def test_single_output_point():
    """Minimum valid case producing exactly one moving average value."""
    data = [0, 1, 2, 3, 0]
    size = 1
    # vrange: first=1, last=3, padded: i=2, j=3
    # i=2: avg(1,2,3) = 2.0
    expected_sma = [2.0]
    expected_vrange = (2, 3)
    sma, vrange = moving_average(data, size)
    assert sma == expected_sma
    assert vrange == expected_vrange


if __name__ == "__main__":
    test_empty_data()
    test_all_zeros()
    test_data_too_short()
    test_simple_sequence()
    test_constant_values()
    test_with_leading_trailing_zeros()
    test_size_zero()
    test_larger_window()
    test_sma_length_matches_range()
    test_non_integer_averages()
    test_zeros_in_middle()
    test_single_output_point()
    print("All tests passed!")
