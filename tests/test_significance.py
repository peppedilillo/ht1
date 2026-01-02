"""
Tests for significance function.
"""

from math import log, sqrt
from ht1.ht1 import significance


def test_x_equal_b():
    """When x equals b, significance should be 0."""
    assert significance(5.0, 5.0) == 0.0
    assert significance(100.0, 100.0) == 0.0


def test_x_less_than_b():
    """When x < b (deficit), significance should be 0."""
    assert significance(3.0, 5.0) == 0.0
    assert significance(0.0, 5.0) == 0.0


def test_b_zero():
    """When b is zero, significance should be 0."""
    assert significance(5.0, 0.0) == 0.0
    assert significance(0.0, 0.0) == 0.0


def test_b_negative():
    """When b is negative, significance should be 0. Note this is allowed."""
    assert significance(5.0, -1.0) == 0.0
    assert significance(-1.0, -2.0) == 0.0


def test_x_negative_b_positive():
    """When x is negative and b is positive, significance should be 0. Note this is allowed."""
    assert significance(-1.0, 5.0) == 0.0


def test_excess_returns_positive():
    """When x > b > 0, significance should be positive."""
    result = significance(10.0, 5.0)
    assert result > 0.0


def test_formula_correctness():
    """Verify the LLR formula: x * log(x/b) - (x - b)."""
    x, b = 10.0, 5.0
    expected = x * log(x / b) - (x - b)
    assert significance(x, b) == expected

def test_formula_approximation():
    b = 1e+6
    assert significance(b + sqrt(b) * 5.01, b) > .5 * 5. ** 2

def test_fractional_values():
    """Test with non-integer values. Note this is allowed."""
    x, b = 7.5, 2.5
    result = significance(x, b)
    expected = x * log(x / b) - (x - b)
    assert abs(result - expected) < 1e-10


if __name__ == "__main__":
    test_x_equal_b()
    test_x_less_than_b()
    test_b_zero()
    test_b_negative()
    test_x_negative_b_positive()
    test_excess_returns_positive()
    test_formula_correctness()

    test_fractional_values()
    print("All tests passed!")
