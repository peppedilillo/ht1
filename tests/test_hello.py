#!/usr/bin/env python3
"""
Test script for hello_test function.
"""

from fm1trig.hello import hello_test


def test_hello():
    """Test that hello_test returns True."""
    result = hello_test()
    assert result is True, f"Expected True, got {result}"
    print("PASS: hello_test returned True")


if __name__ == "__main__":
    test_hello()
    print("\nAll tests passed!")
