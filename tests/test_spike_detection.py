"""
Tests for spike detection and fixing functions.
"""

from fm1trig.search import _find_spikes, _fix_spikes, NQUADRANTS


def test_find_single_quadrant_spike():
    """Spike in a single quadrant should be detected."""
    band_data = (
        [0, 0, 100, 0, 0],  # quad 0 (B): spike at index 2
        [0, 0, 0, 0, 0],    # quad 1 (C): no spike
        [0, 0, 0, 0, 0],    # quad 2 (D): no spike
    )

    spikes = _find_spikes(band_data, thr=5.0)

    assert 2 in spikes, "Spike at index 2 should be detected"
    assert spikes[2] == [0], "Spike should be in quadrant 0 only"


def test_find_multi_quadrant_spike():
    """Spike in 2 quadrants should be detected."""
    band_data = (
        [0, 0, 100, 0, 0],  # quad 0 (B): spike at index 2
        [0, 0, 100, 0, 0],  # quad 1 (C): spike at index 2
        [0, 0, 0, 0, 0],    # quad 2 (D): no spike
    )

    spikes = _find_spikes(band_data, thr=5.0)

    assert 2 in spikes, "Spike at index 2 should be detected"
    assert set(spikes[2]) == {0, 1}, "Spike should be in quadrants 0 and 1"


def test_find_all_quadrant_spike_not_detected():
    """Spike in all 3 quadrants should NOT be detected (real event)."""
    band_data = (
        [0, 0, 100, 0, 0],  # all quadrants have spike
        [0, 0, 100, 0, 0],
        [0, 0, 100, 0, 0],
    )

    spikes = _find_spikes(band_data, thr=5.0)

    assert 2 not in spikes, "All-quadrant spike should not be detected (real event)"
    assert len(spikes) == 0, "No spikes should be returned"


def test_find_spike_at_first_index_not_detected():
    """Spike at index 0 is not detected (boundary)."""
    band_data = (
        [100, 0, 0, 0, 0],  # spike at first index
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    )

    spikes = _find_spikes(band_data, thr=5.0)

    assert 0 not in spikes, "Spike at index 0 is not detected"
    assert len(spikes) == 0, "No spikes should be returned"


def test_find_spike_at_last_index_not_detected():
    """Spike at last index is not detected (boundary)."""
    band_data = (
        [0, 0, 0, 0, 100],  # spike at last index
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    )

    spikes = _find_spikes(band_data, thr=5.0)

    assert 4 not in spikes, "Spike at last index is not detected"
    assert len(spikes) == 0, "No spikes should be returned"


def test_find_multiple_spike():
    """Multiple spikes should all be detected."""
    band_data = (
        [0, 0, 100, 0, 50, 0, 0, 0],  # spikes at indices 2 and 4
        [0, 0, 0, 0, 0, 0, 10, 0],
        [0, 0, 0, 0, 0, 0, 0, 0],
    )

    spikes = _find_spikes(band_data, thr=3.0)

    assert 2 in spikes, "Spike at index 2 should be detected"
    assert 4 in spikes, "Spike at index 4 should be detected"
    assert 6 in spikes, "Spike at index 6 should be detected"
    assert spikes[2] == [0], "Spike at index 2 should be in quadrant 0"
    assert spikes[4] == [0], "Spike at index 4 should be in quadrant 0"
    assert spikes[6] == [1], "Spike at index 6 should be in quadrant 1"


def test_fix_single_spike():
    """Fixing a spike should replace it with average of neighbors."""
    band_data = (
        [1, 2, 100, 4, 5],  # spike at index 2
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    )
    spikes = {2: [0]}

    _fix_spikes(band_data, spikes)

    assert band_data[0][2] == 3.0, "Spike should be replaced with average of neighbors"
    assert band_data[0][1] == 2, "Neighboring values should not change"
    assert band_data[0][3] == 4, "Neighboring values should not change"


def test_fix_only_affects_specified_quadrants():
    """Fixing should only modify quadrants listed in the spike dictionary."""
    band_data = (
        [0, 0, 100, 0, 0],  # quad 0: has spike, will be fixed
        [0, 0, 100, 0, 0],  # quad 1: has spike, but not in fix list
        [0, 0, 0, 0, 0],    # quad 2: no spike
    )
    spikes = {2: [0]}  # Only fix quadrant 0

    _fix_spikes(band_data, spikes)

    assert band_data[0][2] == 0.0, "Quadrant 0 should be fixed"
    assert band_data[1][2] == 100, "Quadrant 1 should NOT be fixed"


def test_fix_multiple_spikes():
    """Fixing multiple spikes should work correctly."""
    band_data = (
        [0, 0, 100, 0, 50, 0, 0],  # spikes at indices 2 and 4
        [0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0],
    )
    spikes = {2: [0], 4: [0]}

    _fix_spikes(band_data, spikes)

    assert band_data[0][2] == 0.0, "Spike at index 2 should be fixed"
    assert band_data[0][4] == 0.0, "Spike at index 4 should be fixed"


def test_fix_spike_across_multiple_quadrants():
    """Fixing spike detected in multiple quadrants should fix all of them."""
    band_data = (
        [0, 0, 100, 0, 0],  # quad 0 (B): spike at index 2
        [0, 0, 100, 0, 0],  # quad 1 (C): spike at index 2
        [0, 0, 0, 0, 0],    # quad 2 (D): no spike
    )
    spikes = {2: [0, 1]}  # Spike in quads 0 and 1

    _fix_spikes(band_data, spikes)

    assert band_data[0][2] == 0.0, "Quadrant 0 should be fixed"
    assert band_data[1][2] == 0.0, "Quadrant 1 should be fixed"
    assert band_data[2][2] == 0, "Quadrant 2 should remain unchanged"


if __name__ == "__main__":
    test_find_single_quadrant_spike()
    test_find_multi_quadrant_spike()
    test_find_all_quadrant_spike_not_detected()
    test_find_spike_at_first_index_not_detected()
    test_find_spike_at_last_index_not_detected()
    test_fix_single_spike()
    test_fix_only_affects_specified_quadrants()
    test_fix_multiple_spikes()
    test_fix_spike_across_multiple_quadrants()

    print("\nAll tests passed!")
