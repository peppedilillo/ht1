"""
Tests for TriggerDyadic class.
"""

from math import sqrt

from ht1.ht1 import TriggerDyadic


def test_threshold_conversion():
    """Threshold sigma is converted to half-squared LLR."""
    t = TriggerDyadic(maxtest=4, threshold=5.0)
    assert t.llr_threshold_halfsq == 0.5 * (5.0**2)


def test_no_trigger_on_quiet_data():
    """Constant counts equal to background should produce no triggers."""
    t = TriggerDyadic(maxtest=4, threshold=3.0)

    # 20 steps of quiet data: counts = background = 10
    xs = [10] * 20
    bs = [10.0] * 20
    vrange = (0, 20)

    hits = t(xs, bs, vrange)
    assert hits == []


def test_trigger_on_single_data():
    """Constant counts equal to background should produce no triggers."""
    t = TriggerDyadic(maxtest=4, threshold=3.0)

    # 20 steps of quiet data: counts = background = 10
    xs = [100]
    bs = [1.0]
    vrange = (0, 1)

    hits = t(xs, bs, vrange)
    assert 0 in [i for i, _ in hits]

def test_obvious_transient_single_bin():
    """A large single-bin spike should trigger."""
    t = TriggerDyadic(maxtest=4, threshold=5.0)

    # Background: 10 counts per bin
    # Spike: 200 counts at index 10
    xs = [10] * 20
    xs[10] = 200
    bs = [10.0] * 20
    vrange = (0, 20)

    hits = t(xs, bs, vrange)

    assert len(hits) > 0
    # we get atn hits at the time of the transient
    assert 10 in [i for i, _ in hits]
    # All hits should be at or near the spike location
    for time_idx, window_len in hits:
        # The spike is at index 10, hit should include it
        interval_start = time_idx - window_len + 1
        assert interval_start <= 10 <= time_idx


def test_threshold():
    """Check that we trigger over a bin just a bit above the threshold."""
    t = TriggerDyadic(maxtest=4, threshold=5.0)

    # Background: 10 counts per bin
    # Spike: 200 counts at index 10
    b = 1e6
    xs = [b] * 20
    xs[10] = b + sqrt(b) * 5.01
    bs = [b] * 20
    vrange = (0, 20)

    hits = t(xs, bs, vrange)

    # check we got the trigger
    assert len(hits) > 0
    assert 10 in [i for i, h in hits]


def test_obvious_transient_multi_bin():
    """A sustained transient over multiple bins should trigger with larger windows."""
    t = TriggerDyadic(maxtest=8, threshold=5.0)

    # Background: 10 counts per bin
    # Sustained transient: 100 counts at indices 12-15 (4 bins)
    xs = [10] * 30
    for i in range(12, 16):
        xs[i] = 100
    bs = [10.0] * 30
    vrange = (0, 30)

    hits = t(xs, bs, vrange)

    assert len(hits) > 0
    # Check that some hits cover the transient region
    transient_hits = [(t_idx, w_len) for t_idx, w_len in hits if t_idx - w_len + 1 <= 15 and t_idx >= 12]
    assert len(transient_hits) > 0


def test_hit_time_indices_in_range():
    """All hit time indices should be within the valid range."""
    t = TriggerDyadic(maxtest=4, threshold=5.0)

    xs = [10] * 20
    xs[12] = 200
    bs = [10.0] * 20
    vrange = (5, 18)

    hits = t(xs, bs, vrange)

    for time_idx, _ in hits:
        assert vrange[0] <= time_idx < vrange[1]


def test_accumulator_correctness():
    """Verify accumulators sum correctly."""
    t = TriggerDyadic(maxtest=4, threshold=0.0)

    xs = [1, 2, 3, 4, 5]
    bs = [0.5, 1.0, 1.5, 2.0, 2.5]

    for x, b in zip(xs, bs):
        t.step(x, b)

    assert t.acc_x == sum(xs)
    assert t.acc_b == sum(bs)


def test_no_false_positives_on_fluctuations():
    """Small fluctuations around background should not trigger at high threshold."""
    t = TriggerDyadic(maxtest=8, threshold=5.0)

    # Background 50, small fluctuations +/- 10
    xs = [50, 55, 45, 52, 48, 53, 47, 51, 49, 54, 46, 52, 48, 55, 45, 50, 51, 49, 53, 47, 52, 48, 54, 46]
    bs = [50.0] * 24
    vrange = (0, 24)

    hits = t(xs, bs, vrange)
    assert hits == []


def test_marginal_excess_below_threshold():
    """An excess just below threshold should not trigger."""
    # threshold=5 means llr_threshold = 12.5
    # For a very mild excess that doesn't cross threshold
    t = TriggerDyadic(maxtest=4, threshold=5.0)

    # Very mild excess: 12 vs background 10
    xs = [10] * 20
    xs[10] = 12
    bs = [10.0] * 20
    vrange = (0, 20)

    hits = t(xs, bs, vrange)
    assert hits == []


def test_call_interface():
    """Test the __call__ interface returns expected structure."""
    t = TriggerDyadic(maxtest=4, threshold=5.0)

    xs = [10] * 15
    xs[10] = 200
    bs = [10.0] * 15
    vrange = (0, 15)

    hits = t(xs, bs, vrange)

    # Hits should be list of tuples (time_index, window_length)
    assert isinstance(hits, list)
    for hit in hits:
        assert isinstance(hit, tuple)
        assert len(hit) == 2
        assert isinstance(hit[0], int)
        assert isinstance(hit[1], int)


if __name__ == "__main__":
    test_threshold_conversion()
    test_no_trigger_on_quiet_data()
    test_obvious_transient_single_bin()
    test_obvious_transient_multi_bin()
    test_hit_time_indices_in_range()
    test_accumulator_correctness()
    test_no_false_positives_on_fluctuations()
    test_marginal_excess_below_threshold()
    test_call_interface()
    print("All tests passed!")
