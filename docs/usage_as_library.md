# Using HT1 as a library

You can install HT1 as a standard Python 3.6+ package for on-ground reproducibility. 
The process is straightforward: git clone this repository, move to your local clone directory and install with `pip install .`.

Once installed, you can run a transient search over a SRA file with custom parameters using `search_filepath`:

```python
from ht1 import search_filepath

# Runs search for triggers over all three quadrants (B, C, D) in MID and HIGH energy bands.
# The LOW energy band is excluded from the search.
# Requires at least 5 out of 6 quadrant/band combinations to trigger simultaneously.
# Parameters: moving average window size 42.1s (size=210), maximum test window size 1.6s (maxtest=16),
# detection threshold 5.0 sigma.
hits = search_filepath("/path/to/srafile.raw", size=210, maxtest=16, threshold=5.0)
```

The variable `hits` contains a list of trigger hits as simple 2-tuples representing (trigger time-index, trigger window length).
Only hits that triggered simultaneously in at least 5 of the 6 quadrant/band combinations are returned, guaranteeing all quadrants and at least two energy bands are involved.

Note that while we validate input when using `HT1` as a script, **minimal or no checks are performed at the individual function level**.
Hence, when using HT1 as a package, it is your responsibility to provide a rational choice for these parameters:
* `size` must be a positive integer. It controls the extension of the moving average window. Moving average is performed over a centered window with length `2 * size + 1`. Recommended value: 210 (42.1s window)
* `maxtest` must be a power of two integer. It controls the length of the largest count interval checked for statistically significant count excess relative to the background. Recommended value: 16 (0.1s, 0.2s, 0.4s, 0.8s, 1.6s windows)
* `threshold` must be positive float. It controls the threshold standard deviation over which an interval test results in a trigger hit. Recommended value: 5.

Other search functions exist, providing more granular control over the algorithm.
The function `search_data` is designed to run the algorithm over data parsed from a SRA file into HT1 format.
In HT1, data are stored as a 3-tuple, one per energy band. Each record of the tuple is itself a 3-tuple, one per quadrant data. Finally, counts for each quadrant-band combinations are stored as lists of integers.

Like `search_filepath`, `search_data` searches all three quadrants in MID and HIGH bands and requires at least 5 out of 6 simultaneous triggers.

```python
from ht1 import search_data, hit_tointerval, sra_parse

data, abts = sra_parse("/path/to/srafile.raw")
hits = search_data(data, size=210, maxtest=16, threshold=5.0)
for interval in map(hit_tointerval, hits):
    ...
```

A function `hit_tointerval` is provided for converting hits into proper time series interval ranges you can use for slicing.

For custom coincidence logic, you can use `search_qbdata` for running on individual quadrant-band combinations.
The following example shows how to search for triggers in the LOW-energy band simultaneously over all quadrants, a custom trigger condition:

```python
from collections import Counter
from ht1 import search_qbdata, sra_parse, EnBand, Quadrant

data, abts = sra_parse("/path/to/srafile.raw")
hits_counter = Counter()
for band, quadrant in (
    (EnBand.LOW, Quadrant.B),
    (EnBand.LOW, Quadrant.C),
    (EnBand.LOW, Quadrant.D),
):
    hits = search_qbdata(data[band][quadrant], size=210, maxtest=16, threshold=5.0)
    for ih in hits:
        hits_counter[ih] += 1
# Require all 3 LOW-band quadrants to trigger
hits = [ih for ih, count in hits_counter.items() if count == 3]
```

For complete control over the algorithm and background estimation, you can instantiate the `TriggerDyadic` class directly and call the `moving_average` function.

```python
from ht1 import TriggerDyadic, moving_average, sra_parse, EnBand, Quadrant

data, abts = sra_parse("/path/to/srafile.raw")
xs = data[EnBand.HIGH][Quadrant.C]
bs, vrange = moving_average(xs, size=210)
trigger = TriggerDyadic(maxtest=16, threshold=5.0)
hits = trigger(xs, bs, vrange)
```

The function `moving_average` returns a 2-tuple:
- The first element contains the moving average estimates (a list of floats)
- The second element is a range (start, end) for aligning the estimates with the data

Because the length of the data and moving average estimates may differ (due to edge effects), `moving_average` provides the valid range over which the estimates are computed.