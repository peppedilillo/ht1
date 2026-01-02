# Scripting with ht1

You can install `ht1` as a standard Python >=3.6 package for on-ground reproducibility. 
The process is straightforward: git clone this repository, move to your local clone directory and install with `pip install .`.

Once installed, you can run a transient search over a SRA file with custom parameters using `search_filepath`:

```python
from ht1 import search_filepath, EnBand, Quadrant

checks = (
    (EnBand.MID, Quadrant.B),
    (EnBand.MID, Quadrant.C),
    (EnBand.MID, Quadrant.D),
)
# runs search for trigger on mid-energy band data with moving average window size 42.0 s, 
# and a maximum test window size of 1.6s, at threshold 4.5 sigma. 
hits = search_filepath("/path/to/srafile", checks=checks, size=210, maxtest=16, threshold=4.5)
```

The variable `hits` contains a list of trigger hits as simple 2-tuple representing trigger time-index and trigger window length.
A function `hit_tointerval` is provided for converting hits into proper time series interval ranges you can use for slicing.
Data are stored as a 3-tuple, one per energy band. Each record of the tuple is itself a 3-tuple, one for quadrant data.
Finally, counts for each quadrant-band combinations are stored as lists of integers. 
Other search functions exists, providing more granular control over the algorithm. 

Note that while we validate input when using `ht1` as a script, *minimal or no checks are performed at the individual function level* when using ht1 as a package.
Hence it is your responsibility to provide a rational choice for these parameters:
* `size` must be a positive integer and controls the extension of the moving average window. Moving average is performed over a centered window with length `2 * size + 1`.
* `maxtest` must be a power of two integer. It controls the lenght of the count interval checked for statistically significant count excess relative to the background.
* `threshold` must be positive float. It controls the threshold standard deviation over which an interval test results in a trigger hit.

The function `search_data` for example is designed to run the algorithm over parsed from a SRA file into the format described above.

```python
from ht1 import search_data, hit_tointerval, sra_parse, EnBand, Quadrant

checks = (
    (EnBand.MID, Quadrant.B),
    (EnBand.MID, Quadrant.C),
)

data, abts = sra_parse("/path/to/srafile")
hits = search_data(data, checks=checks, size=210, maxtest=16, threshold=4.5)
for interval in map(hit_tointerval, hits):
    start, end = interval
    print(data[EnBand.MID][Quadrant.B][start: end])
    # for each trigger hit, prints the over-threshold data interval for quadrant B and high-energy band.
```

Finally, you can use `search_qbdata` for running on data from a particular band-energy combination:

```python
from ht1 import search_qbdata, sra_parse, EnBand, Quadrant

data, abts = sra_parse("/path/to/srafile")
hits = search_qbdata(data[EnBand.HIGH][Quadrant.B], size=210, maxtest=16, threshold=4.5)
```

For complete control over the algorithm and background estimation one can instance the `TriggerDyadic` class.

```python
from ht1 import TriggerDyadic, moving_average, sra_parse, EnBand, Quadrant

data, abts = sra_parse("/path/to/srafile")
xs = data[EnBand.HIGH][Quadrant.C]
bs = moving_average(xs, size=210)
trigger = TriggerDyadic(maxtest=16, threshold=5.)
# note the `*` in  `*bs`. `moving_average` returns a 2-tuple.
# the first field of the tuple contains the moving average estimates.
# since the length of xs and moving averages estimates differs because edges
# do not contain enough data to compute the moving average, `moving_average` 
# also returns an interval in the second tuple field for alignment.
hits = trigger(xs, *bs)
```