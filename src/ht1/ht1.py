"""
A program and module for searching transients in SpIRIT/HERMES-FM1 SRA data.

Author: Giuseppe Dilillo
Date: January 2026
"""

import argparse
from collections import Counter
from collections import deque
from enum import Enum
from enum import IntEnum
from itertools import accumulate
from math import log
from pathlib import Path
import struct
import sys
from typing import List, NamedTuple, Sequence, Tuple, Union


class Quadrant(IntEnum):
    """Quadrant identifier mapping names to indices."""

    B = 0
    C = 1
    D = 2


class EnBand(IntEnum):
    """Energy band identifier mapping names to indices."""

    LOW = 0
    MID = 1
    HIGH = 2


class ErrorCode(IntEnum):
    """Exit codes returned by the main script."""

    OK = 0
    INVALID_FILE = 1
    INVALID_PARAMETERS = 2


BandData = Tuple[List[int], ...]  # one list per quadrant
Data = Tuple[BandData, ...]  # one BandData per band
Hit = Tuple[int, int]  # trigger time-index, window length
Index = Tuple[EnBand, Quadrant]  # band-quadrant combination index
Interval = Tuple[int, int]


_NQUADRANTS: int = 3
_NBANDS: int = 3

_SEARCH_SIZE_DEFAULT = 210
_SEARCH_MAXTEST_DEFAULT = 8
_SEARCH_THRESHOLD_DEFAULT = 5.0
_SEARCH_CHECKS_DEFAULT = (
    (EnBand.MID, Quadrant.B),
    (EnBand.MID, Quadrant.C),
    (EnBand.MID, Quadrant.D),
    (EnBand.HIGH, Quadrant.B),
    (EnBand.HIGH, Quadrant.C),
    (EnBand.HIGH, Quadrant.D),
)

_INVALID_INTERVAL = (-1, -1)


class SRASpecs(NamedTuple):
    """SRA file format specification and header layout constants."""

    BLOCK_SIZE: int = 124
    HEADER_SIZE: int = 33
    FILELEN: Tuple[int, int] = (21, 25)
    FILELEN_FMT: str = "=I"
    ABT_SIZE: int = 4
    ABT_FMT: str = "=I"
    DATA_SIZE: int = 10
    DATA_FMT: str = "10B"


class InvalidSRA(Exception):
    """An error raised when parsing invalid SRA files."""


def sra_parse(filepath: Union[Path, str]) -> Tuple[Data, List[int]]:
    """Parse an SRA file and extract count data and timestamps.

    Data are organized as nested tuples: the outer level indexes energy band,
    the inner level indexes quadrant. For example, `data[2][0]` selects
    quadrant B (index 0) in the high-energy band (index 2).

    Args:
        filepath: Path to the SRA file.

    Returns:
        A tuple containing:
            - data: Nested tuple of count lists indexed by (band, quadrant).
            - abts: List of ABT timestamps, one per data block.

    Raises:
        InvalidSRA: If the file is malformed or header size exceeds actual size.
        FileNotFoundError: If the file does not exist.
    """
    sra_specs = SRASpecs()
    block_fmt = sra_specs.ABT_FMT + sra_specs.DATA_FMT * _NBANDS * (_NQUADRANTS + 1)
    block_struct = struct.Struct(block_fmt)

    fsize = Path(filepath).stat().st_size
    with open(filepath, "rb") as f:
        if fsize < sra_specs.HEADER_SIZE:
            raise InvalidSRA("File size smaller than header size.")
        header = f.read(sra_specs.HEADER_SIZE)
        hfsize, *_ = struct.unpack(sra_specs.FILELEN_FMT, header[slice(*sra_specs.FILELEN)])
        nblocks = (hfsize - sra_specs.HEADER_SIZE) // sra_specs.BLOCK_SIZE
        if nblocks * sra_specs.BLOCK_SIZE > fsize - sra_specs.HEADER_SIZE:
            raise InvalidSRA("Invalid header fsize.")
        abts = [0] * nblocks
        data = tuple(tuple([0] * (nblocks * sra_specs.DATA_SIZE) for _ in range(_NQUADRANTS)) for _ in range(_NBANDS))
        for ib in range(nblocks):
            values = block_struct.unpack(f.read(sra_specs.BLOCK_SIZE))
            abts[ib] = values[0]
            offset = 1  # skip ABT
            for band in range(_NBANDS):
                offset += sra_specs.DATA_SIZE  # skip quadrant A
                for quad in range(_NQUADRANTS):
                    start = ib * sra_specs.DATA_SIZE
                    end = start + sra_specs.DATA_SIZE
                    data[band][quad][start:end] = values[offset : offset + sra_specs.DATA_SIZE]
                    offset += sra_specs.DATA_SIZE
    return data, abts


def ma_range(xs: Sequence[int], size: int) -> Interval:
    """Compute the valid range for moving average calculation.

    Finds the first and last non-zero elements in the sequence and returns
    an interval padded inward by the window half-size to ensure all values
    in the range have complete windows.

    Args:
        xs: Input sequence of integer counts.
        size: Half-window size for the moving average.

    Returns:
        A tuple (start, end) defining the valid range, or (-1, -1) if
        the sequence is empty, all zeros, or too short for the window.
    """
    n = len(xs)
    if n == 0:
        return _INVALID_INTERVAL
    i = 0
    while i < n and xs[i] == 0:
        i += 1
    if i == n:
        return _INVALID_INTERVAL
    j = n - 1
    while xs[j] == 0:
        j -= 1
    i += size
    j -= size - 1
    # `j - i` is the number of moving average values you get
    if j - i < 1:
        return _INVALID_INTERVAL
    # interval closed at left, open to the right: iterate over range(i, j)
    return i, j


def moving_average(xs: Sequence[int], size: int) -> Tuple[List[float], Interval]:
    """Compute a centered simple moving average over the valid data range.

    Uses a window of length `2 * size + 1` centered on each point.

    Args:
        xs: Input sequence of integer counts.
        size: Half-window size (full window is `2 * size + 1`).

    Returns:
        A tuple containing:
            - sma: List of moving average values.
            - vrange: The (start, end) interval over which averages were computed.

    Raises:
        ValueError: If size is not a positive integer.
    """
    if size < 1:
        raise ValueError("Parameter `size` should be a positive integer.")

    vrange = ma_range(xs, size)
    cumsum = [0] + list(accumulate(xs))
    sma = []
    window_len = 2 * size + 1
    for i in range(*vrange):
        window_sum = cumsum[i + size + 1] - cumsum[i - size]
        sma.append(window_sum / window_len)
    return sma, vrange


class TriggerStatus(Enum):
    """Status of the trigger algorithm: ACQUIRING initial data or RUNNING detection."""

    ACQUIRING = 0
    RUNNING = 1


def significance(x: float, b: float) -> float:
    """Compute the half-squared Poisson log-likelihood ratio.

    Measures statistical significance of an observed count excess over
    expected background. Returns zero for deficits or invalid inputs.

    Args:
        x: Observed count.
        b: Expected background count.

    Returns:
        The half-squared log-likelihood ratio, or 0.0 if x <= b or b <= 0.
    """
    if x <= b or b <= 0:
        return 0.0
    return x * log(x / b) - (x - b)


class TriggerDyadic:
    """Dyadic windowed trigger algorithm for detecting excesses in Poisson count series.

    Searches for time segments where significance exceeds a threshold using
    power-of-2 window sizes (1, 2, 4, 8, ...).

    Attributes:
        maxtest: Maximum window size to test (must be power of 2).
        llr_threshold_halfsq: Detection threshold as half-squared LLR.
    """

    def __init__(self, maxtest: int, threshold: float):
        """Initialize the trigger algorithm.

        Args:
            maxtest: Maximum window size to test (must be power of 2).
            threshold: Detection threshold in sigma units (converted internally
                to half-squared log-likelihood ratio).
        """
        self.maxtest = maxtest
        self.llr_threshold_halfsq = 0.5 * (threshold**2)

        self.phase_counter = -1
        self.acc_x = 0.0
        self.acc_b = 0.0
        self.queue_x = deque([0.0], maxlen=maxtest + 1)
        self.queue_b = deque([0.0], maxlen=maxtest + 1)

    def __call__(self, xs: Sequence[int], bs: Sequence[float], vrange: Interval) -> List[Hit]:
        """Run the trigger algorithm on count data with background estimates.

        Args:
            xs: Observed count sequence.
            bs: Background estimate sequence (same length as valid range).
            vrange: The (start, end) interval to process.

        Returns:
            List of hits as (time_index, window_length) tuples.
        """
        hits = []
        for i, j in enumerate(range(*vrange)):
            hits.extend([(j, h) for h in self.step(xs[j], bs[i])])
        return hits

    def status(self) -> TriggerStatus:
        """Return the current trigger status.

        Returns:
            ACQUIRING if still collecting initial data, RUNNING otherwise.
        """
        return TriggerStatus.RUNNING if self.phase_counter >= 0 else TriggerStatus.ACQUIRING

    def maximize(self) -> List[int]:
        """Check dyadic window sizes ending at current time for statistically significant excess.

        Tests windows of size 1, 2, 4, 8, ... up to maxtest. Uses phase_counter
        optimization: a window of length h is tested once every h calls.

        Returns:
            List of window sizes that exceed the significance threshold.
        """
        hs = []
        h = 1
        while h <= self.maxtest:
            # checking `2 * self.phase_counter % h` results in window length h being tested *2* times over h calls.
            if self.phase_counter % h:
                break
            win_x = self.acc_x - self.queue_x[-h]
            win_b = self.acc_b - self.queue_b[-h]
            if significance(win_x, win_b) > self.llr_threshold_halfsq:
                hs.append(h)
            h *= 2
        return hs

    def step(self, x: float, b: float) -> List[int]:
        """Process one time step with observed count and background estimate.

        Args:
            x: Observed count at current time step.
            b: Background estimate at current time step.

        Returns:
            List of window sizes that triggered, or empty list if no trigger
            or still in acquisition phase.
        """
        self.acc_x += x
        self.acc_b += b

        if self.status() == TriggerStatus.RUNNING:
            hs = self.maximize()

            self.queue_x.popleft()
            self.queue_b.popleft()
            self.queue_x.append(self.acc_x)
            self.queue_b.append(self.acc_b)
            self.phase_counter = (self.phase_counter + 1) % self.maxtest
            return hs

        else:
            self.queue_x.append(self.acc_x)
            self.queue_b.append(self.acc_b)
            if len(self.queue_x) == self.queue_x.maxlen:
                self.phase_counter = 0
            return []


def search_qbdata(xs: Sequence[int], size: int, maxtest: int, threshold: float) -> List[Hit]:
    """Search a single count time series for transient events.

    Uses moving average for background estimation and dyadic trigger for detection.

    Args:
        xs: Count time series for one quadrant/band combination.
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples.
    """
    bs, vrange = moving_average(xs, size)
    if vrange == _INVALID_INTERVAL:
        return []
    t = TriggerDyadic(maxtest=maxtest, threshold=threshold)
    return t(xs, bs, vrange=vrange)


def search_data(data: Data, checks: Sequence[Index], size: int, maxtest: int, threshold: float) -> List[Hit]:
    """Search data for transient events across multiple band/quadrant combinations.

    Only returns hits that trigger in all specified combinations simultaneously.

    Args:
        data: Nested tuple of count lists indexed by (band, quadrant).
        checks: Sequence of (band, quadrant) combinations to search.
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples that triggered
        in all specified band/quadrant combinations.
    """
    hits_counter = Counter()
    for band, quadrant in checks:
        hits = search_qbdata(data[band][quadrant], size, maxtest, threshold)
        for ih in hits:
            hits_counter[ih] += 1
    return [ih for ih, count in hits_counter.items() if count == len(checks)]


def search_filepath(
    filepath: Union[Path, str],
    checks: Sequence[Index],
    size: int,
    maxtest: int,
    threshold: float,
) -> List[Hit]:
    """Search an SRA file for transient events.

    Parses the file and searches across specified band/quadrant combinations.

    Args:
        filepath: Path to the SRA file.
        checks: Sequence of (band, quadrant) combinations to search.
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples.
    """
    data, abts = sra_parse(filepath)
    return search_data(data, checks, size, maxtest, threshold)


def hit_tointerval(hit: Hit) -> Interval:
    """Convert a trigger hit to a time interval.

    Args:
        hit: A (time_index, window_length) tuple from the trigger algorithm.

    Returns:
        An (start, end) interval where start is the transient start index
        and end is the trigger time index + 1.
    """
    i, h = hit
    return i - h + 1, i + 1


def summarize(hits: Sequence[Hit]) -> Tuple[int, Interval]:
    """Summarize a list of trigger hits into a count and bounding interval.

    Args:
        hits: Sequence of (time_index, window_length) tuples.

    Returns:
        A tuple containing:
            - nhits: Number of hits.
            - interval: The (earliest_start, latest_end) bounding interval,
              or (-1, -1) if no hits.
    """
    nhits = len(hits)
    if nhits > 0:
        starts, ends = zip(*[hit_tointerval(h) for h in hits])
        return nhits, (min(starts), max(ends))
    return nhits, _INVALID_INTERVAL


def main():
    """Search an SRA file for transient events and write results to output file.

    Parses command-line arguments, runs the trigger algorithm on the specified
    SRA file, and writes a summary to <filename>_trigger.txt if events are found.

    Returns:
        ErrorCode.OK (0) on success, ErrorCode.INVALID_FILE (1) if the input
        file is missing or malformed, ErrorCode.INVALID_PARAMETERS (2) if
        command-line arguments are invalid.
    """

    parser = argparse.ArgumentParser(description="Search for transient events in SpIRIT/HERMES-FM1 SRA data.")
    parser.add_argument(
        "filepath",
        help="Path to SRA file.",
    )
    parser.add_argument(
        "--size",
        type=int,
        default=_SEARCH_SIZE_DEFAULT,
        help=f"moving average half-window size (default: {_SEARCH_SIZE_DEFAULT}).",
    )
    parser.add_argument(
        "--maxtest",
        type=int,
        default=_SEARCH_MAXTEST_DEFAULT,
        help=f"maximum trigger window length, must be power of 2 (default: {_SEARCH_MAXTEST_DEFAULT}).",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=_SEARCH_THRESHOLD_DEFAULT,
        help=f"detection threshold in standard deviation units (default: {_SEARCH_THRESHOLD_DEFAULT}).",
    )
    args = parser.parse_args()

    if args.size < 1:
        print("Error: --size must be a positive integer.")
        return ErrorCode.INVALID_PARAMETERS
    if args.maxtest < 1 or (args.maxtest & (args.maxtest - 1)) != 0:
        print("Error: --maxtest must be a positive power of 2.")
        return ErrorCode.INVALID_PARAMETERS
    if args.threshold <= 0:
        print("Error: --threshold must be a positive number.")
        return ErrorCode.INVALID_PARAMETERS

    filepath = Path(args.filepath)

    try:
        nhits, (start, end) = summarize(
            search_filepath(
                filepath,
                checks=_SEARCH_CHECKS_DEFAULT,
                size=args.size,
                maxtest=args.maxtest,
                threshold=args.threshold,
            )
        )
    except FileNotFoundError:
        print("Input file does not exist. Goodbye.")
        return ErrorCode.INVALID_FILE
    except InvalidSRA as e:
        print(f"Input SRA file is not valid: {e}")
        return ErrorCode.INVALID_FILE

    if nhits > 0:
        out_filepath = filepath.parent / f"{filepath.stem}_trigger.txt"
        with open(out_filepath, "w") as f:
            f.write(f"{nhits} {start} {end}")
    return ErrorCode.OK


if __name__ == "__main__":
    sys.exit(main())
    # ~p26
