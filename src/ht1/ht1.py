"""A program and module for searching transients in SpIRIT/HERMES-FM1 SRA data.

This module provides tools to parse SRA files, compute moving averages for
background estimation, and apply a dyadic windowed trigger algorithm to detect
Poisson count excesses (transients).

Example:
    $ python ht1.py srafile.raw --threshold 5.0 --maxtest 16
"""

import argparse
from collections import Counter
from collections import deque
from enum import IntEnum
from itertools import accumulate
import logging
from math import inf
from math import log
from pathlib import Path
import struct
import sys
from typing import List, NamedTuple, Sequence, Tuple, Union

_logger = logging.getLogger(__name__)
_logger.addHandler(logging.NullHandler())


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
_SEARCH_MAXTEST_DEFAULT = 16
_SEARCH_THRESHOLD_DEFAULT = 5.
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
    """Parses an SRA file and extracts count data and ABT timestamps.

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


def moving_average_range(xs: Sequence[int], size: int) -> Interval:
    """Computes a valid range for a simple moving average.

    The range is based on a centered window of length `2*size + 1`. It finds the
    first and last non-zero elements and returns an interval padded inward by
    the window half-size.

    Args:
        xs: Input sequence of integer counts.
        size: Half-window size for the moving average.

    Returns:
        A tuple (start, end), or (-1, -1) if sequence is empty, all zeros,
        or too short for the window size.
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
    """Computes a simple moving average with a centered window.

    Leading and trailing zeroes are ignored.

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

    vrange = moving_average_range(xs, size)
    cumsum = [0] + list(accumulate(xs))
    sma = []
    window_len = 2 * size + 1
    for i in range(*vrange):
        window_sum = cumsum[i + size + 1] - cumsum[i - size]
        sma.append(window_sum / window_len)
    return sma, vrange


def significance(x: float, b: float) -> float:
    """Compute half the Poisson log-likelihood ratio.

    Args:
        x: Observed count.
        b: Expected background count.

    Returns:
        Half the log-likelihood ratio, or 0.0 if x <= b or b <= 0.
    """
    if x <= b or b <= 0:
        return 0.0
    return x * log(x / b) - (x - b)


class TriggerDyadic:
    """Dyadic windowed trigger algorithm for detecting Poisson count excesses.

    Attributes:
        maxtest: Maximum window size to test (must be power of 2).
        llr_threshold_halfsq: Detection threshold as half-squared LLR.
        phase_counter: Controls the cycle of window checks.
        frequency: Controls how many times a window is checked in a phase cycle.
        acc_x: Accumulated sum of observed counts.
        acc_b: Accumulated sum of background counts.
        queue_x: Deque storing count history.
        queue_b: Deque storing background history.
    """

    def __init__(self, maxtest: int, threshold: float):
        """Initializes the trigger algorithm.

        Args:
            maxtest: Maximum window size to test (must be power of 2).
            threshold: Detection threshold in standard deviation units.
        """
        self.maxtest = maxtest
        # convert to half, squared llr threshold for performances
        self.llr_threshold_halfsq = 0.5 * (threshold**2)

        # the choice of the initial phase value determines the time-index at which
        # a window with length h is tested for the first time
        self.phase_counter = 0
        # frequency=1: checking window length h, 1 time over h calls.
        # frequency=2: checking window length h, 2 times over h calls.
        # frequency=4: ... and so on
        # frequency values must be power of two
        self.frequency = 2
        # using accumulator allows computing counts over a window with just one difference.
        # not using an accumulator would require summing many individual counts.
        self.acc_x = 0
        self.acc_b = 0.0
        # accumulator values are stored in a queue.
        # accumulator queues are pre-filled to avoid index errors.
        # this results in a few wasted significance computation during first cycle.
        # we can live with that.
        maxtest_plus_one = maxtest + 1
        self.queue_x = deque([0] * maxtest_plus_one, maxlen=maxtest_plus_one)
        # the background queue is filled with infinity so that we accept a hit from
        # window size h only if we've processed at least h data points.
        # this prevents large windows from triggering over intervals predating first count,
        # which could happen, if the count time series starts with huge values.
        self.queue_b = deque([inf] * maxtest + [0.], maxlen=maxtest_plus_one)

    def __call__(self, xs: Sequence[int], bs: Sequence[float], vrange: Interval) -> List[Hit]:
        """Runs the trigger algorithm on count data with background estimates.

        Windows are checked with half-length offsets or, in other words, a window
        with length h is checked two times every h iterations.

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

    def maximize(self) -> List[int]:
        """Checks dyadic window sizes for significant excesses.

        Tests windows of size 1, 2, 4, 8, ... up to maxtest ending at the
        current time step.

        Returns:
            List of window sizes that exceed the significance threshold.
        """
        hs = []
        h = 1
        while h <= self.maxtest:
            if self.phase_counter % h:
                break
            win_x = self.acc_x - self.queue_x[-h]
            win_b = self.acc_b - self.queue_b[-h]
            if significance(win_x, win_b) > self.llr_threshold_halfsq:
                hs.append(h)
            h *= 2
        return hs

    def step(self, x: int, b: float) -> List[int]:
        """Processes one time step with observed count and background estimate.

        Args:
            x: Observed count at current time step.
            b: Background estimate at current time step.

        Returns:
            List of window sizes that triggered.
        """
        self.acc_x += x
        self.acc_b += b
        self.phase_counter = (self.phase_counter + self.frequency) % self.maxtest
        hs = self.maximize()
        self.queue_x.append(self.acc_x)
        self.queue_b.append(self.acc_b)
        return hs


def search_qbdata(xs: Sequence[int], size: int, maxtest: int, threshold: float) -> List[Hit]:
    """Searches a single count time series for transient events.

    Args:
        xs: Count time series for one quadrant/band combination.
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples.
    """
    bs, (vmin, vmax) = moving_average(xs, size)
    if (vmin, vmax) == _INVALID_INTERVAL:
        _logger.warning("Invalid MA range: data too short or all zeros")
        return []
    _logger.info(f"Algorithm ran {vmax - vmin} iterations")
    t = TriggerDyadic(maxtest=maxtest, threshold=threshold)
    return t(xs, bs, vrange=(vmin, vmax))


def search_data(data: Data, size: int, maxtest: int, threshold: float) -> List[Hit]:
    """Searches data for transients across multiple band/quadrant combinations.

    Searches six quadrant/band combinations: all three quadrants (B, C, D) in both
    the MID and HIGH energy bands. The LOW energy band is excluded from the search.

    Only returns hits that trigger simultaneously (same time index and window length)
    in at least five of the six combinations. This requirement guarantees that all
    three quadrants triggered and at least two energy bands are over threshold.

    Args:
        data: Nested tuple of count lists indexed by (band, quadrant).
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples that met the
        coincidence requirement.
    """
    hits_counter = Counter()
    for band, quadrant in (
        # we ignore triggers in low-energy band
        (EnBand.MID, Quadrant.B),
        (EnBand.MID, Quadrant.C),
        (EnBand.MID, Quadrant.D),
        (EnBand.HIGH, Quadrant.B),
        (EnBand.HIGH, Quadrant.C),
        (EnBand.HIGH, Quadrant.D),
    ):
        _logger.info(f"Running on quadrant {quadrant.name}, band {band.name}")
        hits = search_qbdata(data[band][quadrant], size, maxtest, threshold)
        for ih in hits:
            hits_counter[ih] += 1
    # why *>4*? i think asking for at least 5 quadrant/band combo to be over threshold is the right
    # number because that's the smallest number of combinations necessary to guarantee:
    #   1. all quadrant to be above threshold;
    #   2. over at least two energy band.
    # we could ask for 6, but doing so the conditions seems to stiff. this, relax it a bit.
    return [ih for ih, count in hits_counter.items() if count > 4]


def search_filepath(
    filepath: Union[Path, str],
    size: int,
    maxtest: int,
    threshold: float,
) -> List[Hit]:
    """Searches an SRA file for transient events across multiple band/quadrant combinations.

    Searches six quadrant/band combinations: all three quadrants (B, C, D) in both
    the MID and HIGH energy bands. The LOW energy band is excluded from the search.

    Only returns hits that trigger simultaneously (same time index and window length)
    in at least five of the six combinations. This requirement guarantees that all
    three quadrants triggered and at least two energy bands are over threshold.

    Args:
        filepath: Path to the SRA file.
        size: Half-window size for moving average background estimate.
        maxtest: Maximum trigger window size (must be power of 2).
        threshold: Detection threshold in sigma units.

    Returns:
        List of hits as (time_index, window_length) tuples that met the
        coincidence requirement.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        InvalidSRA: If the file is malformed or invalid.
    """
    data, abts = sra_parse(filepath)
    return search_data(data, size, maxtest, threshold)


def hit_tointerval(hit: Hit) -> Interval:
    """Converts a trigger hit to a time interval.

    Args:
        hit: A (time_index, window_length) tuple.

    Returns:
        An (start, end) interval where `start` is the window start index
        and `end` is the trigger index + 1.
    """
    i, h = hit
    return i - h + 1, i + 1


def summarize(hits: Sequence[Hit]) -> Tuple[int, Interval]:
    """Summarizes a list of trigger hits into a count and bounding interval.

    Args:
        hits: Sequence of (time_index, window_length) tuples.

    Returns:
        A tuple containing:
            - nhits: Total number of hits.
            - interval: The (earliest_start, latest_end) bounding interval.
    """
    nhits = len(hits)
    if nhits > 0:
        starts, ends = zip(*[hit_tointerval(h) for h in hits])
        return nhits, (min(starts), max(ends))
    return nhits, _INVALID_INTERVAL


def main() -> ErrorCode:
    """Searches an SRA file for transient events and writes results to output file.

    Parses command-line arguments and runs the trigger algorithm on the specified SRA file.
    Searches six quadrant/band combinations (quadrants B, C, D in MID and HIGH energy bands)
    and requires at least five simultaneous triggers to declare a detection. This guarantees
    all three quadrants triggered and at least two energy bands are over threshold.

    If transients are detected, writes a summary to <filename>_trigger.txt containing
    the number of hits and the bounding interval.

    Returns:
        - OK (0) on success
        - INVALID_FILE (1) if the input file is missing or malformed
        - INVALID_PARAMETERS (2) if command-line arguments are invalid.
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
    # not passing --log flag sets logging level to error
    # passing --log flag sets logging level to info
    # passing --log=LEVEL sets logging level to LEVEL
    parser.add_argument(
        "--log",
        nargs="?",
        const="info",
        choices=["info", "warning", "error"],
        default="error",
        help="Set logging level (default: error, --log without value: info).",
    )
    args = parser.parse_args()

    level_map = {"info": logging.INFO, "warning": logging.WARNING, "error": logging.ERROR}
    level = level_map[args.log]
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
    _logger.addHandler(handler)
    _logger.setLevel(level)

    if args.size < 1:
        _logger.error("--size must be a positive integer.")
        return ErrorCode.INVALID_PARAMETERS
    if args.maxtest < 1 or (args.maxtest & (args.maxtest - 1)) != 0:
        _logger.error("--maxtest must be a positive power of 2.")
        return ErrorCode.INVALID_PARAMETERS
    if args.threshold <= 0:
        _logger.error("--threshold must be a positive number.")
        return ErrorCode.INVALID_PARAMETERS

    filepath = Path(args.filepath)

    _logger.info(f"Input: {filepath}")
    _logger.info(f"Parameters: size={args.size}, maxtest={args.maxtest}, threshold={args.threshold}")

    try:
        nhits, (start, end) = summarize(
            search_filepath(
                filepath,
                size=args.size,
                maxtest=args.maxtest,
                threshold=args.threshold,
            )
        )
    except FileNotFoundError:
        _logger.error("Input file does not exist")
        return ErrorCode.INVALID_FILE
    except InvalidSRA as e:
        _logger.error(f"Input SRA file is not valid: {e}")
        return ErrorCode.INVALID_FILE

    if nhits > 0:
        _logger.info(f"Found {nhits} trigger hits between time-index {start} and {end}")
        out_filepath = filepath.parent / f"{filepath.stem}_trigger.txt"
        _logger.info(f"Writing to {out_filepath}")
        with open(out_filepath, "w") as f:
            f.write(f"{nhits} {start} {end}")
    else:
        _logger.info("No trigger hits")
    return ErrorCode.OK


if __name__ == "__main__":
    sys.exit(main())
    # ~p26
