"""
A program for searching transients in SpIRIT/HERMES-FM1 SRA data.
"""

import argparse
import struct
from itertools import accumulate
from math import log
from collections import deque, Counter
from pathlib import Path
from typing import NamedTuple, Tuple, List, Sequence
from enum import Enum, IntEnum


NQUADRANTS: int = 3
NBANDS: int = 3

class SRASpecs(NamedTuple):
    """SRA file specifics and header definition."""
    BLOCK_SIZE: int = 124
    HEADER_SIZE: int = 33
    FILELEN: Tuple[int, int] = (21, 25)
    FILELEN_FMT: str = "=I"
    ABT_SIZE: int = 4
    ABT_FMT: str = "=I"
    DATA_SIZE: int = 10
    DATA_FMT: str = "10B"

class Quadrant(IntEnum):
    """Quadrant name/index"""
    B = 0
    C = 1
    D = 2

class EnBand(IntEnum):
    """Energy band name/index"""
    LOW = 0
    MID = 1
    HIGH = 2

BandData = Tuple[List[int], ...]  # one list per quadrant
Data = Tuple[BandData, ...]  # one BandData per band
Hit = Tuple[int, int]  # trigger time-index, window length
Index = Tuple[EnBand, Quadrant]  # band-quadrant combination index


def sra_parse(filepath: Path) -> Tuple[Data, List[int]]:
    """Parses SRA file and returns data and ABTs as lists.
    Data are organized as a list of lists, the most external layer for the energy band,
    the most internal for the quadrant. According to this convention, the quadrant B (0),
    high-energy band (2) shall is selected with `data[2][0]`

    Raises a ValueError if header file size value is larger than actual file size."""
    sra_specs = SRASpecs()
    block_fmt = sra_specs.ABT_FMT + sra_specs.DATA_FMT * ((NQUADRANTS + 1) * NBANDS)
    block_struct = struct.Struct(block_fmt)

    fsize = filepath.stat().st_size
    with open(filepath, "rb") as f:
        if fsize < sra_specs.HEADER_SIZE:
            raise ValueError("File size smaller than header size.")
        header = f.read(sra_specs.HEADER_SIZE)
        hfsize, = struct.unpack(sra_specs.FILELEN_FMT, header[slice(*sra_specs.FILELEN)])
        nblocks = (hfsize - sra_specs.HEADER_SIZE) // sra_specs.BLOCK_SIZE
        if nblocks * sra_specs.BLOCK_SIZE > fsize - sra_specs.HEADER_SIZE:
            raise ValueError("Invalid header fsize.")
        abts = [0] * nblocks
        data = tuple(tuple([0] * (nblocks * sra_specs.DATA_SIZE) for _ in range(NQUADRANTS)) for _ in range(NBANDS))
        for block in range(nblocks):
            values = block_struct.unpack(f.read(sra_specs.BLOCK_SIZE))
            abts[block] = values[0]
            offset = 1  # skip ABT
            for band in range(NBANDS):
                offset += sra_specs.DATA_SIZE  # skip quadrant A
                for quad in range(NQUADRANTS):
                    start = block * sra_specs.DATA_SIZE
                    end = start + sra_specs.DATA_SIZE
                    data[band][quad][start:end] = values[offset:offset + sra_specs.DATA_SIZE]
                    offset += sra_specs.DATA_SIZE
    return data, abts


def moving_average(data: Sequence[int], size: int) -> List[float]:
    """Computes centered, simple moving average of `data`, at `2 * size + 1` window length.
    At edges, the window size is asymmetrical and shrinks down to `size + 1`."""
    n = len(data)
    cumsum = [0] + list(accumulate(data))

    sma = []
    for i in range(n):
        lower = max(0, i - size)
        upper = min(n, i + size + 1)

        window_sum = cumsum[upper] - cumsum[lower]
        window_len = upper - lower
        sma.append(window_sum / window_len)
    return sma


class TriggerStatus(Enum):
    """Status of the trigger algorithm: ACQUIRING initial data or RUNNING detection."""
    ACQUIRING = 0
    RUNNING = 1


def significance(x: float, b: float) -> float:
    """Returns half-squared Poisson log-likelihood ratio.
    Only returns positive values for excesses (x > b)."""
    if x <= b or b <= 0:
        return 0.0
    return x * log(x / b) - (x - b)


class TriggerDyadic:
    """Dyadic windowed trigger algorithm for detecting excess in Poisson count series.
    Searches for segments with significance exceeding threshold using power-of-2 window sizes."""

    def __init__(self, foreground_len: int, threshold: float):
        """Initializes trigger with test maximum window size and detection threshold.
        The threshold is interpreted as sigma and internally converted to half-squared LLR."""
        self.foreground_len = foreground_len
        self.llr_threshold_halfsq = 0.5 * (threshold ** 2)

        self.phase_counter = -1
        self.acc_x = 0.0
        self.acc_b = 0.0
        self.queue_x = deque([0.0], maxlen=foreground_len + 1)
        self.queue_b = deque([0.0], maxlen=foreground_len + 1)

    def __call__(self, xs: Sequence[int], bs: Sequence[float]) -> List[Hit]:
        """Runs algorithm on count data and associated background estimates."""
        hits = []
        for i, (x, b) in enumerate(zip(xs, bs)):
            hits.extend([(i, h) for h in self.step(x, b)])
        return hits

    def status(self) -> TriggerStatus:
        """Returns current status: ACQUIRING if still collecting initial data, RUNNING otherwise."""
        return TriggerStatus.RUNNING if self.phase_counter >= 0 else TriggerStatus.ACQUIRING

    def maximize(self) -> List[int]:
        """Checks dyadic window sizes (1, 2, 4, 8...) ending at current time.
        Returns list of window sizes that exceed significance threshold.
        Uses phase_counter optimization: a window with length 8 is checked once every 8 calls."""
        hs = []
        h = 1
        while h <= self.foreground_len:
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
        """Processes one time step with observed count x_t and background estimate b_t.
        Returns list of window sizes that triggered, or empty list if no trigger or still acquiring."""
        self.acc_x += x
        self.acc_b += b

        if self.status() == TriggerStatus.RUNNING:
            hs = self.maximize()

            self.queue_x.popleft()
            self.queue_b.popleft()
            self.queue_x.append(self.acc_x)
            self.queue_b.append(self.acc_b)
            self.phase_counter = (self.phase_counter + 1) % self.foreground_len
            return hs

        else:
            self.queue_x.append(self.acc_x)
            self.queue_b.append(self.acc_b)
            if len(self.queue_x) == self.queue_x.maxlen:
                self.phase_counter = 0
            return []


_INVALID_RANGE = (-1, -1)
def data_valid_range(data: Sequence[int], size: int) -> Tuple[int, int]:
    """Return range between first and last non-zero element, padded by the moving average window size.
    In other words, this returns range over which moving average ran over presumably sane data."""
    i = 0
    while data[i] == 0:
        i += 1
    if i == len(data):
        return _INVALID_RANGE
    j = len(data) - 1
    while data[j] == 0:
        j -= 1
    # no need to check here, we have at least one zero
    i = min(i + size, len(data) - 1)
    j = max(0, j - size)
    if j < i:
        return _INVALID_RANGE
    return i, j + 1


def search_data(
        data: Data,
        checks: Sequence[Index],
        size: int,
        foreground_len: int,
        threshold: float
) -> List[Hit]:
    n = len(data[0][0])
    size = size if 2 * size + 1 <= n else n // 2 - 1
    vranges = {}
    hits_counter = Counter()
    for band, quadrant in checks:
        xs = data[band][quadrant]
        bs = moving_average(xs, size)
        t = TriggerDyadic(foreground_len=foreground_len, threshold=threshold)
        for ih in t(xs, bs):
            i, h = ih
            vmin, vmax = vranges.setdefault((band, quadrant), data_valid_range(xs, size))
            if vmin <= i < vmax:
                hits_counter[ih] += 1
    return [ih for ih, count in hits_counter.items() if count == len(checks)]


def search(
        filepath: Path,
        checks: Sequence[Index] = (
            (EnBand.MID, Quadrant.B),
            (EnBand.MID, Quadrant.C),
            (EnBand.MID, Quadrant.D),
            (EnBand.HIGH, Quadrant.B),
            (EnBand.HIGH, Quadrant.C),
            (EnBand.HIGH, Quadrant.D),
        ),
        size: int = 50,
        foreground_len: int = 8,
        threshold: float = 5.0,
) -> List[Hit]:
    data, abts = sra_parse(filepath)
    return search_data(data, checks, size, foreground_len, threshold)


def main():
    """Script interface."""
    parser = argparse.ArgumentParser(description="...")
    parser.add_argument("filepath", help="Path to SRA file.")
    args = parser.parse_args()

    filepath = Path(args.filepath)
    data, abts = sra_parse(filepath)
    ...


if __name__ == "__main__":
    main()