"""
A program for searching transients in SpIRIT/HERMES-FM1 SRA data.
"""

import argparse
import struct
from math import log
from pathlib import Path
from typing import NamedTuple, Tuple, List, Dict


NQUADRANTS: int = 3
NBANDS: int = 3


class SRASpecs(NamedTuple):
    """SRA file specifics and header definition."""
    BLOCK_SIZE: int = 124
    HEADER_SIZE: int = 33
    FILELEN: tuple = (21, 25)
    FILELEN_FMT: str = "=I"
    ABT_SIZE: int = 4
    ABT_FMT: str = "=I"
    DATA_SIZE: int = 10
    DATA_FMT: str = "10B"


BandData = Tuple[List[int], ...] # one list for quadrant
Data = Tuple[BandData, ...]

def sra_parse(filepath: Path) -> Tuple[Data, List[int],]:
    """Parses a SRA file and returns data and ABTs as lists.
    Data are organized as a list of lists, the most external layer for the energy band,
    the most internal for the quadrant. According to this convention, the quadrant B (0),
    high-energy band (2) shall is selected with `data[2][0]`"""
    sra_specs = SRASpecs()
    with open(filepath, "rb") as f:
        header = f.read(sra_specs.HEADER_SIZE)
        fsize, = struct.unpack(sra_specs.FILELEN_FMT, header[slice(*sra_specs.FILELEN)])
        assert (fsize - sra_specs.HEADER_SIZE) % sra_specs.BLOCK_SIZE == 0
        abts = []; data =  tuple(tuple([] for _ in range(NQUADRANTS)) for _ in range(NBANDS))
        for block in range((fsize - sra_specs.HEADER_SIZE) // sra_specs.BLOCK_SIZE):
            abts.extend(struct.unpack(sra_specs.ABT_FMT, f.read(sra_specs.ABT_SIZE)))
            for band in range(NBANDS):
                # the quadrant A of FM1 is not working, we read but not store its data
                _ = struct.unpack(sra_specs.DATA_FMT, f.read(sra_specs.DATA_SIZE))
                for quad in range(NQUADRANTS):
                    data[band][quad].extend(struct.unpack(sra_specs.DATA_FMT, f.read(sra_specs.DATA_SIZE)))
    return data, abts


def _find_spikes(bdata: BandData, thr: float=10.0) -> Dict[int, List[int]]:
    """Finds spike indices and affected quadrants in a single energy band.
    Returns a dictionary mapping spike indices to lists of quadrants where the spike was detected.
    Only returns spikes appearing in fewer than all quadrants (detector noise, not real events)."""
    # TODO: at present, we do not check the first and last bin for spikes.
    spikes = {}
    thr_sq_half = thr ** 2 / 2
    for quad in range(NQUADRANTS):
        qdata = bdata[quad]
        for i in range(1, len(qdata) - 1):
            x, b = qdata[i], qdata[i - 1] + qdata[i + 1] + 1
            if x > 0 and (x * log(x / b) - (x - b)) > thr_sq_half:
                spikes.setdefault(i, []).append(quad)
    return {k: v for k, v in spikes.items() if len(v) < NQUADRANTS}


def _fix_spikes(bdata: BandData, spikes: Dict[int, List[int]]):
    """In-place fixes noise spikes by replacing them with the average of neighboring values.
    Only modifies the specific quadrants where each spike was detected."""
    for i, quads in spikes.items():
        for quad in quads:
            bdata[quad][i] = 0.5 * (bdata[quad][i - 1] + bdata[quad][i + 1])


def find_fix_spikes(data: Data, thr: float=10.0):
    """Detects and fixes single-bin noise spikes across all energy bands.
    Spikes appearing in all quadrants simultaneously are preserved as potential real events."""
    for band in range(NBANDS):
        _fix_spikes(data[band], _find_spikes(data[band], thr))


def main():
    """Script interface."""
    parser = argparse.ArgumentParser(description="...")
    parser.add_argument("filepath", help="Path to SRA file.")
    args = parser.parse_args()

    filepath = Path(args.filepath)
    data, _ = sra_parse(filepath)
    band, quad = 1, 2
    print(data[band][quad], len(data[band][quad]))


if __name__ == "__main__":
    main()