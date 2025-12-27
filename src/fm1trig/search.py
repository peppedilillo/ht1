import argparse
import struct
from pathlib import Path
from typing import NamedTuple, Tuple, List

NQUADRANTS: int = 4
NBANDS: int = 3

class SRASpecs(NamedTuple):
    """SRA file specifics and header definition."""
    BLOCK_SIZE: int = 124
    HEADER_SIZE: int = 33
    # FNAME: tuple = (1, 21)
    # FNAME_FMT: str = "20c"
    FILELEN: tuple = (21, 25)
    FILELEN_FMT: str = "=I"
    # BLOCKSTART: tuple = (25, 27)
    # BLOCKSTART_FMT: str = "=H"
    # BLOCKEND: tuple = (27, 29)
    # BLOCKEND_FMT: str = "=H"
    # ROWSTART: tuple = (29, 31)
    # ROWSTART_FMT: str = "=H"
    # ROWEND: tuple = (31, 33)
    # ROWEND_FMT: str = "=H"
    ABT_SIZE: int = 4
    ABT_FMT: str = "=I"
    DATA_SIZE: int = 10
    DATA_FMT: str = "10B"


def sra_parse(filepath: Path) -> Tuple[List, List]:
    """Parses a SRA file and returns ABTs and data."""
    sra_specs = SRASpecs()
    with open(filepath, "rb") as f:
        header = f.read(sra_specs.HEADER_SIZE)
        fsize, = struct.unpack(sra_specs.FILELEN_FMT, header[slice(*sra_specs.FILELEN)])
        assert (fsize - sra_specs.HEADER_SIZE) % sra_specs.BLOCK_SIZE == 0
        abts = []; data =  [[[] for _ in range(NQUADRANTS)] for _ in range(NBANDS)]
        for block in range((fsize - sra_specs.HEADER_SIZE) // sra_specs.BLOCK_SIZE):
            abts.extend(struct.unpack(sra_specs.ABT_FMT, f.read(sra_specs.ABT_SIZE)))
            for band in range(NBANDS):
                for quad in range(NQUADRANTS):
                    data[band][quad].extend(
                        struct.unpack(sra_specs.DATA_FMT, f.read(sra_specs.DATA_SIZE)))
    return data, abts


def main():
    """Script interface."""
    parser = argparse.ArgumentParser(description="...")
    parser.add_argument("filepath", help="Path to SRA file.")
    args = parser.parse_args()

    filepath = Path(args.filepath)
    data, _ = sra_parse(filepath)
    band, quad = 1, 3
    print(data[band][quad], len(data[band][quad]))


if __name__ == "__main__":
    main()