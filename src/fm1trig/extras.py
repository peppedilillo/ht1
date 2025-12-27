import struct
from pathlib import Path

from fm1trig.search import SRASpecs


def sra_duration(filepath: Path) -> float:
    sra_specs = SRASpecs()
    with open(filepath, "rb") as f:
        header = f.read(sra_specs.HEADER_SIZE)
        fsize, = struct.unpack(sra_specs.FILELEN_FMT, header[slice(*sra_specs.FILELEN)])
    return (fsize - sra_specs.HEADER_SIZE) / sra_specs.BLOCK_SIZE
