import struct
from typing import List
from reedsolo import RSCodec


class Droplet:
    def __init__(
        self,
        data: List[int],
        seed: int,
        num_chunks: List[int] = None,
        rs=0,
        rs_obj: RSCodec = None,
        degree: int = None,
    ):
        self.data: List[int] = data
        self.seed = seed
        self.num_chunks = set(num_chunks)
        self.rs = rs
        self.rs_obj: RSCodec = rs_obj
        self.degree: int = degree
        self.dna: str = None

    def to_dna(self) -> str:
        if self.dna is not None:
            return self.dna
        seed_ord = [c for c in struct.pack("!I", self.seed)]
        message = seed_ord + self.data
        if self.rs > 0:
            # adding RS symbols to the message
            message = self.rs_obj.encode(message)
        # convert to a long sring of binary values
        binary_data = "".join("{0:08b}".format(e) for e in message)
        # convert binary array to a string of 0,1,2,3
        self.dna = "".join(
            str(int(binary_data[i : i + 2], 2)) for i in range(0, len(binary_data), 2)
        )
        return self.dna

    def to_readable_dna(self) -> str:
        return (
            self.to_dna()
            .replace("0", "A")
            .replace("1", "C")
            .replace("2", "G")
            .replace("3", "T")
        )
