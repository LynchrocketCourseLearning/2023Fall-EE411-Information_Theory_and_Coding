from collections import defaultdict
from typing import List, Tuple
from reedsolo import RSCodec
import numpy as np
import operator
from .robust_solution import PRNG
from .droplet import Droplet
from . import scr_rept as sr


class Glass:
    def __init__(
        self,
        num_chunks: int,
        header_size=4,
        rs=0,
        c_dist=0.1,
        delta=0.05,
        gc=0.2,
        max_hamming=100,
        max_homopolymer=4,
    ):
        self.entries = list()
        self.droplets: set[Droplet] = set()
        self.num_chunks = num_chunks
        self.chunks: List[List[int]] = [None] * num_chunks
        self.header_size = header_size
        self.rs = rs
        self.max_hamming = max_hamming
        self.max_homopolymer = max_homopolymer
        self.gc = gc
        self.chunk_to_droplets = defaultdict(set)
        self.done_segments = set()
        self.RSCodec = RSCodec(self.rs)
        self.seen_seeds = set()
        self.PRNG = PRNG(K=self.num_chunks, delta=delta, c=c_dist)

    def _dna_to_int_arr(self, dna_str: str) -> List[int]:
        num_str = (
            dna_str.replace("A", "0")
            .replace("C", "1")
            .replace("G", "2")
            .replace("T", "3")
        )
        s = "".join("{0:02b}".format(int(num_str[i])) for i in range(len(num_str)))
        data = [int(s[t : t + 8], 2) for t in range(0, len(s), 8)]
        return data

    def add_dna(self, dna_str: str) -> Tuple[int, List[int]]:
        data = self._dna_to_int_arr(dna_str)

        try:
            # evaluate the error correcting code
            data_corrected = list(self.RSCodec.decode(data)[0])
        except:
            # could not correct the code
            return -1, None
        # encode the data again to evaluate the correctness of the decoding
        # list() is to convert byte array to int
        data_again = list(self.RSCodec.encode(data_corrected))
        # measuring hamming distance between raw input and expected raw input
        if np.count_nonzero(data != list(data_again)) > self.max_hamming:
            # too many errors to correct in decoding
            return -1, None

        seed_array = data_corrected[: self.header_size]
        seed = sum([int(x) * 256**i for i, x in enumerate(seed_array[::-1])])
        payload = data_corrected[self.header_size :]

        # more error detection (filter seen seeds)
        if seed in self.seen_seeds:
            return -1, None
        self.seen_seeds.add(seed)

        # create droplet from DNA
        self.PRNG.set_seed(seed)
        ix_samples = self.PRNG.get_src_blocks_wrap()[1]
        droplet = Droplet(payload, seed, ix_samples)
        # # more error detection (filter DNA that does not make sense)
        if not sr.screen_repeat(droplet, self.max_homopolymer, self.gc):
            return -1, None
        self.add_droplet(droplet)
        return seed, data

    def add_droplet(self, droplet: Droplet) -> None:
        self.droplets.add(droplet)
        for chunk_num in droplet.num_chunks:
            # document for each chunk all connected droplets
            self.chunk_to_droplets[chunk_num].add(droplet)
        # one round of message passing
        self.update_entry(droplet)

    def update_entry(self, droplet: Droplet) -> None:
        # removing solved segments from droplets
        for chunk_num in droplet.num_chunks & self.done_segments:
            # subtract (ie. xor) the value of the solved segment from the droplet.
            droplet.data = list(map(operator.xor, droplet.data, self.chunks[chunk_num]))
            # cut the edge between droplet and input segment.
            droplet.num_chunks.remove(chunk_num)
            # cut the edge between the input segment to the droplet
            self.chunk_to_droplets[chunk_num].discard(droplet)

        # solving segments when the droplet have exactly 1 segment
        if len(droplet.num_chunks) == 1:
            lone_chunk = droplet.num_chunks.pop()

            # assign the droplet value to the input segment (=entry[0][0])
            self.chunks[lone_chunk] = droplet.data
            # add the lone_chunk to a data structure of done segments.
            self.done_segments.add(lone_chunk)
            # cut the edge between the droplet and input segment
            self.droplets.discard(droplet)
            # cut the edge between the input segment and the droplet
            self.chunk_to_droplets[lone_chunk].discard(droplet)

            # update other droplets
            for other_droplet in self.chunk_to_droplets[lone_chunk].copy():
                self.update_entry(other_droplet)

    def is_done(self) -> bool:
        return self.num_chunks <= len(self.done_segments)

    def chunks_done(self) -> int:
        return len(self.done_segments)

    def flatten_chunks(self) -> str:
        return [ch for chunk in self.chunks for ch in chunk]
