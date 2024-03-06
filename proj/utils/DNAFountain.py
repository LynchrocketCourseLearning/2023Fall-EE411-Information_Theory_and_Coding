from typing import List, Tuple
from .misc import process_raw_input
from .robust_solution import PRNG
from .droplet import Droplet
from . import LFSR
from . import scr_rept as sr
from reedsolo import RSCodec
import operator, math


class DNAFountain:
    def __init__(
        self,
        input_file: str,
        chunk_size: int,
        alpha: float,
        rs=0,
        c_dist=0.1,
        delta=0.5,
        max_homopolymer=3,
        gc=0.05,
        final: int = None,
    ):
        """
        alpha: the redundency level
        final: maximal number of oligos; type = int
        chunk_size: in bytes
        rs: the number of bytes for reed-solomon error correcting code over gf(2^8).
        c_dist: a parameter of the degree distribution
        delta: a parameter of the degree distribution
        max_homopolymer: the largest homopolymer allowed
        gc: the allowable range of gc +- 50%
        """

        # things realted to data:
        data_array, data_len = process_raw_input(input_file, chunk_size)
        self.data_array = data_array
        self.chunk_size = chunk_size
        self.num_chunks = int(math.ceil(data_len / float(chunk_size)))
        self.alpha = alpha
        self.final = (
            final if final is not None else int(self.num_chunks * (1 + self.alpha)) + 1
        )

        # things related to random mnumber generator:
        # starting an lfsr with a certain state and a polynomial for 32bits.
        self.lfsr = LFSR.lfsr_s_p()
        # calculate the length of lsfr in bits
        self.lfsr_l = len("{:b}".format(LFSR.lfsr32p())) - 1
        self.seed = next(self.lfsr)

        # creating the solition distribution object
        self.PRNG = PRNG(K=self.num_chunks, delta=delta, c=c_dist)
        self.PRNG.set_seed(self.seed)

        # things related to error correcting code:
        self.rs = rs  # the number of symbols (bytes) to add
        self.rs_obj = RSCodec(self.rs)  # initalizing an reed solomon object

        # things related to biological screens:
        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.tries: int = 0  # number of times we tried to create a droplet
        self.good: int = 0  # droplets that were screened successfully
        # number of nucleotides in an oligo
        self.oligo_l = (self.chunk_size * 8 + self.lfsr_l + self.rs * 8) / 4

    def _rand_chunk_nums(self) -> Tuple[int, List[int]]:
        """
        This function returns a subset of segments based on the solition distribution.
        It updates the lfsr to generates a new seed.
        This function creates a fresh seed for the droplet and primes the solition inverse cdf sampler
        """
        self.seed = next(self.lfsr)  # deploy one round of lfsr, and read the register.
        self.PRNG.set_seed(self.seed)  # update the seed with the register
        degree, ix_samples = self.PRNG.get_src_blocks_wrap()
        return degree, ix_samples  # return a list of segments.

    def droplet(self) -> Droplet:
        # creating a droplet.
        data = None
        # creating a random list of segments.
        degree, num_chunks = self._rand_chunk_nums()
        for num in num_chunks:  # iterating over each segment
            if data is None:  # first round. data payload is empty.
                data = self.data_array[num]  # just copy the segment to the payload.
            else:  # more rounds. Starting xoring the new segments with the payload.
                data = list(map(operator.xor, data, self.data_array[num]))
        self.tries += 1  # upadte counter.

        return Droplet(
            data=data,
            seed=self.seed,
            rs=self.rs,
            rs_obj=self.rs_obj,
            num_chunks=num_chunks,
            degree=degree,
        )

    def screen(self, droplet) -> bool:
        if sr.screen_repeat(droplet, self.max_homopolymer, self.gc):
            self.good += 1
            return True
        return False
