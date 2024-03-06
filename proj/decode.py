import logging
import os
from collections import defaultdict
from utils.glass import Glass


class Decoder:
    def __init__(
        self,
        input_file: str,
        output_file: str,
        chunk_num=128,
        header_size=4,
        rs=0,
        delta=0.05,
        c_dist=0.1,
        gc=0.5,
        max_homopolymer=4,
        max_hamming=100,
        chunk_size=32,
    ):
        """
        input_file: file to decode; type = str
        output_file: output file; type = str
        chunk_num: the total number of chunks in the file; type = int
        header_size: number of bytes for the header; type = int
        rs: number of bytes for rs codes; type = int
        delta: Degree distribution tuning parameter; type = float
        c_dist: Degree distribution tuning parameter; type = float
        gc: range of gc content; type = float
        max_homopolymer: the largest number of nt in a homopolymer; type = int
        max_hamming: How many differences between sequenced DNA and corrected DNA to tolerate; type = int
        chunk_size: The number of bytes of the data payload in each DNA string; type = int
        """
        logging.basicConfig(level=logging.DEBUG)
        if not os.path.exists(input_file):
            logging.error("{input_file} file not found")
            exit(1)
        self.input_file = input_file
        self.output_file = output_file
        self.chunk_num = chunk_num
        self.header_size = header_size
        self.rs = rs
        self.delta = delta
        self.c_dist = c_dist
        self.gc = gc
        self.max_homopolymer = max_homopolymer
        self.max_hamming = max_hamming
        self.chunk_size = chunk_size
        self.glass = Glass(
            self.chunk_num,
            header_size=self.header_size,
            rs=self.rs,
            c_dist=self.c_dist,
            delta=self.delta,
            gc=self.gc,
            max_homopolymer=max_homopolymer,
            max_hamming=max_hamming,
        )

    def decode(self) -> None:
        line = 0
        errors = 0
        seen_seeds = defaultdict(int)
        with open(self.input_file, "r") as file:
            while True:
                dna = file.readline().rstrip("\n")
                if len(dna) == 0:
                    logging.info("Finished reading input file!")
                    break

                line += 1
                seed, data = self.glass.add_dna(dna)
                # Exclude the sequence with error, which is founded by RS code
                if seed == -1:
                    errors += 1
                else:
                    seen_seeds[seed] += 1

                if line % 1000 == 0:
                    logging.info(
                        "After reading {} lines. {} chunks done. {} rejections.".format(
                            line, self.glass.chunks_done(), errors
                        )
                    )
                if self.glass.is_done():
                    logging.info(
                        "Done! Totally {} lines are read. {} chunks done. {} rejections.".format(
                            line, self.glass.chunks_done(), errors
                        )
                    )
                    break
            if not self.glass.is_done():
                logging.error("Could not decode all file...")
                exit(1)

        out_str = bytes(self.glass.flatten_chunks())
        with open(self.output_file, "wb") as file:
            file.write(out_str)
