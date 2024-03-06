from utils.DNAFountain import DNAFountain
import logging, tqdm


class Encoder:
    def __init__(
        self,
        input_file: str,
        output_file: str,
        chunk_size=128,
        max_homopolymer=4,
        gc=0.2,
        rs=0,
        delta=0.05,
        c_dist=0.1,
        alpha=0.01,
        final: int = None,
    ):
        """
        input_file: file to encode; type = str
        output_file: file to have DNA oligos written in; type = str
        chunk_size: number of information bytes per message; type = int
        max_homopolymer: the largest number of nt in a homopolymer; type = int
        gc: the fraction of gc content above/below 0.5 (example:0.1 means 0.4-0.6); type = restricted_float
        rs: number of bytes for rs codes; type = int
        delta: degree distribution tuning parameter; type = float
        c_dist: degree distribution tuning parameter; type = float
        final: maximal number of oligos; type = int
        alpha: number of more fragments to generate on top of frst k (example: 0.1 will generate 10 percent more fragments); type = float
        """
        logging.basicConfig(level=logging.DEBUG)
        self.output_file = output_file
        if gc < 0.0 or gc > 1.0:
            logging.error("%s not in range [0.0, 1.0]", self.gc)
            exit(1)
        self.dna_fountain = DNAFountain(
            input_file=input_file,
            chunk_size=chunk_size,
            rs=rs,
            max_homopolymer=max_homopolymer,
            gc=gc,
            delta=delta,
            c_dist=c_dist,
            alpha=alpha,
            final=final,
        )

    def encode(self):
        debug_info = self.dna_fountain.PRNG.debug()
        logging.info(
            "Upper bounds on packets for decoding is {} (x{}) with {} probability\n".format(
                debug_info["K_prime"], debug_info["Z"], debug_info["delta"]
            )
        )

        with open(self.output_file, "w") as out, tqdm.tqdm(
            total=self.dna_fountain.final, desc="Valid oligos"
        ) as pbar:
            used_bc: set[int] = set()
            while self.dna_fountain.good < self.dna_fountain.final:
                droplet = self.dna_fountain.droplet()

                if self.dna_fountain.screen(droplet):
                    out.write("{}\n".format(droplet.to_readable_dna()))

                    if droplet.seed in used_bc:
                        logging.error(
                            "Seed %d has been seen before\nDone", droplet.seed
                        )
                        exit(1)

                    used_bc.add(droplet.seed)
                    pbar.update()

            logging.info(
                "Finished. Generated %d packets out of %d tries (%.3f)",
                self.dna_fountain.good,
                self.dna_fountain.tries,
                (self.dna_fountain.good + 0.0) / self.dna_fountain.tries,
            )
