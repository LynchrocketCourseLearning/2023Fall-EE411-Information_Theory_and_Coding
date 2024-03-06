import math, random, numpy
from typing import Dict, Any, List, Tuple, Union


class PRNG:
    def __init__(self, K: float, delta: float, c: float):
        self.K = float(K)
        self.delta = delta
        self.c = c
        self.S = self.c * math.log(self.K / self.delta) * math.sqrt(self.K)
        self.cdf, self.Z = self._gen_rsd_cdf(K, self.S, self.delta)
        self.state: int = 1  # seed

    def _gen_rsd_cdf(self, K, S, delta) -> Tuple[numpy.ndarray, int]:
        pivot = int(math.floor(K / S))
        val1 = [S / K * 1 / d for d in range(1, pivot)]
        val2 = [S / K * math.log(S / delta)]
        val3 = [0 for _ in range(pivot, K)]
        tau = val1 + val2 + val3
        rho = [1.0 / K] + [1.0 / (d * (d - 1)) for d in range(2, K + 1)]
        Z = sum(rho) + sum(tau)
        mu = [(rho[d] + tau[d]) / Z for d in range(K)]
        cdf = numpy.cumsum(mu)
        return cdf, Z

    def set_seed(self, seed: int) -> None:
        self.state = seed

    def get_src_blocks_wrap(self, seed: int = None) -> Tuple[int, List[int]]:
        if seed:
            self.state = seed
        random.seed(self.state)
        p = random.random()
        d = self._sample_d(p)
        nums = random.sample(range(int(self.K)), d)
        return d, nums

    def _sample_d(self, p: float) -> int:
        for ix, v in enumerate(self.cdf):
            if v > p:
                return ix + 1
        return ix + 1

    def debug(self) -> Dict[str, Union[float, int]]:
        return {
            "K": self.K,
            "delta": self.delta,
            "c": self.c,
            "S": self.S,
            "Z": self.Z,
            "K_prime": self.K * self.Z,
        }
