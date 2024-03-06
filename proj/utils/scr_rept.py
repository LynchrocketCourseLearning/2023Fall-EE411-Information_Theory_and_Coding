from utils.droplet import Droplet


def screen_repeat(drop: Droplet, max_homopolymer: int, gc_dev: float) -> bool:
    As = "0" * (max_homopolymer + 1)
    Cs = "1" * (max_homopolymer + 1)
    Gs = "2" * (max_homopolymer + 1)
    Ts = "3" * (max_homopolymer + 1)
    dna = drop.to_dna()
    if As in dna or Cs in dna or Gs in dna or Ts in dna:
        return False
    gc = dna.count("1") + dna.count("2")
    gc = gc / (len(dna) + 0.0)
    if (gc < 0.5 - gc_dev) or (gc > 0.5 + gc_dev):
        return False
    return True
