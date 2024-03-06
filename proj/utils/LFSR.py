def lfsr(state: int, mask: int):
    """
    Galois lfsr:
    """
    result = state
    nbits = mask.bit_length() - 1
    while True:
        result = result << 1
        xor = result >> nbits
        if xor != 0:
            result ^= mask
        yield result


def lfsr32p() -> int:
    """
    This function returns a hard coded polynomial (0b100000000000000000000000011000101).
    The polynomial corresponds to 1 + x^25 + x^26 + x^30 + x^32, which is known
    to repeat only after 32^2-1 tries. Don't change unless you know what you are doing.
    """
    return 0b100000000000000000000000011000101


def lfsr32s() -> int:
    """
    This function returns a hard coded state for the lfsr (0b001010101)
    this state is the inital position in the register. You can change it without a major implication.
    """
    return 0b001010101


def lfsr_s_p():
    return lfsr(lfsr32s(), lfsr32p())
