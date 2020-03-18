import numpy as np


def primes_3_to_n(n: int) -> np.array:
    """Returns a array of primes, 3 <= p < n (Faster & more memory-wise numpy
    code).

    Parameters
    ----------
    n: int
        Given limit for finding all prime numbers

    Returns
    -------
    Numpy array
        All primes from 3 to p.

    Raises
    ------
    ValueError
        If n is equal or less than 2.
    """

    if n <= 2:
        raise ValueError("n must be greater than 2")
    sieving_range = np.ones(n//2, dtype=np.bool)
    for i in range(3, int(n**0.5)+1, 2):
        if sieving_range[i//2]:
            sieving_range[i*i//2::i] = False
    return 2 * np.nonzero(sieving_range)[0][1::] + 1


def prime_polynomials(field_order: int, generator: int = 2,
                      fast_primes: bool = False,
                      first_prime: bool = False) -> list:
    """Compute the list of prime polynomials for the given generator and galois
    field prime exponent. Bruteforce approach.

    Parameters
    ----------
    field_order: int
        Order of the field, i.e 2^p, in which the search for prime po-
        lynomials will be done.

    generator: int, default=2
        Generator number (the "increment" that will be used to walk through the
        field by multiplication, this must be a prime number). This is basical
        ly the base of the logarithm/anti-log tables. Also often noted α (al
        pha) in academic books.

    fast_primes: bool, optional
        Will output less results but will be significantly faster.

    first_prime: bool, optional
        will output the first prime polynomial found, so if all you want is to
        just find one prime polynomial to generate the LUT for Reed-Solomon to
        work, then just use that.

    Returns
    -------
    correct_primes: list
        List of all prime polynomials.
    """

    gf_2 = 2  # we're in GF(2)
    max_field_value = int(gf_2**field_order - 1)
    next_field_max_value = int(gf_2**(field_order+1) - 1)

    prim_candidates = np.empty
    if fast_primes:
        prim_candidates = primes_3_to_n(next_field_max_value)
        prim_candidates = [prim_candidate for prim_candidate
                           in prim_candidates
                           if prim_candidate > max_field_value]
    else:
        prim_candidates = range(max_field_value+2, next_field_max_value, gf_2)

    correct_primes = list()
    for prim_candidate in prim_candidates:
        seen = np.zeros(max_field_value+1)
        conflict = False

        x = 1
        for i in range(max_field_value):
            x = gf_mult_noLUT(x, generator, prim_candidate, max_field_value+1)

            if x > max_field_value or seen[x] == 1:  # not prime
                conflict = True
                break
            else:
                seen[x] = 1

        if not conflict:
            correct_primes.append(prim_candidate)
            if first_prime:
                return prim_candidate

    return correct_primes


def look_up_tables(primitive_polynomial: bin, field_order: int,
                   generator: int = 2):
    """Precompute the logarithm and anti-log tables for faster computation la-
    ter, using the provided primitive polynomial.

    These tables are used for multiplication/division. For each possible value
    in the galois field 2^field_order, we will pre-compute the logarithm and an
    ti-logarithm (exponential) of this value, to do that, we generate the Galo
    is Field GF(2^p) by building a list starting with the element 0 followed by
    the (p-1) succesive powers of the generator a : 1, a, a^1, a^2, ..., a^(p-1).

    Parameters
    ----------
    primitive_polynomial: bin
        Primitive/prime (binary) polynomial and must be irreducible (ie, it
        can't represented as the product of two smaller polynomials). It's a po
        lynomial in the binary sense: each bit is a coefficient,but in fact
        it's an integer between max_field_value+1 and max_field_value*2 and not
        a list of gf values. The prime polynomial will be used to reduce the
        overflows back into the range of the Galois Field without duplicating
        values (all values should be unique).

    generator: int, default=2
        Generator number (the "increment" that will be used to walk through the
        field by multiplication, this must be a prime number). This is basical
        ly the base of the logarithm/anti-log tables. Also often noted α (al
        pha) in academic books.

    field_order: int
        Order of the field, i.e 2^p.

    Returns
    -------
    gf_exp: list
        Anti-log (exponential) table. The first two elements will always be
        [GF(2^p)int(1), generator].

    gf_log: list
        Log table, log[0] is impossible and thus unused.
    """

    global gf_exp, gf_log, max_field_value
    max_field_value = int(2**field_order - 1)
    gf_exp = [0] * (max_field_value * 2)
    gf_log = [0] * (max_field_value + 1)

    x = 1
    for i in range(max_field_value):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult_noLUT(x, generator, max_field_value+1,
                          primitive_polynomial)

    for i in range(max_field_value, max_field_value * 2):
        gf_exp[i] = gf_exp[i - max_field_value]

    return [gf_log, gf_exp]


def gf_add(x, y):
    return x ^ y


def gf_sub(x, y):
    return x ^ y


def gf_neg(x):
    return x


def gf_mul(x, y):
    if x == 0 or y == 0:
        return 0
    return gf_exp[(gf_log[x] + gf_log[y]) % max_field_value]


def gf_div(x, y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return gf_exp[(gf_log[x] + max_field_value - gf_log[y]) % max_field_value]


def gf_pow(x, power):
    return gf_exp[(gf_log[x] * power) % max_field_value]


def gf_inverse(x):
    return gf_exp[max_field_value - gf_log[x]]  # gf_inverse(x) == gf_div(1, x)


def gf_mult_noLUT(x: bin, y: bin, field_order: int, poly_modular=0, carryless=True):
    """Galois Field integer multiplication using Russian Peasant Multiplication
    algorithm (faster than the standard multiplication + modular reduction).

    If poly_modular is 0 and carryless=False, then the function produces the re
    sult for a standard integers multiplication (no carry-less arithmetics nor
    modular reduction).

    Parameters
    ---------
    x: bin
        Polynomial in binary form, must be in the same Galois Field as y.

    y: bin
        Polynomial in binary form, must be in the same Galois Field as x.

    field_order: int
        Order of the field, i.e 2^p.

    poly_modular: bin
        Irreducible primitive polynomial in binary form for modular reduction
        after the multiplication between x and y. To stay in the bound of the
        field, we need to compute the modulo of any value above the Galois Field.

    carryless:
        If true multiplication will be carry-less.

    Returns
    -------
    r: bin
        Multiplication of x times y.

    """
    r = 0
    while y:
        if y & 1:
            r = r ^ x if carryless else r + x
        y = y >> 1
        x = x << 1
        if poly_modular > 0 and x & field_order:
            x = x ^ poly_modular

    return r


if __name__ == "__main__":
    ltables = look_up_tables(0x11D, 8)
    print(bin(gf_mul(0b10001001, 0b00101010)))
