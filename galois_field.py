import numpy as np


def primes_3_to_n(n):
    """Returns a array of primes, 3 <= p < n (Faster & more memory-wise numpy
    code).

    param
    -------
    n: given limit for finding all prime numbers

    returns
    -------
    numpy array of all primes from 3 to p.

    raises
    -------
    ValueError: If n is equal or less than 2.
    """

    if n <= 2:
        raise ValueError("n must be greater than 2")
    sieving_range = np.ones(n//2, dtype=np.bool)
    for i in range(3, int(n**0.5)+1, 2):
        if sieving_range[i//2]:
            sieving_range[i*i//2::i] = False
    return 2 * np.nonzero(sieving_range)[0][1::] + 1


def prime_polynomials(field_order, generator=2, fast_primes=False,
                      first_prime=False):
    """Compute the list of prime polynomials for the given generator and galois
    field prime exponent. Bruteforce approach.

    param
    -------
    field_order: order of the field, i.e 2^p, in which the search for prime po-
    lynomials will be done.

    generator: generator number (the "increment" that will be used to walk
    through the field by multiplication, this must be a prime number). This is
    basically the base of the logarithm/anti-log tables. Also often noted α (al
    pha) in academic books.

    fast_primes: will output less results but will be significantly faster.

    first_prime: will output the first prime polynomial found, so if all you
    want is to just find one prime polynomial to generate the LUT for
    Reed-Solomon to work, then just use that.

    returns
    -------
    correct_primes: list of all prime polynomials

    raises
    -------
    None
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

        # Second loop, build the whole Galois Field
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


def look_up_tables(primitive_polynomial, field_order, generator=2):
    """Precompute the logarithm and anti-log tables for faster computation la-
    ter, using the provided primitive polynomial.These tables are used for mul-
    tiplication/division. For each possible value in the galois field
    2^field_order, we will pre-compute the logarithm and anti-logarithm (expo
    nential) of this value, to do that, we generate the Galois Field GF(2^p)
    by building a list starting with the element 0 followed by the (p-1) succes
    sive powers of the generator a : 1, a, a^1, a^2, ..., a^(p-1).

    param
    -------
    primitive_polynomial: is the primitive/prime (binary) polynomial and must
    be irreducible (ie, it can't represented as the product of two smaller poly
    nomials). It's a polynomial in the binary sense: each bit is a coefficient,
    but in fact it's an integer between max_field_value+1 and max_field_value*2
    and not a list of gf values. The prime polynomial will be used to reduce
    the overflows back into the range of the Galois Field without duplicating
    values (all values should be unique).

    generator: generator number (the "increment" that will be used to walk
    through the field by multiplication, this must be a prime number). This is
    basically the base of the logarithm/anti-log tables. Also often noted α (al
    pha) in academic books.

    field_order: order of the field, i.e 2^p.


    returns
    -------
    gf_exp: anti-log (exponential) table. The first two elements will always be
    [GF256int(1), generator].

    gf_log: log table, log[0] is impossible and thus unused.


    raises
    -------
    ZeroDivisionError

    """

    global gf_exp, gf_log, max_field_value
    max_field_value = int(2**field_order - 1)
    gf_exp = np.zeros(max_field_value * 2)
    gf_log = np.zeros(max_field_value + 1)

    x = 1
    for i in range(max_field_value):
        gf_exp[i] = x
        gf_log[x] = i
        x = gf_mult_noLUT(x, generator, primitive_polynomial, max_field_value+1)

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


def gf_mult_noLUT(x, y, prim=0, field_charac_full=65536, carryless=True):
    """Galois Field integer multiplication using Russian Peasant Multiplication
    algorithm (faster than the standard multiplication + modular reduction). If
    prim is 0 and carryless=False, then the function produces the result for a
    standard integers multiplication (no carry-less arithmetics nor modular re-
    duction).

    param
    -------

    returns
    -------

    raise
    -------

    """
    r = 0
    while y:
        if y & 1:
            r = r ^ x if carryless else r + x
        y = y >> 1
        x = x << 1
        if prim > 0 and x & field_charac_full:
            x = x ^ prim

    return r


if __name__ == "__main__":
    print(prime_polynomials(field_order=8))
