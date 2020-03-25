import numpy as np


class GaloisField():
    """
    Class to represent a Galois Field.

    It's encapsulated all field and polynomial operations on fields.

    """

    def __init__(self, primitive_polynomial: bin, field_order: int,
                 generator: int = 2):
        """
        Parameters
        ----------
        primitive_polynomial: bin
            Primitive/prime (binary) polynomial and must be irreducible (ie, it
            can't represented as the product of two smaller polynomials). It'se
            a polynomial in the binary sense: each bit is a coefficient,but in
            fact it's an integer between max_field_value+1 and max_field_value
            *2 and not a list of gf values. The prime polynomial will be used
            to reduce the overflows back into the range of the Galois Field wi-
            thout duplicating values (all values should be unique).

        generator: int, default=2
            Generator number (the "increment" that will be used to walk through
            the field by multiplication, this must be a prime number). This is
            basically the base of the logarithm/anti-log tables. Also often no-
            ted α (alpha) in academic books.

        field_order: int
            Order of the field, i.e 2^p.

        look_up_tables
            Pre-computes the log and exp table for generated field.
        """

        self.primitive_polynomial = primitive_polynomial
        self.field_order = field_order
        self.generator = generator
        self.look_up_tables = self.look_up_tables(primitive_polynomial,
                                                  field_order,
                                                  generator)

    def primes_3_to_n(self, n: int) -> np.array:
        """Returns a array of primes, 3 <= p < n (Faster & more memory-wise num
        py code).

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

    def prime_polynomials(self, field_order: int, generator: int = 2,
                          fast_primes: bool = False,
                          first_prime: bool = False) -> list:
        """Compute the list of prime polynomials for the given generator and ga
        lois field prime exponent. Bruteforce approach.

        Parameters
        ----------
        field_order: int
            Order of the field, i.e 2^p, in which the search for prime po-
            lynomials will be done.

        generator: int, default=2
           Generator number (the "increment" that will be used to walk through
           the field by multiplication, this must be a prime number). This is
           basically the base of the logarithm/anti-log tables. Also often no-
           ted α (alpha) in academic books.

        fast_primes: bool, optional
            Will output less results but will be significantly faster.

        first_prime: bool, optional
            will output the first prime polynomial found, so if all you want is
            to just find one prime polynomial to generate the LUT for Reed-Solo
            mon to work, then just use that.

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
            prim_candidates = self.primes_3_to_n(next_field_max_value)
            prim_candidates = [prim_candidate for prim_candidate
                               in prim_candidates
                               if prim_candidate > max_field_value]
        else:
            prim_candidates = range(max_field_value+2, next_field_max_value,
                                    gf_2)

        correct_primes = list()
        for prim_candidate in prim_candidates:
            seen = np.zeros(max_field_value+1)
            conflict = False

            x = 1
            for i in range(max_field_value):
                x = self.gf_mult_noLUT(x, generator, prim_candidate,
                                       max_field_value+1)

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

    def look_up_tables(self, primitive_polynomial: bin, field_order: int,
                       generator: int = 2):
        """Precompute the logarithm and anti-log tables for faster computation
        later, using the provided primitive polynomial.

        These tables are used for multiplication/division. For each possible va
        lue in the galois field 2^field_order, we will pre-compute the logari-
        thm and anti-logarithm (exponential) of this value, to do that, we gene
        rate the Galois Field GF(2^p) by building a list starting with the ele
        ment 0 followed by the (p-1) succesive powers of the generator a : 1,
        a, a^1, a^2, ..., a^(p-1).

        Parameters
        ----------
        primitive_polynomial: bin
            Primitive/prime (binary) polynomial and must be irreducible (ie, it
            can't represented as the product of two smaller polynomials). It's
            a polynomial in the binary sense: each bit is a coefficient,but in
            fact it's an integer between max_field_value+1 and max_field_value
            *2 and not a list of gf values. The prime polynomial will be used
            to reduce the overflows back into the range of the Galois Field wi-
            thout duplicating values (all values should be unique).

        generator: int, default=2
            Generator number (the "increment" that will be used to walk through
            the field by multiplication, this must be a prime number). This is
            basically the base of the logarithm/anti-log tables. Also often no-
            ted α (alpha) in academic books.

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
            x = self.gf_mult_noLUT(x, generator, max_field_value+1,
                                   primitive_polynomial)

        for i in range(max_field_value, max_field_value * 2):
            gf_exp[i] = gf_exp[i - max_field_value]

        return [gf_log, gf_exp]

    def gf_add(self, x, y):
        return x ^ y

    def gf_sub(self, x, y):
        return x ^ y

    def gf_neg(self, x):
        return x

    def gf_mul(self, x, y):
        if x == 0 or y == 0:
            return 0
        return gf_exp[(gf_log[x] + gf_log[y]) % max_field_value]

    def gf_div(self, x, y):
        if y == 0:
            raise ZeroDivisionError()
        if x == 0:
            return 0
        return gf_exp[(gf_log[x] + max_field_value - gf_log[y])
                      % max_field_value]

    def gf_pow(self, x, power):
        return gf_exp[(gf_log[x] * power) % max_field_value]

    def gf_inverse(self, x):
        # gf_inverse(x) == gf_div(1, x)
        return gf_exp[max_field_value - gf_log[x]]

    def gf_mult_noLUT(self, x: bin, y: bin, field_order: int, poly_modular=0,
                      carryless=True):
        """Galois Field integer multiplication using Russian Peasant Multiplica
        tion algorithm (faster than the standard multiplication + modular reduc
        tion).

        If poly_modular is 0 and carryless=False, then the function produces
        the result for a standard integers multiplication (no carry-less arith-
        metics nor modular reduction).

        Parameters
        ---------
        x: bin
            Polynomial in binary form, must be in the same Galois Field as y.

        y: bin
            Polynomial in binary form, must be in the same Galois Field as x.

        field_order: int
            Order of the field, i.e 2^p.

        poly_modular: bin
            Irreducible primitive polynomial in binary form for modular reducti
            on after the multiplication between x and y. To stay in the bound
            of the field, we need to compute the modulo of any value above the
            Galois Field.

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

    def gf_poly_scale(self, polynomial: list, scalar: int):
        """Multiplies a polynomial by a scalar.

        In Python, polynomials will be represented by a list of numbers in des-
        cending order of powers of x.

        Parameters
        ----------
        polynomial: list
            Polynomial coefficients, in descending order of powers of x.

        scalar: int
            Degree of the scalar, i.e y ^ x.

        Returns
        -------
        list
            Coefficients in descending order of powers of x times the scalar.

        Raises
        ------
        None

        """
        return [self.gf_mul(polynomial[i], scalar) for i in
                range(len(polynomial))]

    def gf_poly_add(self, p: list, q: list):
        """This function "adds" two polynomials (using exclusive-or).

        Parameters
        ----------
        p: list
            Polynomial coefficients, in descending order of powers of x.

        q: list
            Polynomial coefficients, in descending order of powers of x.

        Returns
        -------
        r: list
            The result of p + q.
        """
        r = np.zeros(max(len(p), len(q)))
        for i in range(len(p)):
            r[i + len(r) - len(p)] = p[i]
        for i in range(len(q)):
            r[i + len(r) - len(q)] ^= q[i]
        return r

    def gf_poly_mul(self, p, q):
        """Multiply two polynomials, inside Galois Field (but the procedure is
        generic). Optimized function by precomputation of log.

        Compute the polynomial multiplication (just like the outer product of
        two vectors, we multiply each coefficients of p with all coefficients
        of q).

        Parameters
        ----------
        p: list
            Polynomial coefficients, in descending order of powers of x.

        q: list
            Polynomial coefficients, in descending order of powers of x.

        Returns
        -------
        r: list
            The result of p * q.
        """
        r = [0] * (len(p) + len(q) - 1)
        lp = [gf_log[p[i]] for i in range(len(p))]
        for j in range(len(q)):
            qj = q[j]  # optimization: load the coefficient once
            if qj != 0:  # log(0) is undefined, we need to check that
                lq = gf_log[qj]  # Optimization: precache the logarithm of the current coefficient of q
                for i in range(len(p)):
                    if p[i] != 0:  # log(0) is undefined, need to check that...
                        r[i + j] ^= gf_exp[lp[i] + lq]
        return r

    def gf_poly_neg(self, poly):
        """Returns the polynomial with all coefficients negated. In GF(2^p),
        negation does not change the coefficient, so we return the polynomial
        as-is."""
        return poly

    def gf_poly_div(self, dividend: list, divisor: list):
        """Fast polynomial division by using Extended Synthetic Division and op
        timized for GF(2^p) computations (doesn't work with standard polynomia
        ls outside of this galois field,see the Wikipedia article for generic
        algorithm).

        CAUTION: this function expects polynomials to follow the opposite con-
        vention at decoding: the terms must go from the biggest to lowest de-
        gree (while most other functions here expect a list from lowest to big-
        gest degree). eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5].

        The resulting msg_out contains both the quotient and the remainder, the
        remainder being the size of the divisor (the remainder has necessarily
        the same degree as the divisor -- not length but degree == length-1 --
        since it's what we couldn't divide from the dividend), so we compute
        the index where this separation is, and return the quotient and remainder.

        Parameters
        ----------
        dividend: list

        divisor: list


        Returns
        -------
        msg_out: list
            quotient, remainder.
        """

        # Copy the dividend list and pad with 0 where the ecc bytes will be computed
        msg_out = list(dividend)
        for i in range(len(dividend) - (len(divisor)-1)):
            coef = msg_out[i]  # precaching
            if coef != 0:
                for j in range(1, len(divisor)):
                    if divisor[j] != 0:
                        # xoring directly is faster): msg_out[i + j] += -divisor[j] * coef
                        msg_out[i + j] ^= self.gf_mul(divisor[j], coef)
        separator = -(len(divisor)-1)
        return msg_out[:separator], msg_out[separator:]

    def gf_poly_eval(self, poly, x):
        """Evaluates a polynomial in GF(2^p) given the value for x. This is ba-
        sed on Horner's scheme for maximum efficiency.

        Parameters
        ----------

        Returns
        -------
        """
        y = poly[0]
        for i in range(1, len(poly)):
            y = self.gf_mul(y, x) ^ poly[i]
        return y


if __name__ == "__main__":
    gf_16 = GaloisField(0x11D, 8)
    print(len(gf_16.look_up_tables[0]))
    print(gf_16.gf_poly_scale([2, 3, 5], 2))
