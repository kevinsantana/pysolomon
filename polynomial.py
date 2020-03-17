################### GALOIS FIELD POLYNOMIALS MATHS ###################

def gf_poly_scale(p, x):
    return [gf_mul(p[i], x) for i in xrange(len(p))]

def gf_poly_add(p,q):
    r = [0] * max(len(p),len(q))
    for i in xrange(len(p)):
        r[i+len(r)-len(p)] = p[i]
    for i in xrange(len(q)):
        r[i+len(r)-len(q)] ^= q[i]
    return r

def gf_poly_mul(p, q):
    '''Multiply two polynomials, inside Galois Field (but the procedure is generic). Optimized function by precomputation of log.'''
    # Pre-allocate the result array
    r = [0] * (len(p) + len(q) - 1)
    # Precompute the logarithm of p
    lp = [gf_log[p[i]] for i in xrange(len(p))]
    # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
    for j in xrange(len(q)):
        qj = q[j] # optimization: load the coefficient once
        if qj != 0: # log(0) is undefined, we need to check that
            lq = gf_log[qj] # Optimization: precache the logarithm of the current coefficient of q
            for i in xrange(len(p)):
                if p[i] != 0: # log(0) is undefined, need to check that...
                    r[i + j] ^= gf_exp[lp[i] + lq] # equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))
    return r

def gf_poly_mul_simple(p, q): # simple equivalent way of multiplying two polynomials without precomputation, but thus it's slower
    '''Multiply two polynomials, inside Galois Field'''
    # Pre-allocate the result array
    r = [0] * (len(p) + len(q) - 1)
    # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
    for j in xrange(len(q)):
        for i in xrange(len(p)):
            r[i + j] ^= gf_mul(p[i], q[j]) # equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j])) -- you can see it's your usual polynomial multiplication
    return r

def gf_poly_neg(poly):
    '''Returns the polynomial with all coefficients negated. In GF(2^p), negation does not change the coefficient, so we return the polynomial as-is.'''
    return poly

def gf_poly_div(dividend, divisor):
    '''Fast polynomial division by using Extended Synthetic Division and optimized for GF(2^p) computations
    (doesn't work with standard polynomials outside of this galois field, see the Wikipedia article for generic algorithm).'''
    # CAUTION: this function expects polynomials to follow the opposite convention at decoding:
    # the terms must go from the biggest to lowest degree (while most other functions here expect
    # a list from lowest to biggest degree). eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]

    msg_out = list(dividend) # Copy the dividend list and pad with 0 where the ecc bytes will be computed
    #normalizer = divisor[0] # precomputing for performance
    for i in xrange(len(dividend) - (len(divisor)-1)):
        #msg_out[i] /= normalizer # for general polynomial division (when polynomials are non-monic), the usual way of using
                                  # synthetic division is to divide the divisor g(x) with its leading coefficient, but not needed here.
        coef = msg_out[i] # precaching
        if coef != 0: # log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization).
            for j in xrange(1, len(divisor)): # in synthetic division, we always skip the first coefficient of the divisior,
                                              # because it's only used to normalize the dividend coefficient
                if divisor[j] != 0: # log(0) is undefined
                    msg_out[i + j] ^= gf_mul(divisor[j], coef) # equivalent to the more mathematically correct
                                                               # (but xoring directly is faster): msg_out[i + j] += -divisor[j] * coef

    # The resulting msg_out contains both the quotient and the remainder, the remainder being the size of the divisor
    # (the remainder has necessarily the same degree as the divisor -- not length but degree == length-1 -- since it's
    # what we couldn't divide from the dividend), so we compute the index where this separation is, and return the quotient and remainder.
    separator = -(len(divisor)-1)
    return msg_out[:separator], msg_out[separator:] # return quotient, remainder.

def gf_poly_eval(poly, x):
    '''Evaluates a polynomial in GF(2^p) given the value for x. This is based on Horner's scheme for maximum efficiency.'''
    y = poly[0]
    for i in xrange(1, len(poly)):
        y = gf_mul(y, x) ^ poly[i]
    return y


if __name__ == "__main__":
    pass