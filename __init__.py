#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################################
### Universal Reed-Solomon Codec
### initially released at Wikiversity
###################################################

################### INIT and stuff ###################


import numpy as np


try:  # compatibility with Python 3+
    xrange
except NameError:
    xrange = range


class ReedSolomonError(Exception):
    pass

# For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple multiplication
# of two numbers can be resolved without calling % 255.


gf_log = np.zeros(65536, dtype=int)
gf_exp = np.zeros(len(gf_log)*2, dtype=int)
max_field_value = int(2**16 - 1)
