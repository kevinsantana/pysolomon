#!/usr/bin/env python
# -*- coding: utf-8 -*-

###################################################
### Universal Reed-Solomon Codec
### initially released at Wikiversity
###################################################

################### INIT and stuff ###################

try: # compatibility with Python 3+
    xrange
except NameError:
    xrange = range

class ReedSolomonError(Exception):
    pass

gf_exp = [0] * 131072 # For efficiency, gf_exp[] has size 2*GF_SIZE, so that a simple multiplication of two numbers can be resolved without calling % 255. For more infos on how to generate this extended exponentiation table, see paper: "Fast software implementation of finite field operations", Cheng Huang and Lihao Xu, Washington University in St. Louis, Tech. Rep (2003).
gf_log = [0] * 65536
field_charac = int(2**16 - 1)