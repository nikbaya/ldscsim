#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: nbaya
"""

import hail as hl

#
#class BetaDistribution(a, b):
#    pass

class InfinitesimalModel():
    def __init__(self):
        pass
    
    def make_random_function(self, mt, h2):
        M = mt.count_rows() # number of variants
        return hl.rand_norm(0, h2/M)
    
    def make_noise(self, mt, h2):
        return hl.rand_norm(0, 1-h2)
    
class InfinitesimalSpikeSlabModel():
    def __init__(self):
        pass
    
    def make_random_function(self, mt, h2, pi): #pi is slab prob
        M = mt.count_rows() # number of variants
        if (hl.rand_unif(0,1) < pi):
            return hl.rand_norm(0,h2/(M*pi))
        else:
            return 0
        
    def make_noise(self, mt, h2):
        return hl.rand_norm(0, 1-h2)

def simulate(mt, beta_distribution, eps_distribution, h2):
    beta_dist = beta_distribution.make_random_function(mt, h2)
    epsilon_dist  = beta_distribution.make_noise(mt, h2)
    
    mt = mt.annotate_rows(beta = beta_dist)
    
    mt = mt.annotate_cols(y_no_noise = hl.agg.sum(mt.beta * mt.GT.n_alt_alleles()))
    mt = mt.annotate_cols(y = epsilon_dist + mt.y_no_noise)
    
    return mt


import pandas as pd
tb = pd.DataFrame([0,1,1,0,0,0,2,1])
tb = hl.Table.from_pandas(tb)
tb = tb.key_by('0')
mt = hl.MatrixTable.from_rows_table(tb)

b = InfinitesimalSpikeSlabModel()
b.make_random_function(mt, 0.1, 0.5)