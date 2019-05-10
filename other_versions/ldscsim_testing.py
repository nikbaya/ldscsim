#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 17:17:25 2019

Testing suite for ldscsim

@author: nbaya
"""

import hail as hl
from hail import dtype
import numpy as np
from ldscsim import *
import scipy.stats as stats
import numpy.random as rand


def get_sim_h2(mt, k):
    if mt.y_no_noise.dtype is dtype('array<float64>'):        
        return list(map(lambda i: mt.aggregate_cols(hl.agg.stats(mt.y_no_noise[i])).stdev**2, range(k)))
    else:
        return [mt.aggregate_cols(hl.agg.stats(mt.y_no_noise)).stdev**2]

def get_sim_rg(mt, k):
    rg = [0]*int(k*(k-1)/2)
    ii = 0
    M = mt.count_rows()
    beta = np.asarray(mt.beta.take(M))
    for i in range(k):
        for j in range(i+1, k):
            rg[ii] = stats.pearsonr(beta[:,i],beta[:,j])[0]
            ii += 1
    return rg
            


def test_suite(mt, genotype, popstrat):
    
    
    
def single_trait(mt, genotype, eps):
    """
    ╔═══════════════╗
    ║ Infinitesimal ║
    ╚═══════════════╝
    """
    # Test type errors
    for h2 in [None, [None], 1+1j, genotype]:
        try:
            sim = simulate_phenotypes(mt, genotype, h2=h2)
            assert False, f'did not raise error when h2={h2}'
        except TypeError:
            pass
    
    # Test assertion errors
    for h2 in [-rand.exponential(), rand.exponential()+1]:
        try:
            sim = simulate_phenotypes(mt, genotype, h2=h2)
            assert False, f'did not raise error when h2={h2}'
        except AssertionError:
            pass
        
    try:
        h2 = rand.uniform(size=(1))
        rg = rand.uniform(size=(2))
        sim = simulate_phenotypes(mt,genotype,h2=h2,rg=rg)

    
    # Test functionality
    for h2 in [0,1, [0], [1]]:
        sim = simulate_phenotypes(mt, genotype, h2=h2)
        assert h2 == get_sim_h2(sim, 1)[0]
        
    for h2 in [rand.uniform().tolist(),
               rand.uniform()]:
        sim = simulate_phenotypes(mt, genotype, h2=h2)
        assert abs(h2 < get_sim_h2(sim)) < eps

    for h2 in [rand.uniform(size=2), rand.uniform(size=3), rand.uniform(size=10),
               rand.uniform(size=50), rand.uniform(size=100)]:
        for rg in [None, [None], ]:
            sim = simulate_phenotypes(mt,genotype,h2=h2,rg=rg)
    
    sim = simulate_phenotypes(mt,genotype,h2=[0.1,0.3,0.43],rg=[0.3, 0.4, 0.1])
    
    
