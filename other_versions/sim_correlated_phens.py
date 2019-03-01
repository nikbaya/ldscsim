#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 11:27:04 2019

Generate correlated betas

@author: nbaya
"""

import hail as hl
import numpy as np
import scipy.stats as stats
import pandas as pd


"""
*******************************
n correlated phenotype solution
*******************************
"""
def sim_correlated_phens(mt, genotype, h2_ls=[], rg_ls=[], cov_array = None, seed = None):
    if type(h2_ls) is not list or type(rg_ls) is not list:
        raise ValueError("Both h2_ls and rg_ls must be lists")
    if (len(h2_ls) == 0 or len(h2_ls) == 0) and cov_array is None:
        raise ValueError("This method requires h2 and rg values")
    elif (len(h2_ls) > 0 or len(h2_ls) > 0) and cov_array is not None:
        raise ValueError("This method only requires h2_ls and rg_ls or cov_array, not both")
    if cov_array is not None: 
        if np.asarray(cov_array).shape[0] != np.asarray(cov_array).shape[1]:
            raise ValueError("cov_array must be a square matrix")
    if seed is None:
        import os
        random_data = os.urandom(4)
        seed= int.from_bytes(random_data, byteorder="big") 
        print(f'Using random seed: {seed}')
    else:
        print(f'Using seed: {seed}')
    if cov_array is None:
        cov_array = create_cov_array(h2_ls, rg_ls)
    mt = mt.annotate_entries(__gt = genotype)
    mt = mt.annotate_globals(__seed = seed)
    mt = normalize_genotypes(mt)
    mt = create_correlated_betas(mt, cov_array, seed)
    return sim_phenotypes(mt, cov_array)

def create_cov_array(h2_ls, rg_ls):
    n_rg = len(rg_ls)
    n_h2 = len(h2_ls)
    exp_n_rg = int((n_h2**2-n_h2)/2)
    if n_rg is not exp_n_rg:
        raise ValueError(f'The number of rg values given is {n_rg}, expected {exp_n_rg}')
    rg_ls = np.asarray(rg_ls)
    cov_array = np.zeros(shape=[n_h2,n_h2])
    cov_array[np.triu_indices(n_h2, k=1)] = rg_ls**2 #sets upper diagonal with covariances
    h2_diag = np.diag(h2_ls)
    cov_array = (np.matmul(np.matmul(cov_array,h2_diag).T,h2_diag)**(1/2))
    cov_array += cov_array.T + h2_diag
    return cov_array

def normalize_genotypes(mt):
    mt1 = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__gt))
    return mt1.annotate_entries(__norm_gt = (mt1.__gt-mt1.__gt_stats.mean)/mt1.__gt_stats.stdev)

def create_correlated_betas(mt, cov_array, seed):
    M = mt.count_rows()
    n_phens = cov_array.shape[0]
    cov_array = (1/M)*cov_array
    randstate = np.random.RandomState(int(seed)) #seed random state for replicability
    if n_phens > 1:
        betas = randstate.multivariate_normal(mean=np.zeros(n_phens),cov=cov_array,size=[M,])
    df = pd.DataFrame(betas,columns=[f'__beta_{i}' for i in range(n_phens)])
    tb = hl.Table.from_pandas(df)
    tb = tb.add_index().key_by('idx')
    mt = mt.add_row_index()
    for i in range(n_phens):
        mt = mt._annotate_all(row_exprs = {f'__beta_{i}': tb[mt.row_idx][f'__beta_{i}']}) 
    return mt

def sim_phenotypes(mt, cov_array):
    h2_ls = np.diag(cov_array)
    n_phens = len(h2_ls)
    for i in range(n_phens):
        mt = mt._annotate_all(col_exprs={f'__y_no_noise_{i}': hl.agg.sum(mt[f'__beta_{i}'] * mt.__norm_gt)})
    for i in range(n_phens):
        mt = mt._annotate_all(col_exprs={f'__y_{i}':mt[f'__y_no_noise_{i}']+hl.rand_norm(0,hl.sqrt(1-h2_ls[i]))})
    return mt
