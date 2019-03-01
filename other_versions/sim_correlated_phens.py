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
def sim_correlated_phens(mt, genotype, h2_ls=[], rg_ls=[], cov_array = None):
    if type(h2_ls) is not list or type(rg_ls) is not list:
        raise ValueError("Both h2_ls and rg_ls must be lists")
    if (len(h2_ls) == 0 or len(h2_ls) == 0) and cov_array is None:
        raise ValueError("This method requires h2 and rg values")
    elif (len(h2_ls) > 0 or len(h2_ls) > 0) and cov_array is not None:
        raise ValueError("This method only requires h2_ls and rg_ls or cov_array, not both")
    if cov_array is not None: 
        if np.asarray(cov_array).shape[0] != np.asarray(cov_array).shape[1]:
            raise ValueError("cov_array must be a square matrix")
    if cov_array is None:
        cov_array = create_cov_array(h2_ls, rg_ls)
    mt = mt.annotate_entries(__gt = genotype)
    mt = normalize_genotypes(mt)
    mt = create_correlated_betas(mt, cov_array)
    return sim_phenotypes(mt, cov_array)

def create_cov_array(h2_ls, rg_ls):
    n_rg = len(rg_ls)
    n_h2 = len(h2_ls)
    exp_n_rg = int((n_h2**2-n_h2)/2) #expected number of rg values passed in rg_ls
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

def create_correlated_betas(mt, cov_array):
    M = mt.count_rows()
    n_phens = cov_array.shape[0]
    cov_array *= (1/M)
    if n_phens > 1:
        betas = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_array,size=[M,])
    df = pd.DataFrame(betas,columns=[f'beta_{i}' for i in range(n_phens)])
    tb = hl.Table.from_pandas(df)
    tb = tb.add_index().key_by('idx')
    mt = mt.add_row_index()
    return mt.annotate_rows(betas = tb[mt.row_idx]) 
   

def sim_phenotypes(mt, cov_array):
    h2_ls = np.diagonal(cov_array)
    n_phens = len(h2_ls)
    for beta_i in range(n_phens):
        mt = mt._annotate_all(col_exprs={f'__y_no_noise_{beta_i}': hl.agg.sum(mt.betas[f'beta_{beta_i}'] * mt.__norm_gt)})
    for beta_i in range(n_phens):
        mt = mt._annotate_all(col_exprs={f'__y_{beta_i}':mt[f'__y_no_noise_{beta_i}']+hl.rand_norm(0,hl.sqrt(1-h2_ls[beta_i]))})
    return mt
    
def check_sim_results(sim_mt, h2_ls, rg_ls):
    n_phen = len(h2_ls)
    obs_h2_ls = [0]*n_phen
    for phen_i in range(n_phen):
        obs_h2_ls[phen_i] =  sim_mt.aggregate_cols(hl.agg.stats(sim_mt['__y_no_noise_{}'.format(phen_i)]).stdev**2)

    betas = sim_mt.key_rows_by().select_rows('betas').rows().to_pandas().values
    corr = np.corrcoef((betas).T)    
    print('Expected h2s: {} \t Observed h2s: {}'.format(h2_ls, obs_h2_ls))
    print('Expected rgs: {} \t Observed rgs: {}'.format(rg_ls, corr[np.triu_indices(n_phen, k=1)]))

if __name__ == '__main__':
    gs_bucket = 'gs://nbaya/ldscsim/'
    mt = hl.read_matrix_table(gs_bucket+'hm3.50_sim_h2_0.08.mt/')
    
    mt.describe()
    
    h2_s = 0.08 #symptom phenotype
    h2_b1 = 0.20 # binary phenotype 1
    h2_b2 = 0.30 # binary phenotype 2
    
    rg_sb1 = 0.90 #symptom-binary phen 1 rg
    rg_sb2 = 0.15 #symptom-binary phen 2 rg
    rg_b1b2 = 0.12 #binary phen 1-binary phen 2 rg
    
    h2_ls = [h2_s, h2_b1, h2_b2]
    rg_ls = [rg_sb1, rg_sb2, rg_b1b2]
    #rg_ls = [0, 0, 0]
    cov_array = create_cov_array(h2_ls, rg_ls)
    
    sim_mt_new = sim_correlated_phens(mt, mt.dosage, cov_array=cov_array)
    
    
    sim_mt_new.cols().write('gs://nbaya/ldscsim/hm3.50_sim_corr_phenotypes_v1_cols.ht',overwrite=True)
    sim_mt_new.rows().write('gs://nbaya/ldscsim/hm3.50_sim_corr_phenotypes_v1_rows.ht',overwrite=True)
    #sim_mt_new.write('gs://nbaya/ldscsim/hm3.50_sim_corr_phenotypes_v1.mt',overwrite=True)
    
    


"""
*******************************
Two phenotype solution
*******************************
"""
##def normalize_genotypes(mt):
##    mt1 = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__gt))
##    return mt1.annotate_entries(__norm_gt = (mt1.__gt-mt1.__gt_stats.mean)/mt1.__gt_stats.stdev)
##
##def create_betas(mt, h2_A, h2_B, rg):
##    M = mt.count_rows()
##    cov_array = np.asarray([[h2_A/M, rg*(h2_A*h2_B)**(1/2)/M],[rg*(h2_A*h2_B)**(1/2)/M, h2_B/M]])    
##    betas = np.random.multivariate_normal(mean=np.zeros(cov_array.shape[0]),cov=cov_array,size=[M,])
##    print('\rExpected h2_A: '+str(h2_A)+'\tObserved h2_A: '+str(np.var(betas[:,0])*M))
##    print('\rExpected h2_B: '+str(h2_B)+'\tObserved h2_B: '+str(np.var(betas[:,1])*M))
##    print('\rExpected rg: '+str(rg)+'\tObserved rg: '+str(stats.pearsonr(betas[:,0],betas[:,1])[0]))
##    df = pd.DataFrame(betas,columns=['beta_A','beta_B'])
##    tb = hl.Table.from_pandas(df)
##    tb = tb.add_index().key_by('idx')
##    mt = mt.add_row_index()
##    return mt.annotate_rows(__beta_A = tb[mt.row_idx].beta_A,
##                            __beta_B = tb[mt.row_idx].beta_B)    
##    
##def sim_phenotypes(mt, h2_A, h2_B):
##    mt1 = mt.annotate_cols(__y_no_noise_A = hl.agg.sum(mt.__beta_A * mt.__norm_gt),
##                           __y_no_noise_B = hl.agg.sum(mt.__beta_B * mt.__norm_gt))
##    return mt1.annotate_cols(__y_A = mt1.__y_no_noise_A + hl.rand_norm(0,hl.sqrt(1-h2_A)),
##                             __y_B = mt1.__y_no_noise_B + hl.rand_norm(0,hl.sqrt(1-h2_B)))
##
##def sim_correlated_phenotypes(mt, genotype, h2_A, h2_B, rg):
##    mt = mt.annotate_entries(__gt = genotype)
##    mt = normalize_genotypes(mt)
##    mt = create_betas(mt, h2_A, h2_B, rg)
##    return  sim_phenotypes(mt, h2_A, h2_B)
#
#
#sim_mt1 = sim_correlated_phens(mt, mt.dosage, h2_A=h2_s, h2_B=h2_b1, rg=rg_sb1)
#sim_mt2 = sim_correlated_phens(mt, mt.dosage, h2_A=h2_s, h2_B=h2_b2, rg=rg_sb2)
#
##check variance of betas
#obs_h2_s1 =  sim_mt1.aggregate_cols(hl.agg.stats(sim_mt1.__y_no_noise_A).stdev**2)
#print('\n'+str(obs_h2_s1))
#obs_h2_b1 = sim_mt1.aggregate_cols(hl.agg.stats(sim_mt1.__y_no_noise_B).stdev**2)
#obs_h2_s2 = sim_mt2.aggregate_cols(hl.agg.stats(sim_mt2.__y_no_noise_A).stdev**2)
#obs_h2_b2 = sim_mt2.aggregate_cols(hl.agg.stats(sim_mt2.__y_no_noise_B).stdev**2)
#print('\n***** sim 1 *****')
#print('Expected h2: '+str(h2_s) +'\tObserved h2: '+str(obs_h2_s1))
#print('Expected h2: '+str(h2_b1)+'\tObserved h2: '+str(obs_h2_b1))
#print('***** sim 2 *****')
#print('Expected h2: '+str(h2_s) +'\tObserved h2: '+str(obs_h2_s2))
#print('Expected h2: '+str(h2_b2)+'\tObserved h2: '+str(obs_h2_b2))
#
##check mean and std of phenotypes
#phen_1A_stats = sim_mt1.aggregate_cols(hl.agg.stats(sim_mt1.__y_A))
#phen_1B_stats = sim_mt1.aggregate_cols(hl.agg.stats(sim_mt1.__y_B))
#phen_2A_stats = sim_mt2.aggregate_cols(hl.agg.stats(sim_mt2.__y_A))
#phen_2B_stats = sim_mt2.aggregate_cols(hl.agg.stats(sim_mt2.__y_B))
#print('***** sim 1 *****')
#print(phen_1A_stats)
#print(phen_1B_stats)
#print('***** sim 2 *****')
#print(phen_2A_stats)
#print(phen_2B_stats)
