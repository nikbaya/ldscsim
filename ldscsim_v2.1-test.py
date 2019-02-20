#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:25:59 2019

Version 2.1 of ldsc simulation framework.

Pros:   Each method can be run independently of the simulate() method.
Cons:   Adds (unnecessary?) code.
        Requires re-annotating field which adds to runtime.

@author: nbaya
"""

import hail as hl
from hail.expr.expressions import *
from hail.typecheck import *
from hail.matrixtable import MatrixTable
import numpy as np
from datetime import datetime, timedelta

@typecheck(h2=oneof(float,int),
           pi=oneof(float,int),
           annot=oneof(nullable(expr_int32),
                       nullable(expr_float64)),
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,int,expr_int32,expr_int64,expr_float32,expr_float64),
           path_to_save=nullable(str))
def print_header(h2, pi, annot, popstrat, popstrat_s2, path_to_save):
    '''Makes the header for the simulation'''
    header =  '\n****************************************\n'
    header += 'Running simulation framework\n'
    header += 'h2 = {}\n'.format(h2)
    header += 'pi = {} (default: 1)\n'.format(pi)
    header += 'Annotation-informed betas?: {}\n'.format('NO' if annot is None else 'YES')
    header += 'Add population stratification?: {}\n'.format('NO' if popstrat is None else 'YES')    
    header += '' if popstrat is None else 'Population stratification scaling constant = {} (default: 1)\n'.format(popstrat_s2)
    header += '' if path_to_save is None else 'Saving to: {}\n'.format(path_to_save)
    header += '****************************************'
    print(header)

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, 
                          expr_float64),
           h2=oneof(float,int),
           pi=oneof(float,int),
           annot=oneof(nullable(expr_int32),
                       nullable(expr_float64)),
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64))
def annotate_w_temp_fields(mt, genotype, h2, pi, annot, popstrat, popstrat_s2):
    '''Annotate mt with temporary fields of simulation parameters'''
    if annot is None and popstrat is None: #Infinitesimal/SpikeSlab
        mt1 = mt._annotate_all(entry_exprs={'__gt_temp':genotype},
                               global_exprs={'__h2_temp':h2, 
                                             '__pi_temp':pi})
    elif annot is not None and popstrat is None: #AnnotationInformed
        mt1 = mt._annotate_all(row_exprs={'__annot_temp':annot},
                               entry_exprs={'__gt_temp':genotype},
                               global_exprs={'__h2_temp':h2})
    elif annot is None and popstrat is not None: #Infinitesimal/SpikeSlab + popstrat
        mt1 = mt._annotate_all(col_exprs={'__popstrat_temp':popstrat},
                               entry_exprs={'__gt_temp':genotype},
                               global_exprs={'__h2_temp':h2,
                                             '__pi_temp':pi,
                                             '__popstrat_s2_temp':popstrat_s2})
    elif annot is not None and popstrat is not None: #AnnotationInformed + popstrat
        mt1 = mt._annotate_all(row_exprs={'__annot_temp':annot},
                               col_exprs={'__popstrat_temp':popstrat},
                               entry_exprs={'__gt_temp':genotype},
                               global_exprs={'__h2_temp':h2,
                                             '__popstrat_s2_temp':popstrat_s2})
    return mt1

@typecheck(mt=MatrixTable, 
           h2=oneof(float,int),
           pi=oneof(float,int),
           annot=oneof(nullable(expr_int32),
                       nullable(expr_float64)))
def make_betas(mt, h2, pi=1, annot=None):
    '''Simulate betas. Options: Infinitesimal model, spike & slab, annotation-informed'''        
    M = mt.count_rows()
    if annot is not None:
        print('\rSimulating annotation-informed betas w/ h2 = {}'.format(h2))
        mt1 = mt._annotate_all(row_exprs={'__annot':annot},
                              global_exprs={'__h2':h2})
        annot_stats = mt1.aggregate_rows(hl.agg.stats(mt1.__annot), _localize=True)
        return mt1.annotate_rows(__beta = (mt1.__annot - annot_stats.mean) / annot_stats.stdev * hl.sqrt(h2 / M))
    else:
        print('Simulating betas using {} model w/ h2 = {}'.format(('infinitesimal' if pi is 1 else 'spike & slab'),h2))
        mt1 = mt.annotate_globals(__h2 = h2, __pi = pi)
        return mt1.annotate_rows(__beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(h2/(M*pi))))

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, 
                          expr_float64))
def normalize_genotypes(mt, genotype):
    '''Normalizes genotypes'''
    print('\rNormalizing genotypes...')
    mt1 = mt.annotate_entries(__gt = genotype)
    mt2 = mt1.annotate_rows(__gt_stats = hl.agg.stats(mt1.__gt))
    return mt2.annotate_entries(__norm_gt = (mt2.__gt-mt2.__gt_stats.mean)/mt2.__gt_stats.stdev)    

@typecheck(mt=MatrixTable, 
           y=oneof(expr_int32,expr_float64),
           popstrat=oneof(nullable(expr_int32),nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64))
def add_pop_strat(mt, y, popstrat, popstrat_s2=1):
    '''Adds popstrat to a phenotype'''
    print('\rAdding population stratification (popstrat_s2 = {})'.format(popstrat_s2))
    mt1 = mt.annotate_cols(__y = y, 
                           __popstrat = popstrat)
    mt1 = mt1.annotate_globals(__popstrat_s2 = popstrat_s2)
    popstrat_stats = mt1.aggregate_cols(hl.agg.stats(mt1.__popstrat), _localize=True)
    mt2 = mt1.annotate_cols(__norm_popstrat = (mt1.__popstrat-popstrat_stats.mean)/popstrat_stats.stdev)
    mt3 = mt2.annotate_cols(__y_w_popstrat = mt2.__y + mt2.__norm_popstrat*mt2.__popstrat_s2)
    y_w_popstrat_stats = mt3.aggregate_cols(hl.agg.stats(mt3.__y_w_popstrat))
    return mt3.annotate_cols(__y_w_popstrat = (mt3.__y_w_popstrat-y_w_popstrat_stats.mean)/y_w_popstrat_stats.stdev)

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32, 
                          expr_int64, 
                          expr_float32, 
                          expr_float64),
           h2=oneof(float,int),
           beta=expr_float64,
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64))
def sim_phenotypes(mt, genotype, h2, beta, popstrat=None, popstrat_s2=1):
    '''Simulate phenotypes given betas and genotypes. Adding population stratification is optional'''
    print('\rCalculating phenotypes{}...'.format('' if popstrat is None else ' w/ population stratification').ljust(81))
    if popstrat is None:
        mt1 = mt._annotate_all(row_exprs={'__beta':beta},
                               entry_exprs={'__gt':genotype})
    else:
        mt1 = mt._annotate_all(row_exprs={'__beta':beta},
                               col_exprs={'__popstrat':popstrat},
                               entry_exprs={'__gt':genotype},
                               global_exprs={'__popstrat_s2':popstrat_s2})

    mt2 = normalize_genotypes(mt1,
                              mt1.__gt)
    mt3 = mt2.annotate_cols(__y_no_noise = hl.agg.sum(mt2.__beta * mt2.__norm_gt))
    mt4 = mt3.annotate_cols(__y = mt3.__y_no_noise + hl.rand_norm(0,hl.sqrt(1-h2)))

    if popstrat is None:
        return mt4        
    else:
        return add_pop_strat(mt4, 
                             y=mt4.__y, 
                             popstrat=mt4.__popstrat, 
                             popstrat_s2=hl.eval(mt4.__popstrat_s2))

@typecheck(mt = MatrixTable, 
           str_expr=str)
def clean_fields(mt, str_expr):
    '''Removes fields with names that have str_expr in them'''
    all_fields = list(mt.col)+list(mt.row)+list(mt.entry)+list(mt.globals)
    for field_to_drop in [x for x in all_fields if str_expr in x]:
        mt = mt.drop(field_to_drop)
    return mt

@typecheck(mt=MatrixTable, 
           h2=oneof(float,int),
           starttime=datetime,
           stoptime=datetime,
           runtime=timedelta,
           pi=oneof(float,int),
           annot=oneof(nullable(expr_int32),
                       nullable(expr_float64)),
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64),
           path_to_save=nullable(str))
def add_sim_description(mt,h2,starttime,stoptime,runtime,pi=1,annot=None,popstrat=None,popstrat_s2=1,path_to_save=None):
    '''Annotates mt with description of simulation'''
    sim_id = 0
    while (str(sim_id) in [x.strip('sim_desc') for x in list(mt.globals) if 'sim_desc' in x]):
        sim_id += 1
    sim_desc = 'h2={h2}\npi={pi}\nis_annot={is_annot}\nis_popstrat={is_popstrat}'.format(
                h2=h2,pi=pi,is_annot=(annot is not None),is_popstrat=(annot is not None))
    sim_desc += '' if popstrat is None else '\npopstrat_s2={}'.format(popstrat_s2)
    sim_desc += '\nstarttime:{}\nstoptime: {}\nruntime:  {}\n'.format(starttime,stoptime,runtime)
    sim_desc += 'Saved to: {}'.format(path_to_save)
    mt = mt._annotate_all(global_exprs={'sim_desc{}'.format(sim_id):sim_desc})
    return mt

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, expr_float64),
           h2=oneof(float,int),
           pi=oneof(float,int),
           annot=oneof(nullable(expr_int32),
                       nullable(expr_float64)),
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,int),
           path_to_save=nullable(str))
def simulate(mt, genotype, h2, pi=1, annot=None, popstrat=None, popstrat_s2 = 1,path_to_save=None):
    ''' Simulates phenotypes. 
        Options: Infinitesimal, spike/slab, annotation-informed, population stratification'''
    print_header(h2=h2, 
                 pi=pi, 
                 annot=annot, 
                 popstrat=popstrat, 
                 popstrat_s2=popstrat_s2, 
                 path_to_save=path_to_save)
    
    starttime = datetime.now()
    
    mt1 = annotate_w_temp_fields(mt=mt, 
                                 genotype=genotype,
                                 h2=h2, 
                                 pi=pi, 
                                 annot=annot, 
                                 popstrat=popstrat, 
                                 popstrat_s2=popstrat_s2)
    
    mt2 = make_betas(mt=mt1, 
                     h2=h2, 
                     pi=pi, 
                     annot=None if annot is None else mt1.__annot_temp)
        
    mt2 = mt2.rename({'__beta':'__beta_temp'})

    if popstrat is None:
        mt3 = sim_phenotypes(mt=mt2, 
                             genotype=mt2.__gt_temp, 
                             h2=h2, 
                             beta=mt2.__beta_temp)
    if popstrat is not None:
        mt3 =  sim_phenotypes(mt=mt2, 
                              genotype=mt2.__gt_temp, 
                              h2=h2, 
                              beta=mt2.__beta_temp,
                              popstrat=mt2.__popstrat_temp, 
                              popstrat_s2=mt2.__popstrat_s2_temp)
        
    mt4 = clean_fields(mt3, '_temp')
    stoptime = datetime.now()
    runtime = stoptime-starttime
    mt5 = add_sim_description(mt4,h2,starttime,stoptime,runtime,pi,annot,popstrat,popstrat_s2,path_to_save)
    print('\rFinished simulation! (runtime={} min)'.format(np.round((runtime.seconds+runtime.microseconds/1e6)/60, 4)).ljust(100))
    if path_to_save is not None:
        print('\rWriting simulation to: {}'.format(path_to_save))
        mt5.write(path_to_save)
    return mt5

def test_suite(mt, genotype, popstrat):
    '''Testing suite for simulation framework'''
    mt = mt._annotate_all(row_exprs={'a1': hl.rand_norm(),
                                     'a2': hl.rand_bool(0.1)},
                          col_exprs={'popstrat':popstrat},
                          entry_exprs={'gt':genotype})
    mt = mt.annotate_rows(annot = mt.a1+mt.a2)
    
    n_sim = 7 #number of simulations
    sim_h2_ls = np.round(np.random.uniform(low=0,high=1,size=n_sim),4)
    obs_h2_ls = []
    sim_mt_ls = []
    
    # Infinitesimal
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[0], genotype = mt.gt))
    # Spike & slab
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[1], pi=0.1, genotype = mt.gt))
    # Annotation-informed
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[1], genotype = mt.gt, annot = mt.annot)) #has same h2 as previous spike and slab to check if sims match
    # Infinitesimal + population stratification, popstrat_s2 = 0.5
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[3], genotype=mt.gt, popstrat=mt.popstrat, popstrat_s2 = 0.5))
    # Infinitesimal + population stratification, popstrat_s2 = 0.25
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[3], genotype=mt.gt, popstrat=mt.popstrat, popstrat_s2 = 0.25))
    # Spike & slab + population stratification
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[5], pi=0.1, genotype=mt.gt, popstrat=mt.popstrat))
    # Annotation-informed + population stratification
    sim_mt_ls.append(simulate(mt=mt, h2=sim_h2_ls[6], genotype=mt.gt, annot = mt.annot, popstrat=mt.popstrat))
    
    for sim_mt in sim_mt_ls:
        print(sim_mt.describe())
    
    for sim_i,sim_mt in enumerate(sim_mt_ls):
        obs_h2_ls.append(np.round(sim_mt.aggregate_cols(hl.agg.stats(sim_mt['__y_no_noise']).stdev**2),4))
    print('\nExpected h2s: {} \nObserved h2s: {}'.format(sim_h2_ls, obs_h2_ls)) 

def gwas(mt, x, y, cov_list=[], with_intercept=True, pass_through=[], path_to_save=None, 
         normalize_x=True, is_std_cov_list=True):
    '''Runs GWAS'''
    
    mt = mt._annotate_all(col_exprs={'y':y},
                           entry_exprs={'x':x})
    if normalize_x:
        mt = normalize_genotypes(mt, mt.x) 
        mt = mt.annotate_entries(x = mt.__norm_gt).drop('__norm_gt')
    
    if is_std_cov_list:
        cov_list = ['isFemale','age','age_squared','age_isFemale',
                    'age_squared_isFemale']+['PC{:}'.format(i) for i in range(1, 21)]
        
    if str in list(map(lambda x: type(x),cov_list)):
        cov_list = list(map(lambda x: mt[x] if type(x) is str else x,cov_list))
        
    cov_list = ([1] if with_intercept else [])+cov_list
    
    print(pass_through)

    gwas_ht = hl.linear_regression_rows(y=mt.y,
                                        x=mt.x,
                                        covariates=cov_list,
                                        pass_through = ['rsid']+pass_through)
    gwas_ht.describe()
    
    gwas_ht = gwas_ht.annotate_globals(with_intercept = with_intercept)
    
    gwas_ht.describe()
    
    gwas_ht = gwas_ht.rename({'rsid':'SNP'}).key_by('SNP')
        
    gwas_ht = gwas_ht.select(Z = gwas_ht.t_stat,
                             N = gwas_ht.n)
    
    sumstats_template = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint64})
    sumstats_template = sumstats_template.key_by('SNP')
    
    gwas_ht.describe()
    
    sumstats = sumstats_template.annotate(Z = gwas_ht[sumstats_template.SNP].Z,
                                          N = gwas_ht[sumstats_template.SNP].N)
    
    if path_to_save is not None:
        sumstats.export(path_to_save)
        
    return gwas_ht
    
#if __name__ == '__main__':    
#    gs_bucket = 'gs://nbaya/ldscsim/'
#    mt0 = hl.read_matrix_table(gs_bucket+'hm3.50_sim_h2_0.08.mt/')
##    sim_mt = simulate(mt=mt0,
##                      h2=0.5,
##                      genotype=mt0.dosage,
##                      popstrat=mt0.PC1,
##                      popstrat_s2=3,
##                      path_to_save=None)
##    sim_mt.cols().write(gs_bucket+'hm3.sim_h2_0.5.w_popstrat_PC1.cols.ht',overwrite=True)
#    cols_ht = hl.read_table(gs_bucket+'hm3.sim_h2_0.5.w_popstrat_PC1.cols.ht')
#    cols_ht.describe()
#    sim_mt = mt0.annotate_cols(__y = cols_ht[mt0.s].__y,
#                               __y_w_popstrat = cols_ht[mt0.s].__y_w_popstrat)
#    gwas(mt=sim_mt,
#         x=sim_mt.dosage,
#         y=sim_mt.__y,
#         pass_through = [],
#         normalize_x=False,
#         path_to_save=gs_bucket+'hm3.sim_h2_0.5.no_popstrat.sumstats.tsv.bgz')
#    gwas(mt=sim_mt,
#         x=sim_mt.dosage,
#         y=sim_mt.__y_w_popstrat,
#         pass_through = [],
#         normalize_x=False,
#         path_to_save=gs_bucket+'hm3.sim_h2_0.5.w_popstrat_PC1.popstrat_s2_3.sumstats.tsv.bgz')

if __name__ == '__main__':    
    gs_bucket = 'gs://nbaya/ldscsim/'
    mt0 = hl.read_matrix_table(gs_bucket+'hm3.50_sim_h2_0.08.mt/')
    ht0 = mt0.cols()
    mt  = hl.read_matrix_table(gs_bucket+'ukbb31063.GT.chr22.MAF_0.05.white_british.mt/')
    mt = mt.annotate_cols(PC1 = ht0[mt.s]['PC1'])
    test_suite(mt=mt, genotype=mt.GT.n_alt_alleles(), popstrat=mt.PC1)
