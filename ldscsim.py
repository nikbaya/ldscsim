#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:25:59 2019

ldsc simulation framework

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
def annotate_w_temp_fields(mt, genotype, h2, pi=1, annot=None, popstrat=None, popstrat_s2=1):
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
        annot_sum = mt1.aggregate_rows(hl.agg.sum(mt1.__annot))
        return mt1.annotate_rows(__beta = hl.rand_norm(0, hl.sqrt(mt1.__annot/annot_sum*h2)))
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
    mt3 = mt2.annotate_cols(__y_w_popstrat = mt2.__y + mt2.__norm_popstrat*hl.sqrt(mt2.__popstrat_s2))
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
                          expr_float32, 
                          expr_float64),
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
