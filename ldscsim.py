#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:25:59 2019

ldsc simulation framework

@author: nbaya
"""

import hail as hl
from hail.expr.expressions import expr_int32,expr_int64,expr_float32,expr_float64
from hail.typecheck import typecheck, oneof, nullable
from hail.matrixtable import MatrixTable
from hail.table import Table
import re
import numpy as np
from datetime import datetime, timedelta

@typecheck(h2=oneof(float,int),
           pi=oneof(float,int),
           is_annot_inf=bool,
           h2_normalize=bool,
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64),
           path_to_save=nullable(str))
def print_header(h2, pi, is_annot_inf, h2_normalize, popstrat, popstrat_s2, path_to_save):
    '''Makes the header for the simulation'''
    header =  '\n****************************************\n'
    header += 'Running simulation framework\n'
    header += 'h2 = {}\n'.format(h2) if h2_normalize else ''
    header += 'pi = {} (default: 1)\n'.format(pi)
    header += 'Annotation-informed betas?: {}\n'.format('YES' if is_annot_inf else 'NO')
    header += 'h2-normalized betas?: {}\n'.format('YES' if h2_normalize else 'NO')
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

def make_tau_ref_dict():
    '''Make tau_ref_dict from tsv?/dataframe?/Hail Table?'''
    pass

@typecheck(mt=MatrixTable,
           annot_list=oneof(list,
                            dict),
           annot_pattern=str)
def add_annot_pattern(mt, annot_list, annot_pattern, prefix=True):
    '''Adds a given pattern to the names of row field annotations listed in a 
    list or from the keysof a dict. Helpful for searching for those annotations 
    in get_tau_dict() or agg_annotations(). If prefix=False, the pattern will be
    added as a suffix.'''
    if type(annot_list) == dict :
        annot_list = list(annot_list.keys)
    row_fields = set(mt.row)
    if set(annot_list) != row_fields:
        annots_in_row_fields= set(annot_list).intersection(row_fields)
        if len(annots_in_row_fields) == 0:
            print('No row fields matched annots in annot_list')
        else:
            print('Ignored fields for adding annot pattern: {}'.format(set(annot_list).difference(annots_in_row_fields)))
    for annot in annot_list:
        new_field = [annot_pattern,annot][prefix,prefix is True].join('')
        if new_field in row_fields:
            print('Row field name collision: {}'.format(new_field))
            print('To avoid name collision, rename field {}'.format(new_field))
        else:
            mt = mt._annotate_all(row_exprs={new_field:mt[annot]})
    return mt
        
@typecheck(tb=oneof(MatrixTable,
                    Table),
           annot_pattern=nullable(str),
           tau_ref_dict=nullable(dict))
def get_tau_dict(tb, annot_pattern=None, tau_ref_dict=None):
    '''Gets annotations matching annot_pattern and pairs with tau reference dict
    Number of annotations returned by annotation search should be less than or 
    equal to number of keys in tau_ref_dict.'''
    assert (annot_pattern != None or tau_ref_dict != None), "annot_pattern and tau_ref_dict cannot both be None"
    if annot_pattern is None: 
        tau_dict = {k: tau_ref_dict[k] for k in tau_ref_dict.keys() if k in tb.row} # take all row fields in mt matching keys in tau_dict
        assert len(tau_dict) > 0, 'None of the keys in tau_ref_dict match any row fields' #if intersect is empty: return error
        return tau_dict #return subset of tau_ref_dict
    else:
        pattern = re.compile(annot_pattern)
        annots = [rf for rf in list(tb.row) if pattern.match(rf)] #annotation search in list of row fields
        assert len(annots) > 0, 'No row fields matched annotation regex pattern: {}'.format(annot_pattern)
        if tau_ref_dict is None:
            print('Assuming tau = 1 for all annotations')
            return {k: 1 for k in annots}
        in_tau_ref_dict = set(annots).intersection(set(tau_ref_dict.keys())) #annots in tau_ref_dict
        if in_tau_ref_dict != set(annots): # if >0 annots returned by search are not in tau_ref_dict
            assert len(in_tau_ref_dict) > 0, 'annotation fields in tau_ref_dict do not match annotation search results' # if none of the annots returned by search are in tau_ref_dict
            print('Ignored fields from annotation search: {}'.format(set(annots).difference(in_tau_ref_dict)))
            print('To include ignored fields, change annot_pattern to regex match desired fields')
            annots = list(in_tau_ref_dict)
        return {k: tau_ref_dict[k] for k in annots}
    
@typecheck(mt=MatrixTable,
           tau_dict=nullable(dict),
           annot_pattern=nullable(str))
def agg_annotations(mt,tau_dict=None,annot_pattern=None):
    '''Aggregates annotations by linear combination. The coefficient are specified
    by tau_dict value, the annotation field name is specified by tau_dict key.'''
    assert (annot_pattern is None and tau_dict is None), "annot_pattern and tau_dict cannot both be None"
    mt = mt._annotate_all(row_exprs={'__annot':0},
                          global_exprs={'__tau_dict':tau_dict if tau_dict is not None else hl.null('dict'),
                                        '__annot_pattern':annot_pattern if annot_pattern is not None else hl.null('str')})
    tau_dict = get_tau_dict(tb=mt,annot_pattern=annot_pattern, tau_ref_dict=tau_dict)
    if str not in map(type,tau_dict.keys()): #if none of the keys are strings (maybe they are row exprs)
        pass
    print('Annotation fields and associated tau values used in annotation aggregation: {}'.format(tau_dict))
    for annot,tau in mt.tau_dict.items():
        mt = mt.annotate_rows(__annot = mt.__annot + tau*mt[annot])

@typecheck(mt=MatrixTable, 
           h2=oneof(float,int),
           pi=oneof(float,int),
           is_annot_inf=bool,
           tau_dict=nullable(dict),
           annot_pattern=nullable(str),
           h2_normalize=bool)
def make_betas(mt, h2, pi=1, is_annot_inf=False, tau_dict=None, annot_pattern=None, h2_normalize=False):
    '''Simulate betas. Options: Infinitesimal model, spike & slab, annotation-informed'''        
    M = mt.count_rows()
    if is_annot_inf:
        print('\rSimulating {} annotation-informed betas {}'.format(
                'h2-normalized' if h2_normalize else '',
                '(default tau: 1)' if tau_dict is None else 'using tau dict'))
        mt1 = agg_annotations(mt,tau_dict,annot_pattern=annot_pattern)
        annot_sum = mt1.aggregate_rows(hl.agg.sum(mt1.__annot))
        return mt1.annotate_rows(__beta = hl.rand_norm(0, hl.sqrt(mt1.__annot/annot_sum*(h2 if h2_normalize else 1)))) # if is_h2_normalized: scale variance of betas to be h2, else: keep unscaled variance
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
           is_annot_inf=bool,
           tau_dict=nullable(dict),
           annot_pattern=nullable(str),
           h2_normalize=bool,
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64),
           path_to_save=nullable(str))
def add_sim_description(mt,h2,starttime,stoptime,runtime,pi=1,is_annot_inf=False,
                        tau_dict=None, annot_pattern=None,h2_normalize=False,
                        popstrat=None,popstrat_s2=1,path_to_save=None):
    '''Annotates mt with description of simulation'''
    sim_id = 0
    while (str(sim_id) in [x.strip('sim_desc') for x in list(mt.globals) if 'sim_desc' in x]):
        sim_id += 1
    for arg in [tau_dict,annot_pattern,popstrat,path_to_save]: #convert None args to hl null expression
        if arg  is None: 
            arg = hl.null('str') # arbitrarily chose 'str' for all args
    sim_desc = hl.struct(h2=h2,starttime=str(starttime),stoptime=str(stoptime),
                         runtime=str(runtime),pi=pi,is_annot_inf=is_annot_inf,
                         annot_pattern=annot_pattern,tau_dict=tau_dict,
                         h2_normalize=h2_normalize, popstrat=popstrat,
                         popstrat_s2=popstrat_s2,path_to_save=path_to_save)
    mt = mt._annotate_all(global_exprs={'sim_desc{}'.format(sim_id):sim_desc})
    return mt

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, 
                          expr_float64),
           h2=oneof(float,int),
           pi=oneof(float,int),
           is_annot_inf=bool,
           tau_dict=nullable(dict),
           annot_pattern=nullable(str),
           h2_normalize=bool,
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)),
           popstrat_s2=oneof(float,
                             int,
                             expr_int32,
                             expr_int64,
                             expr_float32,
                             expr_float64),
           path_to_save=nullable(str))
def simulate(mt, genotype, h2, pi=1, is_annot_inf=False, tau_dict=None, annot_pattern=None,
             h2_normalize=False, popstrat=None, popstrat_s2 = 1,path_to_save=None):
    ''' Simulates phenotypes. 
        Options: Infinitesimal, spike/slab, annotation-informed, population stratification'''
    if is_annot_inf:
        assert (tau_dict != None or annot_pattern != None), 'If using annotation-informed model, both tau_dict and annot_pattern cannot be None'

    starttime = datetime.now()
    
    print_header(h2=h2, 
                 pi=pi, 
                 is_annot_inf=is_annot_inf, 
                 h2_normalize=h2_normalize,
                 popstrat=popstrat, 
                 popstrat_s2=popstrat_s2, 
                 path_to_save=path_to_save)
    
    
    mt1 = annotate_w_temp_fields(mt=mt, 
                                 genotype=genotype,
                                 h2=h2, 
                                 pi=pi, 
                                 is_annot_inf=is_annot_inf,
                                 popstrat=popstrat, 
                                 popstrat_s2=popstrat_s2)
    
    mt2 = make_betas(mt=mt1, 
                     h2=h2, 
                     pi=pi, 
                     is_annot_inf=is_annot_inf,
                     tau_dict=tau_dict,
                     annot_pattern=annot_pattern,
                     h2_normalize=h2_normalize)
        
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
    mt5 = add_sim_description(mt=mt4,h2=h2,starttime=starttime,stoptime=stoptime,
                              runtime=runtime,pi=pi,is_annot_inf=is_annot_inf,
                              tau_dict=tau_dict,annot_pattern=annot_pattern,
                              h2_normalize=h2_normalize,popstrat=popstrat,
                              popstrat_s2=popstrat_s2,path_to_save=path_to_save)
    print('\rFinished simulation! (runtime={} min)'.format(np.round((runtime.seconds+runtime.microseconds/1e6)/60, 4)).ljust(100))
    if path_to_save is not None:
        print('\rWriting simulation to: {}'.format(path_to_save))
        mt5.write(path_to_save)
    return mt5
