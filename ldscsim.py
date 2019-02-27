#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 08:25:59 2019

ldsc simulation framework

@author: nbaya
"""

import hail as hl
from hail.expr.expressions import expr_int32, expr_int64, expr_float32, expr_float64
from hail.typecheck import typecheck, oneof, nullable
from hail.matrixtable import MatrixTable
from hail.table import Table
import re
from datetime import datetime, timedelta

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, 
                          expr_float64),
           h2=oneof(nullable(float),
                    nullable(int)),
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
def simulate(mt, genotype, h2=None, pi=1, is_annot_inf=False, tau_dict=None, annot_pattern=None,
             h2_normalize=True, popstrat=None, popstrat_s2 = 1,path_to_save=None):
    ''' Simulates phenotypes. 
        Options: 
            models for betas: Infinitesimal, spike/slab, annotation-informed
            models for phenotypes: population stratification'''
    check_beta_args(h2=h2,pi=pi,is_annot_inf=is_annot_inf,tau_dict=tau_dict,
                        annot_pattern=annot_pattern,h2_normalize=h2_normalize)
    
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
                                 tau_dict=tau_dict,
                                 annot_pattern=annot_pattern,
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
    print('\rFinished simulation! (runtime={} min)'.format(round((runtime.total_seconds())/60, 4)).ljust(100))
    if path_to_save is not None:
        print('\rWriting simulation to: {}'.format(path_to_save))
        mt5.write(path_to_save)
    return mt5

@typecheck(mt=MatrixTable,
           genotype=oneof(expr_int32,
                          expr_float64),
           beta=nullable(expr_float64),
           popstrat=oneof(nullable(expr_int32),
                          nullable(expr_float64)))
def check_mt_sources(mt,genotype,beta=None,popstrat=None):
    '''checks that mt matches source mt of genotype and popstrat'''
    if beta is not None:
        assert(mt == genotype._indices.source and mt == beta._indices.source), 'mt must match mt source of genotype and beta'
    if popstrat is not None:
        assert(mt == genotype._indices.source and mt == popstrat._indices.source), 'mt must match mt source of genotype and popstrat'
    else:
        assert(mt == genotype._indices.source), 'mt must match mt source of genotype'
        

@typecheck(h2=oneof(nullable(float),
                    nullable(int)),
           pi=oneof(float,int),
           is_annot_inf=bool,
           annot_coef_dict=nullable(dict),
           annot_regex=nullable(str),
           h2_normalize=bool)
def check_beta_args(h2=None, pi=1, is_annot_inf=False, annot_coef_dict=None, 
                    annot_regex=None, h2_normalize=True):
    '''checks beta args for simulate() and make_betas()'''
    if is_annot_inf: #if using the annotation-informed model
        assert (annot_coef_dict != None or annot_regex != None), 'If using annotation-informed model, both tau_dict and annot_pattern cannot be None'
        if h2_normalize:
            assert (h2 != None), 'h2 cannot be None when h2_normalize=True'
            assert (h2 >= 0 and h2 <= 1), 'h2 must be in [0,1]'
        if not h2_normalize and h2 != None and not (h2 >= 0 and h2 <= 0):
            print('Ignoring non-valid h2={} (not in [0,1]) because h2_normalize=False'.format(h2))
    else:
        assert (h2 != None), 'h2 cannot be None, unless running annotation-informed model'
        assert (h2 >= 0 and h2 <= 1), 'h2 must be in [0,1]'
        assert (pi >= 0 and pi <= 1), 'pi must be in [0,1]'
        assert h2_normalize == True, 'h2_normalize cannot be true unless running annotation-informed model'
        if annot_coef_dict != None or annot_regex != None:
            print('Ignoring annotation-informed-related args because is_annot_inf=False')

@typecheck(h2=oneof(nullable(float),
                    nullable(int)),
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
    header += 'h2 = {}\n'.format(h2) if (h2_normalize and h2 != None) else ''
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
           h2=oneof(nullable(float),
                    nullable(int)),
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
                             expr_float64))
def annotate_w_temp_fields(mt, genotype, h2=None, pi=1, is_annot_inf=False, tau_dict=None,
                           annot_pattern=None, h2_normalize=True, popstrat=None, popstrat_s2=1):
    '''Annotate mt with temporary fields of simulation parameters'''
    check_mt_sources(mt=mt,genotype=genotype)
    mt1 = mt._annotate_all(col_exprs={'__popstrat_temp':none_to_null(popstrat)},
                           entry_exprs={'__gt_temp':genotype},
                           global_exprs={'__h2_temp':none_to_null(h2), 
                                         '__pi_temp':pi,
                                         '__is_annot_inf_temp':is_annot_inf,
                                         '__tau_dict_temp':none_to_null(tau_dict),
                                         '__annot_pattern_temp':none_to_null(annot_pattern),
                                         '__h2_normalize_temp':h2_normalize})
    return mt1

@typecheck(arg=oneof(nullable(float),
                     nullable(int),
                     nullable(dict),
                     nullable(str),
                     nullable(expr_int32),
                     nullable(expr_float64)))
def none_to_null(arg):
    '''Converts arg to hl null representation if arg is None'''
    if arg is None:
        return hl.null('str')
    else:
        return arg

@typecheck(mt=MatrixTable, 
           h2=oneof(nullable(float),
                    nullable(int)),
           pi=oneof(float,int),
           is_annot_inf=bool,
           tau_dict=nullable(dict),
           annot_pattern=nullable(str),
           h2_normalize=bool)
def make_betas(mt, h2=None, pi=1, is_annot_inf=False, annot_coef_dict=None, annot_regex=None, h2_normalize=True):
    '''Simulate betas. Options: Infinitesimal model, spike & slab, annotation-informed'''  
    check_beta_args(h2=h2,pi=pi,is_annot_inf=is_annot_inf,tau_dict=annot_coef_dict,
                    annot_pattern=annot_regex,h2_normalize=h2_normalize)
    M = mt.count_rows()
    if is_annot_inf:
        print('\rSimulating {} annotation-informed betas {}'.format(
                'h2-normalized' if h2_normalize else '',
                '(default tau: 1)' if annot_coef_dict is None else 'using tau dict'))
        mt1 = agg_fields(mt=mt,coef_dict=annot_coef_dict,regex=annot_regex)
        annot_sum = mt1.aggregate_rows(hl.agg.sum(mt1.__annot))
        return mt1.annotate_rows(__beta = hl.rand_norm(0, hl.sqrt(mt1.__annot*(h2/annot_sum if h2_normalize else 1)))) # if is_h2_normalized: scale variance of betas to be h2, else: keep unscaled variance
    else:
        print('Simulating betas using {} model w/ h2 = {}'.format(('infinitesimal' if pi is 1 else 'spike & slab'),h2))
        mt1 = mt.annotate_globals(__h2 = none_to_null(h2), __pi = pi)
        return mt1.annotate_rows(__beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(h2/(M*pi))))

@typecheck(mt=MatrixTable,
           coef_dict=nullable(dict),
           regex=nullable(str))
def agg_fields(mt,coef_dict=None,regex=None,axis='rows'):
    '''Aggregates fields by linear combination. The coefficient are specified
    by coef_dict value, the row (or col) field name is specified by coef_dict key.
    By default, it searches row field annotations.'''
    assert (regex != None or coef_dict != None), "regex and coef_dict cannot both be None"
    assert axis is 'rows' or axis is 'cols', "axis must be 'rows' or 'cols'"
    coef_dict = get_coef_dict(tb=mt,annot_pattern=regex, tau_ref_dict=coef_dict,axis=axis)
    axis_field = 'annot' if axis=='rows' else 'cov'
    if axis == 'rows':
        mt = mt._annotate_all(row_exprs={'__agg_annot':0},
                              global_exprs={f'__annot_coef_dict':none_to_null(coef_dict),
                                            f'__annot_regex':none_to_null(regex)})
    else:
        mt = mt._annotate_all(col_exprs={'__agg_cov':0},
                              global_exprs={f'__cov_coef_dict':none_to_null(coef_dict),
                                            f'__cov_regex':none_to_null(regex)})
        mt = mt._annotate_all(col_exprs={'__agg_cov':0})
    print(f'Fields and associated coefficients used in {axis_field} aggregation: {coef_dict}')
    for annot,tau in coef_dict.items():
        mt = mt.annotate_rows(__annot = mt.__annot + tau*mt[annot])
    return mt

@typecheck(mt=oneof(MatrixTable),
           regex=nullable(str),
           coef_ref_dict=nullable(dict),
           axis=str)
def get_coef_dict(mt, regex=None, coef_ref_dict=None,axis='rows'):
    '''Gets annotations matching regex and pairs with coefficient reference dict
    Number of annotations returned by annotation search should be less than or 
    equal to number of keys in coef_ref_dict. By default, this searches row field
    annotations.'''
    assert (regex != None or coef_ref_dict != None), "regex and coef_ref_dict cannot both be None"
    assert axis is 'rows' or axis is 'cols', "axis must be 'rows' or 'cols'"
    fields_to_search = (mt.row if axis=='rows' else mt.col)
    axis_field = 'annotation' if axis=='rows' else 'covariate' #when axis='rows' we're searching for annotations, axis='cols' searching for covariates
    if regex is None: 
        coef_dict = {k: coef_ref_dict[k] for k in coef_ref_dict.keys() if k in fields_to_search} # take all row (or col) fields in mt matching keys in coef_dict
        assert len(coef_dict) > 0, f'None of the keys in coef_ref_dict match any {axis[:-1]} fields' #if intersect is empty: return error
        return coef_dict #return subset of coef_ref_dict
    else:
        pattern = re.compile(regex)
        fields = [rf for rf in list(fields_to_search) if pattern.match(rf)] #regex search in list of row (or col) fields
        assert len(fields) > 0, f'No {axis[:-1]} fields matched regex search: {regex}'
        if coef_ref_dict is None:
            print(f'Assuming coef = 1 for all {axis_field}s')
            return {k: 1 for k in fields}
        in_coef_ref_dict = set(fields).intersection(set(coef_ref_dict.keys())) #fields in coef_ref_dict
        if in_coef_ref_dict != set(fields): # if >0 fields returned by search are not in coef_ref_dict
            assert len(in_coef_ref_dict) > 0, f'None of the {axis_field} fields in coef_ref_dict match search results' # if none of the fields returned by search are in coef_ref_dict
            fields_to_ignore=set(fields).difference(in_coef_ref_dict)
            print(f'Ignored fields from {axis_field} search: {fields_to_ignore}')
            print('To include ignored fields, change regex to match desired fields')
            fields = list(in_coef_ref_dict)
        return {k: coef_ref_dict[k] for k in fields}

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
    check_mt_sources(mt,genotype,beta)
    mt1 = mt._annotate_all(row_exprs={'__beta':beta},
                           col_exprs={'__popstrat':none_to_null(popstrat)},
                           entry_exprs={'__gt':genotype},
                           global_exprs={'__popstrat_s2':popstrat_s2})
    mt2 = normalize_genotypes(mt1.__gt)
    print('\rCalculating phenotypes{}...'.format('' if popstrat is None else ' w/ population stratification').ljust(81))
    mt3 = mt2.annotate_cols(__y_no_noise = hl.agg.sum(mt2.__beta * mt2.__norm_gt))
    mt4 = mt3.annotate_cols(__y = mt3.__y_no_noise + hl.rand_norm(0,hl.sqrt(1-h2)))

    if popstrat is None:
        return mt4        
    else:
        return add_pop_strat(mt4, 
                             y=mt4.__y, 
                             popstrat=mt4.__popstrat, 
                             popstrat_s2=hl.eval(popstrat_s2))

@typecheck(genotypes=oneof(expr_int32,
                          expr_int64, 
                          expr_float32, 
                          expr_float64))
def normalize_genotypes(genotypes):
    '''Normalizes genotypes'''
    print('\rNormalizing genotypes...')
    mt = genotypes._indices.source #get source matrix table of genotypes
    mt1 = mt.annotate_entries(__gt = genotypes)
    mt2 = mt1.annotate_rows(__gt_stats = hl.agg.stats(mt1.__gt))
    return mt2.annotate_entries(__norm_gt = (mt2.__gt-mt2.__gt_stats.mean)/mt2.__gt_stats.stdev)  

@typecheck(mt=MatrixTable, 
           y=oneof(expr_int32,
                   expr_float64),
           popstrat=oneof(expr_int32,
                          expr_float64),
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

@typecheck(mt = MatrixTable, 
           str_expr=str)
def clean_fields(mt, str_expr):
    '''Removes fields with names that have str_expr in them'''
    all_fields = list(mt.col)+list(mt.row)+list(mt.entry)+list(mt.globals)
    for field_to_drop in [x for x in all_fields if str_expr in x]:
        mt = mt.drop(field_to_drop)
    return mt

@typecheck(mt=MatrixTable, 
           h2=oneof(nullable(float),
                    nullable(int)),
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
def add_sim_description(mt,starttime,stoptime,runtime,h2=None,pi=1,is_annot_inf=False,
                        tau_dict=None, annot_pattern=None,h2_normalize=True,
                        popstrat=None,popstrat_s2=1,path_to_save=None):
    '''Annotates mt with description of simulation'''
    sim_id = 0
    while (str(sim_id) in [x.strip('sim_desc') for x in list(mt.globals) if 'sim_desc' in x]):
        sim_id += 1
    sim_desc = hl.struct(h2=none_to_null(h2),pi=pi,starttime=str(starttime),
                         stoptime=str(stoptime),runtime=str(runtime),
                         is_annot_inf=is_annot_inf,tau_dict=none_to_null(tau_dict),
                         annot_pattern=none_to_null(annot_pattern),h2_normalize=h2_normalize, 
                         popstrat=none_to_null(popstrat),popstrat_s2=popstrat_s2,path_to_save=none_to_null(path_to_save))
    mt = mt._annotate_all(global_exprs={'sim_desc{}'.format(sim_id):sim_desc})
    return mt
