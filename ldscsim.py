#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Basic version of ldscsim.

Models for SNP effects:
    - Infinitesimal (can simulate n correlated traits)
    - Spike & slab (can simulate 2 correlated traits)
    - Annotation-informed
    
Field aggregation tools for annotation-informed model and population stratification
using many covariates.

Changes from previous version:
    - annotation-informed model takes only the field of aggregated annotations
    - adding population stratification takes only the field of aggregated covariates 
    - no header printed

@author: nbaya
"""

import hail as hl
from hail import dtype
from hail.typecheck import typecheck, oneof, nullable
from hail.expr.expressions import expr_float64, expr_int32, expr_array
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.utils.java import Env
import numpy as np
import pandas as pd
import random
import string
import scipy.stats as stats

@typecheck(mt=MatrixTable, 
           genotype=oneof(expr_int32,
                          expr_float64),
           h2=nullable(oneof(float,
                             int,
                             list)),
           pi=nullable(oneof(float,
                             int,
                             list)),
           rg=nullable(oneof(float,
                             int,
                             list)),
           annot=nullable(oneof(expr_float64,
                                expr_int32)),
           popstrat=nullable(oneof(expr_int32,
                                   expr_float64)))
def simulate_phenotypes(mt, genotype, h2, pi=None, rg=None, annot=None, popstrat=None):
    """Simulate phenotypes for testing LD score regression.
    
    Simulates betas (SNP effects) under the infinitesimal, spike & slab, or 
    annotation-informed models, depending on parameters passed. Optionally adds
    population stratification.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable containing genotypes to be used. Also should contain 
        variant annotations as row fields if running the annotation-informed
        model or covariates as column fields if adding population stratification.
    genotype : :class:`.Expression`
        Entry field containing genotypes of individuals to be used for the
        simulation.
    h2 : :obj:`float` or :obj:`int` or :obj:`list`
        Heritability of simulated trait. Can only be None if running annotation-
        informed model.
    pi : :obj:`float` or :obj:`int`
        Probability of SNP being causal when simulating under the spike & slab 
        model.
    rg : :obj:`float` or :obj:`int` or :obj:`list`
        Genetic correlation between traits.
    annot : :class:`.Expression`
        Row field to use as our aggregated annotations.
    popstrat: :class:`.Expression`
        Column field to use as our aggregated covariates for adding population
        stratification.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with simulated betas and phenotypes, simulated according
        to specified model.
    """
    tid = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5)) # "temporary id" -- random string to identify temporary intermediate fields generated by this method
    mt = annotate_all(mt=mt,
                      row_exprs={} if annot is None else {'annot_'+tid:annot},
                      col_exprs={} if popstrat is None else {'popstrat_'+tid:popstrat},
                      entry_exprs={'gt_'+tid:genotype},
                      global_exprs={'ldscsim':hl.struct(**{**{'h2':h2},
                                                           **({} if pi is None else {'pi':pi}),
                                                           **({} if rg is None else {'rg':rg}),
                                                           **({} if annot is None else {'is_annot_inf':True}),
                                                           **({} if popstrat is None else {'is_popstrat_inf':True})
                                                           })})
    mt = make_betas(mt=mt, 
                    h2=h2, 
                    pi=1 if pi is None else pi, 
                    annot=None if annot is None else mt['annot_'+tid],
                    rg=rg)
    mt = calculate_phenotypes(mt=mt, 
                              genotype=mt['gt_'+tid], 
                              beta=mt['beta'],
                              h2 = h2,
                              popstrat=None if popstrat is None else mt['popstrat_'+tid])
    mt = _clean_fields(mt, tid)
    return mt
    
@typecheck(mt=MatrixTable, 
           h2=nullable(oneof(float,
                             int,
                             list)),
           pi=oneof(float,
                    int,
                    list),
           annot=nullable(oneof(expr_float64,
                                expr_int32)),
           rg=nullable(oneof(float,
                             int,
                             list)))
def make_betas(mt, h2, pi=1, annot=None, rg=None):
    """Generates betas under different models. 
       
    Simulates betas (SNP effects) under the infinitesimal, spike & slab, or 
    annotation-informed models, depending on parameters passed.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable containing genotypes to be used. Also should contain 
        variant annotations as row fields if running the annotation-informed
        model or covariates as column fields if adding population stratification.
    h2 : :obj:`float` or :obj:`int` or :obj:`list`
        Heritability of simulated trait. Can only be None if running annotation-
        informed model.
    pi : :obj:`float` or :obj:`int` or :obj:`list`
        Probability of SNP being causal when simulating under the spike & slab 
        model. If doing multitrait spike & slab pi is a list of probabilities for
        overlapping causal SNPs
    annot : :class:`.Expression`
        Row field of aggregated annotations for annotation-informed model.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with betas as a row field, simulated according to specified model.
    """
    h2 = [h2] if type(h2) is not list else h2
    pi = [pi] if type(pi) is not list else pi
    rg = [rg] if type(rg) is not list else rg
    assert (all(x >= 0 and x <= 1 for x in h2)), 'h2 values must be between 0 and 1'
    assert (all(x >= 0 and x <= 1 for x in pi)), 'pi values for spike & slab must be between 0 and 1'
    assert (rg==[None] or all(x >= 0 and x <= 1 for x in rg)), 'rg values must be between 0 and 1 or None'
    if annot is not None: #annotation-informed
        M = mt.count_rows()
        annot_sum = mt.aggregate_rows(hl.agg.sum(annot))
        mt = mt.annotate_rows(**{'beta': hl.rand_norm(0, hl.sqrt(annot*h2[0]/annot_sum))}) # if is_h2_normalized: scale variance of betas to be h2, else: keep unscaled variance
        return mt
    elif len(h2)>1 and pi==[1]: #multi-trait infinitesimal
        return multitrait_inf(mt=mt,h2=h2,rg=rg)
    elif len(h2)==2 and len(pi)>1: #multi-trait spike & slab
        return multitrait_ss(mt=mt,h2=h2,rg=0 if rg is [None] else rg[0],pi=pi)
    elif len(h2)==1 and pi==[1]: #infinitesimal/spike & slab
        M = mt.count_rows()
        return mt.annotate_rows(beta = hl.rand_bool(pi[0])*hl.rand_norm(0,hl.sqrt(h2[0]/(M*pi[0]))))
    else:
        raise ValueError('Insufficient parameters')
        
@typecheck(mt=MatrixTable, 
           h2=nullable(oneof(float,
                             int,
                             list)),
           rg=nullable(oneof(float,
                             int,
                             list)),
           cov_matrix=nullable(np.ndarray),
           seed=nullable(int))
def multitrait_inf(mt, h2=None, rg=None, cov_matrix=None, seed=None):
    """Generates correlated betas for multi-trait infinitesimal simulations for 
    any number of phenotypes.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable for simulated phenotype.
    h2 : :obj:`float` or :obj:`int` or :obj:`list`, optional
        Desired heritability of simulated traits. If h2=None, h2 is based on 
        diagonal of cov_matrix.
    rg : :obj:`float` or :obj:`int` or :obj:`list`, optional
        Desired genetic correlation between simulated traits. If simulating more 
        than two correlated traits, rg should be a list of rg values corresponding
        to the upper right triangle of the covariance matrix. If rg=None and 
        cov_matrix=None, rg is assumed to be 0 between traits. If rg and 
        cov_matrix are both not None, rg values from cov_matrix are used instead.
    cov_matrix : :class:`numpy.ndarray`, optional
        Covariance matrix for traits, unscaled by `M`. Overrides h2 and rg if h2
        and rg are not None.
    seed : :obj:`int`, optional
        Seed for random number generator. If seed=None, seed is set randomly.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with simulated SNP effects as a row field of arrays.
    """
    h2 = [h2] if type(h2) is not list else h2
    rg = [rg] if type(rg) is not list else rg
    assert (all(x >= 0 and x <= 1 for x in h2)), 'h2 values must be between 0 and 1'
    assert h2 is not [None] or cov_matrix is not None, 'h2 and cov_matrix cannot both be None'
    seed = seed if seed is not None else int(str(Env.next_seed())[:8])
    M = mt.count_rows()
    if cov_matrix != None:
        n_phens = cov_matrix.shape[0]
    else:
        n_phens = len(h2)
        if rg == [None]:
            print(f'Assuming rg=0 for all {n_phens} traits')
            rg = [0]*int((n_phens**2-n_phens)/2)
        assert (all(x >= -1 and x<= 1 for x in rg)), 'rg values must be between 0 and 1'
        cov_matrix = create_cov_matrix(h2, rg)
    cov_matrix = (1/M)*cov_matrix
    randstate = np.random.RandomState(int(seed)) #seed random state for replicability
    betas = randstate.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
    df = pd.DataFrame([0]*M,columns=['beta'])
    tb = hl.Table.from_pandas(df)
    tb = tb.add_index().key_by('idx')
    tb = tb.annotate(beta = hl.literal(betas.tolist())[hl.int32(tb.idx)])
    mt = mt.add_row_index()
    mt = mt.annotate_rows(beta = tb[mt.row_idx]['beta'])
    return mt

@typecheck(mt=MatrixTable, 
           h2=list,
           pi=list,
           rg=oneof(float,
                    int))
def multitrait_ss(mt, h2, pi, rg=0):
    """Generates spike & slab betas for simulation of two correlated phenotypes.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable for simulated phenotype.
    h2 : :obj:`list`
        Desired heritability of simulated traits.
    pi : :obj:`list`
        List of proportion of SNPs: ptt, ptf, pft
        ptt is the proportion of SNPs that are causal for both traits, ptf is
        the proportion of SNPs that are causal for trait 1 but not trait 2,
        pft is the proportion of SNPs that are causal for trait 2 but not trait 1.
        pff is inferred based on the previous proportions and is the proportion
        of SNPs that are not causal for both traits.
    rg : :obj:`float` or :obj:`int`
        Genetic correlation between traits.
        
    Warning
    -------
    May give inaccurate results if chosen parameters make cov_matrix not positive
    semi-definite.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with simulated SNP effects as a row field of arrays.
    """
    ptt, ptf, pft, pff = pi[0], pi[1], pi[2], 1-sum(pi)
    cov_matrix = np.asarray([[1/(ptt+ptf), rg/ptt],[rg/ptt,1/(ptt+pft)]])
    M = mt.count_cols()
    beta = np.random.multivariate_normal(mean=np.zeros(2),cov=cov_matrix,size=[int(M),])
    zeros = np.zeros(shape=int(M)).T
    beta_matrix = np.stack((np.asarray([zeros,zeros]).T,np.asarray([zeros,beta[:,1]]).T,
                            np.asarray([beta[:,0],zeros]).T,beta),axis=1)
    idx = np.random.choice([0,1,2,3],p=[pff,pft, ptf, ptt],size=int(M))
    betas = beta_matrix[range(int(M)),idx,:]
    betas[:,0] *= (h2[0]/M)**(1/2)
    betas[:,1] *= (h2[1]/M)**(1/2)
    df = pd.DataFrame([0]*M,columns=['beta'])
    tb = hl.Table.from_pandas(df)
    tb = tb.add_index().key_by('idx')
    tb = tb.annotate(beta = hl.literal(betas.tolist())[hl.int32(tb.idx)])
    mt = mt.add_row_index()
    mt = mt.annotate_rows(beta = tb[mt.row_idx]['beta'])
    return mt

@typecheck(h2=list,
           rg=list)
def create_cov_matrix(h2, rg):
    """Creates covariance matrix for simulating correlated betas. 
    
    Parameters
    ----------
    h2 : :obj:`list`
        h2 values for traits. h2 values in list should be ordered by their order 
        in the diagonal of the covariance array, reading from top left to bottom
        right.
    rg : :obj:`list`        
        rg values for traits. rg values should be ordered in the order they appear
        in the upper triangle of the covariance matrix, from left to right, top to
        bottom.
    
    Notes
    -----
    Covariance matrix is not scaled by number of SNPs.
    
    Returns
    -------
    :class:`numpy.ndarray`
    """
    assert (all(x >= 0 and x <= 1 for x in h2)), 'h2 values must be between 0 and 1'
    assert (all(x >= -1 and x<= 1 for x in rg)), 'rg values must be between 0 and 1'
    n_rg = len(rg)
    n_h2 = len(h2)
    exp_n_rg = int((n_h2**2-n_h2)/2) #expected number of rg values, given number of traits
    assert n_rg is exp_n_rg, f'The number of rg values given is {n_rg}, expected {exp_n_rg}'
    rg = np.asarray(rg)
    cov_matrix = np.zeros(shape=[n_h2,n_h2])
    cov_matrix[np.triu_indices(n_h2, k=1)] = rg**2 #sets upper diagonal with covariances
    h2_diag = np.diag(h2)
    cov_matrix = (np.matmul(np.matmul(cov_matrix,h2_diag).T,h2_diag)**(1/2))
    cov_matrix += cov_matrix.T + h2_diag
    return cov_matrix

@typecheck(mt=MatrixTable,
           genotype=expr_int32,
           beta=oneof(expr_float64,
                      expr_array(expr_float64)),
           h2=nullable(oneof(float,
                             int,
                             list)),
           popstrat=nullable(expr_int32))
def calculate_phenotypes(mt, genotype, beta, h2, popstrat=None):
    """Calculates phenotypes by multiplying genotypes and betas.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable with all relevant fields passed as parameters.
    genotype : :class:`.Expression`
        Entry field of genotypes. 
    beta : :class:`.Expression`
        Row field of SNP effects.
    h2 : :obj:`float` or :obj:`int` or :obj:`list`
        Heritability of simulated trait. Can only be None if running annotation-
        informed model.
    popstrat : :class:`.Expression`, optional
        Column field containing population stratification term.
        
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with simulated phenotype as column field.
    """
    tid = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5)) # "temporary id" -- random string to identify temporary intermediate fields generated by this method
    mt = annotate_all(mt=mt,
                      row_exprs={'beta_'+tid:beta},
                      entry_exprs={'gt_'+tid:genotype})
    mt = normalize_genotypes(mt['gt_'+tid])
    if mt['beta_'+tid].dtype == dtype('array<float64>'):
        mt = mt.annotate_cols(y_no_noise = hl.agg.array_agg(lambda beta: hl.agg.sum(beta*mt['norm_gt']),mt['beta_'+tid]))
        mt = mt.annotate_cols(y = mt.y_no_noise + hl.literal(h2).map(lambda x: hl.rand_norm(0,hl.sqrt(1-x))))
    else:
        mt = mt.annotate_cols(y_no_noise = hl.agg.sum(mt['beta_'+tid] * mt['norm_gt']))
        mt = mt.annotate_cols(y = mt.y_no_noise + hl.rand_norm(0, hl.sqrt(1-h2)))
    mt = _clean_fields(mt, tid)
    return mt
    
@typecheck(genotypes=oneof(expr_int32,
                           expr_float64))
def normalize_genotypes(genotypes):
    """Normalizes genotypes to have mean 0 and variance 1
    
    Parameters
    ----------
    genotypes : :class:`.Expression`
        Entry field of genotypes.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with normalized genotypes.
    """
    tid = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5)) # "temporary id" -- random string to identify temporary intermediate fields generated by this method
    mt = genotypes._indices.source
    mt = mt.annotate_entries(**{'gt_'+tid: genotypes})
    mt = mt.annotate_rows(**{'gt_stats_'+tid: hl.agg.stats(mt['gt_'+tid])})
    mt = mt.annotate_entries(norm_gt = (mt['gt_'+tid]-mt['gt_stats_'+tid].mean)/mt['gt_stats_'+tid].stdev)  
    mt = _clean_fields(mt, tid)
    return mt

@typecheck(mt=MatrixTable,
           str_expr=str)
def _clean_fields(mt, str_expr):
    """Removes fields with names that have `str_expr` in them.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable with fields to be removed.
    str_expr : :obj:`str`
        string to filter names of fields to remove.
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable with specified fields removed.
    """
    all_fields = list(mt.col)+list(mt.row)+list(mt.entry)+list(mt.globals)
    return mt.drop(*(x for x in all_fields if str_expr in x))

@typecheck(mt=MatrixTable,
           row_exprs=dict,
           col_exprs=dict,
           entry_exprs=dict,
           global_exprs=dict)
def annotate_all(mt,row_exprs={},col_exprs={},entry_exprs={},global_exprs={}):
    """Equivalent of _annotate_all, but checks source MatrixTable of exprs"""
    exprs = {**row_exprs, **col_exprs, **entry_exprs, **global_exprs}
    for key, value in exprs.items():
        if type(value) == expr_float64 or type(value) == expr_int32:
            assert value._indices.source == mt, 'Cannot combine expressions from different source objects.'
    return mt._annotate_all(row_exprs, col_exprs, entry_exprs, global_exprs)

@typecheck(mt=MatrixTable,
           y=expr_int32,
           P=oneof(int,
                   float))
def ascertainment_bias(mt,y,P):
    """Adds ascertainment bias to a binary phenotype such that it was sample 
    prevalence of `P` = cases/(cases+controls).
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable containing binary phenotype to be used.
    y : :class:`.Expression`
        Column field of binary phenotype.
    P : :obj:`int` or :obj:`float`
        Desired "sample prevalence" of phenotype.
        
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable containing binary phenotype with prevalence of approx. P
    """
    assert P>=0 and P<=1, 'P must be in [0,1]'
    tid = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5)) # "temporary id" -- random string to identify temporary intermediate fields generated by this method
    mt = mt.annotate_cols(y_w_asc_bias = y)
    y_stats = mt.aggregate_cols(hl.agg.stats(mt.y_w_asc_bias))
    K = y_stats.mean
    n = y_stats.n
    assert abs(P-K)<1, 'Specified sample prevalence is incompatible with population prevalence.'
    if P < K:
        p = (1-K)*P/(K*(1-P))
        con = mt.filter_cols(mt.y_w_asc_bias==0) 
        cas = mt.filter_cols(mt.y_w_asc_bias==1).add_col_index(name='col_idx_'+tid)
        keep = round(p*n*K)*[1]+round((1-p)*n*K)*[0]
        cas = cas.annotate_cols(**{'keep_'+tid : hl.literal(keep)[hl.int32(cas['col_idx_'+tid])]})
        cas = cas.filter_cols(cas['keep_'+tid]==1)
        cas = _clean_fields(cas, tid)
        mt = cas.union_cols(con)
    elif P > K:
        p = K*(1-P)/((1-K)*P)
        cas = mt.filter_cols(mt.y_w_asc_bias==1)
        con = mt.filter_cols(mt.y_w_asc_bias==0).add_col_index(name='col_idx_'+tid)
        keep = round(p*n*(1-K))*[1]+round((1-p)*n*(1-K))*[0]
        con = con.annotate_cols(**{'keep_'+tid : hl.literal(keep)[hl.int32(con['col_idx_'+tid])]})
        con = con.filter_cols(con['keep_'+tid]==1)
        con = _clean_fields(con, tid)
        mt = con.union_cols(cas)
    return mt

@typecheck(mt=MatrixTable,
           y=oneof(expr_int32,
                   expr_float64),
           K=oneof(int,
                   float),
           exact=bool)
def binarize(mt,y,K,exact=False):
    """Binarize phenotype `y` such that it has prevalence `K` = cases/(cases+controls)
    Uses inverse CDF of Gaussian to set binarization threshold when exact = False, 
    otherwise uses ranking to determine threshold.
    
    Parameters
    ----------
    mt : :class:`.MatrixTable`
        MatrixTable containing phenotype to be binarized.
    y : :class:`.Expression`
        Column field of phenotype.
    K : :obj:`int` or :obj:`float`
        Desired "population prevalence" of phenotype.
    exact : :obj:`bool`
        Whether to get prevalence as close as possible to K (does not use inverse CDF)
    
    Returns
    -------
    :class:`.MatrixTable`
        MatrixTable containing binary phenotype with prevalence of approx. K
    """
    if exact: 
        key = list(mt.col_key)
        tid = ''.join(random.choices(string.ascii_uppercase + string.ascii_lowercase, k=5)) # "temporary id" -- random string to identify temporary intermediate fields generated by this method
        mt = mt.annotate_cols(**{'y_'+tid: y})
        tb = mt.cols().order_by('y_'+tid)
        tb = tb.add_index('idx_'+tid)
        n = tb.count()
        tb = tb.annotate(y_binarized = tb['idx_'+tid]+1 <= round(n*K)) # "+1" because of zero indexing
        tb, mt = tb.key_by('y_'+tid), mt.key_cols_by('y_'+tid)
        mt = mt.annotate_cols(y_binarized = tb[mt['y_'+tid]].y_binarized)
        mt = mt.key_cols_by(*map(lambda x: mt[x],key))
    else: #use inverse CDF
        y_stats = mt.aggregate_cols(hl.agg.stats(y))
        threshold = stats.norm.ppf(1-K,loc=y_stats.mean,scale=y_stats.stdev)
        mt = mt.annotate_cols(y_binarized = y > threshold)
    return mt

@typecheck(tb=oneof(MatrixTable,
                    Table),
           coef_dict=nullable(dict),
           str_expr=nullable(str),
           axis=str)
def agg_fields(tb,coef_dict=None,str_expr=None,axis='rows'):
    '''Aggregates by linear combination fields matching either keys in `coef_dict`
    or `str_expr`. 
    
    Parameters
    ----------
    tb : :class:`.MatrixTable` or :class:`.Table`
        MatrixTable or Table containing fields to be aggregated.
    coef_dict : :obj:`dict`, optional
        Coefficients to multiply each field. The coefficients are specified by 
        `coef_dict` value, the row (or col) field name is specified by `coef_dict` key.
        If not included, coefficients are assumed to be 1.
    str_expr : :obj:`str`, optional
        String expression to match against row (or col) field names.
    axis : :obj:`str`
        Either 'rows' or 'cols'. If 'rows', this aggregates across row fields. 
        If 'cols', this aggregates across col fields. If tb is a Table, axis = 'rows'.
    
    Returns
    -------
    :class:`.MatrixTable` or :class:`.Table`
        MatrixTable or Table containing aggregated fields.
    '''
    assert (str_expr != None or coef_dict != None), "str_expr and coef_dict cannot both be None"
    assert axis is 'rows' or axis is 'cols', "axis must be 'rows' or 'cols'"
    coef_dict = get_coef_dict(tb=tb,str_expr=str_expr, ref_coef_dict=coef_dict,axis=axis)
    axis_field = 'annot' if axis=='rows' else 'cov'
    annotate_fn = (MatrixTable.annotate_rows if axis=='rows' else MatrixTable.annotate_cols) if type(tb) is MatrixTable else Table.annotate
    tb = annotate_fn(self=tb,**{'agg_'+axis_field: 0})
    print(f'Fields and associated coefficients used in {axis_field} aggregation: {coef_dict}')
    for field,coef in coef_dict.items():
        tb = annotate_fn(self=tb,**{'agg_'+axis_field: tb['agg_'+axis_field]+coef*tb[field]})
    return tb

@typecheck(tb=oneof(MatrixTable,
                    Table),
           str_expr=nullable(str),
           ref_coef_dict=nullable(dict),
           axis=str)
def get_coef_dict(tb, str_expr=None, ref_coef_dict=None,axis='rows'):
    '''Gets either col or row fields matching `str_expr` and take intersection 
    with keys in coefficient reference dict.
    
    Parameters
    ----------
    tb : :class:`.MatrixTable` or :class:`.Table`
        MatrixTable or Table containing row (or col) for `coef_dict`.
    str_expr : :obj:`str`, optional
        String expression pattern to match against row (or col) fields. If left
        unspecified, the intersection of field names is only between existing 
        row (or col) fields in `mt` and keys of `ref_coef_dict`.
    ref_coef_dict : :obj:`dict`, optional
        Reference coefficient dictionary with keys that are row (or col) field 
        names from which to subset. If not included, coefficients are assumed to be 1.
    axis : :obj:`str`
        Field type in which to search for field names.
        
    Returns
    -------
    coef_dict : :obj:`dict`, optional
        Coefficients to multiply each field. The coefficients are specified by 
        `coef_dict` value, the row (or col) field name is specified by `coef_dict` key. 
    '''
    assert (str_expr != None or ref_coef_dict != None), "str_expr and ref_coef_dict cannot both be None"
    assert axis is 'rows' or axis is 'cols', "axis must be 'rows' or 'cols'"
    fields_to_search = (tb.row if axis=='rows' or type(tb) is Table else tb.col)
    axis_field = 'annotation' if axis=='rows' else 'covariate' #when axis='rows' we're searching for annotations, axis='cols' searching for covariates
    if str_expr is None: 
        coef_dict = {k: ref_coef_dict[k] for k in ref_coef_dict.keys() if k in fields_to_search} # take all row (or col) fields in mt matching keys in coef_dict
        assert len(coef_dict) > 0, f'None of the keys in ref_coef_dict match any {axis[:-1]} fields' #if intersect is empty: return error
        return coef_dict #return subset of ref_coef_dict
    else:
        fields = [rf for rf in list(fields_to_search) if str_expr in rf] #str_expr search in list of row (or col) fields
        assert len(fields) > 0, f'No {axis[:-1]} fields matched str_expr search: {str_expr}'
        if ref_coef_dict is None:
            print(f'Assuming coef = 1 for all {axis_field}s')
            return {k: 1 for k in fields}
        in_ref_coef_dict = set(fields).intersection(set(ref_coef_dict.keys())) #fields in ref_coef_dict
        if in_ref_coef_dict != set(fields): # if >0 fields returned by search are not in ref_coef_dict
            assert len(in_ref_coef_dict) > 0, f'None of the {axis_field} fields in ref_coef_dict match search results' # if none of the fields returned by search are in ref_coef_dict
            fields_to_ignore=set(fields).difference(in_ref_coef_dict)
            print(f'Ignored fields from {axis_field} search: {fields_to_ignore}')
            print('To include ignored fields, change str_expr to match desired fields')
            fields = list(in_ref_coef_dict)
        return {k: ref_coef_dict[k] for k in fields}
