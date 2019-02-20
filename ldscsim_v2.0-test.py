#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 08:57:43 2019

Version 2.0 of ldsc simulation framework.

@author: nbaya
"""

import hail as hl

mt = hl.read_matrix_table('/Users/nbaya/Documents/lab/ldscsim/ukbb31063.GT.chr22.MAF_0.05.white_british.mt/')

mt_annot = mt.annotate_entries(genotype = mt.GT.n_alt_alleles())

mt_annot = mt_annot.annotate_rows(a1 = hl.rand_norm(0, 1))
mt_annot = mt_annot.annotate_rows(a2 = hl.rand_bool(0.1))
mt_annot = mt_annot.annotate_rows(annot = mt_annot.a1+mt_annot.a2)


mt_annot.describe()

def normalize_genotypes(mt):
    mt1 = mt.annotate_rows(__gt_stats = hl.agg.stats(mt.__gt))
    return mt1.annotate_entries(__norm_gt = (mt1.__gt-mt1.__gt_stats.mean)/mt1.__gt_stats.stdev)

def make_betas(mt, h2, pi=1, annot=None):
    M = mt.count_rows()
    if annot is not None:
        annot_stats = mt.aggregate_rows(hl.agg.stats(mt.__annot), _localize=True)
        return mt.annotate_rows(__beta = (mt.__annot - annot_stats.mean) / annot_stats.stdev * hl.sqrt(h2 / M))
    else:
        return mt.annotate_rows(__beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(h2/(M*pi))))

def add_pop_strat(mt, popstrat, popstrat_c):
    popstrat_stats = mt.aggregate_cols(hl.agg.stats(mt.__popstrat), _localize=True)
    mt1 = mt.annotate_cols(__popstrat = (mt.__popstrat-popstrat_stats.mean)/popstrat_stats.stdev)
    return mt1.annotate_cols(__y_w_popstrat = mt1.__y + mt1.__popstrat*mt1.__popstrat_c)
    
def sim_phenotypes(mt, h2, popstrat=None, popstrat_c=None):
    mt1 = mt.annotate_cols(__y_no_noise = hl.agg.sum(mt.__beta * mt.__norm_gt))
    mt2 = mt1.annotate_cols(__y = mt1.__y_no_noise + hl.rand_norm(0,hl.sqrt(1-h2)))
    if popstrat is not None:
        return add_pop_strat(mt2, popstrat, popstrat_c)
    else:
        return mt2
  
def simulate(mt, h2, genotype, pi=1, annot=None, popstrat=None, popstrat_c = None):
    if annot is None:
        mt1 = mt._annotate_all(entry_exprs={'__gt':genotype})
    elif popstrat is None:
        mt1 = mt._annotate_all(row_exprs={'__annot':annot},
                               entry_exprs={'__gt':genotype})
    else:
        mt1 = mt._annotate_all(row_exprs={'__annot':annot},
                               col_exprs={'__popstrat':popstrat},
                               entry_exprs={'__gt':genotype},
                               global_exprs={'__popstrat_c':popstrat_c})
    mt2 = normalize_genotypes(mt1)
    mt3 = make_betas(mt2, h2, pi, annot)
    return sim_phenotypes(mt3, h2, popstrat)

# Infinitesimal
sim_mt = simulate(mt=mt, h2=0.5, genotype = mt.GT.n_alt_alleles())
# Spike & slab
sim_mt = simulate(mt=mt, h2=0.4, pi=0.1, genotype = mt.genotype)
# Annotation-informed betas
sim_mt = simulate(mt=mt_annot, h2=0.3, genotype = mt_annot.genotype, annot = mt_annot.annot)




sim_mt.aggregate_entries(hl.agg.stats(sim_mt.norm_dosage))

sim_mt.aggregate_cols(hl.agg.stats(sim_mt.y_no_noise))

y = sim_mt.select_cols(sim_mt.y_no_noise).make_table()

cov = hl.import_table('/Users/nbaya/Documents/lab/ldscsim/ukb31063.gwas_covariates.both_sexes.tsv',impute=True, types={'s': hl.tstr}).key_by('s')

mt0 = sim_mt.annotate_cols(**cov[sim_mt.s])

mt0 = mt0.rename({'__norm_gt__': 'x'})

mt = mt0

cov_list = [ mt['isFemale'], mt['age'], mt['age_squared'], mt['age_isFemale'],
                            mt['age_squared_isFemale'] ]+ [mt['PC{:}'.format(i)] for i in range(1, 21)] 

ht = hl.linear_regression_rows(
                y=mt.y,
                x=mt.x,
                covariates=[1]+cov_list,
                pass_through = ['rsid'])

ht = ht.rename({'rsid':'SNP'}).key_by('SNP')

ht = ht.select(Z = ht.beta/ht.standard_error)

sumstats_template = hl.import_table('gs://nbaya/rg_sex/50_snps_alleles_N.tsv.gz',types={'N': hl.tint64})
sumstats_template = sumstats_template.key_by('SNP')
sumstats_template = sumstats_template.annotate(N = n_samples)
#            sumstats_template.show()

sumstats = sumstats_template.annotate(Z = ht[sumstats_template.SNP]['Z'])