#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 08:59:45 2020

Directly simulate GWAS summary statistics without individual-level genotypes for 
discovery GWAS.

Uses numpy arrays instead of block matrices

@author: nbaya
"""

import hail as hl
import numpy as np
from scipy.sparse import coo_matrix, block_diag, identity
import matplotlib.pyplot as plt

wd = '/Users/nbaya/Documents/lab/risk_gradients/data'

autosomes = range(1,23)

def import_from_plink(ref_panel):

    mt = hl.import_plink(bed=f'{wd}/{ref_panel}.bed',
                         bim=f'{wd}/{ref_panel}.bim',
                         fam=f'{wd}/{ref_panel}.fam')
    X = hl.linalg.BlockMatrix.from_entry_expr(mt.GT.n_alt_alleles())
    X = X.T
    
    X.write(f'{wd}/{ref_panel}.X.bm', overwrite=True)
        
def blockmatrix_to_numpy(ref_panel):
    non_hail_wd = wd.replace('file://','')
    X = hl.linalg.BlockMatrix.read(f'{wd}/{ref_panel}.X.bm')
    np.save(f'{non_hail_wd}/{ref_panel}.X.npy', X.to_numpy(), allow_pickle=False)

def get_X(ref_panel, chrom_list=autosomes, as_list=True):
    r'''
    Returns N_ref x M dim numpy matrix of column-standardized genotypes of LD ref panel
    If `as_list`=True and `chrom_list` is a str or int, `X` is returned as a list 
    with an ndarray as the only element
    '''    
    assert any(map(lambda x: isinstance(chrom_list, x), [int, str, list, range])), '`chrom_list` must be an int, str, list or range'
    if (isinstance(chrom_list, list)|isinstance(chrom_list, range)):
        X = []
        for chrom in chrom_list:
            X_chrom = get_X(ref_panel=ref_panel, chrom_list=chrom, as_list=False)
            X.append(X_chrom)
    else:
        X = np.load(f'{wd}/{ref_panel}.chr{chrom_list}.X.npy')
        invariant_idxs = X.std(axis=0)==0 # indices of invariant SNPs
        if any(invariant_idxs):
            print(f'chr{chrom_list} MAF=0 SNPs: {sum(invariant_idxs)}')
        X = X[:,~invariant_idxs]
        X -= X.mean(axis=0)
        X /= X.std(axis=0)
        if as_list: X = [X]
    
    return X
    
def get_block_idxs(x, transpose=False, as_list=False):
    r'''
    Returns indices for blocks contained in list of arrays `x`
    '''
    if transpose: x = [block.T for block in x]
    block_lengths = [0]+list(map(len, x)) # use 1 instead of 0 because we need to subtract by 1 in next step
    idx =  np.cumsum(block_lengths) # subtract 1 to make this zero indexed
    if as_list: idx=idx.tolist()
    return idx

def get_yg(X : list, beta : np.ndarray):
    r'''
    Returns genetic component of trait
    '''
    X_chrom_idx_list = get_block_idxs(X, transpose=True)
    yg_list = [X_chrom@beta[start:stop] for X_chrom, start, stop in 
               zip(X, X_chrom_idx_list[:-1],X_chrom_idx_list[1:])]
    yg = np.sum(np.asarray(yg_list), axis=0)
    yg = np.squeeze(yg)
    return yg

def get_N(X : list):
    r'''
    Gets length of first dimension of arrays in list of arrays `X`. 
    This checks that the length of the first dimension is the same for all arrays
    in the list.
    '''
    N = X[0].shape[0]
    assert all([x.shape[0]==N for x in X]), 'ERROR: Sample count varies across chromosomes'
    return N

def normalize_beta(X, beta, h2):
    r'''
    Normalize beta to have the right heritability
    '''
    yg = get_yg(X=X, beta=beta)
    N_r = get_N(X)
    s2 = (1/N_r)*(yg.T)@(yg)
    beta *= np.sqrt(h2/s2)
    beta = np.squeeze(beta)
    return beta
    

def get_beta(M, h2, X, pi=None, seed=None):
    r'''
    Returns M-dim vector of true SNP effect sizes
    '''
    assert isinstance(X, list) & all(map(lambda x: isinstance(x, np.ndarray), X)), 'X must be a list of numpy ndarrays'
    
    np.random.seed(seed=seed)
    if pi==None or pi==1: #infinitesimal model
        beta = np.random.normal(scale=np.sqrt(h2/M), size=(M,1))
    if pi!=None: #spike & slab model
        m_causal = round(M*pi)
        causal_beta = np.random.normal(scale=np.sqrt(h2/m_causal), size=(m_causal,1))
        zeros = np.zeros(shape=(M-m_causal, 1))
        beta = np.concatenate((causal_beta,zeros))
        np.random.shuffle(beta)
    beta = normalize_beta(X=X, beta=beta, h2=h2)

    return beta

def create_ld_blocks(chrom_list, peak_radius = 500000, max_ldblk_len=5000000, peak_rec_rate=2e-8):
    r'''
    Creates LD blocks, given provided recombination rates.
    `peak_radius`: minimum radius in base pairs around a local peak in recombination rate
    `peak_rec_rate`: minimum recombination rate for a position to be a "peak", a prospective hotspot
    `max_ldblk_len`: maximum length of an LD block in base pairs
    '''
    # load recombination maps (from https://github.com/nikbaya/risk_gradients/tree/master/data)
    rec_map_df = pd.read_csv('/Users/nbaya/Documents/lab/smiles/data/genetic_map_combined_b37.txt.gz',
                             delim_whitespace=True, compression='gzip', names=['chr','position','rate','cm'])
    assert 'chr' in rec_map_df.columns, 'Chromosome field in `rec_map_df`must be named "chr"'
    
    rec_map_df_dict = {chrom:rec_map_df[rec_map_df.chr==chrom] for chrom in chrom_list}
    first_position_dict = {chrom:[] for chrom in chrom_list} # dict of lists of positions for base pairs at the start of each LD block for each chromosome
    max_position_list = [] # list of maximum base pair positions
    for chrom in chrom_list:
        max_position = rec_map_df.position.max()
        max_position_list.append(max_position)
        rec_map_df_chrom = rec_map_df_dict[chrom]
        first_position = 0 # first position of window defining current LD block (left-most side of window)
        first_position_dict[chrom].append(first_position)
        # initialize with `None` to indicate that we don't have a peak
        peak_idx = None # index of hotspot in recmap dataframe
        peak_position = None # position in base pairs of current peak
        rec_map_positions = rec_map_df_chrom.position.values # list of positions in recmap
        rec_map_rates = rec_map_df_chrom.rate.tolist() # list of rates in recmap
        for idx, position, rate in zip(rec_map_df_chrom.index, rec_map_positions, rec_map_rates):
            if peak_position == None and rate > peak_rec_rate: # if no peak position has been found yet and current position has recombination rate > threshold
                peak_idx = idx    
                peak_position = position
                peak_rate = rate
                continue
            elif peak_position != None and (position > peak_position+peak_radius or position > first_position+max_ldblk_len): # if current position is outside of peak radius and max ld block length
                first_position = position # first position of window defining current LD block (left-most side of window)
                first_position_dict[chrom].append(first_position)
                # reset for new LD block
                peak_idx = None # index of hotspot in recmap dataframe
                peak_position = None # position in base pairs of current peak
                peak_rate = peak_rec_rate # start at baseline of rec rate threshold
            elif rate > peak_rate: # update if still in ld block and at a new maximum rate
                peak_idx = idx    
                peak_position = position
                peak_rate = rate
#    for chrom in chrom_list:
#        print(f'LD blocks in chrom {chrom}: {len(first_position_dict[chrom])}',
#              f'(length: mean={round(np.diff(first_position_dict[chrom]).mean())},',
#              f'std={round(np.diff(first_position_dict[chrom]).std())})')
    
    return first_position_dict

def convert_breakpoints(chrom, X, first_position_dict):
    r'''
    Converts breakpoints from recombination map into breakpoints in bim file.
    '''
    bim = pd.read_csv(f'{wd}/1kg_eur.chr{chrom}.bim',
                      delim_whitespace=True, names=['chr','snp','cm','position','a1','a2'])
    bim = bim.sort_values(by='position').reset_index() # just in case
    
    first_position_chrom = first_position_dict[chrom]
    break_points_chrom = [bim[bim.position>position].index[0] for position in 
                          first_position_chrom]
    
    break_points_chrom = sorted(list(set(break_points_chrom)))
    
    return break_points_chrom
        
def plot_breakpoints_from_recmap():
    for chrom in chrom_list:
        rec_map_df_chrom = rec_map_df_dict[chrom]
        plt.plot(rec_map_df_chrom.position, rec_map_df_chrom.rate)
        plt.plot(first_position_dict[chrom], [rec_map_df_chrom.rate.max()]*len(first_position_dict[chrom]), 'r.')
    

def munge_break_pts_chrom(break_pts_chrom, M):
    r'''
    Prepares list of breakpoints in a given chromosome to ensure that it is 
    ordered by base pair position and includes a zero as the first element and
    `M`, the number of SNPs in the given chromosome, as the final element.
    '''
    break_pts_chrom = sorted(list(set(break_pts_chrom)))
    if 0 not in break_pts_chrom:
        break_pts_chrom.insert(0,0)
    if M not in break_pts_chrom:
        break_pts_chrom.append(M)
    return break_pts_chrom

def get_sparse_R(X, break_pts, as_list=True, decimals=None):
    N_r = get_N(X)
    break_pts = [munge_break_pts_chrom(break_pts_chrom=break_pts_chrom, M=X_chrom.shape[1]) for break_pts_chrom, X_chrom in zip(break_pts, X)]
    
    # TODO: Consider rounding R to the nearest decimal point to remove possible noise
    if decimals==None:
        R = [[X_chrom[:,i:j].T@X_chrom[:,i:j]/N_r for i,j in zip(break_pts_chrom[:-1], break_pts_chrom[1:])] 
             for X_chrom, break_pts_chrom in zip(X, break_pts)] # runtime for 3m variants, 300 samples: 24 sec uncached, 5 sec cached
    else:
        R = [[np.round(X_chrom[:,i:j].T@X_chrom[:,i:j]/N_r, decimals=decimals)
                       for i,j in zip(break_pts_chrom[:-1], break_pts_chrom[1:])] 
         for X_chrom, break_pts_chrom in zip(X, break_pts)] # runtime for 3m variants, 300 samples: 24 sec uncached, 5 sec cached
    # TODO: Consider explicitly setting diagonal entries to 1, this is not guaranteed. Does it affect the results?
    if not as_list:
        R = block_diag([R_block for R_chrom in R for R_block in R_chrom])
        
    return R
    
#    diag_blocks = [coo_matrix(X[:,i:j].T@X[:,i:j]) for i,j in zip(break_pts[:-1], break_pts[1:])]
#    R = block_diag(diag_blocks)
#    R /= N_r
#    return R

def concatenate(x_list: list):
    r'''
    Concatenates elements of list `x`, a list of numpy ndarrays and/or scalars
    '''
    x_list = [x if isinstance(x,np.ndarray) else [x] for x in x_list] # accounts for zero-dimensional elements before np.concatenate
    x_concat = np.concatenate(x_list, axis=0)
    return x_concat
    
def get_alpha(R, beta):
    r'''
    Returns M-dim vector of true marginal SNP effect sizes
    Includes code for R as a list of lists of ndarrays or R as a sparse coo_matrix
    '''
    assert isinstance(R, coo_matrix)|isinstance(R, list), 'R must be coo_matrix or list of list of numpy ndarrays'
    if isinstance(R, list):
        assert all([isinstance(R_block, np.ndarray) for R_chrom in R for R_block in R_chrom]), 'R must be coo_matrix or list of list of numpy ndarrays'
    if isinstance(R, coo_matrix):
        alpha = R@beta
    else:
        R_chrom_idx_list = np.cumsum([0]+[sum(map(len, R_chrom)) for R_chrom in R])
        R_idx_list = [get_block_idxs(R_chrom, as_list=True) for R_chrom in R]
        alpha_list = [R_block@beta[R_chrom_idx+start:R_chrom_idx+stop] 
                      for R_chrom_idx, R_idx, R_chrom in zip(R_chrom_idx_list[:-1], R_idx_list, R) 
                      for start, stop, R_block in zip(R_idx[:-1], R_idx[1:], R_chrom) ]
        alpha = concatenate(x_list=alpha_list)
    return alpha

def get_Z(N_r, seed=None):
    r'''
    Returns `N_r`-dim standard normal random vector
    '''
    np.random.seed(seed=seed)
    Z = np.random.normal(size=(N_r,1)) # N_r-dimensional standard normal random vector
    return Z

def get_alphahat(alpha, N_d, R, seed=None):
    r'''
    Returns M-dim vector of estimated marginal SNP effect sizes
    '''
    
#    if isinstance(Z, type(None)):
#        Z = get_Z(N_r=N_r, seed=seed)
#    alphahat = alpha + 1/np.sqrt(N_d*N_r)*X.T@Z
    
    np.random.seed(seed=seed)
    R_chrom_idx_list = np.cumsum([0]+[sum(map(len, R_chrom)) for R_chrom in R])
    
    R_idx_list = [get_block_idxs(R_chrom) for R_chrom in R]
    alphahat_list = [np.random.multivariate_normal(mean=np.squeeze(alpha)[R_chrom_idx+start:R_chrom_idx+stop],
                                                   cov=1/np.sqrt(N_d)*R_block) 
                    if len(R_block)>1 else np.random.normal(loc=np.squeeze(alpha)[R_chrom_idx+start:R_chrom_idx+stop],
                                                            scale=1/np.sqrt(N_d)*np.squeeze(R_block)) 
                    for R_chrom_idx, R_idx, R_chrom in zip(R_chrom_idx_list[:-1], R_idx_list, R) 
                    for start, stop, R_block in zip(R_idx[:-1], R_idx[1:], R_chrom) ]
    alphahat = concatenate(x_list=alphahat_list)
    return alphahat

def get_corr(R, beta, betahat):
    r'''
    Returns correlation coefficient between a standard normal phenotype with 
    true SNP effect sizes `beta` and polygenic scores calculated with `betahat`
    as the polygenic score coefficients
    '''
#    corr = (beta.T@R@betahat)/np.sqrt(betahat.T@R@betahat)
    corr = np.squeeze((beta.T@R@betahat)**2/(betahat.T@R@betahat))
    return corr
    
def get_y(yg, h2):
    r'''
    Returns standardized phenotype using genetic component `yg` and heritability `h2`
    '''
    e = np.random.normal(size=yg.shape) # wait to scale by 1-h2
    e *= np.sqrt(1-h2)/e.std()
    y = yg + e
    y -= y.mean()
    y /= y.std()
    return y

def plot_comparison(x,y, xlabel='', ylabel='', title=None, save=True, logscale=False):
    plt.figure(figsize=(6,4))
    plt.plot(x, y, '.')
    r = np.corrcoef(x,y)[0,1]
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title((title if title is not None else f'{xlabel} v {ylabel}') + f'\n(r={round(r,3)}, r^2={round(r**2,3)})')
    if save:
        plt.savefig(f'/Users/nbaya/Downloads/{xlabel}_vs_{ylabel}.{title}.png',dpi=300)

def plot_break_pt_length(break_pts, chrom_list):
    for break_pts_chrom, chrom in zip(break_pts, chrom_list):
        break_pts_diff = np.diff(break_pts_chrom)
        mean = break_pts_diff.mean()
        plt.figure()
        n, _, _ = plt.hist(break_pts_diff, 50)
        plt.plot([mean]*2, [0, n.max()], 'k--')
        plt.xlim(0)
        plt.xlabel('LD block length (measured in consec. variants)')
        plt.ylabel('Density')
        plt.title(f'Distribution of LD block lengths (chr{chrom})\nmean: {round(mean,3)}')

def main():

    seed = None
    ref_panel = '1kg_eur'
    make_test_cohort = False
    chrom_list = [20,21,22]
    X = get_X(ref_panel=ref_panel, chrom_list=chrom_list)
    
#    chrom_break_pts = get_block_idxs(X, transpose=True)

    if make_test_cohort:
        X_all = X
        test_frac = 0.2
        test_idx = np.random.choice(X.shape[0], round(X.shape[0]*test_frac), replace=False)
        discovery_idx = [x for x in range(X.shape[0]) if x not in test_idx]
        X_test = [X_chrom[test_idx, :] for X_chrom in X_all]
        X = [X_chrom[discovery_idx, :] for X_chrom in X_all]
        
    M = sum([X_chrom.shape[1] for X_chrom in X])
    
    corr_dict = {}
    
    diag_block_size_list = [50,10,5]
    
    reps = 1
    h2=0.5
    model='inf'
    
    if model=='inf':
        pi=None
    elif model=='spike':
        pi=0.001
#    get_rep_id = lambda: f'corr_Xbeta_Xalpha.h2_{h2}.{model}{f".pi_{pi}" if pi!=None else ""}.blocksize_{diag_block_size}.reps_{reps}'
#    get_rep_id = lambda: f'corr_y_Xtestalphahat.ntest_{X_test.shape[0]}.Nd_{N_d}.h2_{h2}.{model}{f".pi_{pi}" if pi!=None else ""}.blocksize_{diag_block_size}.reps_{reps}'
#    get_rep_id = lambda: f'corr_Xtestalpha_Xtestalphahat.ntest_{X_test.shape[0]}.Nd_{N_d}.h2_{h2}.{model}{f".pi_{pi}" if pi!=None else ""}.blocksize_{diag_block_size}.reps_{reps}'
    for diag_block_size in diag_block_size_list:            
#        if get_rep_id() in corr_dict.keys():
#            pass
#        else:
#            print(f'diag block size: {diag_block_size}')
#        diag_block_size = 1
#        break_pts = [[x*diag_block_size for x in range(round(chrom.shape[1]/diag_block_size))] for chrom in X] # uniform breakpoints along bim file
        peak_radius = 500000 # decrease to get smaller LD blocks, increase for larger LD blocks (default: 500000 base pairs)
        max_ldblk_len = peak_radius*10 # decrease to get smaller LD blocks, increase for larger LD blocks (default: 10*peak_radius)
        first_position_dict = create_ld_blocks(chrom_list=chrom_list, 
                                               peak_radius=peak_radius,
                                               max_ldblk_len=max_ldblk_len)
        break_pts = [convert_breakpoints(chrom=chrom, X=X, first_position_dict=first_position_dict) for chrom in chrom_list]
        plot_break_pt_length(break_pts, chrom_list)
        
#        %%timeit
        R = get_sparse_R(X=X, break_pts=break_pts) # block_size=5, chr=21,22: {as_list=False: 32 s, as_list=True: 160 ms} -> 200x faster when as_list=True
        
#        R = [[np.ones(shape=1)]*X_chrom.shape[1] for X_chrom in X] # for testing only, this LD matrix is the identity
        
        corr_list = []
        for i in range(reps):
        
#            beta = get_beta(M=M, h2=h2, X=X, pi=pi, seed=seed)
            # sim 1
            b1 = get_beta(M=M, h2=0.5, X=X, pi=0.05, seed=1)
            b2 = get_beta(M=M, h2=0.5, X=X, pi=0.05, seed=2)
            b3 = normalize_beta(X=X, beta=b1+b2, h2=0.5)#get_beta(M=M, h2=0.5, X=X, pi=0.05, seed=seed)
            
            # sim 2
            betas = [get_beta(M=M, h2=0.5, X=X, pi=0.001, seed=seed) for seed in range(3,50+3)]
            betas = np.array(betas)
            b1 = normalize_beta(X=X, beta=betas[:10,:].sum(axis=0), h2=0.9)
            b2 = normalize_beta(X=X, beta=betas.sum(axis=0), h2=0.2)
            np.random.seed(seed=54)
            b3_weights = np.random.uniform(size=b2.shape)
            b3 = normalize_beta(X=X, beta=b3_weights*betas.sum(axis=0), h2=0.4)
            beta = b1

                        
            
#            %%timeit
            alpha = get_alpha(R=R, beta=beta) # block_size=5, chr=21,22: {R as_list=False: 936 µs, R as_list=True: 59.3 ms} -> 63x faster when R is coo_matrix
            alpha2 = get_alpha(R=R, beta=b2)
            alpha3 = get_alpha(R=R, beta=b3)
            
            plot_comparison(x=beta, 
                            y=alpha, 
                            xlabel='beta',
                            ylabel='alpha',
                            title=f'peak_radius={peak_radius}, max_ldblk_len={max_ldblk_len}, {model}{f", pi={pi}" if pi is not None else ""}')
#                            title=f'block_size={diag_block_size}, {model}{f", pi={pi}" if pi is not None else ""}')
            
#            Z = get_Z(N_r=N_r, seed=seed)
            
#            N_d = 1000000
            
#            alphahat = get_alphahat(alpha=alpha, N_d=N_d, N_r=N_r, X=X, Z=Z)
            
            alphahat = get_alphahat(alpha=alpha, N_d=100000, R=R)
            alphahat2 = get_alphahat(alpha=alpha2, N_d=100000, R=R)
            alphahat3 = get_alphahat(alpha=alpha3, N_d=100000, R=R)
            
#            np.savetxt('/Users/nbaya/Downloads/beta1.h2_0.5.pi_0.05.Nd_100000.seed_1.tsv',alphahat, delimiter='\t')
#            np.savetxt('/Users/nbaya/Downloads/beta2.h2_0.5.pi_0.05.Nd_100000.seed_2.tsv',alphahat, delimiter='\t')
#            np.savetxt('/Users/nbaya/Downloads/beta3.h2_0.5.Nd_100000.tsv',alphahat, delimiter='\t')
            
            np.savetxt('/Users/nbaya/Downloads/sim2.v2.beta1.h2_0.9.Nd_100000.tsv',alphahat, delimiter='\t')
            np.savetxt('/Users/nbaya/Downloads/sim2.v2.beta2.h2_0.2.Nd_100000.tsv',alphahat2, delimiter='\t')
            np.savetxt('/Users/nbaya/Downloads/sim2.v2.beta3.h2_0.4.Nd_100000.tsv',alphahat3, delimiter='\t')
            
            for a_hat in [alphahat, alphahat2, alphahat3]:
                plt.plot(a_hat[:10000], '.', alpha=0.1)
                plt.yscale('symlog',linthreshy=1e-4)
            
            for N_d in map(int, [1e4, 1e5, 1e6, 1e7]):
                alphahat = get_alphahat(alpha=alpha, N_d=N_d, R=R)
                
                plot_comparison(x=alpha, y=alphahat, xlabel='alpha', ylabel='alphahat',
                                title=(f'{model}{f", pi={pi}" if pi is not None else ""}, N_d={N_d}'+
                                       f', peak_radius={peak_radius}, max_ldblk_len={max_ldblk_len}'))
#                                       f', block_size={diag_block_size}'))
                
            
            
            # sim 3
            betas = [get_beta(M=M, h2=0.5, X=X, pi=0.001, seed=seed) for seed in range(3,50+3)]
            betas = np.array(betas)
            b1 = normalize_beta(X=X, beta=betas[:10,:].sum(axis=0), h2=0.9)
            b2 = normalize_beta(X=X, beta=betas.sum(axis=0), h2=0.2)
            np.random.seed(seed=54)
            b3_weights = np.random.uniform(size=b2.shape)
            b3 = normalize_beta(X=X, beta=b3_weights*betas.sum(axis=0), h2=0.4)

            for seed, beta in enumerate(betas,3):
                print(seed)
                alpha = get_alpha(R=R, beta=beta)
                alphahat = get_alphahat(alpha=alpha, N_d=100000, R=R)
                np.savetxt(f'/Users/nbaya/Downloads/sim3.beta_X_Y1.h2_0.5.pi_0.001.Nd_100000.seed_{seed}.tsv',alphahat, delimiter='\t')
                
            alpha1 = get_alpha(R=R, beta=b1) # block_size=5, chr=21,22: {R as_list=False: 936 µs, R as_list=True: 59.3 ms} -> 63x faster when R is coo_matrix
            alpha2 = get_alpha(R=R, beta=b2)
            alpha3 = get_alpha(R=R, beta=b3)
            
            alphahat1 = get_alphahat(alpha=alpha1, N_d=100000, R=R)
            alphahat2 = get_alphahat(alpha=alpha2, N_d=100000, R=R)
            alphahat3 = get_alphahat(alpha=alpha3, N_d=100000, R=R)
            
            np.savetxt('/Users/nbaya/Downloads/sim3.beta1.h2_0.9.Nd_100000.tsv',alphahat1, delimiter='\t')
            np.savetxt('/Users/nbaya/Downloads/sim3.beta2.h2_0.2.Nd_100000.tsv',alphahat2, delimiter='\t')
            np.savetxt('/Users/nbaya/Downloads/sim3.beta3.h2_0.4.Nd_100000.tsv',alphahat3, delimiter='\t')
            
            
#            yg = get_yg(X=X, beta=beta)
#            y = get_y(yg=yg, h2=h2)
##            
#            yhat = get_yg(X=X, beta=alpha)
##            
##            corr = get_corr(R=R, beta=beta, betahat=alphahat)
#            corr = get_corr(R=R, beta=beta, betahat=alpha)
##            
#
#            corr_list+=[corr]
            
            # out of sample prediction
#            yg = get_yg(X=X_test, beta=beta)
            
#            y = get_y(yg=yg, h2=h2)
            yg = get_yg(X=X, beta=alpha)
            yhat = get_yg(X=X, beta=alphahat)
            
#            corr_list+=[np.corrcoef(y, yhat)[0,1]]
            corr_list+=[np.corrcoef(yg, yhat)[0,1]]

        corr_dict.update({get_rep_id():corr_list})
        
#    for diag_block_size in diag_block_size_list:
#        corr_list = corr_dict[f'corr_Xbeta_Xalpha.h2_{h2}.{model}{f".pi_{pi}" if pi!=None else ""}.blocksize_{diag_block_size}.reps_{reps}']
#        plt.figure(figsize=(6,4))
#        plt.hist(np.asarray(corr_list)**2, 50)
#        plt.title(f'h2={h2}, model={model}{f", pi={pi}" if pi!=None else ""}, block_size={diag_block_size}\n(mean={round(np.mean(corr_list),3)}, reps={reps})')
#        plt.xlabel('Corr(X@beta, X@alpha)^2')
#        plt.ylabel('density')
#        plt.savefig(f'/Users/nbaya/Downloads/corr_Xbeta_Xalpha.h2_{h2}.{model}{f".pi_{pi}" if pi!=None else ""}.blocksize_{diag_block_size}.reps_{reps}.png',dpi=300)
        
    plt.figure(figsize=(6,4))
    bins = np.linspace(0,1,101)
    for diag_block_size in diag_block_size_list:
        corr_list = corr_dict[get_rep_id()]
        plt.hist(np.asarray(corr_list)**2, bins=bins, alpha=0.8)
        plt.title(f'h2={h2}, model={model}{f", pi={pi}" if pi!=None else ""}, N_d={N_d}\n(reps={reps})')
#        plt.xlabel('Corr(X@beta, X@alpha)^2')
#        plt.xlabel('Corr(X_test@beta, X_test@alphahat)^2')
        plt.xlabel('Corr(X_test@alpha, X_test@alphahat)^2')
        plt.ylabel('density')
        plt.xlim([0,1])
        plt.legend(diag_block_size_list)
    plt.savefig('/Users/nbaya/Downloads/'+get_rep_id().replace(f'.blocksize_{diag_block_size}','')+'.png',dpi=300)
    
    
    

if __name__=="__main__":

    main()