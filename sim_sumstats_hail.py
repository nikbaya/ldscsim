#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 11:41:43 2020

Simulates GWAS summary statistics for a standardized quantitative trait.

Input:
- Genotype data for reference samples
- Hail Block Matrices of LD matrix of reference samplesin block matrix form 
  (a list of Block Matrices, one for each LD block) 

@author: nbaya
"""
#
import hail as hl
from hail.linalg import BlockMatrix
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

bucket= 'gs://nbaya/risk_gradients'

ref_panel = '1kg_eur'

wd = '/Users/nbaya/Documents/lab/risk_gradients/data'

autosomes = range(1,23)

def import_from_plink(ref_panel):

    mt = hl.import_plink(bed=f'{wd}/{ref_panel}.bed',
                         bim=f'{wd}/{ref_panel}.bim',
                         fam=f'{wd}/{ref_panel}.fam')
    X = hl.linalg.BlockMatrix.from_entry_expr(mt.GT.n_alt_alleles())
    X = X.T
    
    X.write(f'{wd}/{ref_panel}.X.bm', overwrite=True)

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

def get_N(X : list):
    r'''
    Gets length of first dimension of arrays in list of arrays `X`. 
    This checks that the length of the first dimension is the same for all arrays
    in the list.
    '''
    N = X[0].shape[0]
    assert all([x.shape[0]==N for x in X]), 'ERROR: Sample count varies across chromosomes'
    return N

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

def get_toy_R(M, n_blocks, identity=False):
    r'''
    Creates "toy" LD matrix as a list of Hail Block Matrices for testing purposes.
    The list has length=`n_blocks`.
    '''
    R = []
    block_snp_idxs = np.array_split(range(M),n_blocks)
    block_sizes = [len(block) for block in block_snp_idxs]
    for block_size in block_sizes:
        if identity:
            R_block = BlockMatrix.from_numpy(np.identity(n=block_size))
        else:
            A = np.random.uniform(low=-1, high=1, size=(block_size,1))**11 # exponentiate to number (odd to preserve negative sign) to avoid highly correlated SNPs
            cov = A@A.T
            np.fill_diagonal(cov, 1)
            R_block = BlockMatrix.from_numpy(cov)
        R.append(R_block)
    return R

def _get_len_R(R):
    r'''
    Returns number of variants in `R`, a list of Block Matrices representing
    the LD matrix in block matrix format.
    '''
    return sum([R_block.shape[0] for R_block in R])

    
def _cast(array, astype):
    r'''
    Casts `array` to specificed type
    '''
    assert type(array) in {np.ndarray, BlockMatrix}, 'type of array is not supported'
    if astype in {np.ndarray, 'numpy', 'np'}:
        if isinstance(array, np.ndarray):
            return array
        elif isinstance(array, BlockMatrix):
            return array.to_numpy()
    elif astype in {BlockMatrix, 'BlockMatrix', 'bm'}:
        if isinstance(array, BlockMatrix):
            return array
        elif isinstance(array, np.ndarray):
            return BlockMatrix.from_numpy(array)

def get_ref_X(ref_panel, overwrite=False):
    r'''
    Returns N_ref x M dim matrix of column-standardized genotypes of LD ref panel
    '''
    X_bm_path = f'{bucket}/{ref_panel}.X.bm'
    
    if overwrite or not hl.hadoop_is_file(f'{X_bm_path}/_SUCCESS'):
        mt = hl.import_plink(bed=f'{bucket}/{ref_panel}.bed',
                             bim=f'{bucket}/{ref_panel}.bim',
                             fam=f'{bucket}/{ref_panel}.fam')
        
        mt = mt.annotate_rows(stats = hl.agg.stats(mt.GT.n_alt_alleles()))
        mt = mt.annotate_entries(X = (mt.GT.n_alt_alleles()-mt.stats.mean)/mt.stats.stdev)
        
        X = BlockMatrix.from_entry_expr(mt.X)
        X = X.T
        
        X.write(f'{bucket}/{ref_panel}.X.bm', overwrite=True)
        
    X = BlockMatrix.read(X_bm_path)
    
    return X

def get_beta(M, h2, pi=1, seed=None, astype=BlockMatrix):
    r'''
    Returns M-dim vector of true SNP effect sizes
    '''
    assert pi>=0 and pi<=1, '`pi` (proportion of causal variants) must be in the interval [0,1]'
    kwargs = {'seed':seed} if type(seed)!=type(None) else {}
    if pi==1: #infinitesimal model
        beta = np.sqrt(h2/M)*BlockMatrix.random(n_rows=M, 
                                                n_cols=1, 
                                                gaussian=True,
                                                **kwargs)
    else:
        np.random.seed(seed)
        M_causal= round(M*pi) # number of causal variants
        causal_beta = np.random.normal(loc=0,scale=np.sqrt(h2/M_causal),size=(M_causal,1))
        causal_idx = np.random.choice(M, replace=False, size=M_causal)
        causal_idx.sort()
        beta = np.zeros(shape=(M,1))
        beta[causal_idx] = causal_beta

    return _cast(beta, astype=astype)

def get_Z(N_r):
    r'''
    Returns `N_r`-dim standard normal random vector
    '''
    Z = BlockMatrix.random(n_rows=N_r,n_cols=1,gaussian=True) # N_r-dimensional standard normal random vector
    return Z

def get_yg(X, beta):
    r'''
    Returns genetic component of trait
    '''
    yg = X@beta
    return yg

#def get_rectangles(aldi_break_pts, M):
#    r'''
#    Get rectangles for `sparsify_rectangles`
#    '''
#    if not 0 in aldi_break_pts:
#        aldi_break_pts.insert(0,0)
#    if not M in aldi_break_pts:
#        aldi_break_pts.append(M)
#    rectangles = [[x1,x2]*2 for x1, x2 in zip(aldi_break_pts[:-1], aldi_break_pts[1:])]    
#    rectangles.append([M]*4)
#    return rectangles
#
#
#def get_R(X, aldi_break_pts):
#    r'''
#    Returns M x M sparse LD matrix
#    '''
#    N_r, M = X.shape
#    rectangles = get_rectangles(aldi_break_pts=aldi_break_pts, M=M)  
#    R = (X.T@X).sparsify_rectangles(rectangles)
#    R = R/N_r
#    return R

def get_alpha(R, beta, astype=BlockMatrix):
    r'''
    Returns M-dim vector of true marginal SNP effect sizes
    '''
    assert type(R) in {list}
    alpha_list = []

    if isinstance(R,list):
        assert _get_len_R(R)==beta.shape[0], 'Number of variants in `R` and `beta` do not match'
        start_idx = 0
        for R_block in R:
            stop_idx = start_idx+R_block.shape[0]
            alpha_block = R_block@beta[start_idx:stop_idx,:]
            
            alpha_list.append(alpha_block.to_numpy())
            start_idx = stop_idx
        alpha = np.concatenate(alpha_list)
        return _cast(alpha, astype=astype)

    elif isinstance(R,BlockMatrix):
        return R@beta

def old_get_alphahat(alpha, N_d, N_r, X, Z):
    r'''
    Returns M-dim vector of estimated marginal SNP effect sizes
    '''
    alphahat = alpha + 1/np.sqrt(N_d*N_r)*X.T@Z
    return alphahat


def get_alphahat(alpha, N_d, R, astype=BlockMatrix, seed=None):
    r'''
    Returns M-dim vector of estimated marginal SNP effect sizes
    '''
    np.random.seed(seed)
    alphahat_list = []
    start_idx=0
    for R_block in R:
        stop_idx = start_idx+R_block.shape[0]
        alpha_block = np.squeeze(alpha[start_idx:stop_idx,:].to_numpy())
        alphahat_block = np.random.multivariate_normal(mean=alpha_block, 
                                                       cov=(1/N_d)*R_block.to_numpy())
        alphahat_list.append(alphahat_block)
        start_idx = stop_idx
    alphahat = np.concatenate(alphahat_list)
    alphahat = alphahat.reshape((alpha.shape[0],1))
    return _cast(alphahat, astype=astype)

def checkpoint_bm(bm, path, read_if_exists=True):
    if not hl.hadoop_is_file(f'{path}/_SUCCESS'):
        bm.write(path)
    bm = BlockMatrix.read(path)
    return bm

def get_ld_scores(R, block_matrix=True):
    ld_scores = []
    for R_block in R:
        ld_scores_block = ((R_block.to_numpy() if block_matrix else R_block)**2).sum(axis=0)-1
        ld_scores.append(ld_scores_block)
    ld_scores = np.concatenate(ld_scores, axis=0)
    return ld_scores
#    
#    
#def E():
#    Z = np.random.normal(0,1,size=N)
#    Z = BlockMatrix.from_numpy(Z)
#    Z = Z.T
#    return (1/np.sqrt(N))*(X.T)@Z
#
#for h2 in [0.5]: #np.linspace(0.1, 0.5, 4):
#
#    # infinitesimal
#    beta = get_beta(M=M, )
#
#    yg = get_yg(X=X, beta=beta)
#    
#    s2 = (1/N)*(yg.T)@(yg)
#    
#    alpha *= np.sqrt(h2/s2) # comment out later?
#    
#    beta = (1/N)*X.T@yg
#                    
#    print(f'var(alpha)*M/h2 = {alpha.var()*M/h2}')
#    print(f'var(beta)*M/h2 = {beta.var()*M/h2}')
#    
#    print(f'corr(alpha, beta) = {np.corrcoef(alpha, beta)[0,1]}')
#    

def main():
    

    ref_panel = '1kg_eur'
    make_test_cohort = True
    chrom_list = [20,21,22]
    X = get_X(ref_panel=ref_panel, chrom_list=chrom_list)

    N = get_N(X)

    if make_test_cohort:
        X_all = X
        test_frac = 0.2
        test_idx = np.random.choice(N, round(N*test_frac), replace=False)
        discovery_idx = [x for x in range(N) if x not in test_idx]
        X_test = [X_chrom[test_idx, :] for X_chrom in X_all]
        X = [X_chrom[discovery_idx, :] for X_chrom in X_all]
        
    M = sum([X_chrom.shape[1] for X_chrom in X])
#    M = 10000
    
    h2 = 0.5 # SNP heritability of trait
    pi = 1 # 1: infinitesimal model, <1 : spike & slab
    seed = 1
    
    beta = get_beta(M=M, h2=h2, pi=pi, seed=seed) # true SNP effect sizes
    beta_np = get_beta(M=M, h2=h2, pi=pi, seed=seed, astype=np.ndarray)
    
#    n_blocks = 1600
#    ld_type = 'random'
#    R = get_toy_R(M=M, n_blocks=n_blocks, identity=(ld_type=='identity')) #R_random R_identity #R_identity

    peak_radius = 5000 # decrease to get smaller LD blocks, increase for larger LD blocks (default: 500000 base pairs)
    max_ldblk_len = peak_radius*10 # decrease to get smaller LD blocks, increase for larger LD blocks (default: 10*peak_radius)
    first_position_dict = create_ld_blocks(chrom_list=chrom_list, 
                                           peak_radius=peak_radius,
                                           max_ldblk_len=max_ldblk_len)
    break_pts = [convert_breakpoints(chrom=chrom, X=X, first_position_dict=first_position_dict) for chrom in chrom_list]
    R = get_sparse_R(X=X, break_pts=break_pts)
    R_np = [R_block for R_chrom in R for R_block in R_chrom ]
    ld_scores = get_ld_scores(R_np, block_matrix=False)
    R = [BlockMatrix.from_numpy(R_block) for R_block in R_np]
    
    plt.figure(figsize=(6*1.5,4*1.5))
    plt.plot(ld_scores)
    plt.title(f'peak radius: {peak_radius}')
#    R_path = f'{bucket}/{ref_panel}.R.h2_{h2}.{model}.seed_{seed}.bm'
#    R = checkpoint_bm(bm=R, 
#                      path=R_path, 
#                      read_if_exists=True)
    
    alpha = get_alpha(R, beta) # M-dimensional vector of true marginal SNP effect sizes
    alpha_np = alpha.to_numpy() # M-dimensional vector of true marginal SNP effect sizes
    
    r_beta_alpha = np.corrcoef(beta_np.T,alpha_np.T)[0,1]
    
    plt.figure(figsize=(1.5*6, 1.5*4))
#    plt.plot(beta_np, alpha_np,'.')
    plt.scatter(beta_np, alpha_np, s= ld_scores/ld_scores.max()*20, alpha=0.1)
    plt.plot(*[[beta_np.min(), beta_np.max()]]*2, 'k--')
    plt.xlabel('beta')
    plt.ylabel('alpha')
    plt.title(f'r2={round(beta_alpha_r**2,3)}\nh2: {h2}, pi: {pi}, LD type: {ld_type}, {n_blocks} blocks (seed: {seed})')
    
    N_d = 100000
    
    alphahat = get_alphahat(alpha=alpha, N_d=N_d, R=R, seed=seed)
    alphahat_np = alphahat.to_numpy()
    
    r_alpha_alphahat = np.corrcoef(alpha_np.T,alphahat_np.T)[0,1]
    
    plt.figure(figsize=(1.5*6,1.5*4))    
#    plt.plot(alpha_np, alphahat_np, '.')
    plt.scatter(alpha_np, alphahat_np, s = ld_scores/ld_scores.max()*20, alpha=0.1)
    plt.plot(*[[alpha_np.min(), alpha_np.max()]]*2, 'k--')
    plt.xlabel('alpha')
    plt.ylabel('alphahat')
    plt.title(f'r2={round(r_alpha_alphahat**2,3)}\nh2: {h2}, pi: {pi}, N_d: {N_d}, LD type: {ld_type}, {n_blocks} n_blocks (seed: {seed})')
    
#    alphahat_path = f'{bucket}/{ref_panel}.alphahat.h2_{h2}.{model}.seed_{seed}.bm'
#    alphahat = checkpoint_bm(bm=alphahat, 
#                             path=alphahat_path, 
#                             read_if_exists=True)
    
#    alphahat_np = alphahat.to_numpy()
    
    yg = get_yg(X=X, beta=beta)    
    yg_np =np.squeeze(yg.to_numpy())
    
    yhat = get_yg(X=X, beta=alphahat)
    yhat_np = np.squeeze(yhat.to_numpy())
    
    print(yg_np[:10])
    print(yhat_np[:10])
    
    yg_yhat_corr = np.corrcoef(yg_np, yhat_np)[0,1]
    print(f'yg-yhat correlation: {yg_yhat_corr}')
    
    
    
if __name__=='__main__':

    
    main()
    
    
#    R = (1/N)*X.T@X
#    R_bm_path = f'{gs_bucket}/{ref_panel}.R.bm'
#    if not hl.hadoop_is_file(f'{R_bm_path}/_SUCCESS'):
#        R.write(R_bm_path)
#    R = BlockMatrix.read(R_bm_path)
    
    
    
#    E_bm_path = f'{gs_bucket}/{ref_panel}.E.bm'
#    if not hl.hadoop_is_file(f'{E_bm_path}/_SUCCESS'):
#        E.write(E_bm_path)
#    E = BlockMatrix.read(E_bm_path)
#        
#    seed = 1
#    np.random.seed(seed=seed)
#    h2 = 0.2
#    alpha = np.random.normal(loc=0, scale=np.sqrt(h2/M), size=M)
#    alpha = BlockMatrix.from_numpy(alpha).T
#    
#    beta = R@alpha
#    
#    beta_bm_path = f'{gs_bucket}/{ref_panel}.beta.bm'
#    if not hl.hadoop_is_file(f'{beta_bm_path}/_SUCCESS'):
#        beta.write(beta_bm_path, overwrite=True)
#    beta = BlockMatrix.read(beta_bm_path)
#    
#    N_d = 10000
#    betahat = beta + (1/np.sqrt(N_d))*E
#    
#    betahat_bm_path = f'{gs_bucket}/{ref_panel}.betahat.bm'
#    if not hl.hadoop_is_file(f'{betahat_bm_path}/_SUCCESS'):
#        betahat.write(betahat_bm_path, overwrite=True)
#    betahat = BlockMatrix.read(betahat_bm_path)
#        
#    
#    alphahat = betahat
#    r_prs = h2*(alphahat.T@R@alpha)/np.sqrt((alphahat.T@R@alphahat)*(alpha.T@R@alpha))
#    
#    
#    tb_r_prs = r_prs.to_table_row_major()
#    tb_r_prs.show()
#    
##    r_prs = r_prs.checkpoint(f'{gs_bucket}/tmp.1kg_eur.r_prs.bm',
##                             overwrite=True)
#    
#
#else:
#    
    
    
    
    
    
    
    
    
    
    