#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 12:08:42 2019

@author: nbaya
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from ldscsim_basic import create_cov_matrix
import seaborn as sns
from ldscsim_basic import simulate_phenotypes



fig, ax = plt.subplots(figsize=(12,9))
rgs = np.linspace(0,1,3)
for rg in rgs:
    cov_matrix = create_cov_matrix(h2=[0.01,0.9],rg=[rg])
    M = int(1e6)
    n_phens =2
    cov_matrix = cov_matrix/M
    beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])    
#    beta0 = simulate_phenotypes(mt,mt.gt,h2=[0.5,0.5],rg=[rg]).y_no_noise.take(int(1e3))
#    beta0 = np.asarray(beta0)
    r=[]
    thresholds = np.linspace(min(beta0[:,1])*0.9999,max(beta0[:,1])*0.9999,51)
    for i in thresholds:
        beta = beta0.copy()
        beta[:,1] = np.where(beta[:,1]>i,1,0)
        r.append(stats.pearsonr(beta[:,0],beta[:,1])[0])
    plt.plot(thresholds, r, color=[0, 0, rg])
plt.legend([f'rg = {round(x,3)}' for x in rgs])





# Simulate two correlated spike & slab traits
# randomly set non-causal SNPs

for pi_trait2 in np.logspace(-2,-1,11):
    fig, ax = plt.subplots(figsize=(12,9))
    rgs = np.linspace(0,1,11)
    for rg in rgs:
        cov_matrix = create_cov_matrix(h2=[0.5,0.5],rg=[rg])
        M = int(1e5)
        n_phens =2
        cov_matrix = cov_matrix/M
    
        r=[]
        pi_ls = np.linspace(1/50,1,50)
        for pi in pi_ls:
            beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])    
            beta0[np.random.choice(M,size=int((1-pi)*M),replace=False),0] = 0
            beta0[np.random.choice(M,size=int((1-pi_trait2)*M),replace=False),1] = 0
            r.append(stats.pearsonr(beta0[:,0],beta0[:,1])[0])
        ax.plot(pi_ls, r, color=[0, 0, rg])
    plt.legend([f'rg = {round(x,3)}' for x in rgs])
    plt.xlabel('pi for trait 1')
    plt.ylabel('rg')
    plt.title(f'pi trait 2 = {pi_trait2}')
    plt.ylim([0,1])
    
    
# Simulate two correlated spike & slab traits
# draw from multivariate normal, set a certain number of overlapping SNPs to be 
# non-causal

r_diff = []
rg_ls_ls = []
for start_rg in np.linspace(0.1,1,10)[::-1]:
    start_rg = round(start_rg,5)
    print(f'start rg: {start_rg}')
    cov_matrix = create_cov_matrix(h2=[0.1,0.1],rg=[start_rg])
    M = int(1e5)
    n_phens =2
    cov_matrix *= 1/M
    beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
    pi1 = 0.01
    pi2 = 0.2
    beta0[:,0] = (1/pi1)*beta0[:,0]
    beta0[:,1] = (1/pi2)*beta0[:,1]
    rg = 0.0
    beta = beta0.copy()
    beta[0:int(M*min(pi1,pi2)),:] = 0
    i=0
    rg_ls =[]
    prev_r = stats.pearsonr(beta[:,0],beta[:,1])[0]-rg
    beta = beta0.copy()
    if pi2>pi1:
        i=int(M*(pi2-pi1))
    beta[0:int(M*(1-pi1)),0] = 0
    beta[i:int(M*(1-pi2))+i,1] = 0
    r = stats.pearsonr(beta[:,0],beta[:,1])[0]
    i += 1
    while abs(r-rg)>0.01 and r > rg:
        beta = beta0.copy()
        beta[0:int(M*(1-pi1)),0] = 0
        if int(M*(1-pi2))+i>M:
            beta[0:(int(M*(1-pi2))+i)%M,1] = 0
            beta[i:M,1] = 0
        else:
            beta[i:int(M*(1-pi2))+i,1] = 0
        r = stats.pearsonr(beta[:,0],beta[:,1])[0]
#        print(r)
        rg_ls.append(r)
        if i>=M:#*(1-pi1) and r-rg>0.01:
#                print(f'converged before reaching rg = {rg}')
            print(f'***{start_rg}: {start_rg-r}***')
            r_diff.append(start_rg-r)
            break
#        prev_r = r
        i += 1
    rg_ls_ls.append(rg_ls)
#plt.plot(range(len(rg_ls)),rg_ls,'.-')
#plt.plot(range(len(rg_ls)),[start_rg]*len(rg_ls),'k--',alpha=0.5)
#plt.title(f'start rg = {start_rg}')

fig,ax=plt.subplots(figsize=(8,6))
for i,ls in enumerate(rg_ls_ls):
    ax.plot(range(len(ls)),ls,'-',color=[(i+1)/10,0,0],lw=5)
plt.legend([f'start_rg={round(x,2)}' for x in np.linspace(0.1,1,10)[::-1]])
plt.title(f'[pi1,pi2]={[pi1,pi2]}')



plt.plot(np.linspace(0.1,1,10)[::-1],r_diff,'.-')
plt.plot(np.linspace(0.1,1,10)[::-1],(max(r_diff)-min(r_diff))/0.9*np.linspace(0.1,1,10)[::-1],'.-')
print(f'slope for {[pi1,pi2]}: {(max(r_diff)-min(r_diff))/0.9}')

# plot max correlation given two pi values
n_h2 = 51
rmax = np.zeros(shape=(n_h2,n_h2))
pi1_array = np.logspace(-3,0,num=n_h2)
pi2_array = np.logspace(-3,0,num=n_h2)
#pi1_array = np.linspace(0.01,1,num=n_h2)
#pi2_array = np.linspace(0.01,1,num=n_h2)
for i1, pi1 in enumerate(pi1_array):
    for i2, pi2 in enumerate(pi2_array):   
        cov_matrix = create_cov_matrix(h2=[0.1,0.1],rg=[1])
        M = int(1e5)
        n_phens =2
        cov_matrix *= 1/M
        beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
        beta0[:,0] = (1/pi1)*beta0[:,0]
        beta0[:,1] = (1/pi2)*beta0[:,1]
        beta = beta0.copy()
        i=0
        if pi2>pi1:
            i=int(M*(pi2-pi1))
        beta[0:int(M*(1-pi1)),0] = 0
        beta[i:int(M*(1-pi2))+i,1] = 0
        r = stats.pearsonr(beta[:,0],beta[:,1])[0]
        rmax[i1,i2] = r
sns.heatmap(rmax)
plt.xticks(range(n_h2)[::10],[round(x,3) for x in pi2_array[::10]],rotation=20)
plt.yticks(range(n_h2)[::10],[round(x,3) for x in pi1_array[::10]],rotation=20)
plt.xlabel('pi2')
plt.ylabel('pi1')



# plot rmax as a function of max pi ratio (rmax = start_rg*ratio^(-1/2))
rmax = []
pi_ratio = []
#h2_ls =[]
start_rg = 0.1
pi1_array = 10**(-3*np.random.uniform(size=(n_h2))) #np.logspace(-3,0,num=n_h2)
pi2_array = 10**(-3*np.random.uniform(size=(n_h2))) #np.logspace(-3,0,num=n_h2)
for i1, pi1 in enumerate(pi1_array):
    for i2, pi2 in enumerate(pi2_array): 
        h2=np.random.uniform(size=2).tolist()
        cov_matrix = create_cov_matrix(h2=h2,rg=[start_rg])
        M = int(1e5)
        n_phens =2
        cov_matrix *= 1/M
        beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
        beta0[:,0] = (1/pi1)*beta0[:,0]
        beta0[:,1] = (1/pi2)*beta0[:,1]
        beta = beta0.copy()
        i=0
        if pi2>pi1:
            i=int(M*(pi2-pi1))
        beta[0:int(M*(1-pi1)),0] = 0
        beta[i:int(M*(1-pi2))+i,1] = 0
        r = stats.pearsonr(beta[:,0],beta[:,1])[0]
        pi_ratio.append(max(pi1/pi2,pi2/pi1))
        rmax.append(r)
#        h2_ls.append(max(h2))
plt.scatter(pi_ratio,rmax)
plt.plot(np.logspace(0,3,2),start_rg*np.logspace(0,3,2)**(-1/2),'r-',alpha=0.5)
plt.xscale('log')
plt.yscale('log')
#plt.xlim([1,max(pi_ratio)])
#plt.ylim([10**(-2),10**0])
plt.xlabel('max pi ratio')
plt.ylabel('max r')
    

#get number of iterations required to get correlation given target rg and start rg
#high number of iterations => less overlap in non-causal SNPs
#low number of iterations  => more overlap in non-causal SNPs
target_rg_ls = []
n_iter = []
pi_max= []
miss_ct = 0
for target_rg in np.linspace(0,1,10)[::-1]:
    for i1, pi1 in enumerate(np.logspace(-2,-1,10)):
        for i2, pi2 in enumerate(np.logspace(-2,-1,10)): 
            h2=np.random.uniform(size=2).tolist()
            rmax = target_rg*max(pi1/pi2,pi2/pi1)**(-1/2)
            if target_rg > rmax:
                miss_ct += 1
                next
            cov_matrix = create_cov_matrix(h2=h2,rg=[min(rmax+.1,1)])
            M = int(1e5)
            n_phens =2
            cov_matrix *= 1/M
            beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
            beta0[:,0] = (1/pi1)*beta0[:,0]
            beta0[:,1] = (1/pi2)*beta0[:,1]
            beta = beta0.copy()
            i=0
#            if pi2>pi1:
#                i=int(M*(pi2-pi1))
            beta[0:int(M*(1-pi1)),0] = 0
            beta[i:int(M*(1-pi2))+i,1] = 0
            r = stats.pearsonr(beta[:,0],beta[:,1])[0]
#            print(r)
            while abs(r-target_rg)>0.01 and r > target_rg:
                i += 1
                beta = beta0.copy()
                beta[0:int(M*(1-pi1)),0] = 0
                if int(M*(1-pi2))+i>M:
                    beta[0:(int(M*(1-pi2))+i)%M,1] = 0
                    beta[i:M,1] = 0
                else:
                    beta[i:int(M*(1-pi2))+i,1] = 0
                r = stats.pearsonr(beta[:,0],beta[:,1])[0]
                if i>=M:#*(1-pi1) and r-rg>0.01:
        #                print(f'converged before reaching rg = {rg}')
                    print(f'***{rmax+.1}: {start_rg-r}***')
                    print(f'pi={[pi1,pi2]}')
    #                r_diff.append(start_rg-r)
                    next
#                    break
        #        prev_r = r
            target_rg_ls.append(target_rg)
            n_iter.append(i)
            pi_max.append(max(pi1/pi2,pi2/pi1))
#plt.plot(target_rg_ls, n_iter,'.')
plt.scatter(target_rg_ls, n_iter,c=pi_max)    
print(miss_ct/(16*10*10))


#draft of method for correlated spike & slab

def multitrait_spikeslab(h2,pi,target_rg,M = int(1e5)):
    d = 1/M
#    print(f'attempting target_rg={target_rg} with error {d}, pi={pi}')
    pi1 = pi[0]
    pi2 = pi[1]
    rmax = max(pi1/pi2,pi2/pi1)**(-1/2) #rmax assuming rg between causal SNPs is 1
    if  rmax < target_rg:
        print(f'chosen pi={pi} cannot achieve rg = {target_rg}')
    elif target_rg==0: #independent spike & slab
        cov_matrix = create_cov_matrix(h2=h2,rg=[target_rg])
        cov_matrix *= 1/M
        beta0 = np.random.multivariate_normal(mean=np.zeros(cov_matrix.shape[0]),cov=cov_matrix,size=[M,])
        beta0 = beta0*([x**(-1) for x in pi])
        for i,p in enumerate(pi):
            beta0[0:int(M*(1-p)),i] = 0
        np.random.shuffle(beta0)
        return beta0
    elif target_rg==1: #perfect correlation
        cov_matrix = create_cov_matrix(h2=h2,rg=[target_rg])
        cov_matrix *= 1/M
        beta0 = np.random.multivariate_normal(mean=np.zeros(cov_matrix.shape[0]),cov=cov_matrix,size=[M,])
        beta0 = beta0*([x**(-1) for x in pi])
        beta0[0:int(M*(1-pi1)),0] = 0
        beta0[0:int(M*(1-pi2)),1] = 0
        return beta0
    elif max(pi)==1 and len(set(pi))==1: #multitrait infinitesimal
        cov_matrix = create_cov_matrix(h2=h2,rg=[target_rg])
        cov_matrix *= 1/M
        beta0 = np.random.multivariate_normal(mean=np.zeros(cov_matrix.shape[0]),cov=cov_matrix,size=[M,])
        return beta0
    else:
#        start_rg_ls = []
        mix=False
        if max(pi)==1:
            mix = True
            start_rg_prev = target_rg+(1-target_rg)*0.1
        else:
            start_rg = 1
            start_rg_prev = start_rg
            while rmax > target_rg+0.1:
    #            start_rg_ls.append(start_rg)
                start_rg_prev = start_rg
                start_rg = target_rg+(start_rg-target_rg)*0.99 #approach target rg exponentially
                rmax = (start_rg)*max(pi1/pi2,pi2/pi1)**(-1/2)
#        plt.plot(range(len(start_rg_ls)),start_rg_ls,'.')
#        plt.plot(range(len(start_rg_ls)),[target_rg]*len(start_rg_ls),'k--',alpha=0.5)
        success=False
        iteration=1
        flip = False
        while success is False and iteration <= 3:
            cov_matrix = create_cov_matrix(h2=h2,rg=[start_rg_prev])
            n_phens =2
            cov_matrix *= 1/M
            beta0 = np.random.multivariate_normal(mean=np.zeros(n_phens),cov=cov_matrix,size=[M,])
            beta0[:,0] = (1/pi1)*beta0[:,0]
            beta0[:,1] = (1/pi2)*beta0[:,1]
            beta = beta0.copy()
            if pi2>pi1:
                flip = (flip==False)
                pi1, pi2 = pi2, pi1
                beta = beta[:,::-1]
            i=0 #int(M*(pi2-pi1))
            beta[0:int(M*(1-pi1)),0] = 0
            beta[i:int(M*(1-pi2))+i,1] = 0
            r = stats.pearsonr(beta[:,0],beta[:,1])[0]
            i += 1
            if mix:
                print(f'{i}: {r}')
            while r-target_rg > d and i < M:
                beta = beta0.copy()
                beta[0:int(M*(1-pi1)),0] = 0
                if int(M*(1-pi2))+i>M: #if the window of zeros reaches the end of the array, wrap around
                    beta[0:(int(M*(1-pi2))+i)%M,1] = 0
                    beta[i:M,1] = 0
                else:
                    beta[i:int(M*(1-pi2))+i,1] = 0
                r = stats.pearsonr(beta[:,0],beta[:,1])[0]
                if i%1000 ==0:
                    print(r)
                i += 1
            if i>=M:#*(1-pi1) and r-rg>0.01:
                print(f'iteration {iteration} failed, retrying with updated initial conditions')
                print(f'start_rg_prev: {start_rg_prev},pi={[pi1, pi2]},rg={target_rg}')
                iteration += 1
                if r < target_rg:
                    start_rg_prev = start_rg_prev+(1-start_rg_prev)*0.9
                else:
                    start_rg_prev = start_rg_prev*0.1
            else:
                success=True
                if iteration > 1:
                    print(f'iteration {iteration} was successful!')
                return beta
        if iteration == 3:
            print('failed too many times')

exp_rg=[]
obs_rg=[]
fail_ct = 0
for pi1 in np.logspace(-2,0,5)[:-1]:
    for pi2 in [1]: #np.logspace(-2,0,5):
        for target_rg in np.linspace(0,0.5,6):
            print('\n')
            beta = multitrait_spikeslab(h2=[0.3,0.5],pi=[pi1,pi2],target_rg=round(target_rg,3))
            if beta is not None:
                print(f'pi={[pi1,pi2]}')
                print(f'target rg={round(target_rg,3)}, rg={stats.pearsonr(beta[:,0],beta[:,1])[0]}')
                exp_rg.append(round(target_rg,3))
                obs_rg.append(stats.pearsonr(beta[:,0],beta[:,1])[0])
            else:
                fail_ct += 1
print(f'success rate: {1-fail_ct/(10*10*11)}')

plt.plot(exp_rg,obs_rg,'.')
plt.plot([min(exp_rg),max(exp_rg)],[min(exp_rg),max(exp_rg)],'k--',alpha=0.5)

beta = multitrait_spikeslab(h2=[0.3,0.5],pi=[0.03162277660168379, 1.0],target_rg=0)
