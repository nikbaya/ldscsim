#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 16:09:22 2018

@author: nbaya
"""
import hail as hl

class SimDistribution():
    def __init__(self, h2):
        self.h2 = h2
        
    def make_noise(self, mt):
        return hl.rand_norm(0, 1-self.h2)

class InfinitesimalModel(SimDistribution): 
    def make_random_function(self, mt):
        M = mt.count_rows() # number of variants
        return hl.rand_norm(0, hl.sqrt(self.h2/M))                              #SQUARE ROOT?

class InfinitesimalSpikeSlabModel(SimDistribution):
    def __init__(self, h2, pi):
        super().__init__(h2)
        self.pi = pi
    
    def make_random_function(self, mt): #pi is slab prob
        M = mt.count_rows() # number of variants
#        return hl.cond(hl.rand_unif(0,1) < self.pi, hl.rand_norm(0,self.h2/(M*self.pi)), 0)
#        return hl.cond(hl.rand_bool(self.pi), hl.rand_norm(0,hl.sqrt(self.h2/(M*self.pi))), 0)
        return hl.rand_bool(self.pi) * hl.rand_norm(0,self.h2/(M*self.pi))
    
class AnnotationInformedModel(SimDistribution):
    def __init__(self, h2, f, annotations): 
        self.h2 = h2
        self.f = f #function for combining annotations
        self.a_ht = annotations #hail table holding all annotations
        
    def make_random_function(self, mt): 
        from functools import reduce
        #check that row key of annotations matches row key of mt
        mt = mt.add_row_index()
        rows = [rf for rf in self.a_ht.row]
        self.a_ht = self.a_ht.annotate(__a__ = reduce(self.f, map(lambda x: self.a_ht[rows[x]],range(len(rows)))))
        std=self.a_ht.aggregate(hl.agg.stats(self.a_ht.__a__)).stdev
        self.a_ht = self.a_ht.annotate(__a__ = self.a_ht.__a__*hl.sqrt(self.h2/std))
        return mt.annotate_rows(beta = hl.literal(self.a_ht.__a__.take(mt.count_rows()))[hl.int32(mt.row_idx)])
    
class AnnotationInformedModel_test(SimDistribution):
    def __init__(self, h2, annotations): 
        self.h2 = h2
#        self.f = f #function for combining annotations
#        self.a = annotations #hail table holding final annotations
        
    def make_random_function(self, mt): 
        mt = mt.annotate(beta = get_agg_annotations)
        #check that row key of annotations matches row key of mt
        return 
    
    def get_agg_annotations(self, f, anno_ht):
        from functools import reduce
        rows = [rf for rf in anno_ht.row]
        anno_ht= anno_ht.annotate(__a__ = reduce(f, map(lambda x: anno_ht[rows[x]],range(len(rows)))))
        return hl.literal(anno_ht.__a__.take(anno_ht.count()))
    
    def normalize_for_h2(self):
        std=self.a_ht.aggregate(hl.agg.stats(self.a_ht.__a__)).stdev
        self.a_ht = self.a_ht.annotate(__a__ = self.a_ht.__a__*hl.sqrt(self.h2/std))

def simulate(mt, beta_distribution, eps_distribution):
    beta_dist = beta_distribution.make_random_function(mt)
    epsilon_dist  = beta_distribution.make_noise(mt)
    
    if beta_distribution==AnnotationInformedModel:
        mt = beta_dist
    
    mt = mt.annotate_rows(beta = beta_dist)
    
    mt = mt.annotate_cols(y_no_noise = hl.agg.sum(mt.beta * mt.dosage))
    mt = mt.annotate_cols(y = epsilon_dist + mt.y_no_noise)
    
    return mt




#from gs://nbaya/split/ukb31063.dosage.autosomes.hm3_variants.head100.mt 
#mt = hl.read_matrix_table('/Users/nbaya/Documents/lab/ldscsim/ukb31063.dosage.autosomes.hm3_variants.head100.mt')
mt = hl.read_matrix_table('/Users/nbaya/Documents/lab/ldscsim/ukb31063.hm3_variants.gwas_samples_50.batch_1.downsample.mt')

import pandas as pd
import numpy as np


anno = pd.DataFrame({'a1':[2]*100, 'a2':np.random.normal(0,1,100), 'a3':np.random.normal(0,1,100)<0,'a4':np.random.normal(0,1,100)})
ht = hl.Table.from_pandas(anno)

#a_ht = ht.add_index().rename({'idx':'idx1'}).add_index()
#a_mt = a_ht.to_matrix_table(row_key='idx',col_key='idx1')

f_sum = lambda x,y: x+y 
f_mul = lambda x,y: x*y

from functools import reduce
h2 = 0.1
rows = [rf for rf in ht.row]
a_ht = ht.annotate(__a__ = reduce(f_mul, map(lambda x: ht[rows[x]], range(len(rows))))) #
std =a_ht.aggregate(hl.agg.stats(a_ht.__a__)).stdev
a_ht.__a__*hl.sqrt(h2/std)
a_ht = a_ht.annotate(__a__ = a_ht.__a__*hl.sqrt(h2/std))

mt = mt.add_row_index()


mt = mt.annotate_rows(beta = hl.literal(a_ht.__a__.take(mt.count_rows()))[hl.int32(mt.row_idx)])
mt.row.show()


import datetime
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

model1 = InfinitesimalModel(h2=0)
model1.make_random_function(mt)
#
#model2 = InfinitesimalSpikeSlabModel(h2=0.1,pi=0.01)

model3 = AnnotationInformedModel(h2=0.1,f=f_sum, annotations = ht)
beta_dist = model3.make_random_function(mt)

beta_dist.row.show()
#starttime = datetime.datetime.now()
sim_mt = simulate(mt, model1, model1)

beta = np.asarray(sim_mt.beta.take(sim_mt.count_rows()))

np.sum(beta**2)

plt.subplot(1,2,1)
sns.kdeplot(beta)
plt.title('beta')
plt.subplot(1,2,2)
y = sim_mt.y.take(60100)
print(np.std(y))
sns.kdeplot(y)
plt.title('y')
fig = plt.gcf()
fig.set_size_inches(12, 4)
#endtime = datetime.datetime.now()
#elapsed = endtime-starttime
#print('Iteration time: '+str(round(elapsed.seconds/60, 2))+' minutes')

###############################################################################
"""
Gamma distribution phenotype
"""
mt = hl.read_matrix_table('/Users/nbaya/Documents/lab/ldscsim/ukb31063.hm3_variants.gwas_samples_50.batch_1.downsample.mt')

mt1=mt.filter_cols(mt.group_id<1).sample_rows(0.01)
M = mt1.count()[0] #number of rows
N = mt1.count()[1] #number of cols

X_vec = mt1.dosage.take(M*N) #get dosages
X = np.reshape(X_vec,(M,N)).T #reshape dosages into matrix

# Test different betas (DO NOT EDIT) ##########################################
beta =np.asarray([np.asarray(i) for i in 1/(np.mean(X,0,keepdims=True))**2]).T
noise = np.asarray(np.random.normal(0,0.5,size=[M,1]).tolist())
beta = (beta)/100
y = np.matmul(X, beta)
plt.subplot(2,1,1)
plt.hist(beta.flatten(),50)
plt.title('beta distribution')
plt.subplot(2,1,2)
plt.hist(y.flatten(),50)
plt.title('phenotype distribution')
fig=plt.gcf()
fig.set_size_inches(6,8)
fig.savefig('/Users/nbaya/Desktop/beta_phenotype_dist.png')

# Test different betas (DO NOT EDIT) ##########################################
beta_gamma = np.random.gamma(0.05,1,size=[M,1])*10
beta_gamma = beta_gamma[np.argsort(beta_gamma,axis=None)]
#beta_gamma = beta_gamma - beta_gamma%1
beta_from_dosage =np.asarray([np.asarray(i) for i in -(np.mean(X,0,keepdims=True))]).T
beta = beta_gamma[np.argsort(np.argsort(beta_from_dosage,axis=None),axis=None)]
y = np.matmul(X, beta)

ax1 = plt.subplot(2,2,1)
ax1.hist(beta.flatten(),50)
plt.title('beta distribution')
ax2 = plt.subplot(2,2,2)
ax2.hist(y.flatten(),50)
plt.title('phenotype distribution')
ax3 = plt.subplot(2,1,2)
ax3.set_xlabel('index')
ax3.set_ylabel('beta', color='tab:blue')
ax3.plot(range(len(beta_gamma)),beta_gamma,'.', color='tab:blue')
ax3.tick_params(axis='y', labelcolor='tab:blue')
ax4 = ax3.twinx()  # instantiate a second axes that shares the same x-axis
ax4.set_ylabel('dosage', color='tab:red')  # we already handled the x-label with ax1
ax4.plot(range(len(beta_from_dosage)),-beta_from_dosage[np.argsort(beta_from_dosage,axis=None)],'.', color='tab:red')
ax4.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig=plt.gcf()
fig.set_size_inches(8,8)
fig.savefig('/Users/nbaya/Desktop/gammadist_beta_gammadist_phenotype.png')

# Test different betas (DO NOT EDIT) ##########################################
# Doesn't always give good results
beta_gamma = np.random.normal(0,0.2,size=[M,1])**3 
pi = 0.5
beta_gamma = [[i] for i in (np.random.rand(M)<0.1)]* np.random.normal(0,1/(M*pi),size=[M,1]) 
beta_gamma = beta_gamma[np.argsort(beta_gamma,axis=None)]
beta_from_dosage =np.asarray([np.asarray(i) for i in 1/(np.mean(X,0,keepdims=True))]).T
beta = beta_gamma[np.argsort(np.argsort(beta_from_dosage,axis=None),axis=None)]+ np.random.normal(0,0.001,size=[M,1])+0.0011
y = np.matmul(X, beta)
plt.subplot(2,1,1)
plt.hist(beta.flatten(),20)
plt.title('beta distribution')
plt.subplot(2,1,2)
plt.hist(y.flatten(),20)
plt.title('phenotype distribution')
fig=plt.gcf()
fig.set_size_inches(6,8)

# Test different betas ########################################################
#beta_init = [i/abs(i) for i in (np.random.normal(0,1,size=[M,1]))]*np.random.exponential(0.1,size=[M,1])
beta_init = np.random.normal(0,0.2,size=[M,1])
beta_init = beta_init[np.argsort(beta_init,axis=None)]
beta_from_dosage =np.asarray([np.asarray(i) for i in 1/(np.mean(X,0,keepdims=True))]).T
beta = beta_init[np.argsort(np.argsort(beta_from_dosage,axis=None),axis=None)]+ np.random.normal(0,0.01,size=[M,1])+0.023
y = np.matmul(X, beta)
plt.subplot(2,1,1)
plt.hist(beta.flatten(),50)
plt.title('beta distribution')
plt.subplot(2,1,2)
plt.hist(y.flatten(),50)
plt.title('phenotype distribution')
fig=plt.gcf()
fig.set_size_inches(6,8)



# Test lstsq for beta with normal distribution (Done)
beta = np.random.normal(0.3,1,size=[M,1])
y = np.dot(X, beta)
y.shape
plt.hist(y.flatten())
beta_pred = np.linalg.lstsq(X, y)[0]
np.linalg.norm(beta-beta_pred)



#beta = np.random.beta(1,3,size=[M,1])
    
beta_pred = np.linalg.lstsq(X, y)[0]
np.linalg.norm(beta-beta_pred)


# Test lstsq for y with gamma dist

y = np.random.gamma(1,3,size=[1202,1])
y = np.random.normal(0,3,size=[1202,1])
y_0 = np.random.normal(0,2,size=[1202,1])
y= np.matmul(np.diag(y_0.flatten(),0),y_0)
plt.hist(y.flatten(),20)
#X = np.random.randint(0,3,size=[1000,1000])
#X = np.random.uniform(0,2,size=[1000,1000])
plt.hist(X.flatten(),20)
beta = np.matmul(np.linalg.pinv(X),y)
beta = np.linalg.lstsq(X, y)[0]
y_pred = np.matmul(X, beta)


gamma_dist = beta.flatten()
norm_dist = beta.flatten()
sns.kdeplot([norm_dist,gamma_dist])
plt.hist([norm_dist, gamma_dist],20)
plt.legend(['normal dist','gamma dist'])
plt.hist(np.matmul(X,beta).flatten())

plt.hist(y,20)
sns.kdeplot(y.flatten(), y_pred.flatten(), shade=False)
plt.ylabel('y_pred')
plt.xlabel('y')
plt.xlim([-2, 10])
plt.ylim([-2, 10])
fig=plt.gcf()
fig.set_size_inches(12,8)
plt.legend(['orig','matmul'])



###############################################################################
#try with simulated genotypes
import scipy.stats as stats

X_large = np.concatenate((X,np.random.permutation(X),np.random.permutation(X),np.random.permutation(X),
                np.random.permutation(X),np.random.permutation(X),np.random.permutation(X),
                np.random.permutation(X),np.random.permutation(X),np.random.permutation(X)),axis=1)
                          
beta_dist = np.random.laplace(0,1/12,size=[X_large.shape[1],1])+ np.random.normal(0,0.005,size=[X_large.shape[1],1])

beta_init = beta_dist[np.argsort(beta_dist,axis=None)].flatten()
#   beta_init = np.concatenate((beta_init[300:1230],beta_init[330:500],beta_init[0:330]))
#        beta_init = beta_init[::-1]
#        beta_init = np.concatenate((beta_init[0:30],beta_init[400:1230],beta_init[30:400][::-1]))
#            beta_init = beta_init[::-1]
#            beta_init = np.concatenate((beta_init[0:80],beta_init[1000:1230][::1],beta_init[80:1000]))
#beta_init = np.concatenate((np.random.permutation(beta_init[40:1220][::-1]),np.random.permutation(beta_init[0:40][::-1]),beta_init[1220:]))
#                i = 97 #79, 95, 237, 546, 7
#                j = 1163 #1198, 1129, 966, 1026, 1115
#                beta_init = np.concatenate((np.random.permutation(beta_init[i:j][::-1]),np.random.permutation(beta_init[0:i][::-1]),beta_init[j:]))
#i = 546
#j = 1026
#beta_init = np.concatenate((np.random.permutation(beta_init[i:j][::-1]),np.random.permutation(beta_init[0:i][::-1]),beta_init[j:]))
#beta_init = np.concatenate((np.random.permutation(beta_init[1040:1200]),np.random.permutation(beta_init[0:30]),
#                            np.random.permutation(beta_init[30:650]),beta_init[1200:],
#                            np.random.permutation(beta_init[650:1040])))
#                beta_init = np.concatenate((np.random.permutation(np.append(beta_init[0:15],beta_init[1200:])),
#                                            np.random.permutation(beta_init[15:1150]).flatten(),
#                                            np.random.permutation(beta_init[1150:1200]).flatten()))
#    i = 546
#    j = 1026
#    beta_init = np.concatenate((np.random.permutation(beta_init[i:j]).flatten(),
#                                        np.random.permutation(beta_init[j:1200]).flatten(),
#                                        np.random.permutation(np.append(beta_init[0:i],beta_init[1200:]))
#                                        ))  
#                            beta_init = np.random.permutation(beta_init.reshape((82,15))).flatten()  
#                            beta_init = np.asarray([np.random.permutation(x) for x in beta_init.reshape((41,30))]).flatten()
#for i in range(10):
#    beta_init = np.concatenate((beta_init.reshape(82*3,5)[0::2],beta_init.reshape(82*3,5)[1::2])).flatten()
#
#i = 40
#j = 1134
#beta_init = beta_init[::-1]
#                #def chunks(l, n):
#                #    for i in range(0, len(l), n):
#                #        yield l[i:i + n]
#                #beta_init = np.random.permutation(np.asarray(list(chunks(beta_init,50))).flatten()).flatten()
#                #beta_init = np.concatenate(beta_init)
#mid = np.random.permutation(beta_init[50:1180])
#beta_init = np.concatenate((mid[0:1130-i:2],
#                            beta_init[:50],
#                            mid[1130-i:],
#                            beta_init[1180:],
#                            mid[1:1130-i:2]))
beta_from_dosage =np.asarray([np.asarray(i) for i in -(np.mean(X_large,0,keepdims=True))]).T
beta_ind = np.argsort(beta_from_dosage,axis=None)
beta = beta_init[np.argsort(beta_ind,axis=None)]
y = np.matmul(X_large, beta)

ax1 = plt.subplot(2,2,1)
ax1.hist(beta.flatten(),50)
plt.title('beta distribution')
ax2 = plt.subplot(2,2,2)
ax2.hist(y.flatten(),50)
plt.title('phenotype distribution')
ax3 = plt.subplot(2,1,2)
ax3.set_xlabel('index')
ax3.set_ylabel('beta', color='tab:blue')
ax3.plot(range(len(beta_init)),beta[np.argsort(np.argsort(beta_ind),axis=None)],'.', color='tab:blue')
ax3.tick_params(axis='y', labelcolor='tab:blue')
ax4 = ax3.twinx()  # instantiate a second axes that shares the same x-axis
ax4.set_ylabel('dosage', color='tab:red')  # we already handled the x-label with ax1
ax4.plot(range(len(beta_from_dosage)),-beta_from_dosage[np.argsort(beta_from_dosage,axis=None)],'.', color='tab:red')
ax4.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig=plt.gcf()
fig.set_size_inches(8,8)




ii=0
ij=0
i_vec = [0]*1000
j_vec = [0]*1000
beta_all = [beta_init]*100*50*20
beta_sel = [beta_init]*1000
fig, ax = plt.subplots()
for i in range(0,100):
    for j in range(0,50):
        for k in range(20):
            beta_init = beta_dist[np.argsort(beta_dist,axis=None)]
            #            beta_init = np.concatenate((np.random.permutation(np.append(beta_init[0:i],beta_init[1200:])),
            #                                        np.random.permutation(beta_init[i:j]).flatten(),
            #                                        np.random.permutation(beta_init[j:1200]).flatten()))        
            beta_init = beta_init[::-1]
            mid = np.random.permutation(beta_init[j:1230-j])
            beta_init = np.concatenate((mid[0:len(mid)-i:2],
                                        beta_init[:j],
                                        mid[len(mid)-i:],
                                        beta_init[1230-j:],
                                        mid[1:len(mid)-i:2]))        
            beta_from_dosage =np.asarray([np.asarray(i) for i in -(np.mean(X_large,0,keepdims=True))]).T
            beta_ind = np.argsort(beta_from_dosage,axis=None)
            beta = beta_init[np.argsort(beta_ind,axis=None)]
            y = np.matmul(X_large, beta)
            beta_all[ij] = beta_init
            ij += 1
            if (stats.skew(y.flatten())>0.25) and (np.mean(y.flatten())>0):
                beta_sel[ii] = beta_init
                i_vec[ii] = i
                j_vec[ii] = j
                ii += 1
                print(stats.skew(y.flatten()))
                print(np.mean(y.flatten()))
                print('i = '+str(i)+', j = '+str(j))

                ax.hist(y.flatten(),50,alpha=0.1,color='tab:blue')
                plt.title('phenotype distribution')
            
beta_all = np.asarray(beta_all[:ij-1])
beta_sel = np.asarray(beta_sel[:ii-1])
i_vec = np.asarray(i_vec[:ii-1])
j_vec = np.asarray(j_vec[:ii-1])

mean_beta = np.mean(beta_all,axis=0).flatten()
sd_beta = np.std(beta_all,axis=0).flatten()
plt.subplot(211)
plt.plot(mean_beta,'.',color='k',alpha=0.8)
plt.fill_between(range(len(mean_beta)),mean_beta-2*sd_beta, mean_beta+2*sd_beta,
                     color='k',alpha=0.3)
mean_beta = np.mean(beta_sel,axis=0).flatten()
sd_beta = np.std(beta_sel,axis=0).flatten()
plt.subplot(212)
plt.plot(mean_beta,'.',color='tab:blue',alpha=0.8)
plt.fill_between(range(len(mean_beta)),mean_beta-2*sd_beta, mean_beta+2*sd_beta,
                     color='tab:blue',alpha=0.3)
fig=plt.gcf()
fig.set_size_inches(8,6)            

import seaborn as sns
sns.kdeplot(i_vec, j_vec,n_levels=30,shade=True, shade_lowest=False)
plt.xlabel('i index')
plt.ylabel('j index')
            
            
            
max_i = 500
beta_all = [beta_init]*max_i
skew = [0]*max_i
i = 0     
fig, ax = plt.subplots()      
while i < max_i:
    beta_init = beta_dist[np.argsort(beta_dist,axis=None)].flatten()
    beta_init = np.random.permutation(np.asarray([np.random.permutation(x) for x in beta_init.reshape((41,30))])).flatten()
    beta_from_dosage =np.asarray([np.asarray(i) for i in -(np.mean(X_large,0,keepdims=True))]).T
    beta_ind = np.argsort(beta_from_dosage,axis=None)
    beta = beta_init[np.argsort(beta_ind,axis=None)]
    y = np.matmul(X_large, beta)
#    if (stats.skew(y.flatten())>0.15) and (abs(np.mean(y.flatten())-15)<12):
    if (stats.skew(y.flatten())>0.2):

        beta_all[i] = beta_init
#        print(np.array2string(beta,separator=','))
        print(stats.skew(y.flatten()))
        skew[i] = stats.skew(y.flatten())
        print(np.mean(y.flatten()))
#        print(np.array2string(beta_init.flatten(),separator=','))
        print('\n')
#        plt.subplot(5,2,2*i+1)
#        plt.plot(range(len(beta_init)),beta_init,'.')
#        plt.title(i)
#        plt.subplot(5,2,2*(i+1))
#        plt.hist(y.flatten(),50)
#        plt.title('phenotype distribution')
#        fig=plt.gcf()
#        fig.set_size_inches(12,15)
        ax.hist(y.flatten(),50,alpha=0.1,color='tab:blue')
        plt.title('phenotype distribution')
        i += 1
fig=plt.gcf()
fig.set_size_inches(8,6)

beta_array = np.asarray(beta_all[:i])
weight_array = np.asarray(skew[:i])-0.20
mean_beta = np.mean(beta_array,axis=0)
weighted_mean = np.average(beta_array,axis=0, weights=weight_array)
sd_beta = np.std(beta_array,axis=0)

fig, ax = plt.subplots()
ax.set_xlabel('index')
ax.set_ylabel('beta', color='tab:blue')
plt.plot(mean_beta,'.',color='tab:blue',alpha=0.1)
plt.plot(weighted_mean,'.',color='tab:blue',alpha=0.5)
ax.tick_params(axis='y', labelcolor='tab:blue')
ax1 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax1.set_ylabel('dosage', color='tab:red')  # we already handled the x-label with ax1
ax1.plot(range(len(beta_from_dosage)),-beta_from_dosage[np.argsort(beta_from_dosage,axis=None)],'.', color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig=plt.gcf()
fig.set_size_inches(8,4)




chunks = 41
max_i = 500
beta_sel = [beta_init]*max_i
skew = [0]*max_i
idx_vec = [np.arange(chunks)]*max_i
i = 0     
fig, ax = plt.subplots()      
while i < max_i:
    beta_init = beta_dist[np.argsort(beta_dist,axis=None)].flatten()
    rand_idx = np.random.permutation(np.arange(chunks,dtype=int))
#    beta_init = (np.asarray([np.random.permutation(x) for x in beta_init.reshape((chunks,int(1230/chunks)))]))[rand_idx].flatten()
    beta_init = (np.asarray([x[::1-2*np.random.randint(2)] for x in beta_init.reshape((chunks,int(1230/chunks)))]))[rand_idx].flatten()
    beta_from_dosage =np.asarray([np.asarray(i) for i in -(np.mean(X_large,0,keepdims=True))]).T
    beta_ind = np.argsort(beta_from_dosage,axis=None)
    beta = beta_init[np.argsort(beta_ind,axis=None)]
    y = np.matmul(X_large, beta)
    if (stats.skew(y.flatten())>0.20) and (abs(np.mean(y.flatten())-15)<12):
#    if (stats.skew(y.flatten())>0.2):
        beta_sel[i] = beta_init
        idx_vec[i] = rand_idx
        print(stats.skew(y.flatten()))
        skew[i] = stats.skew(y.flatten())
        print(np.mean(y.flatten()))
        print('\n')
        ax.hist(y.flatten(),50,alpha=0.1,color='tab:blue')
        plt.title('phenotype distribution')
        i += 1
fig=plt.gcf()
fig.set_size_inches(8,6)

beta_array = np.asarray(beta_sel[:i])
weight_array = np.asarray(skew[:i])-0.15
mean_beta = np.mean(beta_array,axis=0)
weighted_mean = np.average(beta_array,axis=0, weights=weight_array)
sd_beta = np.std(beta_array,axis=0)

fig, ax = plt.subplots()
ax.set_xlabel('index')
ax.set_ylabel('beta', color='tab:blue')
plt.plot(mean_beta,'.',color='k',alpha=0.1)
plt.plot(weighted_mean,'.',color='tab:blue',alpha=0.5)
ax.tick_params(axis='y', labelcolor='tab:blue')
ax1 = ax.twinx()  # instantiate a second axes that shares the same x-axis
ax1.set_ylabel('dosage', color='tab:red')  # we already handled the x-label with ax1
ax1.plot(range(len(beta_from_dosage)),-beta_from_dosage[np.argsort(beta_from_dosage,axis=None)],'.', color='tab:red')
ax1.tick_params(axis='y', labelcolor='tab:red')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig=plt.gcf()
fig.set_size_inches(8,4)


idx_vec = np.asarray(idx_vec[:i])
mean_idx = np.mean(idx_vec, axis=0)
int_mean_idx= [np.argsort(np.argsort(mean_idx))]
beta_init = (np.asarray([x[::1-2*np.random.randint(2)] for x in beta_init.reshape((chunks,int(1230/chunks)))]))[int_mean_idx].flatten()
std_idx = np.std(idx_vec, axis=0)
plt.plot(mean_idx, '.')
plt.fill_between(range(len(mean_idx)),mean_idx-2*std_idx, mean_idx+2*std_idx,
                     color='tab:blue',alpha=0.3)
fig=plt.gcf()
fig.set_size_inches(8,4)