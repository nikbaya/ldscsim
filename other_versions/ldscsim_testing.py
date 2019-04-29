import hail as hl
import numpy as np
import matplotlib.pyplot as plt


def check_sim_h2(mt):
    return mt.aggregate_cols(hl.agg.stats(mt.__y_no_noise)).stdev**2
    
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


# Check expected vs. observed h2
n_reps = 3
h2_exp = np.asarray([np.arange(0,1.1,0.1)]*n_reps)
h2_obs1 = h2_exp*0
h2_obs2 = h2_exp*0
pi=0.1
for i in range(n_reps):
    for j in range(1,h2_exp.shape[1]):
        sim_mt1 = simulate(mt,mt.genotype, h2=h2_exp[i,j],pi=pi)
        h2_obs1_temp = check_sim_h2(sim_mt1)
        print('Sim h2: {}'.format(h2_obs1_temp))
        h2_obs1[i,j] = h2_obs1_temp
        
fig,ax = plt.subplots(figsize=(6,4))
h2_mean = np.mean(h2_obs1,axis=0)
h2_std = np.std(h2_obs1,axis=0)
ax.plot([0,1],[0,1],'k--',lw=2)
ax.plot(h2_exp[0,:],h2_obs_mean,c=[1,0,0],ls='.-')
ax.fill_between(h2_exp[0,:],h2_mean+2*h2_std,h2_mean-2*h2_std,color=[1,0,0],alpha=0.2)
plt.title('Spike and slab w/ pi = {}\n(reps per h2 value = {})'.format(pi, reps))
plt.xlabel('Expected h2')
plt.ylabel('Observed h2')
fig=plt.gcf()
fig.set_size_inches(6,4)
