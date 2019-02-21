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
