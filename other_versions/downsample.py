def downsample(mt, frac, phen, for_cases=None, seed = None):
    start = dt.datetime.now()
    assert frac <= 1 and frac >= 0, "frac must be in [0,1]"
    phen_name = phen._ir.name
    n = mt.count_cols()
    n_cas = mt.filter_cols(mt[phen_name]==1).count_cols()
    if frac == 1:
        return mt, n, n_cas
    else:
        seed = seed if seed is not None else int(str(Env.next_seed())[:8])
        header = '\n############\n'
        header += 'Downsampling '+('all' if for_cases is None else ('cases'*for_cases+'controls'*(for_cases==0)))+f' by frac = {frac}\n'
        header += f'n: {n}\n'
        header += f'n_cas: {n_cas}\nn_con: {n-n_cas}\nprevalence: {round(n_cas/n,6)}\n' if for_cases != None else ''
        header += f'seed: {seed}\n'
        header += 'Time: {:%H:%M:%S (%Y-%b-%d)}\n'.format(dt.datetime.now())
        header += '############'
        print(header)        
        randstate = np.random.RandomState(int(seed)) #seed random state for replicability
        col_key  = mt.col_key
        mt = mt.add_col_index('tmp_col_idx')
        mt = mt.key_cols_by('tmp_col_idx')
        for_cases = bool(for_cases) if for_cases != None else None
        filter_arg = (mt[phen_name] == (for_cases==0)) if for_cases != None else (~hl.is_defined(mt[phen_name]))
#        mtA = mt.filter_cols(filter_arg) #keep all individuals in mtA
        mtB = mt.filter_cols(filter_arg , keep=False) #downsample individuals in this mt
        mtB = mtB.add_col_index('col_idx_tmpB')
        mtB = mtB.key_cols_by('col_idx_tmpB')
        nB = n_cas*for_cases + (n-n_cas)*(for_cases==0) if for_cases is not None else n
        n_keep = int(nB*frac)
        labels = ['A']*(n_keep)+['B']*(nB-n_keep)
        randstate.shuffle(labels)
        mtB = mtB.annotate_cols(label = hl.literal(labels)[hl.int32(mtB.col_idx_tmpB)])
        mtB = mtB.filter_cols(mtB.label == 'B') #filter to samples we wish to discard from original mt
        mtB = mtB.key_cols_by('tmp_col_idx')
        mt1 = mt.anti_join_cols(mtB.cols())
        mt1 = mt1.key_cols_by(*col_key)
        n_new = mt1.count_cols()
        n_cas_new = mt1.filter_cols(mt1[phen_name]==1).count_cols()
        elapsed = dt.datetime.now()-start
        print('\n############')
        print('Finished downsampling '+('all' if for_cases is None else ('cases'*for_cases+'controls'*(for_cases==0)))+f' by frac = {frac}')
        print(f'n: {n} -> {n_new} ({round(100*n_new/n,3)}% of original)')
        if n_cas != 0 and n_new != 0 :
            print(f'n_cas: {n_cas} -> {n_cas_new} ({round(100*n_cas_new/n_cas,3)}% of original)')
            print(f'n_con: {n-n_cas} -> {n_new-n_cas_new} ({round(100*(n_new-n_cas_new)/(n-n_cas),3)}% of original)')
            print(f'prevalence: {round(n_cas/n,6)} -> {round(n_cas_new/n_new,6)} ({round(100*(n_cas_new/n_new)/(n_cas/n),6)}% of original)')
        print(f'Time for downsampling: '+str(round(elapsed.seconds/60, 2))+' minutes')
        print('############')
        return mt1, n_new, n_cas_new
