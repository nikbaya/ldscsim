###############################################################################
# Version 2

import hail as hl

#class BetaDistribution(a, b):
#    pass

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
    
###############################################################################################
# Version 3
class SimulateBetas():
	def __init__(self, h2, mt):
		self.h2 = h2
        self.mt = mt
        self.n_rows = self.mt.count_rows()
        if 'beta' in [rf for rf in self.mt.row]:
            self.mt = self.mt.drop('beta')
            
    def infinitesimal(self):
        self.mt = self.mt.annotate_rows(beta = hl.rand_norm(0,hl.sqrt(self.h2/self.n_rows)))

    def infinitesimal_spike_slab(self,pi):
        self.mt = self.mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(self.n_rows*pi))))
        
    def annotation_informed(self,a_ht):
        var = a_ht.aggregate(hl.agg.stats(a_ht.a)).stdev ** 2
        self.mt = self.mt.annotate_rows(anno = a_ht[mt.rsid]['a']/var*self.h2/self.n_rows) #requires both mt and a_ht to have rsid row fields and a_ht to have field 'a' with annotations reduced to a single number. assumes mt does not contain row field 'anno'
        self.mt = self.mt.annotate_rows(beta = mt.anno/var)

    def simulate(self,with_noise):
        self.mt = self.mt.annotate_rows(y_no_noise = hl.agg.sum(self.beta*self.dosage)) #assumes that mt has 'dosage' entry field
        self.mt = self.mt.annotate_rows(y = self.beta*self.dosage + with_noise*hl.rand_norm(0,1-self.h2))
        return self.mt




################################################################################################################
# Version 4
class SimulateBetas():
	def __init__(self, h2, mt):
		self.h2 = h2
        self.mt = mt
        self.n_rows = self.mt.count_rows()
        if 'beta' in [rf for rf in self.mt.row]:
            self.mt = self.mt.drop('beta')
            
    def infinitesimal(self):
        return self.mt.annotate_rows(beta = hl.rand_norm(0,hl.sqrt(self.h2/self.n_rows)))

    def infinitesimal_spike_slab(self,pi):
        return self.mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(self.n_rows*pi))))
        
    def annotation_informed(self,a_ht):
        var = a_ht.aggregate(hl.agg.stats(a_ht.a)).stdev ** 2
        return self.mt.annotate_rows(beta = a_ht[mt.rsid]['a']/var*self.h2/self.n_rows) #requires both mt and a_ht to have rsid row fields and a_ht to have field 'a' with annotations reduced to a single number. assumes mt does not contain row field 'anno'

def simulate(mt, h2, model,with_noise,pi = 1,a_ht = 0):
    simbetas = SimulateBetas(h2, mt)
    if model is 'infinitesimal':
        mt = simbetas.infinitesimal()
    elif model is 'initesimal_spike_slab':
        mt = simbetas.infinitesimal_spike_slab(pi)
    elif model is 'annotation_informed':
        mt = simbetas.annotation_informed(a_ht)
        
    self.mt = self.mt.annotate_cols(y_no_noise = hl.agg.sum(self.beta*self.dosage)) #assumes that mt has 'dosage' entry field
    self.mt = self.mt.annotate_cols(y = self.beta*self.dosage + with_noise*hl.rand_norm(0,1-self.h2))
    return self.mt
        
        
###############################################################################################################
# Version 5
def simulate(h2, mt, pi=1, annot=0):
    
    if annot is 0:
        self.mt = self.mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(self.n_rows*pi))))
    else:
        var = annot.aggregate(hl.agg.stats(annot.a)).stdev ** 2
        self.mt =  self.mt.annotate_rows(beta = annot[mt.rsid]['a']/var*self.h2/self.n_rows) #requires both mt and annot to have rsid row fields and annot to have field 'a' with annotations reduced to a single number. assumes mt does not contain row field 'anno'
    self.mt = self.mt.annotate_cols(y_no_noise = hl.agg.sum(self.beta*self.dosage)) #assumes that mt has 'dosage' entry field
    self.mt = self.mt.annotate_cols(y = self.beta*self.dosage + hl.rand_norm(0,1-self.h2))
    return self.mt

###############################################################################################################
# Version 6
class SimDistribution():
    def __init__(self):
        pass
    
    def make_random_function(self):
        pass   
        
    def add_noise(self, mt, h2):
        return mt.annotate_cols(y = mt.y_no_noise + hl.rand_norm(0,1-self.h2))
        
    def simulate(self, mt, h2, pi, annot):
        self.make_random_function(mt, h2, pi, annot)
        
class InfinitesimalModel(SimDisribution):
    def make_random_function(self, mt, h2, pi=1, annot=0):
        M = mt.count_rows()
        return mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(M*pi)))

class AnnotationInformedModel(SimDistribution):
    def make_random_function(self, mt, h2, pi=0, annot)
        var = annot.aggregate(hl.agg.stats(annot.a)).stdev ** 2
        return self.mt.annotate_rows(beta = annot[mt.rsid]['a']/var*self.h2/self.n_rows) #requires both mt and annot to have rsid row fields and annot to have field 'a' with annotations reduced to a single number. assumes mt does not contain row field 'anno'
        
###############################################################################################################
# Version 7.0
class SimDistribution():
  def __init__(self):
      pass

  def make_betas(self, mt, h2, pi=1, annot=0):
      M = mt.count_rows()
      if annot is not 0:
          var = annot.aggregate(hl.agg.stats(annot.a)).stdev ** 2
          return self.mt.annotate_rows(beta = annot[mt.rsid]['a'] / var * h2 / M)
      else:
          return mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(M*pi))))
          
  def add_noise(self, mt, h2):
      return mt.annotate_cols(y = mt.y_no_noise + hl.rand_norm(0,1-self.h2))

  def simulate(self, mt, h2, pi, annot):
      mt = self.make_betas(mt, h2, pi, annot)
      mt = mt.annotate_cols(y_no_noise = hl.agg.sum(mt.beta * mt.dosage))
      return self.add_noise(mt, h2)
      
###############################################################################################################
# Version 7.1
      
def make_betas(self, mt, h2, annot=None, pi=1):
      M = mt.count_rows()
      if annot is not None:
          var = annot.aggregate(hl.agg.stats(annot).stdev ** 2, _localize=False)
          return self.mt.annotate_rows(beta = annot / var * h2 / M)
      else :
          return mt.annotate_rows(beta = hl.rand_bool(pi)*hl.rand_norm(0,hl.sqrt(self.h2/(M*pi))))
      
def add_noise(self, mt, h2):
    return mt.annotate_cols(y = mt.y_no_noise + hl.rand_norm(0,1-self.h2))
  
def simulate(self, mt, h2, pi, annot):
    mt = self.make_betas(mt, h2, pi, annot)
    mt = mt.annotate_cols(y_no_noise = hl.agg.sum(mt.beta * mt.dosage))
    return self.add_noise(mt, h2)
    
    
