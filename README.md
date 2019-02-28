# ldscsim
`ldscsim` is a module to simulate phenotypes. It was originally designed to generate phenotypes for testing [`ldsc`](https://github.com/bulik/ldsc), but is extensible to general use.
<br>

## Model descriptions
#### Models for SNP effects or "betas"
* **Infinitesimal**: All SNP effects are drawn from a normal distribution with mean=0, variance=`h2/M`, where `h2` is the desired heritability of the simulated trait and `M` is the number of SNPs in the genotype matrix given given by the user.
* **Spike & Slab**: SNPs have probability `pi` of being causal. If causal, the SNP effect is drawn from a normal distribution with variance `h2/(M*pi)`. If not causal, the SNP effect is zero.
* **Annotation-Informed**: The effect for SNP `j` are drawn from a normal distribution with mean=0, variance=`a[j]*h2`, where `a[j]` is the relative heritability contributed by SNP `j` and `a` is a vector of the relative heritability contributed by each SNP. `a[j]` is calculated by taking the linear combination across all annotations of SNP `j`, scaling each annotation by coefficients specified by the user, or assumed to be 1.

#### Modeling population stratification
After calculating the matrix product of genotypes and SNP effects, it is possible to add population stratification. Population stratification is a term added to the phenotype, which is the linear combination of covariates scaled by given coefficients.
<br>

## Getting started
`simulate()` is the method which wraps other methods in the package. However, all methods are self-contained. For instance, if you just want to calculate simulated phenotypes using betas (SNP effects) generated outside of the `ldscsim` framework, you can use `sim_phenotypes()` independent of `simulate()`.

Assume for all examples that we have the following MatrixTable `mt`:

```python
>>> mt.describe()
```

```
––––––––––––––––––––––––––––––––––––––––
Global fields:
    None
––––––––––––––––––––––––––––––––––––––––
Column fields:
    's': str 
    'PC1': float64 
    'PC2': float64 
    'PC3': float64 
––––––––––––––––––––––––––––––––––––––––
Row fields:
    'rsid': str 
    'a1': bool 
    'a2': bool 
    'a3': bool 
––––––––––––––––––––––––––––––––––––––––
Entry fields:
    'gt': int32 
––––––––––––––––––––––––––––––––––––––––
Column key: ['s']
Row key: ['rsid']
––––––––––––––––––––––––––––––––––––––––
```

* `mt.s` : Sample IDs
* `mt.PC1`,`mt.PC2`,`mt.PC3` : Covariates we want to use for population stratification
* `mt.rsid` : SNP IDs
* `mt.a1`,`mt.a2`,`mt.a3` : Annotations we want to use for an annotation-informed model
* `mt.gt` : Genotypes of individuals

### Infinitesimal Model
Simulate a phenotype under the infinitesimal model with heritability = 0.1

```python
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,h2=0.1)
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 1 (default: 1)
Annotation-informed betas?: NO
h2-normalized betas?: YES
Add population stratification?: NO
****************************************
Simulating betas using infinitesimal model w/ h2 = 0.1
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
Finished simulation! (runtime=0.1041 min)
```

### Spike & Slab Model
Simulate a phenotype with heritability = 0.1 and probability of a SNP being causal = 0.01

```python
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,h2=0.1,pi=0.01)
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 0.01 (default: 1)
Annotation-informed betas?: NO
h2-normalized betas?: YES
Add population stratification?: NO
****************************************
Simulating betas using spike & slab model w/ h2 = 0.1
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
Finished simulation! (runtime=0.1118 min)  
```

### Annotation-Informed Model
Simulate a phenotype with heritability = 0.1 using row fields that have names that match the regex expression `'a[0-9]'` as our annotations, without scaling the annotations:

```python
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,h2=0.1,is_annot_inf=True,annot_regex='a[0-9]')
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 1 (default: 1)
Annotation-informed betas?: YES
h2-normalized betas?: YES
Add population stratification?: NO
****************************************
Simulating h2-normalized annotation-informed betas (default coef: 1)
Assuming coef = 1 for all annotations
Fields and associated coefficients used in annot aggregation: {'a1': 1, 'a2': 1, 'a3': 1}
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
Finished simulation! (runtime=0.4891 min)
```

Simulate a phenotype with heritability = 0.1 using row fields that match keys in the `annot_coef_dict` as our annotations, scaling annotations by coefficients in `annot_coef_dict`:

```python
>>> annot_coef_dict={'a1':0.2,'a2':0.4,'a3':0.1}
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,h2=0.1,is_annot_inf=True,annot_coef_dict=annot_coef_dict)
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 1 (default: 1)
Annotation-informed betas?: YES
h2-normalized betas?: YES
Add population stratification?: NO
****************************************
Simulating h2-normalized annotation-informed betas using annot_coef_dict
Fields and associated coefficients used in annot aggregation: {'a1': 0.2, 'a2': 0.4, 'a3': 0.1}
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
Finished simulation! (runtime=0.5021 min)
```

Simulate a phenotype with heritability = 0.1 using row fields that match keys in the `annot_coef_dict` as our annotations, scaling annotations by coefficients in `annot_coef_dict` , without scaling the variance of the betas to match a desired heritability:

```python
>>> annot_coef_dict={'a1':0.0001,'a2':0.0004,'a3':0.0002}
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,is_annot_inf=True,annot_coef_dict=annot_coef_dict,h2_normalize=False)
```

```
****************************************
Running simulation framework
pi = 1 (default: 1)
Annotation-informed betas?: YES
h2-normalized betas?: NO
Add population stratification?: NO
****************************************
Simulating  annotation-informed betas using annot_coef_dict
Fields and associated coefficients used in annot aggregation: {'a1': 0.0001, 'a2': 0.0004, 'a3': 0.0002}
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
Finished simulation! (runtime=1.444 min)     
```

It is possible to have SNP-based h2 > 1 if the coefficients in `annot_coef_dict` are too large. In this case, no environmental noise is added:

```
****************************************
Running simulation framework
pi = 1 (default: 1)
Annotation-informed betas?: YES
h2-normalized betas?: NO
Add population stratification?: NO
****************************************
Simulating  annotation-informed betas using annot_coef_dict
Fields and associated coefficients used in annot aggregation: {'a1': 0.001, 'a2': 0.004, 'a3': 0.0002}
Normalizing genotypes...                                                        
Calculating phenotypes...                                                       
WARNING: Total SNP-based h2 = 2.696333128695569 (>1)
Not adding environmental noise
Finished simulation! (runtime=1.4135 min)      
```

### Population Stratification
Simulate a phenotype with heritability = 0.1 (under the infinitesimal model) using covariates that match the regex expression `PC[0-9]` as our covariates for the population stratification, without scaling the covariates.

```python
>>> sim_mt = simulate(mt=mt,genotype=mt.gt,h2=0.1,is_popstrat=True,cov_regex='PC[0-9]')
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 1 (default: 1)
Annotation-informed betas?: NO
h2-normalized betas?: YES
Add population stratification?: YES
****************************************
Simulating betas using infinitesimal model w/ h2 = 0.1
Normalizing genotypes...                                                        
Calculating phenotypes w/ population stratification...                          
Adding population stratification...                                             
Assuming coef = 1 for all covariates
Fields and associated coefficients used in cov aggregation: {'PC1': 1, 'PC2': 1, 'PC3': 1}
Finished simulation! (runtime=0.1212 min) 
```

Simulate a phenotype with heritability = 0.1 using covariates that match keys in `cov_coef_dict` as our covariates for the population stratification, scaling covariates by coefficients in `cov_coef_dict`.

```python
>>> cov_coef_dict={'b1':0.4,'b2':0.1,'b3':0.7}
>>> simulate(mt,mt.gt,h2=0.1,is_popstrat=True,cov_coef_dict=cov_coef_dict)
```

```
****************************************
Running simulation framework
h2 = 0.1
pi = 1 (default: 1)
Annotation-informed betas?: NO
h2-normalized betas?: YES
Add population stratification?: YES
****************************************
Simulating betas using infinitesimal model w/ h2 = 0.1
Normalizing genotypes...                                                        
Calculating phenotypes w/ population stratification...                          
Adding population stratification...                                             
Fields and associated coefficients used in cov aggregation: {'b1': 0.4, 'b2': 0.1, 'b3': 0.7}
Finished simulation! (runtime=0.136 min)
```

## Methods
###### `simulate(mt, genotype, h2=None, pi=1, is_annot_inf=False, annot_coef_dict=None, annot_regex=None, h2_normalize=True,  is_popstrat=False, cov_coef_dict=None, cov_regex=None, path_to_save=None)`
Simulates phenotypes under various models.
<br>

*Parameters*: 
* **`mt`** :  Hail MatrixTable holding all the matrix fields passed as parameters to the model
* **`genotype`** : `mt` entry field of genotypes. Accepted data types: `expr_int32`,`expr_int64`,`expr_float32`,`expr_float64`
* **`h2`** (optional): Desired heritability of simulated trait. Accepted data types: `float`, `int`
* **`pi`** (optional) : Probability of SNP being causal. Used for spike & slab model. Accepted data types: `float`,`int`
* **`is_annot_inf`** (optional): Whether the desired simulation uses annotation-informed betas. Accepted data types: `bool`
* **`annot_coef_dict`** (optional): Dictionary for annotation coefficients. (keys=name of annotation row field, values=coefficients). Accepted data types: `dict`
* **`annot_regex`** (optional): Regular expression for searching for row fields of desired annotations. Accepted data types: `str`
* **`h2_normalize`** (optional): Whether to set the heritability of the SNP effects to h2 when using the annotation-informed model. Accepted data types: `bool`
* **`is_popstrat`** (optional): Whether the desired simulation adds population stratification. Accepted data types: `bool`
* **`cov_coef_dict`** (optional): Dictionary for covariate coefficients. (keys=name of covariate column field, values=coefficients). Accepted data types: `dict`
* **`cov_regex`** (optional): Regular expression for searching for column fields of desired covariates. Accepted data types: `str`
* **`path_to_save`** (optional) : Path to save MatrixTable containing simulation output. Accepted data type: `str`
<br>

###### `check_beta_args`
Checks beta args for simulate() and make_betas()
<br>

###### `check_popstrat_args`
Checks popstrat args for simulate() and add_popstrat()
<br>

###### `print_header(h2, pi, is_annot_inf, h2_normalize, is_popstrat, path_to_save):`
Prints header describing simulation parameters
<br>

###### `annotate_w_temp_fields(mt, genotype, h2, pi=1, annot=None, popstrat=None, popstrat_s2=1)`
Annotate mt with temporary fields of simulation parameters
<br>

###### `make_betas(mt, h2, pi=1, annot=None)`
Simulate betas. Options: Infinitesimal model, spike & slab, annotation-informed
<br>

###### `normalize genotypes(mt, genotype)`
Normalizes genotypes
<br>

###### `add_pop_strat(mt, y, popstrat, popstrat_s2=1)`
Adds population stratification `popstrat` scaled by `popstrat_s2` to phenotype `y`
<br>

###### `sim_phenotypes(mt, genotype, h2, beta, popstrat=None, popstrat_s2=1)`
Simulate phenotypes given SNP effects (`beta`) and genotypes. Adding population stratification is optional.
<br>

###### `clean_fields(mt, str_expr)`
Removes fields with names that have str_expr in them
<br>

###### `add_sim_description(mt,h2,starttime,stoptime,runtime,pi=1,annot=None,popstrat=None,popstrat_s2=1,path_to_save=None)`
Annotates `mt` with description of simulation

