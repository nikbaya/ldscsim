# ldscsim
`ldscsim` is a module written with [Hail](https://hail.is/) to simulate phenotypes. It is designed to test [LD score regression](https://github.com/bulik/ldsc), but is extensible to general use.
<br>
### Outline
* **[Model descriptions](#model-descriptions)**
    - [Phenotype model](#phenotype-model)
    - [Models for SNP effects](#models-for-snp-effects)
        * *[Infinitesimal Model](#infinitesimal-model)*
        * *[Spike and Slab](#spike-and-slab)*
        * *[Annotation-Informed Betas](#annotation-informed-betas)*
    - [Model for Population Stratification](#model-for-population-stratification)
    - [Models for Correlated Traits](#models-for-correlated-traits)
        * *[Multi-Trait Infinitesimal Model](#multi-trait-infinitesimal-model)*
        * *[Two-Trait Spike and Slab](#two-trait-spike-and-slab)*
    - [Models for Ascertainment Bias](#models-for-ascertainment-bias)
* **[Examples](#examples)**
    - [Infinitesimal Model](#infinitesimal-model)
        * *[Simulate a phenotype under the infinitesimal model](#simulate-a-phenotype-under-the-infinitesimal-model)*
        * *[Simulate two phenotypes under the infinitesimal model](#simulate-two-phenotypes-under-the-infinitesimal-model)*
        * *[Simulate three phenotypes under the infinitesimal model](#simulate-three-phenotypes-under-the-infinitesimal-model)*
    - [Spike and Slab Model](#spike-and-slab-model)
        * *[Simulate a single spike and slab phenotype](#simulate-a-single-spike-and-slab-phenotype)*
        * *[Simulate two correlated spike and slab phenotypes](#simulate-two-correlated-spike-and-slab-phenotypes)*
    - [Annotation-Informed Model](#annotation-informed-model)
    - [Population Stratification](#population-stratification)
<br>

## Model descriptions
### Phenotype model
y<sub>i</sub> = &sum;<sub>j</sub> X<sub>ij</sub>&beta;<sub>j</sub> + &epsilon;<sub>i</sub>
* y<sub>i</sub> : Phenotype of individual i
* X<sub>ij</sub> : Genotype of individual i at SNP j
* &beta;<sub>j</sub> : Effect size of SNP j
* &epsilon;<sub>i</sub> : Environmental noise for individual i

The phenotype of each individual is calculating the dot product of that individual's genotypes with the SNP effects and then adding random environmental noise.
<br>
<br>
*Note: Genotypes as described above are normalized such that genotypes at a given SNP have mean=0, var=1.*
<br>
<br>
Relevant functions: 
* `simulate_phenotypes()`
* `calculate_phenotypes()`

### Models for SNP effects
#### Infinitesimal Model
&beta; ~ N(0, h<sup>2</sup>/M)
* h<sup>2</sup> : SNP-based heritability of phenotype
* M : Number of SNPs in simulation

All SNP effects are drawn from a normal distribution with mean=0, variance= h2/M, where h2 is the desired heritability of the simulated trait and M is the number of SNPs in the genotype matrix given by the user. 
<br>
<br>
Relevant functions: 
* `make_betas()`
* `multitrait_inf()`

#### Spike and Slab
&beta; = N(0, h<sup>2</sup>/(&pi;M)), w/ probability &pi;
<br>
&beta; = 0, otherwise
* h<sup>2</sup> : SNP-based heritability of phenotype
* M : Number of SNPs in simulation
* &pi; : Probability of a SNP being causal

SNPs are assigned to be causal with probability &pi;. If a SNP is causal, its effect size is drawn from the normal distribution with mean=0 and variance=h2/(&pi;M). If a SNP is not causal it has an effect size of 0.
<br>
<br>
Relevant functions: 
* `make_betas()`
* `multitrait_ss()`

#### Annotation-Informed Betas
&beta;<sub>j</sub> = N(0, a<sub>j</sub>h<sup>2</sup>/(Var(a<sub>j</sub>)&middot;M)
<br>
a<sub>j</sub> = &sum;<sub>C</sub> &tau;<sub>C</sub>a<sub>Cj</sub>
* h<sup>2</sup> : SNP-based heritability of phenotype
* a<sub>j</sub> : Annotation for SNP j summed across categories
* &tau;<sub>C</sub> : Heritability contributed by annotation category C
* a<sub>jC</sub> : Annotation for SNP j in category C
* M : Number of SNPs in simulation

The effect for SNP j are drawn from a normal distribution with mean=0, variance=a[j]&middot;h2/(Var(a[j])&middot;M), where a[j] is the relative heritability contributed by SNP j and a is a vector of the relative heritability contributed by each SNP. a[j] is calculated by taking the linear combination across all annotation categories for SNP j, scaling each annotation by coefficient &tau;<sub>C</sub>, the heritability contributed by category C.
<br>
<br>
Relevant functions: 
* `make_betas()`
* `multitrait_inf()`
* `agg_fields()`
* `get_coef_dict()`

### Model for Population Stratification
y<sub>i</sub> = &sum;<sub>j</sub> X<sub>ij</sub>&beta;<sub>j</sub> + &varepsilon;<sub>i</sub> + s<sub>i</sub>
<br>
s<sub>i</sub> = &sum;<sub>k</sub> &sigma;<sub>k</sub>s<sub>ik</sub>
* y<sub>i</sub> : Phenotype of individual i
* X<sub>ij</sub> : Genotype of individual i at SNP j
* &beta;<sub>j</sub> : Effect size of SNP j
* &varepsilon;<sub>i</sub> : Environmental noise for individual i
* s<sub>i</sub> : Population stratification term for individual i
* &sigma;<sub>k</sub> : Square root of variance contributed by covariate k
* s<sub>ik</sub> : Covariate k for individual i

After calculating the matrix product of genotypes and SNP effects, it is possible to add population stratification. Population stratification is a term added to the phenotype, which is the linear combination of covariates scaled by coefficients provided by the user.
<br>
<br>
Relevant functions: 
* `calculate_phenotypes()`
* `agg_fields()`
* `get_coef_dict()`


### Models for Correlated Traits
#### Multi-Trait Infinitesimal Model
&Beta; = (&beta;<sup>1</sup>,...,&beta;<sup>k</sup>)<sup>T</sup>
<br>
&Beta; ~ N<sub>k</sub>(0, &Omega;)
* &beta;<sup>k</sup> : SNP effect sizes for trait k
* k : Number of correlated traits
* &Omega; : Variance-covarance matrix, k x k


SNP effects are drawn from a multivariate normal distribution with mean=0 and variance-covariance matrix defined by heritability and genetic correlation values. If the specified heritability and genetic correlation values result in a covariance matrix that is not positive semi-definite the framework will adjust the genetic correlation values to make the covariance matrix positive semi-definite.
<br>
<br>
Relevant functions: 
* `make_betas()`
* `multitrait_inf()`
* `get_cov_matrix()`
* `_nearpsd()`


#### Two-Trait Spike and Slab
(&beta;<sub>j</sub><sup>A<sub>0</sub></sup>, &beta;<sub>j</sub><sup>B<sub>0</sub></sup>) ~ N(0, &Omega;<sub>SS</sub>)
<br>
&Omega;<sub>SS</sub> = [[<sup>1</sup>/<sub>(p<sub>TT</sub>+p<sub>TF</sub>)</sub>, <sup>r<sub>g</sub></sup>/<sub>p<sub>TT</sub></sub>],
<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[<sup>r<sub>g</sub></sup>/<sub>p<sub>TT</sub></sub>, <sup>1</sup>/<sub>(p<sub>TT</sub>+p<sub>FT</sub>)</sub>]]
<br>
&beta;<sub>j</sub><sup>A</sup> = &beta;<sub>j</sub><sup>A<sub>0</sub></sup> and &beta;<sub>j</sub><sup>B</sup> = &beta;<sub>j</sub><sup>B<sub>0</sub></sup>, w/ probability p<sub>TT</sub>
<br>
&beta;<sub>j</sub><sup>A</sup> = &beta;<sub>j</sub><sup>A<sub>0</sub></sup> and &beta;<sub>j</sub><sup>B</sup> = 0, w/ probability p<sub>TF</sub>
<br>
&beta;<sub>j</sub><sup>A</sup> = 0 and &beta;<sub>j</sub><sup>B</sup> = &beta;<sub>j</sub><sup>B<sub>0</sub></sup>, w/ probability p<sub>FT</sub>
* &beta;<sub>j</sub><sup>A<sub>0</sub></sup> : Initial effect size of SNP j for trait A
* &beta;<sub>j</sub><sup>B<sub>0</sub></sup> : Initial effect size of SNP j for trait B
* &Omega;<sub>SS</sub> : Spike & slab variance-covariance matrix
* p<sub>TT</sub> : Probability of a SNP being causal for both traits
* p<sub>TF</sub> : Probability of a SNP being causal for trait A but not trait B
* p<sub>FT</sub> : Probability of a SNP being causal for trait B but not trait A
* r<sub>g</sub> : Genetic correlation between traits A and B
* &beta;<sub>j</sub><sup>A</sup> : Final effect size of SNP j for trait A
* &beta;<sub>j</sub><sup>B</sup> : Final effect size of SNP j for trait B

SNP effects are drawn from bivariate normal distribution and then set to be causal or non-causal based on parameters defining probability of being causal.
<br>
<br>
*Note: Correlated two-trait spike & slab model is from page 30 of the supplement of [Turley et al. 2018](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-017-0009-4/MediaObjects/41588_2017_9_MOESM1_ESM.pdf)*
<br>
<br>
Relevant functions: 
* `make_betas()`
* `multitrait_ss()`

### Models for Ascertainment Bias
[see docstrings for more information]
<br>
<br>
Relevant functions: 
* `binarize()`
* `ascertainment_bias()`


## Examples
`simulate_phenotypes()` is the main method wrapping other methods in the package. However, all methods are self-contained and thus can be run independently. 

Assume for the examples of the infinitesimal and spike & slab models that we have the following MatrixTable `mt`:

```
>>> mt.describe()
----------------------------------------
Global fields:
    None
----------------------------------------
Column fields:
    's': str 
----------------------------------------
Row fields:
    'rsid': str 
----------------------------------------
Entry fields:
    'gt': int32 
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```

* `mt.s` : Sample IDs
* `mt.rsid` : SNP IDs
* `mt.gt` : Genotypes of individuals

### Infinitesimal Model
#### Simulate a phenotype under the infinitesimal model
Simulation specifications:
* Single trait
* Infinitesimal model
* Heritability = 0.1

```
>>> sim = simulate_phenotypes(mt=mt,genotype=mt.gt,h2=0.1)
>>> sim.describe()
----------------------------------------
Global fields:
    'ldscsim': struct {
        h2: float64
    }
----------------------------------------
Column fields:
    's': str
    'y_no_noise': float64
    'y': float64
----------------------------------------
Row fields:
    'rsid': str
    'beta': float64
----------------------------------------
Entry fields:
    'gt': int32
    'norm_gt': float64
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```
* `sim.beta` are the SNP effects for the simulated trait and have the approximate distribution N(0,`h2`/M), where M is the number of SNPs.
* `sim.y_no_noise` is equivalent to the dot product of the genotype matrix by the vector of SNP effects and is modeled to be distributed as N(0,`h2`). 
* `sim.y` is the phenotype after adding random environmental noise. Random environmental noise is scaled to be distributed as N(0,1-`h2`). The `y_no_noise` and the random environmental noise are independent, therefore their sum will have variance equal to the sum of their variances: `h2` + (1-`h2`) = 1. Both phenotype and the random noise have mean equal to zero, therefore the mean of `y` is also 0. Thus, `y` is modeled to have distribution N(0,1), just like any real phenotype that has been standardized. 
* `sim.norm_gt` is the normalized genotype used to calculate the phenotype. Genotypes are normalized such that at a given SNP, the distribution of genotypes across individuals is approximately N(0,1). 
* `sim.ldscsim.h2` is the `h2` parameter passed to the simulation.

To check heritability of simulated trait:
```
>>> hl.eval(sim.ldscsim.h2) # the expected h2 passed as a parameter
0.1 
>>> sim.aggregate_cols(hl.agg.stats(sim.y_no_noise)).stdev**2 # calculating observed h2
0.09413353719961155 
```

#### Simulate two phenotypes under the infinitesimal model
Simulation specifications:
* Multi-trait
* Infinitesimal model
* Trait 1 heritability = 0.3
* Trait 2 heritability = 0.4
* Genetic correlation = 0.6

```
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.3, 0.4], rg=0.6)
>>> sim.describe()
----------------------------------------
Global fields:
    'ldscsim': struct {
        h2: array<float64>, 
        rg: float64
    }
----------------------------------------
Column fields:
    's': str
    'y_no_noise': array<float64>
    'y': array<float64>
----------------------------------------
Row fields:
    'rsid': str
    'beta': array<float64>
----------------------------------------
Entry fields:
    'gt': int32
    'norm_gt': float64
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```
Compared to the previous simulation, the main changes are that `beta`,`y_no_noise`, and `y` are all arrays because they hold values for multiple traits.
```
>>> sim.cols().show()
+-----------+-----------------------+-----------------------+
| s         | y_no_noise            | y                     |
+-----------+-----------------------+-----------------------+
| str       | array<float64>        | array<float64>        |
+-----------+-----------------------+-----------------------+
| "1002230" | [3.74e-02,9.02e-01]   | [1.99e+00,1.44e-01]   |
| "1006422" | [-1.52e+00,-1.07e+00] | [-1.06e+00,-9.58e-01] |
| "1011645" | [-4.12e-01,-7.80e-01] | [-2.21e-01,-2.42e-01] |
| "1014723" | [7.33e-01,6.17e-01]   | [1.53e+00,-7.43e-01]  |
| "1021993" | [2.58e-01,-2.40e-01]  | [-6.44e-02,5.31e-01]  |
| "1022263" | [7.90e-01,8.05e-01]   | [2.34e+00,3.04e-01]   |
| "1029991" | [8.99e-02,1.87e-01]   | [-7.68e-01,5.34e-01]  |
| "1034880" | [-8.72e-02,-1.08e-01] | [-5.53e-01,-1.39e+00] |
| "1040108" | [1.03e-02,-8.61e-01]  | [2.39e-01,-7.03e-01]  |
| "1046474" | [-3.16e-01,1.34e-01]  | [-1.54e+00,9.61e-01]  |
+-----------+-----------------------+-----------------------+
```
Each field of arrays is indexed by trait. For instance, `mt.y_no_noise[0]` is the `y_no_noise` field for the first trait and `mt.y[0]` is the `y` field corresponding to the first trait, and `mt.y_no_noise[1]` is the `y_no_noise` field for the second trait and `mt.y[1]` is the `y` field corresponding to the second trait. The same rule applied for `mt.beta[0]` and `mt.beta[1]`.

To check the heritabilities and genetic correlations between traits:
```
>>> hl.eval(sim.ldscsim.h2) # the expected h2 values passed as parameters
[0.3, 0.4]
>>> [sim.aggregate_cols(hl.agg.stats(sim.y_no_noise[x])).stdev**2 for x in range(2)] # calculating observed h2
[0.30562739116200704, 0.43268586014716076]
```

#### Simulate three phenotypes under the infinitesimal model
Simulation specifications:
* Multi-trait
* Infinitesimal model
* Trait 1 heritability = 0.1
* Trait 2 heritability = 0.2
* Trait 3 heritability = 0.7
* Genetic correlation trait 1 & trait 2 = 0.8
* Genetic correlation trait 1 & trait 3 = 0.5
* Genetic correlation trait 2 & trait 3 = 0.4

```
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.1, 0.2, 0.7], rg=[0.8, 0.5, 0.4])
>>> sim.describe()
----------------------------------------
Global fields:
    'ldscsim': struct {
        h2: array<float64>, 
        rg: array<float64>
    }
----------------------------------------
Column fields:
    's': str
    'y_no_noise': array<float64>
    'y': array<float64>
----------------------------------------
Row fields:
    'rsid': str
    'beta': array<float64>
----------------------------------------
Entry fields:
    'gt': int32
    'norm_gt': float64
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```
This produces a similar MatrixTable to the previous simulation. However, all array fields have three elements rather than two.
```
>>> sim.cols().show()
+-----------+---------------------------------+---------------------------------+
| s         | y_no_noise                      | y                               |
+-----------+---------------------------------+---------------------------------+
| str       | array<float64>                  | array<float64>                  |
+-----------+---------------------------------+---------------------------------+
| "1002230" | [-3.47e-01,-4.53e-01,1.09e-01]  | [-1.50e+00,2.08e-01,7.32e-01]   |
| "1006422" | [-3.44e-01,-4.23e-01,3.10e-01]  | [-8.73e-01,7.69e-01,5.73e-01]   |
| "1011645" | [-2.09e-01,-4.78e-01,-5.55e-01] | [-6.01e-01,-5.94e-01,-6.69e-01] |
| "1014723" | [5.52e-01,1.72e-01,-3.88e-01]   | [1.44e+00,3.07e-01,-1.52e-01]   |
| "1021993" | [2.89e-02,-1.43e-01,-1.12e-01]  | [-4.30e-01,-1.93e+00,-6.49e-03] |
| "1022263" | [6.97e-03,2.57e-01,2.13e-01]    | [-1.25e+00,1.12e+00,1.99e-01]   |
| "1029991" | [-2.84e-01,-3.34e-01,5.57e-01]  | [6.09e-01,-1.91e+00,2.34e-01]   |
| "1034880" | [-1.87e-01,3.98e-02,-5.90e-01]  | [-1.28e+00,-2.08e+00,-3.99e-01] |
| "1040108" | [4.09e-01,3.11e-01,-5.32e-01]   | [1.70e+00,8.05e-01,-1.07e+00]   |
| "1046474" | [1.39e-01,-1.63e-01,2.44e-01]   | [1.26e-01,2.37e-01,1.22e-01]    |
+-----------+---------------------------------+---------------------------------+
```

### Spike and Slab Model
#### Simulate a single spike and slab phenotype
Simulation specifications:
* Single trait
* Spike & slab
* Heritability = 0.1 
* Probability of a SNP being causal (`pi`) = 0.01

```
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=0.1, pi=0.01)
>>> sim.describe()
----------------------------------------
Global fields:
    'ldscsim': struct {
        h2: float64, 
        pi: float64
    }
----------------------------------------
Column fields:
    's': str
    'y_no_noise': float64
    'y': float64
----------------------------------------
Row fields:
    'rsid': str
    'beta': float64
----------------------------------------
Entry fields:
    'gt': int32
    'norm_gt': float64
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```

#### Simulate two correlated spike and slab phenotypes
Simulation specifications:
* Two traits
* Spike & slab
* Trait 1 heritability = 0.8
* Trait 2 heritability = 0.9
* Genetic correlation = 0.5
* Probability a SNP is causal for both traits = 0.3
* Probability SNP is causal for trait 1 but not trait 2 = 0.1
* Probability SNP is causal for trait 2 but not trait 1 = 0.2 

```
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.8, 0.9], pi=[0.3, 0.1, 0.2], rg =0.5)
```

### Annotation-Informed Model
Assume for this example we begin with the original MatrixTable shown above but we add new fields that represent annotations:

```
>>> mt1 = mt.annotate_rows(a1 = hl.rand_bool(p=0.01), 
                          a2 = hl.rand_bool(p=0.05),
                          a3 = hl.rand_norm())
>>> mt1.describe()
----------------------------------------
Global fields:
    None
----------------------------------------
Column fields:
    's': str
----------------------------------------
Row fields:
    'rsid': str
    'a1': bool
    'a2': bool
    'a3': float64
----------------------------------------
Entry fields:
    'gt': int32
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```
* `a1`,`a2`,`a3` : Annotations we want to use for an annotation-informed model

We can aggregate those annotations as a linear combination using default coefficients to scale each field before summing across fields.

```
>>> mt2 = agg_fields(tb=mt1, str_expr='a', axis='rows') # aggregate across row fields using default coefficients (all equal to 1)
Assuming coef = 1 for all annotations
Fields and associated coefficients used in annot aggregation: {'a1': 1, 'a2': 1, 'a3': 1}
```
We can also use a dictionary to specify coefficients to use for each row field when aggregating as a linear combination.
```
>>> coef_dict = {'a1':0.3, 'a2':0.1, 'a3':0.4} # coefficients to multiply each field by before summing
>>> mt2 = agg_fields(tb=mt1, coef_dict=coef_dict, axis='rows') # aggregate across row fields, specifying coefficients
Fields and associated coefficients used in annot aggregation: {'a1': 0.3, 'a2': 0.1, 'a3': 0.4}
```
In either case, the output MatrixTable will be:
```
>>> mt2.describe()
----------------------------------------
Global fields:
    None
----------------------------------------
Column fields:
    's': str
----------------------------------------
Row fields:
    'rsid': str
    'a1': bool
    'a2': bool
    'a3': float64
    'agg_annot': float64
----------------------------------------
Entry fields:
    'gt': int32
----------------------------------------
Column key: ['s']
Row key: ['rsid']
----------------------------------------
```
Note that when using `agg_fields()`, specifying `str_expr` but not `coef_dict` can unintentionally include fields in the aggregation with names that match `str_expr`. For instance, if we try using `mt2` as our input MatrixTable:
```
>>> mt3 = agg_fields(tb=mt2, str_expr='a', axis='rows')
Assuming coef = 1 for all annotations
Fields and associated coefficients used in annot aggregation: {'a1': 1, 'a2': 1, 'a3': 1, 'agg_annot': 1}
```
Because our `str_expr` matches the field `agg_annot` it also includes it in the fields to be aggregated. To avoid this, either rename the desired target fields to have a more unique pattern to search for, or drop `agg_annot`.

*See docstring of* `make_betas()` *for more information*

### Population Stratification
[*In development*]
* `mt.PC1`,`mt.PC2`,`mt.PC3` : Covariates we want to use for population stratification
