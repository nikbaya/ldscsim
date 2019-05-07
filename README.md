# ldscsim
`ldscsim` is a module to simulate phenotypes. It is designed to test [`ldsc`](https://github.com/bulik/ldsc), but is extensible to general use.
<br>

## Model descriptions
#### Models for single trait SNP effects
* **Infinitesimal**: All SNP effects are drawn from a normal distribution with mean=0, variance=`h2/M`, where `h2` is the desired heritability of the simulated trait and `M` is the number of SNPs in the genotype matrix given given by the user. 
* **Spike & Slab**: SNPs have probability `pi` of being causal. If causal, the SNP effect is drawn from a normal distribution with variance `h2/(M*pi)`. If not causal, the SNP effect is zero.
* **Annotation-Informed**: The effect for SNP `j` are drawn from a normal distribution with mean=0, variance=`a[j]`, where `a[j]` is the relative heritability contributed by SNP `j` and `a` is a vector of the relative heritability contributed by each SNP. `a[j]` is calculated by taking the linear combination across all annotations of SNP `j`, scaling each annotation by coefficients (often written as tau) specified by the user, or assumed to be 1.

#### Models for correlated multi-trait/two-trait SNP effects
* **Infinitesimal** (multi-trait): SNP effects are drawn from a multivariate normal distribution with mean=0 and variance-covariance matrix defined by heritability and genetic correlation values.
* **Spike & Slab** (two-trait): SNP effects are drawn from bivariate normal distribution and then set to be causal or non-causal based on parameters defining probability of being causal.


#### Modeling population stratification
After calculating the matrix product of genotypes and SNP effects, it is possible to add population stratification. Population stratification is a term added to the phenotype, which is the linear combination of covariates scaled by coefficients provided by the user (if not provided, coefficients are assumed to be 1).
<br>

## Getting started
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
Simulate a phenotype under the infinitesimal model with heritability = 0.1

```python
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
In `sim`, `y` is the simulated phenotype with approximate distribution N(mean=0,var=1). `y_no_noise` is the simulated phenotype without noise added and is approximately distributed as Normal(mean=0,var=`h2`). `beta` are the SNP effects for the simulated trait and have the approximate distribution N(mean=0,var=`h2`/M), where M is the number of SNPs. `norm_gt` is the normalized genotype used to calculate the phenotype. Genotypes are normalized such that at a given SNP, the distribution of genotypes across individuals is approximately N(mean=0,var=1). In the global variables, `mt.ldscsim.h2` is the `h2` parameter passed to the simulation.

To check heritability of simulated trait:
```
>>> hl.eval(sim.ldscsim.h2) # the expected h2 passed as a parameter
0.1 
>>> sim.aggregate_cols(hl.agg.stats(sim.y_no_noise)).stdev**2 # calculating observed h2
0.09413353719961155 
```

Simulate two phenotypes under the infinitesimal model with heritabilities 0.3 and 0.4 and genetic correlation of 0.6

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
Compared to the previous simulation, the main changes are that `beta`,`y_no_noise`, and `y` are all arrays because they hold values for all traits.
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
```python
>>> hl.eval(sim.ldscsim.h2) # the expected h2 values passed as parameters
[0.3, 0.4]
>>> [sim.aggregate_cols(hl.agg.stats(sim.y_no_noise[x])).stdev**2 for x in range(2)] # calculating observed h2
[0.30562739116200704, 0.43268586014716076]
```

Simulate three phenotypes under the infinitesimal model with heritabilities 0.1, 0.2, 0.7 and the following genetic correlations: trait 1 & trait 2 = 0.8, trait 1 & trait 3 = 0.5, trait 2 & trait 3 = 0.4

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
```python
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

### Spike & Slab Model
Simulate a phenotype with heritability = 0.1 and probability of a SNP being causal = 0.01

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

Simulate two correlated phenotypes with heritabilities 0.8 and 0.9, genetic correlation of 0.5, and the following probabilities of SNPs being causal: probability a SNP is causal for both traits = 0.3, probability SNP is causal for trait 1 but not trait 2 = 0.1, probability SNP is causal for trait 2 but not trait 1 = 0.2. Expected proportion of SNPs causal for trait 1: 0.1 + 0.3 = 0.4, expected proportion of SNPs causal for trait 2: 0.1 + 0.2 = 0.3

```python
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.8, 0.9], pi=[0.3, 0.1, 0.2], rg =0.5)
```

### Annotation-Informed Model
[see method docstring for more information]

Assume for this example we have the following MatrixTable `mt`:

```python
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
* `mt.a1`,`mt.a2`,`mt.a3` : Annotations we want to use for an annotation-informed model

### Population Stratification
[see method docstring for more information]
* `mt.PC1`,`mt.PC2`,`mt.PC3` : Covariates we want to use for population stratification
