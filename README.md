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

Assume for all examples that we have the following MatrixTable `mt`:

```python
>>> mt.describe()
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
––––––––––––––––––––––––––––––––––––––––"
```

* `mt.s` : Sample IDs
* `mt.PC1`,`mt.PC2`,`mt.PC3` : Covariates we want to use for population stratification
* `mt.rsid` : SNP IDs
* `mt.a1`,`mt.a2`,`mt.a3` : Annotations we want to use for an annotation-informed model
* `mt.gt` : Genotypes of individuals

### Infinitesimal Model
Simulate a phenotype under the infinitesimal model with heritability = 0.1

```python
>>> sim = simulate_phenotypes(mt=mt,genotype=mt.gt,h2=0.1)
```

Simulate two phenotypes under the infinitesimal model with heritabilities 0.3 and 0.4 and genetic correlation of 0.6

```python
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.3, 0.4], rg=0.6)
```

Simulate three phenotypes under the infinitesimal model with heritabilities 0.1, 0.2, 0.7 and the following genetic correlations: trait 1 & trait 2 = 0.8, trait 1 & trait 3 = 0.5, trait 2 & trait 3 = 0.4

```python
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.1, 0.2, 0.7], rg=[0.8, 0.5, 0.4])
```

### Spike & Slab Model
Simulate a phenotype with heritability = 0.1 and probability of a SNP being causal = 0.01

```python
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=0.1, pi=0.01)
```

Simulate two correlated phenotypes with heritabilities 0.8 and 0.9, genetic correlation of 0.5, and the following probabilities of SNPs being causal: probability a SNP is causal for both traits = 0.3, probability SNP is causal for trait 1 but not trait 2 = 0.1, probability SNP is causal for trait 2 but not trait 1 = 0.2. Expected proportion of SNPs causal for trait 1: 0.1 + 0.3 = 0.4, expected proportion of SNPs causal for trait 2: 0.1 + 0.2 = 0.3

```python
>>> sim = simulate_phenotypes(mt=mt, genotype=mt.gt, h2=[0.8, 0.9], pi=[0.3, 0.1, 0.2], rg =0.5)
```

### Annotation-Informed Model
[see method docstring for more information]

### Population Stratification
[see method docstring for more information]
