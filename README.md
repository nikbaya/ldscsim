# ldscsim
`ldscsim` is a package containing methods designed to generate simulated phenotypes using different models. 
<br>
## Getting started
`simulate()` is the method which wraps all other methods in the package. However, all methods are self-contained. For instance, if you just want to calculate simulated phenotypes using betas (SNP effects) generated outside of the `ldscsim` framework, you can use `sim_phenotypes()` independent of `simulate()`.
<br>

## Model options
#### Models for SNP effects or "betas"
- **Infinitesimal Model (stable)**: All SNP effects are drawn from a normal distribution with mean=0, variance=`h2/M`, where `h2` is the desired heritability of the simulated trait and `M` is the number of SNPs in the genotype matrix given given by the user.
- **Spike & Slab Model (stable)**: SNPs have probability `pi` of being causal. If causal, the SNP effect is drawn from a normal distribution with variance `h2/(M*pi)`. If not causal, the SNP effect is zero.
- **Annotation-Informed Betas (in progress)**: The effect for SNP `j` are drawn from a normal distribution with mean=0, variance=`a[j]*h2`, where `a[j]` is the relative heritability contributed by SNP `j` and `a` is a vector of the relative heritability contributed by each SNP. The sum of `a` is 1.

#### Models for population stratification (in progress)
After calculating the phenotype by multiplying genotypes by betas and then adding environmental noise, a term is added with population stratification. This term is `popstrat_s2*popstrat`. `popstrat` is a column field containing the desired population stratification, normalized in the code to have mean=0, variance=1. `popstrat_s2` is the desired relative amount of variance contributed by the population stratification. The simulated phenotype has variance=1 before adding population stratification so if `popstrat_s2`=2 then population stratification will have variance twice that of the phenotype before adding population stratification.
<br>


## Methods
###### `simulate(mt, genotype, h2, pi=1, annot=None, popstrat=None, popstrat_s2=1, path_to_save=None)`
Simulates phenotypes under various models.
- **`mt`** (required):  Hail MatrixTable holding all the matrix fields passed as parameters to the model
- **`genotype`** (required): `mt` entry field of genotypes. Accepted data types: `expr_int32`,`expr_int64`,`expr_float32`,`expr_float64`
- **`h2`** (required): Desired heritability of simulated trait. Accepted data types: `float`, `int`
- **`pi`** (default: 1) : Probability of SNP being causal. Used for spike & slab model. Accepted data types: `float`,`int`
- **`annot`** (default: `None`) : `mt` row field of annotations for annotation-informed model. Accepted data types: `expr_int32`,`expr_int64` 
- **`popstrat`** (default: `None`) : `mt` column field of population stratification. Accepted data types: `expr_int32`,`float_64`
- **`popstrat_s2`** (default: 1) : Amount of variance contributed by population stratification. Accepted data types: `float`,`int`
- **`path_to_save`** (default: `None`) : Path to save MatrixTable containing simulation output. Accepted data type: `str`
<br>

###### `print_header(h2, pi, annot, popstrat, popstrat_s2, path_to_save)`
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

<!-- google-site-verification: google9796b225e0522c44.html -->
