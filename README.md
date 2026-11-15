# Optimal Sampling for Bayesian Quadrature

This repository contains the R implementation of the numerical experiments associated with the article Rate-Optimal Randomised Kernel Quadrature. The experiments compare sampling distributions for one-dimensional Bayesian quadrature and investigate how the choice of evaluation nodes affects integration accuracy and posterior uncertainty on bounded and unbounded domains.

## Experiments

The code covers three integration measures:

- the standard Gaussian distribution;
- the Student-$t$ distribution with 8 degrees of freedom;
- the Beta(3,3) distribution on $[0,1]$.

Each measure is tested with Matérn $3/2$ and Matérn $7/2$ kernels, giving
six experiments in total. Every experiment compares:

- **Tempered target:** the normalized density $\pi^\gamma$;
- **Optimized sampling:** Bach's regularized leverage-score distribution
  $q_\lambda^\ast$;
- **Target sampling:** direct sampling from the integration measure $\pi$.

The integrand, fitted Gaussian-process hyperparameters, kernel mean, and prior
integral variance are held fixed within each comparison.

## Repository structure

- `R/` contains the reusable implementation:
  - `kernels_measures.R`: kernels, kernel means, and prior integral variances;
  - `hyperparameters.R`: Metropolis-within-Gibbs hyperparameter estimation;
  - `bq_core.R`: sequential BQ with Cholesky updates;
  - `bach_sampling.R`: numerical construction of $q_\lambda^\ast$;
  - `experiments.R`: shared experiment workflow;
  - `postprocess.R`: replication summaries and plotting functions.
- `scripts/` contains the executable scripts for hyperparameter estimation,
  prior-variance calculations, BQ experiments, and combined plots.
- `results/` is intentionally empty and receives generated `.rds` files.
- `figures/` is intentionally empty and receives generated PDF figures.

## Requirements

- R
- `ggplot2`
- `cubature`

Install the required packages with:

```r
install.packages(c("ggplot2", "cubature"))
```

## Running the experiments

Run all commands from the repository root.

First estimate the kernel hyperparameters and compute the required prior
variances:

```r
source("scripts/run_hyperparams_matern_gaussian.R")
source("scripts/run_hyperparams_matern72_gaussian.R")

source("scripts/run_hyperparams_matern_student.R")
source("scripts/run_hyperparams_matern72_student.R")
source("scripts/compute_prior_variance_matern_student.R")

source("scripts/run_hyperparams_matern_beta33.R")
source("scripts/run_hyperparams_matern72_beta33.R")
source("scripts/compute_prior_variance_matern_beta33.R")
```

Then run the BQ comparisons and combined variance plots:

```r
source("scripts/run_bq_matern_gaussian.R")
source("scripts/run_bq_matern72_gaussian.R")
source("scripts/plot_gaussian_variance_matern_combined.R")

source("scripts/run_bq_matern_student.R")
source("scripts/run_bq_matern72_student.R")
source("scripts/plot_student_variance_matern_combined.R")

source("scripts/run_bq_matern_beta33.R")
source("scripts/run_bq_matern72_beta33.R")
source("scripts/plot_beta33_variance_matern_combined.R")
```

## Outputs

Each BQ script writes an `.rds` file to `results/`. It contains the experiment
settings, diagnostics for Bach's distribution, a representative run for each
sampling method, and summaries across replications.

The plotting routines write PDF files to `figures/`, including posterior-mean
bands, averaged confidence intervals, and conditional posterior-variance
comparisons. The reported pooled variance follows the law of total variance.

## Reduced test runs

The main scripts use the full experiment budgets and 100 replications. A small
run can be used to check the installation:

```r
Sys.setenv(BQ_N = 20, BQ_N_REP = 2, BACH_GRID_N = 101)
source("scripts/run_bq_matern_gaussian.R")
```

For the Student-$t$ and Beta experiments, numerical integration uses a
relative tolerance of `1e-11`. Student-$t$ calculations are transformed to
the unit interval before calling `integrate()` or `cubature::hcubature()`.
