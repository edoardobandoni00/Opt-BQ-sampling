# Bayesian Quadrature on Unbounded Domains

This repository contains reproducible code for Bayesian quadrature on unbounded domains, with optimized sampling rules and repeated experiment summaries.

Main goals:

- compare standard and optimized sampling distributions;
- fit GP hyperparameters with MCMC for each kernel and measure setup;
- run one-shot and repeated BQ experiments;
- report averaged posterior mean, pooled posterior variance, and 95 percent confidence intervals.

## Repository layout

- `R/`: shared functions for kernels, measures, sampling, hyperparameters, BQ core, and plotting.
- `scripts/`: experiment entry points.
- `results/`: saved `.rds` outputs.
- `figures/`: exported plots.
- `BQ_UD_myversion/`: reference scripts used for style and checks.

## Experiments implemented

1. RBF kernel with Gaussian measure  
   Script: `scripts/run_bq_rbf_gaussian.R`

2. Matern 3/2 kernel with Gaussian measure  
   Script: `scripts/run_bq_matern_gaussian.R`

3. Matern 3/2 kernel with Student t measure and change of measure  
   Script: `scripts/run_bq_matern_student.R`

## Hyperparameter scripts

- `scripts/run_hyperparams_rbf_gaussian.R`
- `scripts/run_hyperparams_matern_gaussian.R`
- `scripts/run_hyperparams_matern_student.R`

Hyperparameters are loaded once from `results/*.rds` in each BQ script and reused for all replications.

## Quick start

1. Open R in the project root.
2. Install required package:

```r
install.packages("ggplot2")
```

3. Fit hyperparameters:

```r
source("scripts/run_hyperparams_rbf_gaussian.R")
source("scripts/run_hyperparams_matern_gaussian.R")
source("scripts/run_hyperparams_matern_student.R")
```

4. Run BQ experiments:

```r
source("scripts/run_bq_rbf_gaussian.R")
source("scripts/run_bq_matern_gaussian.R")
source("scripts/run_bq_matern_student.R")
```

## Plot settings used now

- Mean trend plots can start after burn-in style iterations:
  - RBF Gaussian: `mean_plot_start <- 70`
  - Matern Gaussian: `mean_plot_start <- 100`
  - Matern Student: `mean_plot_start <- 100`
- This start is used in both:
  - single run mean plus/minus 2 standard deviation plot;
  - averaged mean with 95 percent confidence interval plot.
- Variance plots still use all iterations.

## Repeated experiment summary

For each experiment, repeated runs are combined with:

- posterior mean: average of posterior means across runs;
- pooled posterior variance: law of total variance;
- 95 percent confidence interval: simulated draws from Gaussian posteriors, then empirical quantiles.

## Reproducibility conventions

- fixed random seed per script;
- configuration variables grouped at the top of scripts;
- math and plotting helpers isolated in `R/`;
- no package installation inside analysis scripts.
