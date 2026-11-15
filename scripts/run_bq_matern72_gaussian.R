source("scripts/_setup.R")

integrand <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

true_integral <- 1
kernel_nu <- 7 / 2
dimension <- 1
sobolev_alpha <- kernel_nu + dimension / 2
gamma_temper <- 2 * dimension / (2 * sobolev_alpha + dimension)
power_sd <- sqrt(1 / gamma_temper)

hyperparameters <- readRDS("results/hyperparams_matern72_gaussian.rds")
ell <- hyperparameters$posterior$ell_mean
sigma_f2 <- hyperparameters$posterior$sigmaf2_mean

n <- env_integer("BQ_N", 300L)
n_rep <- env_integer("BQ_N_REP", 1000L)
bach_grid_n <- env_integer("BACH_GRID_N", 1500L, minimum = 51L)

sample_target <- function(iteration, n_total) {
  rnorm(1, mean = 0, sd = 1)
}
sample_power <- function(iteration, n_total) {
  rnorm(1, mean = 0, sd = power_sd)
}

run_bq_comparison(
  experiment_id = "bq_matern72_gaussian",
  n = n,
  integrand = integrand,
  kernel_fun = matern72_kernel,
  kernel_mean_fun = matern72_kernel_mean_gaussian,
  prior_var = matern72_prior_var_gaussian(ell, sigma_f2),
  ell = ell,
  sigma_f2 = sigma_f2,
  sample_target = sample_target,
  sample_power = sample_power,
  target_quantile = qnorm,
  alpha = sobolev_alpha,
  dimension = dimension,
  measure_name = "N(0,1)",
  scheme_labels = c(
    power = "pi^gamma",
    bach = "Bach q*_lambda",
    target = "Target: N(0,1)"
  ),
  settings = list(
    true_integral = true_integral,
    gamma_temper = gamma_temper,
    power_sd = power_sd,
    kernel_nu = kernel_nu,
    sobolev_alpha = sobolev_alpha,
    kernel = "matern72"
  ),
  n_rep = n_rep,
  replication_seed = 6000L,
  single_seed = 1L,
  ci_draws = 100L,
  mean_plot_start = min(100L, n),
  mean_plot_end = n,
  plots = c("mean_bands", "mean_ci_avg"),
  bach_grid_n = bach_grid_n
)

