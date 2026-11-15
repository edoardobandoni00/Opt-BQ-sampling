source("scripts/_setup.R")

integrand <- function(x) {
  1 + sin(2*pi*x)
}

target_shape <- 3
true_integral <- integrate(
  function(x) integrand(x) * dbeta(x, target_shape, target_shape),
  lower = 0,
  upper = 1,
  rel.tol = 1e-12
)$value
kernel_nu <- 7 / 2
dimension <- 1
sobolev_alpha <- kernel_nu + dimension / 2
gamma_temper <- 2 * dimension / (2 * sobolev_alpha + dimension)
power_shape <- 1 + gamma_temper * (target_shape - 1)

hyperparameters <- readRDS(
  "results/hyperparams_matern72_beta33.rds"
)
ell <- hyperparameters$posterior$ell_mean
sigma_f2 <- hyperparameters$posterior$sigmaf2_mean
prior_variances <- readRDS(
  "results/prior_variance_matern_beta33_quadrature.rds"
)

n <- env_integer("BQ_N", 80L)
n_rep <- env_integer("BQ_N_REP", 1000L)
bach_grid_n <- env_integer("BACH_GRID_N", 1500L, minimum = 51L)

kernel_mean <- function(x, ell, sigma_f2 = 1) {
  matern72_kernel_mean_beta(
    x,
    ell,
    sigma_f2,
    shape1 = target_shape,
    shape2 = target_shape
  )
}
sample_target <- function(iteration, n_total) {
  rbeta(1, shape1 = target_shape, shape2 = target_shape)
}
sample_power <- function(iteration, n_total) {
  rbeta(1, shape1 = power_shape, shape2 = power_shape)
}
target_quantile <- function(probability) {
  qbeta(probability, shape1 = target_shape, shape2 = target_shape)
}

run_bq_comparison(
  experiment_id = "bq_matern72_beta33",
  n = n,
  integrand = integrand,
  kernel_fun = matern72_kernel,
  kernel_mean_fun = kernel_mean,
  prior_var = prior_variances$matern72$value,
  ell = ell,
  sigma_f2 = sigma_f2,
  sample_target = sample_target,
  sample_power = sample_power,
  target_quantile = target_quantile,
  alpha = sobolev_alpha,
  dimension = dimension,
  measure_name = "Beta(3,3)",
  scheme_labels = c(
    power = "pi^gamma",
    bach = "Bach q*_lambda",
    target = "Target: Beta(3,3)"
  ),
  settings = list(
    true_integral = true_integral,
    gamma_temper = gamma_temper,
    target_shape = target_shape,
    power_shape = power_shape,
    kernel_nu = kernel_nu,
    sobolev_alpha = sobolev_alpha,
    kernel = "matern72",
    kernel_mean_relative_tolerance = 1e-11,
    prior_variance_relative_tolerance =
      prior_variances$relative_tolerance,
    prior_variance_absolute_error =
      prior_variances$matern72$abs.error
  ),
  n_rep = n_rep,
  replication_seed = 8000L,
  single_seed = 1L,
  ci_draws = 100L,
  mean_plot_start = min(40L, n),
  mean_plot_end = n,
  plots = c("mean_bands", "mean_ci_avg"),
  bach_grid_n = bach_grid_n
)

