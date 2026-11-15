source("scripts/_setup.R")

df <- 8
integrand <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

true_integral <- integrate(
  function(x) integrand(x) * dt(x, df = df),
  lower = -Inf,
  upper = Inf,
  rel.tol = 1e-12
)$value

kernel_nu <- 7 / 2
dimension <- 1
sobolev_alpha <- kernel_nu + dimension / 2
gamma_temper <- 2 * dimension / (2 * sobolev_alpha + dimension)
power_df <- gamma_temper * (df + 1) - 1
power_scale <- sqrt(df / power_df)

hyperparameters <- readRDS("results/hyperparams_matern72_student.rds")
ell <- hyperparameters$posterior$ell_mean
sigma_f2 <- hyperparameters$posterior$sigmaf2_mean
prior_variances <- readRDS(
  "results/prior_variance_matern_student_quadrature.rds"
)

n <- env_integer("BQ_N", 250L)
n_rep <- env_integer("BQ_N_REP", 100L)
bach_grid_n <- env_integer("BACH_GRID_N", 1500L, minimum = 51L)

kernel_mean <- function(x, ell, sigma_f2 = 1) {
  matern72_kernel_mean_student(x, ell, sigma_f2, df = df)
}
sample_target <- function(iteration, n_total) {
  rt(1, df = df)
}
sample_power <- function(iteration, n_total) {
  power_scale * rt(1, df = power_df)
}
target_quantile <- function(probability) {
  qt(probability, df = df)
}

run_bq_comparison(
  experiment_id = "bq_matern72_student",
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
  measure_name = sprintf("t_%s(0,1)", df),
  scheme_labels = c(
    power = "pi^gamma",
    bach = "Bach q*_lambda",
    target = sprintf("Target: t_%s(0,1)", df)
  ),
  settings = list(
    true_integral = true_integral,
    df = df,
    gamma_temper = gamma_temper,
    power_df = power_df,
    power_scale = power_scale,
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
  replication_seed = 10000L,
  single_seed = 1L,
  ci_draws = 100L,
  mean_plot_start = min(40L, n),
  mean_plot_end = n,
  plots = c("mean_bands", "mean_ci_avg"),
  bach_grid_n = bach_grid_n
)
