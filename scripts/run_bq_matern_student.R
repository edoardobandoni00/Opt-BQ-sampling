source("scripts/_setup.R")

set.seed(1)

# Target integral: integral of (1 + sin(2*pi*x)) with respect to t_5.
# BQ is run under t_4.49 as dominating measure for the GP surrogate.
# The integrand includes the change of measure ratio t_5 / t_4.49.
df_measure <- 4.49
df_target <- 5
d <- 1
# Matern 3/2 kernel smoothness parameter.
alpha_matern <- 3 / 2

f <- function(x) {
  (1 + sin(2 * pi * x)) * dt(x, df = df_target) / dt(x, df = df_measure)
}

# Integral of (1 + sin(2*pi*x)) with respect to the target measure t_5.
# The sine term integrates to 0 by symmetry, so the true value is 1.
true_integral <- 1

# Hyperparameters are fixed from this RDS only (never refit inside replication loops).
hyp <- readRDS("results/hyperparams_matern_student.rds")
ell <- hyp$posterior$ell_mean
sigma_f2 <- hyp$posterior$sigmaf2_mean
n <- 500
mean_plot_start <- 100L
n_max_mean_plot <- 250L
tot_n <- 100L
mc_seed <- 7000L
s_draws_ci <- 100L

exp_opt <- 2 * alpha_matern / (alpha_matern + df_measure + d / 2)
label_opt <- sprintf("t_4.49(0, n^{%s})", format(round(exp_opt, 4), nsmall = 4, trim = TRUE))
label_std <- "t_4.49(0,1)"

sample_standard <- function(i, n_total) {
  sample_point("standard", n_total, measure = "student", df = df_measure)
}
sample_optimal <- function(i, n_total) {
  sample_point(
    "optimal",
    n_total,
    measure = "student",
    df = df_measure,
    matern_smoothness = alpha_matern,
    dimension = d
  )
}

prior_var <- matern32_prior_var_student_mc(ell = ell, sigma_f2 = sigma_f2, df = df_measure, n_mc = 1e8)

kernel_mean_student <- function(x, ell_arg, sigma_f2_arg) {
  matern32_kernel_mean_student(x, ell = ell_arg, sigma_f2 = sigma_f2_arg, df = df_measure)
}

run_std <- function() {
  run_bq(
    n = n,
    integrand = f,
    sample_fun = sample_standard,
    kernel_fun = matern32_kernel,
    kernel_mean_fun = kernel_mean_student,
    prior_var = prior_var,
    ell = ell,
    sigma_f2 = sigma_f2
  )
}

run_opt <- function() {
  run_bq(
    n = n,
    integrand = f,
    sample_fun = sample_optimal,
    kernel_fun = matern32_kernel,
    kernel_mean_fun = kernel_mean_student,
    prior_var = prior_var,
    ell = ell,
    sigma_f2 = sigma_f2
  )
}

res_standard <- run_std()
res_optimal <- run_opt()

rep <- replicate_bq_two_schemes(tot_n, n, mc_seed, run_opt, run_std)

summ_optimal <- summarize_repeated_runs(
  rep$optimal$mean_mat,
  rep$optimal$var_mat,
  s_draws = s_draws_ci
)
summ_standard <- summarize_repeated_runs(
  rep$standard$mean_mat,
  rep$standard$var_mat,
  s_draws = s_draws_ci
)

saveRDS(
  list(
    settings = list(
      n = n,
      ell = ell,
      sigma_f2 = sigma_f2,
      df_measure = df_measure,
      df_target = df_target,
      d = d,
      alpha_matern = alpha_matern,
      true_integral = true_integral,
      optimal_sampling_exponent = exp_opt,
      n_rep = tot_n,
      mc_seed = mc_seed,
      s_draws_ci = s_draws_ci
    ),
    standard = res_standard,
    optimal = res_optimal,
    replicated = list(
      summ_optimal = summ_optimal,
      summ_standard = summ_standard
    )
  ),
  file = "results/bq_matern_student.rds"
)

print(c(standard_var_n = tail(res_standard$var, 1), optimal_var_n = tail(res_optimal$var, 1)))

plot_bq_posterior_mean_bands(
  res_optimal,
  res_standard,
  true_value = true_integral,
  label_optimal = label_opt,
  label_standard = label_std,
  iter_start = mean_plot_start,
  n_max = n_max_mean_plot,
  outfile = "figures/bq_matern_student_mean_bands.pdf"
)

plot_bq_posterior_variance_log(
  res_optimal,
  res_standard,
  label_optimal = label_opt,
  label_standard = label_std,
  outfile = "figures/bq_matern_student_variance.pdf"
)

plot_bq_averaged_mean_ci(
  summ_optimal,
  summ_standard,
  true_value = true_integral,
  label_optimal = label_opt,
  label_standard = label_std,
  iter_start = mean_plot_start,
  n_max = n_max_mean_plot,
  outfile = "figures/bq_matern_student_mean_ci_avg.pdf"
)

plot_bq_averaged_variance_log(
  summ_optimal,
  summ_standard,
  label_optimal = label_opt,
  label_standard = label_std,
  outfile = "figures/bq_matern_student_variance_avg.pdf"
)
