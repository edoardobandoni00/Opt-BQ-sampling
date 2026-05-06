source("scripts/_setup.R")

set.seed(1)

f <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

true_integral <- 1

# Hyperparameters are fixed from this RDS file and reused in all replications.
hyp <- readRDS("results/hyperparams_matern_gaussian.rds")
ell <- hyp$posterior$ell_mean
sigma_f2 <- hyp$posterior$sigmaf2_mean
n <- 500
mean_plot_start <- 100L
n_max_mean_plot <- 500L

tot_n <- 10L
mc_seed <- 5000L
s_draws_ci <- 100L

sample_standard <- function(i, n_total) {
  sample_point("standard", n_total, measure = "gaussian")
}
sample_optimal <- function(i, n_total) {
  sample_point("optimal", n_total, measure = "gaussian")
}

prior_var <- matern32_prior_var_gaussian(ell = ell, sigma_f2 = sigma_f2)

run_std <- function() {
  run_bq(
    n = n,
    integrand = f,
    sample_fun = sample_standard,
    kernel_fun = matern32_kernel,
    kernel_mean_fun = matern32_kernel_mean_gaussian,
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
    kernel_mean_fun = matern32_kernel_mean_gaussian,
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
      true_integral = true_integral,
      kernel = "matern32",
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
  file = "results/bq_matern_gaussian.rds"
)

print(c(standard_var_n = tail(res_standard$var, 1), optimal_var_n = tail(res_optimal$var, 1)))

label_std <- "N(0,1)"
label_opt <- "N(0,log(n))"

plot_bq_posterior_mean_bands(
  res_optimal,
  res_standard,
  true_value = true_integral,
  label_optimal = label_opt,
  label_standard = label_std,
  iter_start = mean_plot_start,
  n_max = n_max_mean_plot,
  outfile = "figures/bq_matern_gaussian_mean_bands.pdf"
)

plot_bq_posterior_variance_log(
  res_optimal,
  res_standard,
  label_optimal = label_opt,
  label_standard = label_std,
  outfile = "figures/bq_matern_gaussian_variance.pdf"
)

plot_bq_averaged_mean_ci(
  summ_optimal,
  summ_standard,
  true_value = true_integral,
  label_optimal = label_opt,
  label_standard = label_std,
  iter_start = mean_plot_start,
  n_max = n_max_mean_plot,
  outfile = "figures/bq_matern_gaussian_mean_ci_avg.pdf"
)

plot_bq_averaged_variance_log(
  summ_optimal,
  summ_standard,
  label_optimal = label_opt,
  label_standard = label_std,
  outfile = "figures/bq_matern_gaussian_variance_avg.pdf"
)
