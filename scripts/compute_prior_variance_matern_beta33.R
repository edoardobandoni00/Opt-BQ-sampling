source("scripts/_setup.R")

relative_tolerance <- 1e-11
target_shape <- 3

hyperparameters32 <- readRDS(
  "results/hyperparams_matern_beta33.rds"
)
hyperparameters72 <- readRDS(
  "results/hyperparams_matern72_beta33.rds"
)

result <- list(
  matern32 = matern_prior_var_beta(
    kernel_fun = matern32_kernel,
    ell = hyperparameters32$posterior$ell_mean,
    sigma_f2 = hyperparameters32$posterior$sigmaf2_mean,
    shape1 = target_shape,
    shape2 = target_shape,
    rel_tol = relative_tolerance
  ),
  matern72 = matern_prior_var_beta(
    kernel_fun = matern72_kernel,
    ell = hyperparameters72$posterior$ell_mean,
    sigma_f2 = hyperparameters72$posterior$sigmaf2_mean,
    shape1 = target_shape,
    shape2 = target_shape,
    rel_tol = relative_tolerance
  ),
  relative_tolerance = relative_tolerance,
  target_shape = target_shape
)

saveRDS(
  result,
  file = "results/prior_variance_matern_beta33_quadrature.rds"
)
print(result)

