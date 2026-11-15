source("scripts/_setup.R")

relative_tolerance <- 1e-11
df <- 8

hyperparameters32 <- readRDS("results/hyperparams_matern_student.rds")
hyperparameters72 <- readRDS("results/hyperparams_matern72_student.rds")

result <- list(
  matern32 = matern_prior_var_student(
    kernel_fun = matern32_kernel,
    ell = hyperparameters32$posterior$ell_mean,
    sigma_f2 = hyperparameters32$posterior$sigmaf2_mean,
    df = df,
    rel_tol = relative_tolerance
  ),
  matern72 = matern_prior_var_student(
    kernel_fun = matern72_kernel,
    ell = hyperparameters72$posterior$ell_mean,
    sigma_f2 = hyperparameters72$posterior$sigmaf2_mean,
    df = df,
    rel_tol = relative_tolerance
  ),
  relative_tolerance = relative_tolerance,
  df = df
)

saveRDS(
  result,
  file = "results/prior_variance_matern_student_quadrature.rds"
)
print(result)
