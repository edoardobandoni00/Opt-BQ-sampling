source("scripts/_setup.R")

set.seed(123)

# Same dominating measure and integrand as BQ under Student-t (Matern-3/2 throughout).
df_measure <- 4.49
df_target <- 5

f <- function(x) {
  (1 + sin(2 * pi * x)) * dt(x, df = df_target) / dt(x, df = df_measure)
}

n <- 100
x <- rt(n, df = df_measure)
y <- f(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "matern32",
  n_iter = 1000,
  burn = 200
)

saveRDS(
  fit,
  file = "results/hyperparams_matern_student.rds"
)
print(fit$posterior)
