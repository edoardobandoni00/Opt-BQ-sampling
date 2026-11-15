source("scripts/_setup.R")

set.seed(123)

f <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

# Keep this training design identical to the Matérn 3/2 fit.
n <- 100
training_sd <- sqrt(9 / 2)
x <- rnorm(n, mean = 0, sd = training_sd)
y <- f(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "matern72",
  n_iter = 1000,
  burn = 200
)

saveRDS(fit, file = "results/hyperparams_matern72_gaussian.rds")
print(fit$posterior)

