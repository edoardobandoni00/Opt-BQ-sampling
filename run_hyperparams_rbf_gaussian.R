source("scripts/_setup.R")

set.seed(123)

f <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

n <- 100
x <- rnorm(n, mean = 0, sd = sqrt(log(n)))
y <- f(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "rbf",
  n_iter = 1000,
  burn = 200
)

saveRDS(fit, file = "results/hyperparams_rbf_gaussian.rds")
print(fit$posterior)
