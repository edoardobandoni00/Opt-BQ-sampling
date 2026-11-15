source("scripts/_setup.R")

set.seed(321)

f <- function(x) {
  1+sin(2*pi*x)#exp(-((x - 0.5) / 0.08)^2) + 0.25 * sin(4 * pi * x)
}

# Keep this training design identical to the Matérn 3/2 fit.
n <- 100L
x <- rbeta(n, shape1 = 1.8, shape2 = 1.8)
y <- f(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "matern72",
  n_iter = 1000,
  burn = 200
)

saveRDS(fit, file = "results/hyperparams_matern72_beta33.rds")
print(fit$posterior)
