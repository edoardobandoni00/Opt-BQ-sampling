source("scripts/_setup.R")

set.seed(321)

f <- function(x) {
  1+sin(2*pi*x)
}

# Common training design for both Matérn kernels. Draw from the flattened
# proposal so the endpoints are seen.
n <- 100L
x <- rbeta(n, shape1 = 1.8, shape2 = 1.8)
y <- f(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "matern32",
  n_iter = 1000,
  burn = 200
)

saveRDS(fit, file = "results/hyperparams_matern_beta33.rds")
print(fit$posterior)
