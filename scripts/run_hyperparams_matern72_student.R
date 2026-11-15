source("scripts/_setup.R")

set.seed(123)

df <- 8
integrand <- function(x) {
  sqrt(3) * exp(-x^2) + sin(2 * pi * x) / (1 + x^2)
}

# Keep this training design identical to the Matérn 3/2 fit.
dimension <- 1
kernel_nu <- 7 / 2
sobolev_alpha <- kernel_nu + dimension / 2
gamma_temper <- 2 * dimension / (2 * sobolev_alpha + dimension)
power_df <- gamma_temper * (df + 1) - 1
power_scale <- sqrt(df / power_df)

n <- 100L
x <- power_scale * rt(n, df = power_df)
y <- integrand(x)

fit <- gp_mwg_hyperparams(
  x = x,
  y = y,
  kernel_name = "matern72",
  n_iter = 1000,
  burn = 200
)

saveRDS(fit, file = "results/hyperparams_matern72_student.rds")
print(fit$posterior)
