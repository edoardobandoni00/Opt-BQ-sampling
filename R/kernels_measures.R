rbf_kernel <- function(x, y, ell, sigma_f2 = 1) {
  sigma_f2 * exp(-(outer(x, y, "-")^2) / (2 * ell^2))
}

# RBF / squared-exponential kernel against standard Gaussian N(0, 1) (1D).
rbf_kernel_mean_gaussian <- function(x, ell, sigma_f2 = 1) {
  sf <- sqrt(ell^2 / (1 + ell^2))
  sigma_f2 * sf * exp(-x^2 / (2 * (1 + ell^2)))
}

rbf_prior_var_gaussian <- function(ell, sigma_f2 = 1) {
  sigma_f2 * ell / sqrt(2 + ell^2)
}

matern32_kernel <- function(x, y, ell, sigma_f2 = 1) {
  d <- abs(outer(x, y, "-"))
  lambda <- sqrt(3) / ell
  sigma_f2 * (1 + lambda * d) * exp(-lambda * d)
}

erfc <- function(z) {
  2 * pnorm(-z * sqrt(2))
}

matern32_kernel_mean_gaussian <- function(x, ell, sigma_f2 = 1) {
  lambda <- sqrt(3) / ell
  beta1 <- lambda - x
  beta2 <- lambda + x
  exp_pref <- exp(-x^2 / 2)

  term <- function(beta) {
    i0 <- sqrt(pi / 2) * exp(beta^2 / 2) * erfc(beta / sqrt(2))
    i1 <- 1 - beta * i0
    sigma_f2 * (i0 + lambda * i1)
  }

  exp_pref / sqrt(2 * pi) * (term(beta1) + term(beta2))
}

matern32_prior_var_gaussian <- function(ell, sigma_f2 = 1) {
  lambda <- sqrt(3) / ell
  pref <- 2 / sqrt(4 * pi)
  i0 <- sqrt(pi) * exp(lambda^2) * erfc(lambda)
  i1 <- 2 - 2 * lambda * i0
  sigma_f2 * pref * (i0 + lambda * i1)
}

matern32_kernel_mean_student <- function(x, ell, sigma_f2 = 1, df = 8) {
  lambda <- sqrt(3) / ell
  vapply(
    x,
    function(xi) {
      integrand <- function(u) {
        weight <- (1 + lambda * u) * exp(-lambda * u)
        density <- dt(xi - u, df = df) + dt(xi + u, df = df)
        weight * density
      }
      sigma_f2 * integrate(integrand, lower = 0, upper = Inf, rel.tol = 1e-10)$value
    },
    numeric(1)
  )
}

matern32_prior_var_student_mc <- function(ell, sigma_f2 = 1, df = 8, n_mc = 100000) {
  lambda <- sqrt(3) / ell
  x <- rt(n_mc, df = df)
  y <- rt(n_mc, df = df)
  r <- abs(x - y)
  mean(sigma_f2 * (1 + lambda * r) * exp(-lambda * r))
}
