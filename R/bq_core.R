run_bq <- function(
    n,
    integrand,
    sample_fun,
    kernel_fun,
    kernel_mean_fun,
    prior_var,
    ell,
    sigma_f2 = 1,
    jitter = 1e-12) {
  mu <- numeric(n)
  var <- numeric(n)

  x_nodes <- sample_fun(1, n)
  y_nodes <- integrand(x_nodes)

  k11 <- as.numeric(kernel_fun(x_nodes, x_nodes, ell, sigma_f2)) + jitter
  l_chol <- matrix(sqrt(k11), 1, 1)
  rho <- kernel_mean_fun(x_nodes, ell, sigma_f2)

  alpha <- backsolve(t(l_chol), forwardsolve(l_chol, rho))
  beta <- forwardsolve(l_chol, rho)

  mu[1] <- sum(alpha * y_nodes)
  var[1] <- max(prior_var - sum(beta^2), jitter)

  if (n == 1) {
    return(list(mu = mu, var = var, x = x_nodes, y = y_nodes))
  }

  for (i in 2:n) {
    x_new <- sample_fun(i, n)
    y_new <- integrand(x_new)

    k_old_new <- as.numeric(kernel_fun(x_nodes, x_new, ell, sigma_f2))
    v <- forwardsolve(l_chol, k_old_new)
    s_new <- sqrt(max(sigma_f2 - sum(v^2) + jitter, jitter))

    n_old <- length(x_nodes)
    l_chol <- rbind(
      cbind(l_chol, rep(0, n_old)),
      c(v, s_new)
    )

    x_nodes <- c(x_nodes, x_new)
    y_nodes <- c(y_nodes, y_new)
    rho <- c(rho, kernel_mean_fun(x_new, ell, sigma_f2))

    alpha <- backsolve(t(l_chol), forwardsolve(l_chol, rho))
    beta <- forwardsolve(l_chol, rho)
    mu[i] <- sum(alpha * y_nodes)
    var[i] <- max(prior_var - sum(beta^2), jitter)
  }

  list(mu = mu, var = var, x = x_nodes, y = y_nodes)
}
