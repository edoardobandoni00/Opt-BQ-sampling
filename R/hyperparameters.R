gp_mwg_hyperparams <- function(
    x,
    y,
    kernel_name = "rbf",
    n_iter = 1000,
    burn = 200,
    sigma2_jitter = 1e-8,
    ell_init = 1,
    sigmaf2_init = 2,
    proposal_sd_log_ell = 0.2,
    a_f = 2,
    b_f = 2) {
  kernel_fun <- switch(
    kernel_name,
    rbf = function(a, b, ell) exp(-(outer(a, b, "-")^2) / (2 * ell^2)),
    matern32 = function(a, b, ell) {
      d <- abs(outer(a, b, "-"))
      lambda <- sqrt(3) / ell
      (1 + lambda * d) * exp(-lambda * d)
    },
    stop("kernel_name must be 'rbf' or 'matern32'.")
  )

  log_prior_ell <- function(ell, mu = 0, sd = 100) {
    if (ell <= 0) {
      return(-Inf)
    }
    dlnorm(ell, meanlog = mu, sdlog = sd, log = TRUE)
  }

  ell_samples <- numeric(n_iter)
  sigmaf2_samples <- numeric(n_iter)
  ell_samples[1] <- ell_init
  sigmaf2_samples[1] <- sigmaf2_init
  accept_ell <- 0

  n <- length(y)
  eye <- diag(n)

  for (t in 2:n_iter) {
    ell_curr <- ell_samples[t - 1]
    u_curr <- log(ell_curr)
    u_prop <- u_curr + rnorm(1, 0, proposal_sd_log_ell)
    ell_prop <- exp(u_prop)

    k_curr <- kernel_fun(x, x, ell_curr) + sigma2_jitter * eye
    k_prop <- kernel_fun(x, x, ell_prop) + sigma2_jitter * eye

    quad_curr <- as.numeric(t(y) %*% solve(k_curr, y))
    quad_prop <- as.numeric(t(y) %*% solve(k_prop, y))

    logpost_curr <- -0.5 * n * log(sigmaf2_samples[t - 1]) -
      0.5 * quad_curr / sigmaf2_samples[t - 1] +
      log_prior_ell(ell_curr) + u_curr
    logpost_prop <- -0.5 * n * log(sigmaf2_samples[t - 1]) -
      0.5 * quad_prop / sigmaf2_samples[t - 1] +
      log_prior_ell(ell_prop) + u_prop

    if (log(runif(1)) < (logpost_prop - logpost_curr)) {
      ell_samples[t] <- ell_prop
      accept_ell <- accept_ell + 1
    } else {
      ell_samples[t] <- ell_curr
    }

    k_ell <- kernel_fun(x, x, ell_samples[t]) + sigma2_jitter * eye
    shape <- a_f + n / 2
    rate <- b_f + 0.5 * as.numeric(t(y) %*% solve(k_ell, y))
    sigmaf2_samples[t] <- 1 / rgamma(1, shape = shape, rate = rate)
  }

  idx <- burn:n_iter
  list(
    samples = data.frame(iteration = seq_len(n_iter), ell = ell_samples, sigmaf2 = sigmaf2_samples),
    posterior = list(
      ell_mean = mean(ell_samples[idx]),
      sigmaf2_mean = mean(sigmaf2_samples[idx]),
      ell_ci95 = quantile(ell_samples[idx], c(0.025, 0.975)),
      sigmaf2_ci95 = quantile(sigmaf2_samples[idx], c(0.025, 0.975)),
      acceptance_ell = accept_ell / n_iter
    )
  )
}
