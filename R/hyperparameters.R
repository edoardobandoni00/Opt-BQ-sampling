# Metropolis-within-Gibbs estimation of Matérn length scale and amplitude.

gp_mwg_hyperparams <- function(
    x,
    y,
    kernel_name = c("matern32", "matern72"),
    n_iter = 1000,
    burn = 200,
    sigma2_jitter = 1e-8,
    ell_init = 1,
    sigmaf2_init = 2,
    proposal_sd_log_ell = 0.2,
    a_f = 2,
    b_f = 2) {
  kernel_name <- match.arg(kernel_name)
  if (burn < 0L || burn >= n_iter) {
    stop("burn must be between 0 and n_iter - 1.")
  }
  if (length(x) != length(y) || length(y) < 1L ||
      any(!is.finite(x)) || any(!is.finite(y))) {
    stop("x and y must be finite, non-empty vectors of equal length.")
  }

  unit_kernel <- switch(
    kernel_name,
    matern32 = function(a, b, ell) matern32_kernel(a, b, ell, sigma_f2 = 1),
    matern72 = function(a, b, ell) matern72_kernel(a, b, ell, sigma_f2 = 1)
  )

  log_prior_ell <- function(ell, meanlog = 0, sdlog = 100) {
    if (ell <= 0) {
      return(-Inf)
    }
    dlnorm(ell, meanlog = meanlog, sdlog = sdlog, log = TRUE)
  }

  ell_samples <- numeric(n_iter)
  sigmaf2_samples <- numeric(n_iter)
  ell_samples[1] <- ell_init
  sigmaf2_samples[1] <- sigmaf2_init
  accepted_ell <- 0L

  n <- length(y)
  identity_matrix <- diag(n)

  kernel_statistics <- function(ell) {
    kernel_matrix <- unit_kernel(x, x, ell) +
      sigma2_jitter * identity_matrix
    chol_factor <- chol(kernel_matrix)
    solved_y <- backsolve(
      chol_factor,
      forwardsolve(t(chol_factor), y)
    )
    list(
      log_determinant = 2 * sum(log(diag(chol_factor))),
      quadratic = as.numeric(crossprod(y, solved_y))
    )
  }

  for (iteration in 2:n_iter) {
    ell_current <- ell_samples[iteration - 1L]
    log_ell_current <- log(ell_current)
    log_ell_proposal <- log_ell_current +
      rnorm(1, mean = 0, sd = proposal_sd_log_ell)
    ell_proposal <- exp(log_ell_proposal)

    statistics_current <- kernel_statistics(ell_current)
    statistics_proposal <- kernel_statistics(ell_proposal)

    current_amplitude <- sigmaf2_samples[iteration - 1L]
    log_posterior_current <- -0.5 * n * log(current_amplitude) -
      0.5 * statistics_current$log_determinant -
      0.5 * statistics_current$quadratic / current_amplitude +
      log_prior_ell(ell_current) + log_ell_current
    log_posterior_proposal <- -0.5 * n * log(current_amplitude) -
      0.5 * statistics_proposal$log_determinant -
      0.5 * statistics_proposal$quadratic / current_amplitude +
      log_prior_ell(ell_proposal) + log_ell_proposal

    if (log(runif(1)) < log_posterior_proposal - log_posterior_current) {
      ell_samples[iteration] <- ell_proposal
      selected_statistics <- statistics_proposal
      accepted_ell <- accepted_ell + 1L
    } else {
      ell_samples[iteration] <- ell_current
      selected_statistics <- statistics_current
    }

    shape <- a_f + n / 2
    rate <- b_f + 0.5 * selected_statistics$quadratic
    sigmaf2_samples[iteration] <- 1 / rgamma(1, shape = shape, rate = rate)
  }

  retained <- seq.int(burn + 1L, n_iter)
  list(
    samples = data.frame(
      iteration = seq_len(n_iter),
      ell = ell_samples,
      sigmaf2 = sigmaf2_samples
    ),
    posterior = list(
      ell_mean = mean(ell_samples[retained]),
      sigmaf2_mean = mean(sigmaf2_samples[retained]),
      ell_ci95 = quantile(ell_samples[retained], c(0.025, 0.975)),
      sigmaf2_ci95 = quantile(
        sigmaf2_samples[retained],
        c(0.025, 0.975)
      ),
      acceptance_ell = accepted_ell / (n_iter - 1L)
    )
  )
}
