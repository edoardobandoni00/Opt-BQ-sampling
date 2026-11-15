# Numerical sampler for Bach's optimized leverage-score distribution
#
#   q*_lambda(x) = pi(x) / d(lambda) *
#     sum_m [mu_m / (mu_m + lambda)] e_m(x)^2.
#
# The Mercer operator is approximated on equal-probability cells of pi. This
# avoids choosing an arbitrary finite x-domain and treats Gaussian and
# heavy-tailed Student-t measures with the same construction.

build_bach_distribution <- function(
    kernel_fun,
    ell,
    sigma_f2,
    quantile_fun,
    n_total,
    alpha,
    dimension = 1,
    measure_name = "target",
    n_grid = 801L,
    lambda = NULL,
    eigen_tolerance = 1e-12) {
  n_grid <- as.integer(n_grid)
  if (n_grid < 51L) {
    stop("n_grid must be at least 51.")
  }
  if (n_total < 1 || alpha <= 0 || dimension <= 0) {
    stop("n_total, alpha, and dimension must be positive.")
  }
  if (is.null(lambda)) {
    lambda <- n_total^(-2 * alpha / dimension)
  }
  if (!is.finite(lambda) || lambda <= 0) {
    stop("lambda must be finite and positive.")
  }
  requested_lambda <- lambda

  # Midpoint quadrature in probability space: int h(x) dpi(x)
  # is approximated by mean_j h(Q_pi((j - 1/2) / M)).
  probabilities <- (seq_len(n_grid) - 0.5) / n_grid
  nodes <- as.numeric(quantile_fun(probabilities))
  if (length(nodes) != n_grid || any(!is.finite(nodes))) {
    stop("quantile_fun must return one finite value per probability.")
  }
  quadrature_weights <- rep(1 / n_grid, n_grid)

  kernel_matrix <- kernel_fun(
    nodes,
    nodes,
    ell = ell,
    sigma_f2 = sigma_f2
  )
  operator_matrix <- kernel_matrix / n_grid
  eigensystem <- eigen(operator_matrix, symmetric = TRUE)
  all_eigenvalues <- pmax(eigensystem$values, 0)
  largest_eigenvalue <- max(all_eigenvalues)
  if (!is.finite(largest_eigenvalue) || largest_eigenvalue <= 0) {
    stop("The discretized Mercer operator has no positive eigenvalues.")
  }

  numerical_eigenvalue_floor <- largest_eigenvalue * max(
    eigen_tolerance,
    10 * n_grid * .Machine$double.eps
  )
  keep <- all_eigenvalues > numerical_eigenvalue_floor
  eigenvalues <- all_eigenvalues[keep]
  eigenvectors <- eigensystem$vectors[, keep, drop = FALSE]
  # Truncate unresolved numerical modes, but preserve the requested
  # regularization parameter on every resolved eigenpair.
  lambda <- requested_lambda
  spectral_weights <- eigenvalues / (eigenvalues + requested_lambda)
  effective_dimension <- sum(spectral_weights)
  if (!is.finite(effective_dimension) || effective_dimension <= 0) {
    stop("Bach effective dimension is not finite and positive.")
  }

  # Nyström eigenfunctions satisfy e_m(x_j) = U_jm / sqrt(w_j).
  leverage <- as.vector((eigenvectors^2) %*% spectral_weights) /
    quadrature_weights
  cell_mass <- quadrature_weights * leverage / effective_dimension
  cell_mass <- pmax(cell_mass, 0)
  cell_mass <- cell_mass / sum(cell_mass)

  if (any(!is.finite(cell_mass)) || abs(sum(cell_mass) - 1) > 1e-10) {
    stop("Bach cell probabilities failed to normalize.")
  }

  structure(
    list(
      nodes = nodes,
      probability_edges = seq(0, 1, length.out = n_grid + 1L),
      cell_mass = cell_mass,
      quantile_fun = quantile_fun,
      lambda = lambda,
      requested_lambda = requested_lambda,
      numerical_eigenvalue_floor = numerical_eigenvalue_floor,
      effective_dimension = effective_dimension,
      eigenvalues = eigenvalues,
      spectral_weights = spectral_weights,
      leverage = leverage,
      measure_name = measure_name,
      n_grid = n_grid,
      n_total = n_total,
      alpha = alpha,
      dimension = dimension,
      diagnostics = list(
        probability_mass = sum(cell_mass),
        minimum_cell_mass = min(cell_mass),
        maximum_cell_mass = max(cell_mass),
        retained_eigenpairs = length(eigenvalues),
        discarded_trace = sum(all_eigenvalues[!keep]),
        operator_trace = sum(all_eigenvalues)
      )
    ),
    class = "bach_distribution"
  )
}

sample_bach <- function(distribution, n = 1L) {
  if (!inherits(distribution, "bach_distribution")) {
    stop("distribution must be created by build_bach_distribution().")
  }
  n <- as.integer(n)
  if (n < 1L) {
    return(numeric())
  }

  cell <- sample.int(
    distribution$n_grid,
    size = n,
    replace = TRUE,
    prob = distribution$cell_mass
  )
  lower <- distribution$probability_edges[cell]
  upper <- distribution$probability_edges[cell + 1L]
  probability <- lower + runif(n) * (upper - lower)
  probability <- pmin(pmax(probability, .Machine$double.eps), 1 - .Machine$double.eps)
  as.numeric(distribution$quantile_fun(probability))
}

bach_metadata <- function(distribution) {
  if (!inherits(distribution, "bach_distribution")) {
    stop("distribution must be created by build_bach_distribution().")
  }
  list(
    lambda = distribution$lambda,
    requested_lambda = distribution$requested_lambda,
    numerical_eigenvalue_floor = distribution$numerical_eigenvalue_floor,
    effective_dimension = distribution$effective_dimension,
    measure = distribution$measure_name,
    n_grid = distribution$n_grid,
    n_total = distribution$n_total,
    alpha = distribution$alpha,
    dimension = distribution$dimension,
    diagnostics = distribution$diagnostics
  )
}

