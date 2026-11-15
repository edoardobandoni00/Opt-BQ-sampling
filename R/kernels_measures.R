# Kernels and kernel integrals used by the retained experiments.

matern32_kernel <- function(x, y, ell, sigma_f2 = 1) {
  distance <- abs(outer(x, y, "-"))
  scaled_distance <- sqrt(3) * distance / ell
  value <- sigma_f2 * (1 + scaled_distance) * exp(-scaled_distance)
  value[is.infinite(scaled_distance)] <- 0
  value
}

# Kernel mean of Matérn 3/2 against Lebesgue/uniform on [a, b]:
#   ∫_a^b (1 + λ|x-y|) exp(-λ|x-y|) dy,  λ = √3 / ℓ.
# The antiderivative of (1 + λu) exp(-λu) is -exp(-λu)(u + 2/λ).
matern32_kernel_mean_uniform <- function(x, ell, sigma_f2 = 1, a = 0, b = 1) {
  lambda <- sqrt(3) / ell
  two_over_lambda <- 2 / lambda
  antideriv <- function(u) {
    -exp(-lambda * u) * (u + two_over_lambda)
  }
  x <- as.numeric(x)
  vapply(
    x,
    function(xi) {
      if (!is.finite(xi)) {
        return(0)
      }
      if (xi <= a) {
        val <- antideriv(b - xi) - antideriv(a - xi)
      } else if (xi >= b) {
        val <- antideriv(xi - a) - antideriv(xi - b)
      } else {
        left <- xi - a
        right <- b - xi
        val <- (two_over_lambda - exp(-lambda * left) * (left + two_over_lambda)) +
          (two_over_lambda - exp(-lambda * right) * (right + two_over_lambda))
      }
      sigma_f2 * val
    },
    numeric(1)
  )
}

matern32_prior_var_uniform <- function(ell, sigma_f2 = 1, a = 0, b = 1) {
  lambda <- sqrt(3) / ell
  length <- b - a
  if (!is.finite(length) || length <= 0) {
    stop("Uniform kernel-mean interval must have positive length.")
  }
  exp_term <- exp(-lambda * length)
  sigma_f2 * (
    4 * length / lambda -
      6 / lambda^2 +
      2 * exp_term * (length / lambda + 3 / lambda^2)
  )
}

matern32_kernel_mean_gaussian <- function(x, ell, sigma_f2 = 1) {
  matern_kernel_mean_gaussian_poly(
    x = x,
    ell = ell,
    sigma_f2 = sigma_f2,
    lambda = sqrt(3) / ell,
    coefs = c(1, 1)
  )
}

matern32_prior_var_gaussian <- function(ell, sigma_f2 = 1) {
  matern_prior_var_gaussian_poly(
    sigma_f2 = sigma_f2,
    lambda = sqrt(3) / ell,
    coefs = c(1, 1)
  )
}

# Half-integer Matérn kernel means against N(0, 1). For
# P(t) = sum_k coefs[k + 1] t^k and lambda = sqrt(2 nu) / ell, the required
# half-line integrals satisfy I_k = (k - 1) I_{k-2} - beta I_{k-1}.
matern_kernel_mean_gaussian_poly <- function(x, ell, sigma_f2, lambda, coefs) {
  degree <- length(coefs) - 1L

  half_line_contribution <- function(beta, x_value) {
    integrals <- numeric(degree + 1L)
    # This is phi(x) I_0(beta), evaluated jointly to avoid Inf * 0
    # when x is far in a Gaussian tail.
    integrals[1] <- exp(
      (beta^2 - x_value^2) / 2 +
        pnorm(-beta, log.p = TRUE)
    )
    if (degree >= 1L) {
      integrals[2] <- dnorm(x_value) - beta * integrals[1]
    }
    if (degree >= 2L) {
      for (k in 2:degree) {
        integrals[k + 1L] <- (k - 1) * integrals[k - 1L] -
          beta * integrals[k]
      }
    }
    sigma_f2 * sum(coefs * lambda^(0:degree) * integrals)
  }

  vapply(
    x,
    function(x_value) {
      half_line_contribution(lambda - x_value, x_value) +
        half_line_contribution(lambda + x_value, x_value)
    },
    numeric(1)
  )
}

matern_prior_var_gaussian_poly <- function(sigma_f2, lambda, coefs) {
  degree <- length(coefs) - 1L
  integrals <- numeric(degree + 1L)
  integrals[1] <- exp(
    log(2 * sqrt(pi)) + lambda^2 +
      pnorm(-sqrt(2) * lambda, log.p = TRUE)
  )
  if (degree >= 1L) {
    integrals[2] <- 2 - 2 * lambda * integrals[1]
  }
  if (degree >= 2L) {
    for (k in 2:degree) {
      integrals[k + 1L] <- 2 * (k - 1) * integrals[k - 1L] -
        2 * lambda * integrals[k]
    }
  }
  sigma_f2 * sum(coefs * lambda^(0:degree) * integrals) / sqrt(pi)
}

matern72_kernel <- function(x, y, ell, sigma_f2 = 1) {
  distance <- abs(outer(x, y, "-"))
  scaled_distance <- sqrt(7) * distance / ell
  value <- sigma_f2 * (
    1 + scaled_distance + 2 * scaled_distance^2 / 5 +
      scaled_distance^3 / 15
  ) * exp(-scaled_distance)
  value[is.infinite(scaled_distance)] <- 0
  value
}

matern72_kernel_mean_gaussian <- function(x, ell, sigma_f2 = 1) {
  matern_kernel_mean_gaussian_poly(
    x = x,
    ell = ell,
    sigma_f2 = sigma_f2,
    lambda = sqrt(7) / ell,
    coefs = c(1, 1, 2 / 5, 1 / 15)
  )
}

matern72_prior_var_gaussian <- function(ell, sigma_f2 = 1) {
  matern_prior_var_gaussian_poly(
    sigma_f2 = sigma_f2,
    lambda = sqrt(7) / ell,
    coefs = c(1, 1, 2 / 5, 1 / 15)
  )
}

matern_kernel_mean_beta <- function(
    x,
    ell,
    sigma_f2,
    kernel_fun,
    shape1,
    shape2,
    rel_tol = 1e-11) {
  if (shape1 <= 0 || shape2 <= 0) {
    stop("Beta shapes must be positive.")
  }
  if (any(!is.finite(x)) || any(x < 0) || any(x > 1)) {
    stop("x must contain finite values in [0, 1].")
  }
  log_normaliser <- lbeta(shape1, shape2)
  vapply(
    x,
    function(x_value) {
      integrate(
        function(theta) {
          y <- sin(theta / 2)^2
          weight <- exp(
            (2 * shape1 - 1) * log(pmax(sin(theta / 2), .Machine$double.eps)) +
              (2 * shape2 - 1) * log(pmax(cos(theta / 2), .Machine$double.eps)) -
              log_normaliser
          )
          as.numeric(kernel_fun(
            x_value,
            y,
            ell = ell,
            sigma_f2 = sigma_f2
          )) * weight
        },
        lower = 0,
        upper = pi,
        rel.tol = rel_tol,
        subdivisions = 200L
      )$value
    },
    numeric(1)
  )
}

matern32_kernel_mean_beta <- function(
    x,
    ell,
    sigma_f2 = 1,
    shape1,
    shape2) {
  matern_kernel_mean_beta(
    x,
    ell,
    sigma_f2,
    matern32_kernel,
    shape1,
    shape2
  )
}

matern72_kernel_mean_beta <- function(
    x,
    ell,
    sigma_f2 = 1,
    shape1,
    shape2) {
  matern_kernel_mean_beta(
    x,
    ell,
    sigma_f2,
    matern72_kernel,
    shape1,
    shape2
  )
}

matern_prior_var_beta <- function(
    kernel_fun,
    ell,
    sigma_f2 = 1,
    shape1,
    shape2,
    rel_tol = 1e-11) {
  if (ell <= 0 || sigma_f2 <= 0 || shape1 <= 0 || shape2 <= 0 ||
      rel_tol <= 0) {
    stop("ell, sigma_f2, Beta shapes, and rel_tol must be positive.")
  }
  if (!requireNamespace("cubature", quietly = TRUE)) {
    stop("Install the 'cubature' package to compute Beta prior variances.")
  }

  result <- cubature::hcubature(
    f = function(point) {
      as.numeric(kernel_fun(
        point[1],
        point[2],
        ell = ell,
        sigma_f2 = sigma_f2
      )) *
        dbeta(point[1], shape1 = shape1, shape2 = shape2) *
        dbeta(point[2], shape1 = shape1, shape2 = shape2)
    },
    lowerLimit = c(0, 0),
    upperLimit = c(1, 1),
    tol = rel_tol
  )
  list(
    value = result$integral,
    abs.error = result$error,
    evaluations = result$functionEvaluations,
    return.code = result$returnCode
  )
}

matern_kernel_mean_student <- function(
    x,
    ell,
    sigma_f2,
    kernel_fun,
    df = 8,
    rel_tol = 1e-11) {
  if (ell <= 0 || sigma_f2 <= 0 || df <= 0 || rel_tol <= 0) {
    stop("ell, sigma_f2, df, and rel_tol must be positive.")
  }
  if (any(!is.finite(x))) {
    stop("x must contain only finite values.")
  }
  values <- numeric(length(x))
  for (index in seq_along(x)) {
    x_value <- x[[index]]
    integrand <- function(probability) {
      probability <- pmin(
        pmax(probability, .Machine$double.eps),
        1 - .Machine$double.eps
      )
      as.numeric(kernel_fun(
        x_value,
        qt(probability, df = df),
        ell = ell,
        sigma_f2 = sigma_f2
      ))
    }
    # Split at F(x) so the Matérn peak sits at an endpoint of each piece.
    # An absolute tolerance is required in the tails, where the integral is
    # near zero and a pure relative tolerance is unattainable.
    peak <- min(
      max(pt(x_value, df = df), .Machine$double.eps),
      1 - .Machine$double.eps
    )
    integrate_piece <- function(lower, upper) {
      if (upper <= lower) {
        return(0)
      }
      result <- integrate(
        integrand,
        lower = lower,
        upper = upper,
        rel.tol = rel_tol,
        abs.tol = rel_tol,
        subdivisions = 200L,
        stop.on.error = FALSE
      )
      value <- result$value
      if (length(value) == 1L && is.finite(value)) value else 0
    }
    values[index] <- integrate_piece(0, peak) + integrate_piece(peak, 1)
  }
  values
}

matern32_kernel_mean_student <- function(
    x,
    ell,
    sigma_f2 = 1,
    df = 8,
    rel_tol = 1e-11) {
  matern_kernel_mean_student(
    x,
    ell,
    sigma_f2,
    matern32_kernel,
    df,
    rel_tol
  )
}

matern72_kernel_mean_student <- function(
    x,
    ell,
    sigma_f2 = 1,
    df = 8,
    rel_tol = 1e-11) {
  matern_kernel_mean_student(
    x,
    ell,
    sigma_f2,
    matern72_kernel,
    df,
    rel_tol
  )
}

matern_prior_var_student <- function(
    kernel_fun,
    ell,
    sigma_f2 = 1,
    df = 8,
    rel_tol = 1e-11) {
  if (ell <= 0 || sigma_f2 <= 0 || df <= 0 || rel_tol <= 0) {
    stop("ell, sigma_f2, df, and rel_tol must be positive.")
  }
  if (!requireNamespace("cubature", quietly = TRUE)) {
    stop("Install the 'cubature' package to compute Student prior variances.")
  }

  result <- cubature::hcubature(
    f = function(probability) {
      probability <- pmin(
        pmax(probability, .Machine$double.eps),
        1 - .Machine$double.eps
      )
      as.numeric(kernel_fun(
        qt(probability[1], df = df),
        qt(probability[2], df = df),
        ell = ell,
        sigma_f2 = sigma_f2
      ))
    },
    lowerLimit = c(0, 0),
    upperLimit = c(1, 1),
    tol = rel_tol
  )
  list(
    value = result$integral,
    abs.error = result$error,
    evaluations = result$functionEvaluations,
    return.code = result$returnCode
  )
}
