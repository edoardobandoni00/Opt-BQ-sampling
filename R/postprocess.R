# Plotting helpers aligned with BQ_UD_myversion (theme_bw, ribbons, legend placement).

theme_bq_myversion <- function() {
  ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = "white", color = "black")
    )
}

# Posterior mean with plus/minus 2 posterior standard deviations
# for two proposal schemes, with optional true integral line.
# res_optimal and res_standard are lists returned by run_bq()
# and must contain vectors mu and var.
plot_bq_posterior_mean_bands <- function(
    res_optimal,
    res_standard,
    true_value,
    label_optimal,
    label_standard,
    iter_start = 1L,
    n_max = NULL,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  n <- length(res_optimal$mu)
  if (length(res_standard$mu) != n) {
    stop("res_optimal and res_standard must have the same length.")
  }
  if (is.null(n_max)) {
    n_max <- n
  }
  n_max <- min(as.integer(n_max), n)
  idx <- iter_start:n_max

  df1 <- data.frame(
    t = idx,
    mu = res_optimal$mu[idx],
    sd = sqrt(res_optimal$var[idx]),
    Proposals = label_optimal
  )
  df2 <- data.frame(
    t = idx,
    mu = res_standard$mu[idx],
    sd = sqrt(res_standard$var[idx]),
    Proposals = label_standard
  )
  df_all <- rbind(df1, df2)

  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = t, y = mu, color = Proposals, fill = Proposals)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = mu - 2 * sd, ymax = mu + 2 * sd),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = true_value, color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::labs(x = "Iteration n", y = "Posterior Mean +/- 2 SD") +
    theme_bq_myversion() +
    ggplot2::theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0)
    )

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  invisible(p)
}

# Posterior variance of the integral on log y scale for two proposal schemes.
plot_bq_posterior_variance_log <- function(
    res_optimal,
    res_standard,
    label_optimal,
    label_standard,
    iter_start = 1L,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150,
    x_label = "Iteration n",
    y_label = "Posterior Variance",
    legend_title = "Proposal") {
  n <- length(res_optimal$var)
  if (length(res_standard$var) != n) {
    stop("res_optimal and res_standard must have the same length.")
  }
  idx <- iter_start:n

  df_var_all <- rbind(
    data.frame(iterations = idx, sigma_n_squared = res_optimal$var[idx], series = label_optimal),
    data.frame(iterations = idx, sigma_n_squared = res_standard$var[idx], series = label_standard)
  )

  p <- ggplot2::ggplot(df_var_all, ggplot2::aes(x = iterations, y = sigma_n_squared, color = series)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = x_label, y = y_label, color = legend_title) +
    theme_bq_myversion() +
    ggplot2::theme(
      legend.position = c(0.95, 0.75),
      legend.justification = c(1, 0)
    )

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  invisible(p)
}

# Replicate BQ for two sampling schemes (for example optimal and standard).
# Each replication uses independent RNG streams for the two runs:
# seed_start + 2*r - 2 and seed_start + 2*r - 1.
replicate_bq_two_schemes <- function(n_rep, n, seed_start, run_optimal, run_standard) {
  means_o <- matrix(NA_real_, n_rep, n)
  vars_o <- matrix(NA_real_, n_rep, n)
  means_s <- matrix(NA_real_, n_rep, n)
  vars_s <- matrix(NA_real_, n_rep, n)
  for (r in seq_len(n_rep)) {
    set.seed(seed_start + 2L * r - 2L)
    res_o <- run_optimal()
    set.seed(seed_start + 2L * r - 1L)
    res_s <- run_standard()
    if (length(res_o$mu) != n || length(res_s$mu) != n) {
      stop("Each run_bq result must have length n.")
    }
    means_o[r, ] <- res_o$mu
    vars_o[r, ] <- res_o$var
    means_s[r, ] <- res_s$mu
    vars_s[r, ] <- res_s$var
  }
  list(
    optimal = list(mean_mat = means_o, var_mat = vars_o),
    standard = list(mean_mat = means_s, var_mat = vars_s)
  )
}

# Averaged posterior mean with 95 percent confidence interval ribbon.
# This matches the repeated experiment plot style used in BQ_UD_myversion.
# summ_optimal and summ_standard are outputs of summarize_repeated_runs().
plot_bq_averaged_mean_ci <- function(
    summ_optimal,
    summ_standard,
    true_value,
    label_optimal,
    label_standard,
    iter_start = 1L,
    n_max = NULL,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  n <- length(summ_optimal$mean)
  if (length(summ_standard$mean) != n) {
    stop("summ_optimal and summ_standard must have the same length.")
  }
  if (is.null(n_max)) {
    n_max <- n
  }
  n_max <- min(as.integer(n_max), n)
  idx <- iter_start:n_max

  df1 <- data.frame(
    t = idx,
    mu = summ_optimal$mean[idx],
    lower = summ_optimal$ci_lower[idx],
    upper = summ_optimal$ci_upper[idx],
    Proposals = label_optimal
  )
  df2 <- data.frame(
    t = idx,
    mu = summ_standard$mean[idx],
    lower = summ_standard$ci_lower[idx],
    upper = summ_standard$ci_upper[idx],
    Proposals = label_standard
  )
  df_all <- rbind(df1, df2)

  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = t, y = mu, color = Proposals, fill = Proposals)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(yintercept = true_value, color = "red", linewidth = 0.5, linetype = "dashed") +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_bq_myversion() +
    ggplot2::theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0)
    )

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  invisible(p)
}

# Pooled posterior variance using the law of total variance on log scale.
# This matches the style used by the averaged plot in myversion.
plot_bq_averaged_variance_log <- function(
    summ_optimal,
    summ_standard,
    label_optimal,
    label_standard,
    iter_start = 1L,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  n <- length(summ_optimal$var)
  if (length(summ_standard$var) != n) {
    stop("summ_optimal and summ_standard must have the same length.")
  }
  idx <- iter_start:n

  df_var_all <- rbind(
    data.frame(iterations = idx, sigma_n_squared = summ_optimal$var[idx], series = label_optimal),
    data.frame(iterations = idx, sigma_n_squared = summ_standard$var[idx], series = label_standard)
  )

  p <- ggplot2::ggplot(df_var_all, ggplot2::aes(x = iterations, y = sigma_n_squared, color = series)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = NULL, y = NULL, color = "Proposal") +
    theme_bq_myversion() +
    ggplot2::theme(
      legend.position = c(0.95, 0.75),
      legend.justification = c(1, 0)
    )

  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  }
  invisible(p)
}

summarize_repeated_runs <- function(mean_mat, var_mat, s_draws = 50) {
  post_mean <- colMeans(mean_mat)
  post_var <- colMeans(var_mat) + apply(mean_mat, 2, var)

  r <- nrow(mean_mat)
  k <- ncol(mean_mat)
  means_rep <- mean_mat[rep(seq_len(r), each = s_draws), , drop = FALSE]
  vars_rep <- var_mat[rep(seq_len(r), each = s_draws), , drop = FALSE]
  draws <- means_rep + sqrt(vars_rep) * matrix(rnorm(r * s_draws * k), nrow = r * s_draws, ncol = k)

  list(
    mean = post_mean,
    var = post_var,
    ci_lower = apply(draws, 2, quantile, probs = 0.025),
    ci_upper = apply(draws, 2, quantile, probs = 0.975)
  )
}
