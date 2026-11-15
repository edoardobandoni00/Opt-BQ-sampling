# Shared replication summaries and plotting style.

bq_scheme_order <- c("target", "bach", "power")
bq_scheme_colors <- c(
  target = "#0072B2",
  bach = "#D55E00",
  power = "#009E73"
)
variance_log_plot_floor <- 1e-20

scale_iteration_log10 <- function() {
  ggplot2::scale_x_log10(
    breaks = function(limits) {
      lower <- max(limits[1], .Machine$double.eps)
      candidates <- c(1, 2, 5, 10, 20, 40, 50, 100, 200, 300, 500, 1000)
      candidates[candidates >= lower & candidates <= limits[2]]
    },
    labels = function(breaks) {
      format(breaks, trim = TRUE, scientific = FALSE)
    }
  )
}

theme_bq <- function() {
  ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(
        fill = "white",
        color = "black"
      )
    )
}

validate_schemes <- function(objects, labels) {
  if (is.null(names(objects)) || is.null(names(labels))) {
    stop("Results and labels must be named by scheme ID.")
  }
  if (!identical(sort(names(objects)), sort(bq_scheme_order)) ||
      !identical(sort(names(labels)), sort(bq_scheme_order))) {
    stop("Results and labels must contain power, bach, and target.")
  }
  bq_scheme_order
}

scheme_legend_labels <- function(labels) {
  if (!all(bq_scheme_order %in% names(labels))) {
    stop("Labels must be named with power, bach, and target.")
  }
  list(
    target = expression(pi),
    bach = expression(q[lambda]^{"*"}),
    power = expression(pi^gamma)
  )[bq_scheme_order]
}

combined_matern_labels <- function() {
  list(
    target_32 = expression(pi * ", Mat\u00e9rn 3/2"),
    target_72 = expression(pi * ", Mat\u00e9rn 7/2"),
    bach_32 = expression(q[lambda]^{"*"} * ", Mat\u00e9rn 3/2"),
    bach_72 = expression(q[lambda]^{"*"} * ", Mat\u00e9rn 7/2"),
    power_32 = expression(pi^gamma * ", Mat\u00e9rn 3/2"),
    power_72 = expression(pi^gamma * ", Mat\u00e9rn 7/2")
  )
}

scheme_scales <- function(labels) {
  legend_labels <- scheme_legend_labels(labels)
  list(
    ggplot2::scale_color_manual(
      values = bq_scheme_colors,
      breaks = bq_scheme_order,
      labels = legend_labels,
      drop = FALSE
    ),
    ggplot2::scale_fill_manual(
      values = bq_scheme_colors,
      breaks = bq_scheme_order,
      labels = legend_labels,
      drop = FALSE
    )
  )
}

replicate_bq_schemes <- function(n_rep, n, seed_start, run_funs) {
  if (is.null(names(run_funs)) || any(!nzchar(names(run_funs)))) {
    stop("run_funs must be a named list.")
  }
  scheme_names <- names(run_funs)
  scheme_count <- length(run_funs)
  output <- setNames(
    lapply(scheme_names, function(name) {
      list(
        mean_mat = matrix(NA_real_, nrow = n_rep, ncol = n),
        var_mat = matrix(NA_real_, nrow = n_rep, ncol = n)
      )
    }),
    scheme_names
  )

  for (replication in seq_len(n_rep)) {
    for (scheme_index in seq_along(run_funs)) {
      set.seed(
        seed_start +
          scheme_count * (replication - 1L) +
          scheme_index - 1L
      )
      result <- run_funs[[scheme_index]]()
      if (length(result$mu) != n || length(result$var) != n) {
        stop(sprintf(
          "run_bq result for '%s' must have length n.",
          scheme_names[scheme_index]
        ))
      }
      output[[scheme_index]]$mean_mat[replication, ] <- result$mu
      output[[scheme_index]]$var_mat[replication, ] <- result$var
    }
  }
  output
}

summarize_repeated_runs <- function(mean_mat, var_mat, s_draws = 100L) {
  if (!identical(dim(mean_mat), dim(var_mat))) {
    stop("mean_mat and var_mat must have identical dimensions.")
  }
  replication_count <- nrow(mean_mat)
  posterior_mean <- colMeans(mean_mat)
  conditional_variance <- colMeans(var_mat)
  design_variance <- if (replication_count < 2L) {
    rep(0, ncol(mean_mat))
  } else {
    colMeans(
      (mean_mat - matrix(
        posterior_mean,
        nrow = replication_count,
        ncol = ncol(mean_mat),
        byrow = TRUE
      ))^2
    )
  }
  posterior_variance <- pmax(conditional_variance + design_variance, 0)

  repeated_means <- mean_mat[
    rep(seq_len(replication_count), each = s_draws),
    ,
    drop = FALSE
  ]
  repeated_variances <- var_mat[
    rep(seq_len(replication_count), each = s_draws),
    ,
    drop = FALSE
  ]
  draws <- repeated_means + sqrt(pmax(repeated_variances, 0)) *
    matrix(
      rnorm(replication_count * s_draws * ncol(mean_mat)),
      nrow = replication_count * s_draws
    )

  list(
    mean = posterior_mean,
    var = posterior_variance,
    conditional_var = conditional_variance,
    design_var = design_variance,
    ci_lower = apply(draws, 2, quantile, probs = 0.025),
    ci_upper = apply(draws, 2, quantile, probs = 0.975)
  )
}

make_mean_plot_data <- function(objects, labels, iter_start, n_max, averaged) {
  scheme_ids <- validate_schemes(objects, labels)
  n <- if (averaged) length(objects[[1]]$mean) else length(objects[[1]]$mu)
  n_max <- min(as.integer(n_max), n)
  iter_start <- max(1L, min(as.integer(iter_start), n_max))
  iterations <- seq.int(iter_start, n_max)

  do.call(
    rbind,
    lapply(scheme_ids, function(scheme_id) {
      object <- objects[[scheme_id]]
      if (averaged) {
        data.frame(
          iteration = iterations,
          center = object$mean[iterations],
          lower = object$ci_lower[iterations],
          upper = object$ci_upper[iterations],
          scheme = factor(scheme_id, levels = bq_scheme_order)
        )
      } else {
        standard_deviation <- sqrt(pmax(object$var[iterations], 0))
        data.frame(
          iteration = iterations,
          center = object$mu[iterations],
          lower = object$mu[iterations] - 2 * standard_deviation,
          upper = object$mu[iterations] + 2 * standard_deviation,
          scheme = factor(scheme_id, levels = bq_scheme_order)
        )
      }
    })
  )
}

make_variance_plot_data <- function(
    objects,
    labels,
    iter_start,
    averaged,
    variance_field = "var") {
  scheme_ids <- validate_schemes(objects, labels)
  if (any(vapply(
    objects[scheme_ids],
    function(object) is.null(object[[variance_field]]),
    logical(1)
  ))) {
    stop(sprintf("Variance field '%s' is missing.", variance_field))
  }
  n <- length(objects[[1]][[variance_field]])
  iter_start <- max(1L, min(as.integer(iter_start), n))
  iterations <- seq.int(iter_start, n)

  do.call(
    rbind,
    lapply(scheme_ids, function(scheme_id) {
      data.frame(
        iteration = iterations,
        variance = pmax(
          objects[[scheme_id]][[variance_field]][iterations],
          variance_log_plot_floor
        ),
        scheme = factor(scheme_id, levels = bq_scheme_order)
      )
    })
  )
}

save_bq_plot <- function(plot, outfile, width, height, dpi) {
  if (!is.null(outfile)) {
    ggplot2::ggsave(
      filename = outfile,
      plot = plot,
      width = width,
      height = height,
      dpi = dpi
    )
  }
  invisible(plot)
}

plot_bq_posterior_mean_bands_multi <- function(
    results,
    labels,
    true_value,
    iter_start = 1L,
    n_max = NULL,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  if (is.null(n_max)) {
    n_max <- length(results[[1]]$mu)
  }
  plot_data <- make_mean_plot_data(
    results,
    labels,
    iter_start,
    n_max,
    averaged = FALSE
  )
  scales <- scheme_scales(labels)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = iteration,
      y = center,
      color = scheme,
      fill = scheme
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(
      yintercept = true_value,
      color = "grey35",
      linewidth = 0.5,
      linetype = "dotted"
    ) +
    scales +
    ggplot2::labs(x = NULL, y = NULL, color = "Proposal", fill = "Proposal") +
    theme_bq() +
    ggplot2::theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0)
    )

  save_bq_plot(plot, outfile, width, height, dpi)
}

plot_bq_posterior_variance_log_multi <- function(
    results,
    labels,
    iter_start = 10L,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  plot_data <- make_variance_plot_data(
    results,
    labels,
    iter_start,
    averaged = FALSE
  )
  scales <- scheme_scales(labels)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = iteration, y = variance, color = scheme)
  ) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_iteration_log10() +
    ggplot2::scale_y_log10() +
    scales[[1]] +
    ggplot2::labs(x = NULL, y = NULL, color = "Proposal") +
    theme_bq() +
    ggplot2::theme(
      legend.position = c(0.05, 0.05),
      legend.justification = c(0, 0)
    )

  save_bq_plot(plot, outfile, width, height, dpi)
}

plot_bq_averaged_mean_ci_multi <- function(
    summaries,
    labels,
    true_value,
    iter_start = 1L,
    n_max = NULL,
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  if (is.null(n_max)) {
    n_max <- length(summaries[[1]]$mean)
  }
  plot_data <- make_mean_plot_data(
    summaries,
    labels,
    iter_start,
    n_max,
    averaged = TRUE
  )
  scales <- scheme_scales(labels)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = iteration,
      y = center,
      color = scheme,
      fill = scheme
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(
      yintercept = true_value,
      color = "grey35",
      linewidth = 0.5,
      linetype = "dotted"
    ) +
    scales +
    ggplot2::labs(x = NULL, y = NULL, color = "Proposal", fill = "Proposal") +
    theme_bq() +
    ggplot2::theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0)
    )

  save_bq_plot(plot, outfile, width, height, dpi)
}

plot_bq_averaged_variance_log_multi <- function(
    summaries,
    labels,
    iter_start = 10L,
    variance_field = "var",
    outfile = NULL,
    width = 8,
    height = 5,
    dpi = 150) {
  plot_data <- make_variance_plot_data(
    summaries,
    labels,
    iter_start,
    averaged = TRUE,
    variance_field = variance_field
  )
  scales <- scheme_scales(labels)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = iteration, y = variance, color = scheme)
  ) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_iteration_log10() +
    ggplot2::scale_y_log10() +
    scales[[1]] +
    ggplot2::labs(x = NULL, y = NULL, color = "Proposal") +
    theme_bq() +
    ggplot2::theme(
      legend.position = c(0.05, 0.05),
      legend.justification = c(0, 0)
    )

  save_bq_plot(plot, outfile, width, height, dpi)
}
