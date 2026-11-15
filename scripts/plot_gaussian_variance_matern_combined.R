# Combined average posterior-variance plot for the Gaussian-target
# experiments: Matérn 3/2 and Matérn 7/2, three proposals each.
# Reads the saved RDS files and does not rerun the experiments.

source("scripts/_setup.R")

matern32 <- readRDS(file.path(BQ_UD_ROOT, "results", "bq_matern_gaussian.rds"))
matern72 <- readRDS(file.path(BQ_UD_ROOT, "results", "bq_matern72_gaussian.rds"))

if (is.null(matern32$replicated) || is.null(matern72$replicated)) {
  stop("Both RDS files must contain the replicated summaries.")
}

# Two shades of each family colour used in postprocess.R.
# Darker / solid = Matérn 7/2; lighter / dashed = Matérn 3/2.
series_colors <- c(
  target_32 = "#56B4E9",
  target_72 = "#0072B2",
  bach_32 = "#E69F00",
  bach_72 = "#D55E00",
  power_32 = "#66C2A5",
  power_72 = "#009E73"
)
series_linetypes <- c(
  target_32 = "dashed",
  target_72 = "solid",
  bach_32 = "dashed",
  bach_72 = "solid",
  power_32 = "dashed",
  power_72 = "solid"
)
series_labels <- combined_matern_labels()
series_order <- names(series_colors)

extract_averaged_variance <- function(
    experiment,
    kernel_tag,
    variance_field,
    normalizer = 1) {
  summaries <- experiment$replicated
  missing <- setdiff(bq_scheme_order, names(summaries))
  if (length(missing) > 0L) {
    stop("Missing replicated schemes: ", paste(missing, collapse = ", "))
  }
  do.call(
    rbind,
    lapply(bq_scheme_order, function(scheme_id) {
      variance <- summaries[[scheme_id]][[variance_field]] / normalizer
      variance <- pmax(variance, variance_log_plot_floor)
      data.frame(
        iteration = seq_along(variance),
        variance = variance,
        series = paste0(scheme_id, "_", kernel_tag),
        stringsAsFactors = FALSE
      )
    })
  )
}

make_combined_plot <- function(
    variance_field,
    normalizers = c(`32` = 1, `72` = 1),
    y_label = NULL) {
  plot_data <- rbind(
    extract_averaged_variance(
      matern32,
      "32",
      variance_field,
      normalizers[["32"]]
    ),
    extract_averaged_variance(
      matern72,
      "72",
      variance_field,
      normalizers[["72"]]
    )
  )
  plot_data <- plot_data[
    plot_data$iteration >= 40L & plot_data$iteration <= 300L,
    ,
    drop = FALSE
  ]
  plot_data$series <- factor(plot_data$series, levels = series_order)

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = iteration,
      y = variance,
      color = series,
      linetype = series
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    scale_iteration_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(
      values = series_colors,
      breaks = series_order,
      labels = unname(series_labels),
      drop = FALSE
    ) +
    ggplot2::scale_linetype_manual(
      values = series_linetypes,
      breaks = series_order,
      labels = unname(series_labels),
      drop = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_label,
      color = "Proposal",
      linetype = "Proposal"
    ) +
    theme_bq() +
    ggplot2::theme(
      legend.position = c(0.05, 0.05),
      legend.justification = c(0, 0)
    )
}

save_combined_plot <- function(plot, stem) {
  save_bq_plot(
    plot,
    file.path(BQ_UD_ROOT, "figures", paste0(stem, ".pdf")),
    width = 8,
    height = 5,
    dpi = 150
  )
}

conditional_plot <- make_combined_plot("conditional_var")
save_combined_plot(
  conditional_plot,
  "bq_gaussian_matern32_matern72_conditional_variance_avg"
)

message("Wrote combined Gaussian variance plots to ", file.path(
  BQ_UD_ROOT,
  "figures"
))
invisible(conditional_plot)
