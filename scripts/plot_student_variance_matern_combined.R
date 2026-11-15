# Combined average conditional posterior-variance plot for the
# Student-t experiments with Matérn 3/2 and Matérn 7/2 kernels.

source("scripts/_setup.R")

matern32 <- readRDS(
  file.path(BQ_UD_ROOT, "results", "bq_matern_student.rds")
)
matern72 <- readRDS(
  file.path(BQ_UD_ROOT, "results", "bq_matern72_student.rds")
)

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

extract_conditional_variance <- function(experiment, kernel_tag) {
  summaries <- experiment$replicated
  missing <- setdiff(bq_scheme_order, names(summaries))
  if (length(missing) > 0L) {
    stop("Missing replicated schemes: ", paste(missing, collapse = ", "))
  }
  do.call(
    rbind,
    lapply(bq_scheme_order, function(scheme_id) {
      variance <- pmax(
        summaries[[scheme_id]]$conditional_var,
        variance_log_plot_floor
      )
      data.frame(
        iteration = seq_along(variance),
        variance = variance,
        series = paste0(scheme_id, "_", kernel_tag),
        stringsAsFactors = FALSE
      )
    })
  )
}

plot_data <- rbind(
  extract_conditional_variance(matern32, "32"),
  extract_conditional_variance(matern72, "72")
)
plot_data <- plot_data[
  plot_data$iteration >= 40L & plot_data$iteration <= 250L,
  ,
  drop = FALSE
]
plot_data$series <- factor(plot_data$series, levels = series_order)

plot <- ggplot2::ggplot(
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
  ggplot2::labs(x = NULL, y = NULL, color = "Proposal", linetype = "Proposal") +
  theme_bq() +
  ggplot2::theme(
    legend.position = c(0.05, 0.05),
    legend.justification = c(0, 0)
  )

save_bq_plot(
  plot,
  file.path(
    BQ_UD_ROOT,
    "figures",
    "bq_student_matern32_matern72_conditional_variance_avg.pdf"
  ),
  width = 8,
  height = 5,
  dpi = 150
)

message("Wrote combined Student-t variance plots to ", file.path(
  BQ_UD_ROOT,
  "figures"
))
invisible(plot)
