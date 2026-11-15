# Shared orchestration for the retained BQ experiments.

run_bq_comparison <- function(
    experiment_id,
    n,
    integrand,
    kernel_fun,
    kernel_mean_fun,
    prior_var,
    ell,
    sigma_f2,
    sample_target,
    sample_power,
    target_quantile,
    alpha,
    dimension,
    measure_name,
    scheme_labels,
    settings = list(),
    n_rep = 100L,
    replication_seed = 5000L,
    single_seed = 1L,
    ci_draws = 100L,
    mean_plot_start = 1L,
    mean_plot_end = n,
    plots = c("mean_bands", "mean_ci_avg"),
    bach_grid_n = 801L,
    bach_lambda = NULL,
    results_dir = "results",
    figures_dir = "figures") {
  scheme_ids <- bq_scheme_order
  if (!identical(sort(names(scheme_labels)), sort(scheme_ids))) {
    stop("scheme_labels must be named with: power, bach, target.")
  }

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  bach_distribution <- build_bach_distribution(
    kernel_fun = kernel_fun,
    ell = ell,
    sigma_f2 = sigma_f2,
    quantile_fun = target_quantile,
    n_total = n,
    alpha = alpha,
    dimension = dimension,
    measure_name = measure_name,
    n_grid = bach_grid_n,
    lambda = bach_lambda
  )
  sample_bach_scheme <- function(iteration, n_total) {
    sample_bach(bach_distribution, 1L)
  }

  samplers <- list(
    power = sample_power,
    bach = sample_bach_scheme,
    target = sample_target
  )
  run_functions <- lapply(samplers, function(sampler) {
    force(sampler)
    function() {
      run_bq(
        n = n,
        integrand = integrand,
        sample_fun = sampler,
        kernel_fun = kernel_fun,
        kernel_mean_fun = kernel_mean_fun,
        prior_var = prior_var,
        ell = ell,
        sigma_f2 = sigma_f2
      )
    }
  })

  single_results <- setNames(vector("list", length(scheme_ids)), scheme_ids)
  for (scheme_index in seq_along(scheme_ids)) {
    scheme_id <- scheme_ids[scheme_index]
    set.seed(single_seed + scheme_index - 1L)
    single_results[[scheme_id]] <- run_functions[[scheme_id]]()
  }

  replicated_raw <- replicate_bq_schemes(
    n_rep = n_rep,
    n = n,
    seed_start = replication_seed,
    run_funs = run_functions
  )
  summaries <- lapply(replicated_raw, function(replicated_scheme) {
    summarize_repeated_runs(
      replicated_scheme$mean_mat,
      replicated_scheme$var_mat,
      s_draws = ci_draws
    )
  })

  standard_settings <- list(
    experiment_id = experiment_id,
    n = n,
    ell = ell,
    sigma_f2 = sigma_f2,
    sobolev_alpha = alpha,
    dimension = dimension,
    measure = measure_name,
    n_rep = n_rep,
    replication_seed = replication_seed,
    single_seed = single_seed,
    ci_draws = ci_draws,
    scheme_order = scheme_ids,
    scheme_labels = scheme_labels
  )
  output <- c(
    list(
      settings = utils::modifyList(settings, standard_settings),
      bach_metadata = bach_metadata(bach_distribution)
    ),
    single_results,
    list(replicated = summaries)
  )

  # Save the numerical results before plotting so a graphics-device failure
  # does not discard a completed experiment.
  saveRDS(
    output,
    file = file.path(results_dir, paste0(experiment_id, ".rds"))
  )

  if ("mean_bands" %in% plots) {
    plot_bq_posterior_mean_bands_multi(
      results = single_results,
      labels = scheme_labels,
      true_value = settings$true_integral,
      iter_start = mean_plot_start,
      n_max = mean_plot_end,
      outfile = file.path(
        figures_dir,
        paste0(experiment_id, "_mean_bands.pdf")
      )
    )
  }
  if ("mean_ci_avg" %in% plots) {
    plot_bq_averaged_mean_ci_multi(
      summaries = summaries,
      labels = scheme_labels,
      true_value = settings$true_integral,
      iter_start = mean_plot_start,
      n_max = mean_plot_end,
      outfile = file.path(
        figures_dir,
        paste0(experiment_id, "_mean_ci_avg.pdf")
      )
    )
  }
  final_summary <- data.frame(
    scheme = scheme_ids,
    label = unname(scheme_labels[scheme_ids]),
    single_absolute_error = vapply(
      single_results,
      function(x) abs(tail(x$mu, 1) - settings$true_integral),
      numeric(1)
    ),
    single_variance = vapply(single_results, function(x) tail(x$var, 1), numeric(1)),
    mean_absolute_error = vapply(
      replicated_raw[scheme_ids],
      function(x) {
        mean(abs(x$mean_mat[, n] - settings$true_integral))
      },
      numeric(1)
    ),
    averaged_mean_absolute_error = vapply(
      replicated_raw[scheme_ids],
      function(x) {
        abs(mean(x$mean_mat[, n]) - settings$true_integral)
      },
      numeric(1)
    ),
    averaged_variance = vapply(
      summaries[scheme_ids],
      function(x) tail(x$var, 1),
      numeric(1)
    ),
    averaged_conditional_variance = vapply(
      summaries[scheme_ids],
      function(x) tail(x$conditional_var, 1),
      numeric(1)
    ),
    design_variance = vapply(
      summaries[scheme_ids],
      function(x) tail(x$design_var, 1),
      numeric(1)
    ),
    row.names = NULL
  )
  output$final_summary <- final_summary
  saveRDS(
    output,
    file = file.path(results_dir, paste0(experiment_id, ".rds"))
  )
  print(final_summary)
  invisible(output)
}
