# Locate the BQ_UD root independently of the current R working directory.
# The first candidate works when this file is sourced; the remaining candidates
# cover running it interactively from BQ_UD, its parent, or scripts/.
setup_source_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(error) NULL
)
root_candidates <- c(
  if (!is.null(setup_source_file)) dirname(dirname(setup_source_file)),
  getwd(),
  file.path(getwd(), "BQ_UD"),
  dirname(getwd())
)
root_candidates <- unique(normalizePath(
  root_candidates,
  winslash = "/",
  mustWork = FALSE
))
is_bq_root <- vapply(
  root_candidates,
  function(path) {
    file.exists(file.path(path, "R", "bq_core.R")) &&
      file.exists(file.path(path, "R", "bach_sampling.R")) &&
      file.exists(file.path(path, "R", "experiments.R"))
  },
  logical(1)
)
if (!any(is_bq_root)) {
  stop(
    "Cannot locate the BQ_UD repository root. Open BQ_UD.Rproj or source ",
    "BQ_UD/scripts/_setup.R from the repository or its parent directory."
  )
}

BQ_UD_ROOT <- root_candidates[which(is_bq_root)[1L]]
source(file.path(BQ_UD_ROOT, "R", "kernels_measures.R"))
source(file.path(BQ_UD_ROOT, "R", "hyperparameters.R"))
source(file.path(BQ_UD_ROOT, "R", "bq_core.R"))
source(file.path(BQ_UD_ROOT, "R", "bach_sampling.R"))
source(file.path(BQ_UD_ROOT, "R", "postprocess.R"))
source(file.path(BQ_UD_ROOT, "R", "experiments.R"))

# The experiment scripts use project-relative input and output paths.
if (!identical(
    normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    BQ_UD_ROOT
)) {
  setwd(BQ_UD_ROOT)
}

dir.create(file.path(BQ_UD_ROOT, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BQ_UD_ROOT, "figures"), showWarnings = FALSE, recursive = TRUE)

env_integer <- function(name, default, minimum = 1L) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < minimum) {
    stop(sprintf("%s must be an integer >= %d.", name, minimum))
  }
  parsed
}

