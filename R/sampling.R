sample_point <- function(
    strategy,
    n,
    measure = "gaussian",
    df = 8,
    exponent = 3 / 20,
    matern_smoothness = NULL,
    dimension = 1) {
  if (measure == "gaussian") {
    if (strategy == "standard") {
      return(rnorm(1, mean = 0, sd = 1))
    }
    if (strategy == "optimal") {
      return(rnorm(1, mean = 0, sd = sqrt(log(n))))
    }
  }

  if (measure == "student") {
    if (strategy == "standard") {
      return(rt(1, df = df))
    }
    if (strategy == "optimal") {
      if (!is.null(matern_smoothness)) {
        exp_opt <- matern_smoothness / (matern_smoothness + df + dimension / 2)
        return((n^exp_opt) * rt(1, df = df))
      }
      return((n^exponent) * rt(1, df = df))
    }
  }

  stop("Unsupported sampling strategy/measure combination.")
}
