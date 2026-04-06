#' @importFrom stats dchisq pchisq rchisq uniroot
#' @importFrom graphics points segments
NULL

#' Validate and normalize mixture inputs
#'
#' Internal utility used by all exported functions.
#'
#' @param df Numeric vector of degrees of freedom. The value 0 is allowed and
#'   represents a point mass at zero.
#' @param w Numeric vector of mixture weights.
#'
#' @return A list containing normalized weights and validated degrees of freedom.
#' @keywords internal
.validate_mixchisq_inputs <- function(df, w) {
  if (missing(df) || missing(w)) {
    stop("Both 'df' and 'w' must be supplied.", call. = FALSE)
  }
  if (!is.numeric(df) || !is.numeric(w)) {
    stop("'df' and 'w' must both be numeric vectors.", call. = FALSE)
  }
  if (length(df) != length(w)) {
    stop("'df' and 'w' must have the same length.", call. = FALSE)
  }
  if (length(df) == 0L) {
    stop("'df' and 'w' must be non-empty.", call. = FALSE)
  }
  if (any(!is.finite(df)) || any(df < 0)) {
    stop("All entries of 'df' must be finite and nonnegative.", call. = FALSE)
  }
  if (any(!is.finite(w)) || any(w < 0)) {
    stop("All entries of 'w' must be finite and nonnegative.", call. = FALSE)
  }
  sw <- sum(w)
  if (sw <= 0) {
    stop("The weights must sum to a positive value.", call. = FALSE)
  }
  list(df = df, w = w / sw)
}

#' Density for a finite mixture of chi-square distributions
#'
#' Computes the density of a finite mixture of chi-square distributions. The
#' special case `df = 0` is interpreted as a point mass at zero, so the returned
#' value describes only the continuous density. At `x = 0`, the function returns
#' `Inf` when there is atom at zero and `log = FALSE`, and `Inf` on the log scale
#' when `log = TRUE`.
#'
#' @param x Vector of quantiles.
#' @param df Numeric vector of degrees of freedom. The value `0` is allowed and
#'   is treated as a point mass at zero.
#' @param w Numeric vector of mixture weights.
#' @param log Logical; if `TRUE`, returns the log-density.
#'
#' @return Numeric vector of density values.
#' @export
#'
#' @examples
#' dmixchisq(c(0, 1, 2), df = c(0, 1, 3), w = c(0.5, 0.25, 0.25))
dmixchisq <- function(x, df, w, log = FALSE) {
  pars <- .validate_mixchisq_inputs(df, w)
  df <- pars$df
  w <- pars$w
  atom0 <- sum(w[df == 0])
  pos_df <- df[df > 0]
  pos_w <- w[df > 0]

  out <- vapply(x, function(xx) {
    if (xx < 0 || is.na(xx)) return(NaN)
    if (xx == 0 && atom0 > 0) {
      if (log) return(Inf)
      return(Inf)
    }
    dens <- if (length(pos_df)) sum(pos_w * dchisq(xx, df = pos_df)) else 0
    if (log) log(dens) else dens
  }, numeric(1))

  out
}

#' CDF for a finite mixture of chi-square distributions
#'
#' Computes the cumulative distribution function for a finite mixture of
#' chi-square distributions. Components with `df = 0` contribute a point mass at
#' zero.
#'
#' @param q Vector of quantiles.
#' @param df Numeric vector of degrees of freedom. The value `0` is allowed and
#'   is treated as a point mass at zero.
#' @param w Numeric vector of mixture weights.
#' @param lower.tail Logical; if `TRUE`, probabilities are `P(X \le q)`,
#'   otherwise `P(X > q)`.
#'
#' @return Numeric vector of probabilities.
#' @export
#'
#' @examples
#' pmixchisq(c(0, 1, 2), df = c(0, 1, 3), w = c(0.5, 0.25, 0.25))
pmixchisq <- function(q, df, w, lower.tail = TRUE) {
  pars <- .validate_mixchisq_inputs(df, w)
  df <- pars$df
  w <- pars$w
  atom0 <- sum(w[df == 0])
  pos_df <- df[df > 0]
  pos_w <- w[df > 0]

  cdf <- vapply(q, function(qq) {
    if (is.na(qq)) return(NaN)
    if (qq < 0) return(0)
    if (qq == 0) return(atom0)
    atom0 + if (length(pos_df)) sum(pos_w * pchisq(qq, df = pos_df)) else 0
  }, numeric(1))

  if (lower.tail) cdf else 1 - cdf
}

#' Quantile function for a finite mixture of chi-square distributions
#'
#' Computes quantiles numerically by inverting the mixture CDF with `uniroot()`.
#' If the requested probability is less than or equal to the atom at zero, the
#' quantile is returned as zero.
#'
#' @param p Vector of probabilities in `[0, 1]`.
#' @param df Numeric vector of degrees of freedom. The value `0` is allowed and
#'   is treated as a point mass at zero.
#' @param w Numeric vector of mixture weights.
#' @param lower.tail Logical; if `TRUE`, probabilities are `P(X \le x)`,
#'   otherwise `P(X > x)`.
#' @param interval Initial search interval for the numerical inversion.
#' @param tol Numerical tolerance passed to `uniroot()`.
#'
#' @return Numeric vector of quantiles.
#' @export
#'
#' @examples
#' qmixchisq(c(0.5, 0.9, 0.95), df = c(0, 1, 3), w = c(0.5, 0.25, 0.25))
qmixchisq <- function(p, df, w, lower.tail = TRUE,
                      interval = c(0, 1e4), tol = 1e-10) {
  pars <- .validate_mixchisq_inputs(df, w)
  df <- pars$df
  w <- pars$w
  atom0 <- sum(w[df == 0])

  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("All probabilities in 'p' must lie in [0, 1].", call. = FALSE)
  }
  if (!lower.tail) p <- 1 - p

  out <- vapply(p, function(pp) {
    if (pp == 0) return(0)
    if (pp <= atom0) return(0)
    if (pp == 1) return(Inf)

    f <- function(x) pmixchisq(x, df = df, w = w, lower.tail = TRUE) - pp

    upper <- interval[2]
    while (f(upper) < 0) {
      upper <- upper * 2
      if (upper > 1e12) {
        stop("Could not bracket the requested quantile.", call. = FALSE)
      }
    }

    uniroot(f, interval = c(interval[1], upper), tol = tol)$root
  }, numeric(1))

  out
}

#' Random generation from a finite mixture of chi-square distributions
#'
#' Generates random values from a finite mixture of chi-square distributions by
#' first sampling a component label and then sampling from the corresponding
#' component. Components with `df = 0` generate exact zeros.
#'
#' @param n Number of draws.
#' @param df Numeric vector of degrees of freedom. The value `0` is allowed and
#'   is treated as a point mass at zero.
#' @param w Numeric vector of mixture weights.
#'
#' @return Numeric vector of length `n`.
#' @export
#'
#' @examples
#' set.seed(1)
#' rmixchisq(5, df = c(0, 1, 3), w = c(0.5, 0.25, 0.25))
rmixchisq <- function(n, df, w) {
  pars <- .validate_mixchisq_inputs(df, w)
  df <- pars$df
  w <- pars$w

  if (length(n) != 1L || !is.finite(n) || n < 0 || n != as.integer(n)) {
    stop("'n' must be a single nonnegative integer.", call. = FALSE)
  }

  z <- sample.int(length(w), size = n, replace = TRUE, prob = w)
  out <- numeric(n)
  zero_idx <- df[z] == 0
  out[zero_idx] <- 0
  out[!zero_idx] <- rchisq(sum(!zero_idx), df = df[z][!zero_idx])
  out
}

#' Plot the density or CDF of a chi-square mixture
#'
#' Convenience plotting wrapper for quick inspection of a mixture reference law.
#'
#' @param df Numeric vector of degrees of freedom. The value `0` is allowed and
#'   is treated as a point mass at zero.
#' @param w Numeric vector of mixture weights.
#' @param xlim Numeric vector of length two giving the plotting range.
#' @param n_grid Number of grid points.
#' @param type Either `"cdf"` or `"density"`.
#' @param ... Additional graphical arguments passed to `plot()`.
#'
#' @return Invisibly returns a data frame with the plotted grid.
#' @export
#'
#' @examples
#' plot_mixchisq(df = c(2, 3), w = c(0.5, 0.5))
plot_mixchisq <- function(df, w, xlim = c(0, 10), n_grid = 500,
                          type = c("cdf", "density"), ...) {
  type <- match.arg(type)
  x <- seq(xlim[1], xlim[2], length.out = n_grid)
  y <- switch(type,
    cdf = pmixchisq(x, df = df, w = w),
    density = dmixchisq(x, df = df, w = w)
  )
  y_plot <- y
  y_plot[!is.finite(y_plot)] <- NA_real_
  plot(x, y_plot, type = "l", xlab = "x", ylab = type, ...)
  if (type == "cdf" && any(df == 0)) {
    atom0 <- sum(.validate_mixchisq_inputs(df, w)$w[df == 0])
    points(0, atom0, pch = 19)
    segments(0, 0, 0, atom0, lty = 2)
  }
  invisible(data.frame(x = x, y = y))
}
