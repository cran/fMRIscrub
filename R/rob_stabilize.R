#' Robust linear model on DCT bases
#'
#' Fit a linear model regressing an input vector on DCT bases, robustly.
#'
#' @param x The input vector to regress on DCT bases
#' @param nDCT The number of DCT bases to use. Default: \code{4}
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' @param seed Set a seed right before the call to \code{robustbase::lmrob}?
#'  If \code{NULL}, do not set a seed. If numeric (default: \code{0}), will use
#'  as the seed.
#'
#' @return The output of \code{robustbase::lmrob}
#'
#' @importFrom fMRItools is_integer dct_bases
#' @keywords internal
rob_trend <- function(x, nDCT=4, lmrob_method="MM", seed=0) {
  x <- as.vector(x)
  T_ <- length(x)

  nDCT <- as.numeric(nDCT)
  stopifnot(fMRItools::is_integer(nDCT, nneg=TRUE))
  if (nDCT == 0) {
    mat <- data.frame(rep(1, T_))
    colnames(mat) <- "x_int"
  } else {
    mat <- data.frame(cbind(1, dct_bases(T_, nDCT)))
    colnames(mat) <- c("x_int", paste0("x_dct", seq(nDCT)))
  }
  mat$y <- x

  if (!is.null(seed)) {
    with(
      set.seed(seed),
      robustbase::lmrob(y~., mat, method=lmrob_method, setting="KS2014")
    )
  } else {
    robustbase::lmrob(y~., mat, method=lmrob_method, setting="KS2014")
  }
}

#' Stabilize the center and scale of a timeseries robustly 
#' 
#' Stabilize the center and scale of a timeseries using robust regression of
#'  DCT bases on the first and second moments.
#' 
#' @param x The timeseries to stabilize. 
#' @param center,scale Center and scale? Default: \code{TRUE} for both. If
#'  scaling but not centering, the data must already be centered; otherwise,
#'  the results will be invalid. Can also be the number of DCT bases to use for
#'  robust stabilization of center/scale; \code{TRUE} will use \code{4}.
#' @param lmrob_method The \code{lmrob_method} argument to \code{robustbase::lmrob}.
#' @param rescale After stabilizing \code{x}, re-center and re-scale
#'  to the original mean and variance? Default: \code{TRUE}.
#'
#' @return the timeseries with its center and scale stabilized
#' 
#' @importFrom fMRItools is_integer
#' @export
#' 
rob_stabilize <- function(x, center=TRUE, scale=TRUE, lmrob_method="MM", rescale=TRUE) {
  EPS <- 1e-6

  if (length(x) < 5) { warning("Timeseries too short to variance stabilize."); return(x) }
  if (any(is.na(x) | is.nan(x))) { stop("NA/NaN values in `x` are not supported.") }
  x_mean <- mean(x); x_var <- var(x)
  x <- as.numeric(scale(x))

  if (isTRUE(center)) {
    center <- 4
  } else if (isFALSE(center)) { 
    center <- 0
  } else {
    center <- as.numeric(center)
    stopifnot(fMRItools::is_integer(center, nneg=TRUE))
  }

  if (isTRUE(scale)) {
    scale <- 4
  } else if (isFALSE(scale)) { 
    scale <- 0
  } else {
    scale <- as.numeric(scale)
    stopifnot(fMRItools::is_integer(scale, nneg=TRUE))
  }

  if (center > 0) {
    m <- as.numeric(rob_trend(x, nDCT=center, lmrob_method)$fitted.values)
    x <- scale(x - m)
  }

  const_mask <- abs(x) < EPS

  if (scale > 0) {
    x2 <- ifelse(const_mask, NA, x)
    s <- as.numeric(rob_trend(log(x2^2), nDCT=scale, lmrob_method)$fitted.values)
    s <- sqrt(exp(s))
    if (any(s < EPS)) { stop("Error: near-constant variance detected.") } # TEMPORARY
    x[!const_mask] <- x[!const_mask] / s
    x <- scale(x)
  }

  if (rescale) { x <- (x * sqrt(x_var)) + x_mean }
  x
}