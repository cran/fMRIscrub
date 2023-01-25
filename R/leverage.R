#' Leverage
#' 
#' Computes the leverage of each observation in the PC score (U) or 
#'  IC mixing (M) matrix for projection scrubbing. Can threshold to 
#'  flag potential outliers.
#'
#' @param Comps The \eqn{n} by \eqn{Q} PC score matrix/IC mixing matrix.
#' @param are_norm Assume the columns of \code{Comps} are orthogonal
#'  and have 2-norms equal to 1? Speeds up the computation.
#' @param median_cutoff The outlier cutoff, in multiples of the median leverage.
#'  Default: \code{NULL} (do not compute outliers).
#' 
#' @return A list with entries \code{"meas"} (the leverage values), 
#'  \code{"cut"} (the leverage cutoff value) and 
#'  \code{"flag"} (logical vector indicating the outliers). If 
#'  \code{!is.null(median_cutoff)}, \code{"cut"} and \code{"flag"} are omitted.
#'  
#' @importFrom stats median
#' @importFrom fMRItools hat_matrix
#' 
#' @export
leverage <- function(Comps, are_norm=FALSE, median_cutoff=NULL){
  lev <- if (are_norm) {
    colSums(Comps^2)
  } else {
    diag(hat_matrix(Comps))
  }

  out <- list(meas=lev)
  if (!is.null(median_cutoff)){
    out$cut <- median_cutoff * median(lev)
    out$flag <- out$meas > out$cut
  }
  out
}