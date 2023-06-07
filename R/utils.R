#' Get internal name of this projection
#' 
#' Get name of projection given PESEL & kurtosis.
#' 
#' @param projection,PESEL,kurt_quantile See \code{\link{pscrub}}.
#' 
#' @return The name of the projection 
#'
#' @keywords internal
projection_name <- function(projection, PESEL, kurt_quantile){
  if (!PESEL) { projection <- paste0(projection, "2") }
  if (kurt_quantile > 0) { 
    projection <- paste0(projection, "_kurt")
  }
  projection
}