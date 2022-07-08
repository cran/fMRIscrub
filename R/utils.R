#' Is this an integer?
#' 
#' @param x The putative integer
#' @param nneg Require \code{x>=0} (non-negative) too?
#' @return Logical indicating whether \code{x} is an integer
#' 
#' @keywords internal
is_integer <- function(x, nneg=FALSE){
  y <- FALSE
  if (length(x)==1 && is.numeric(x)) {
    if (x%%1==0) {
      if (x>=0 || !nneg) { y <- TRUE }
    }
  } 
  y
}

#' Is this numeric vector constant?
#' 
#' @param x The numeric vector
#' @param TOL minimum range of \code{x} to be considered non-constant.
#'  Default: \code{1e-8}
#' 
#' @return Is \code{x} constant? 
#' 
#' @keywords internal
is_constant <- function(x, TOL=1e-8) {
  abs(max(x) - min(x)) < TOL
}

#' Check design matrix
#' 
#' @param design The design matrix
#' 
#' @return The (modified) design matrix
#' 
#' @keywords internal
check_design_matrix <- function(design, T_=nrow(design)) {
  class(design) <- "numeric"
  if (identical(design, 1)) { design <- matrix(1, nrow=T_) }
  design <- as.matrix(design)
  stopifnot(nrow(design) == T_)
  # Set constant columns (intercept regressor) to 1, and scale the other columns.
  design_const_mask <- apply(design, 2, is_constant)
  design[,design_const_mask] <- 1
  design[,!design_const_mask] <- scale(design[,!design_const_mask])
  design
}

#' Robust scaling
#' 
#' Centers and scales the columns of a matrix robustly
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). Constant-valued columns are set to \code{NA}
#' (or removed if \code{drop_const}) and a warning is raised. If all 
#' MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#' @param TOL minimum MAD to consider a column non-constant.
#'  Default: \code{1e-8}
#' @param drop_const Drop
#'
#' @return The input matrix with its columns centered and scaled.
#'
#' @importFrom robustbase rowMedians
scale_med <- function(mat, TOL=1e-8, drop_const=TRUE){
  # Transpose.
  mat <- t(mat)

  #	Center.
  mat <- mat - c(rowMedians(mat, na.rm=TRUE))

  # Scale.
  mad <- 1.4826 * rowMedians(abs(mat), na.rm=TRUE)
  mad <- as.numeric(mad)
  const_mask <- mad < TOL
  if (any(const_mask)) {
    if (all(const_mask)) {
    stop("All columns are zero-variance.\n")
    } else {
      warning(paste0(
        "Warning: ", sum(const_mask),
        " constant columns (out of ", length(const_mask),
        " ). These will be removed.\n"
      ))
    }
  }
  mad <- mad[!const_mask]
  mat[const_mask,] <- NA
  mat[!const_mask,] <- mat[!const_mask,] / mad

  # Revert transpose.
  mat <- t(mat)

  if (drop_const) { mat <- mat[!const_mask,] }

  mat
}

#' Infer fMRI data format
#'
#' @param BOLD The fMRI data
#' @param verbose Print the format? Default: \code{FALSE}.
#' @return The format: \code{"CIFTI"} file path, \code{"xifti"} object,
#'  \code{"NIFTI"} file path, \code{"nifti"} object, or \code{"data"}.
#' @keywords internal
infer_BOLD_format <- function(BOLD, verbose=FALSE){

  # Character vector: CIFTI or NIFTI
  if (is.character(BOLD)) {
    stopifnot(length(BOLD)==1)
    if (endsWith(BOLD, ".dtseries.nii") | endsWith(BOLD, ".dscalar.nii")) {
      format <- "CIFTI"
    } else if (endsWith(BOLD, "gii")) {
      format <- "GIFTI"
    } else {
      format <- "NIFTI"
    }

  } else if (inherits(BOLD, "xifti")) {
    format <- "xifti"
  } else if (inherits(BOLD, "niftiExtension") | inherits(BOLD, "niftiImage") | inherits(BOLD, "nifti")) {
    format <- "nifti"
  } else if (inherits(BOLD, "gifti")) {
    format <- "gifti"

  } else {
    if (is.list(BOLD)) { stop("Unknown `BOLD` format.") }
    if (is.matrix(BOLD)) { 
      format <- "data"
    } else if (length(dim(BOLD))==4) {
      format <- "nifti"
    } else {
      stop("Unknown `BOLD` format.")
    }
  }
  if (verbose) { cat("Inferred input format:", format, "\n") }
  format
}

#' Wrapper to common functions for reading NIFTIs
#' 
#' Tries \code{RNifti::readNifti}, then \code{oro.nifti::readNIfTI}. If
#'  neither package is available an error is raised.
#' 
#' @param nifti_fname The file name of the NIFTI.
#' @return The NIFTI
#' @keywords internal
read_nifti <- function(nifti_fname){
  if (requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::readNifti(nifti_fname))
  } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
    return(oro.nifti::readNIfTI(nifti_fname, reorient=FALSE))
  } else {
    stop(
      "Package \"RNifti\" or \"oro.nifti\" needed to read", nifti_fname, 
      ". Please install at least one", call. = FALSE
    )
  }
}

#' Convert CIFTI, NIFTI, or GIFTI file input to \eqn{T} by \eqn{V} matrix
#' 
#' Convert input data to \eqn{T} by \eqn{V} matrix. Reads in a CIFTI, NIFTI, 
#'  or GIFTI file. Transposes CIFTI and xifti object input. Masks NIFTI and 
#'  nifti object input.
#' 
#' @param x The object to coerce to a matrix
#' @param sortSub Sort subcortex by labels? Default: \code{FALSE}
#' @param verbose Print updates? Default: \code{FALSE}
#' @return x as a matrix.
#' @keywords internal
as.matrix2 <- function(x, sortSub=FALSE, verbose=FALSE) {

  format <- infer_BOLD_format(x)
  if (format %in% c("CIFTI", "xifti")) {
    if (format == "CIFTI") {
      if (!requireNamespace("ciftiTools", quietly = TRUE)) {
        stop("Package \"ciftiTools\" needed to read input data. Please install it", call. = FALSE)
      }
      x <- ciftiTools::read_cifti(x, brainstructures=ciftiTools::info_cifti(x)$cifti$brainstructures)
    }
    if (sortSub && !is.null(x$data$subcort)) {
      x$data$subcort <- x$data$subcort[order(x$meta$subcort$labels),]
    }
    x <- t(as.matrix(x))
  } else if (format %in% c("GIFTI", "gifti")) {
    if (format == "GIFTI") {
      if (!requireNamespace("gifti", quietly = TRUE)) {
        stop("Package \"gifti\" needed to read `X`. Please install it", call. = FALSE)
      }
      x <- gifti::read_gifti(x)
    }
    x <- t(do.call(cbind, x$data))
    if (verbose) {cat("GIFTI dimensions:\n"); print(dim(x))}
  } else if (format %in% c("NIFTI", "nifti")) {
    if (format == "NIFTI") {
      x <- read_nifti(x)
    }
    if (verbose) {cat("Masking NIFTI by removing locations with constant zero, NA, or NaN.\n")}
    z <- array(x %in% c(0, NA, NaN), dim=dim(x))
    mask <- !apply(z, seq(3), all)
    x <- matrix(x[rep(mask, dim(x)[4])], ncol=dim(x)[4])
    x <- t(x)
    if (verbose) {cat("NIFTI dimensions:\n"); print(dim(x))}
  }

  x
}