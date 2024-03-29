% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/robdist.R
\name{RD_univOut}
\alias{RD_univOut}
\title{Univariate outlier detection for robust distance}
\usage{
RD_univOut(data, cutoff = 4, trans = c("none", "robust-YJ", "SHASH"))
}
\arguments{
\item{data}{The data}

\item{cutoff}{Default: \code{4}.}

\item{trans}{Transform the data? Default: \code{"none"}. The other choice is
\code{"robust-YJ"}. The \code{"SHASH"} method has not been implemented yet.}
}
\value{
The univariate outliers.
}
\description{
Identify the univariate outliers with robust distance.
}
\keyword{internal}
