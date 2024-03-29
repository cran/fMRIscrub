% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{scrub_plot}
\alias{scrub_plot}
\title{\code{"scrub"} plot sub-function}
\usage{
scrub_plot(
  meas,
  cut = NULL,
  flag_intersect = FALSE,
  colors = NULL,
  log_y = FALSE,
  geom = "point",
  ylim_min = 0,
  ylim_max = max(meas$measure),
  ...
)
}
\arguments{
\item{meas}{A \eqn{T} by \code{m} numeric data.frame with each column being the timecourse for an
outlyingness measure. The names of the columns will be used to label the plot.}

\item{cut}{A length \eqn{m} numeric vector with each value being the cutoff for an
outlyingness measure (each column in \code{meas}).}

\item{flag_intersect}{Only flag timepoints at which all measures are outliers?
Default: \code{FALSE}.}

\item{colors}{A length \eqn{m} character vector giving the colors of each
measure (each column in \code{meas})}

\item{log_y}{Use log scale for y-axis? Default: \code{FALSE}}

\item{geom}{"point" (default) or "line"}

\item{ylim_min, ylim_max}{The range of the y-axis.}

\item{...}{Additional arguments to ggplot: main, sub, xlab, ...}
}
\value{
A ggplot
}
\description{
Plot outlyingness measure(s) with the corresponding threshold(s). Requires
the \code{cowplot} and \code{ggplot2} packages
}
\keyword{internal}
