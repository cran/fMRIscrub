% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/high_kurtosis.R
\name{high_kurtosis}
\alias{high_kurtosis}
\title{Which components have high kurtosis?}
\usage{
high_kurtosis(Comps, kurt_quantile = 0.99, n_sim = 5000, min_1 = FALSE)
}
\arguments{
\item{Comps}{A matrix; each column is a component. For PCA, this is the U
matrix. For ICA, this is the M matrix.}

\item{kurt_quantile}{components with kurtosis of at least this quantile are kept.}

\item{n_sim}{The number of simulation data to use for estimating the sampling
distribution of kurtosis. Only used if a new simulation is performed. (If
\eqn{n<1000} and the quantile is 90\%, a pre-computed value is used instead.
If \eqn{n>1000}, the theoretical asymptotic distribution is used instead.)}

\item{min_1}{Require at least one component to be selected? In other words, if
no components meet the quantile cutoff, should the component with the highest
kurtosis be returned? Default: \code{FALSE}.}
}
\value{
A logical vector indicating whether each component has high kurtosis.
}
\description{
The kurtosis cutoff is a high quantile (default 0.99) of the sampling distribution
of kurtosis for Normal iid data of the same length as the components; it is
estimated by simulation or calculated from the theoretical asymptotic
distribution if the components are long enough.
}
\details{
The components should not have any strong low-frequency trends, because trends
can affect kurtosis in unpredictable ways unrelated to outlier presence.
}
