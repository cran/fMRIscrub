
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fMRIscrub

<!-- badges: start -->

[![R-CMD-check](https://github.com/mandymejia/fMRIscrub/workflows/R-CMD-check/badge.svg)](https://github.com/mandymejia/fMRIscrub/actions)
[![Codecov test
coverage](https://codecov.io/gh/mandymejia/fMRIscrub/branch/master/graph/badge.svg)](https://app.codecov.io/gh/mandymejia/fMRIscrub?branch=master)
<!-- badges: end -->

`fMRIscrub` is a collection of routines for data-driven scrubbing
(projection scrubbing and DVARS), motion scrubbing, and other fMRI
denoising strategies such as anatomical CompCor, detrending, and
nuisance regression. Projection scrubbing is also applicable to other
outlier detection tasks involving high-dimensional data.

## Installation

You can install the development version of fMRIscrub from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mandymejia/fMRIscrub")
```

## Quick start guide

``` r
s_Dat1 <- scrub(Dat1)
plot(s_Dat1)
Dat1_cleaned <- Dat1[!s_Dat1$outlier_flag,]
```

## Data

Two scans from the [ABIDE
I](http://fcon_1000.projects.nitrc.org/indi/abide/abide_I.html) are
included in `fMRIscrub`: `Dat1` has many artifacts whereas `Dat2` has
few visible artifacts. Both are vectorized sagittal slices stored as
numeric matrices. They are loaded into the environment upon loading the
package.

We acknowledge the corresponding funding for the ABIDE I data:

> Primary support for the work by Adriana Di Martino was provided by the
> (NIMH K23MH087770) and the Leon Levy Foundation. Primary support for
> the work by Michael P. Milham and the INDI team was provided by gifts
> from Joseph P. Healy and the Stavros Niarchos Foundation to the Child
> Mind Institute, as well as by an NIMH award to MPM ( NIMH
> R03MH096321).

## Vignette

See [this
link](https://github.com/mandymejia/fMRIscrub/blob/master/vignettes/projection_scrubbing.rmd)
to view the tutorial vignette.

## Citation

If using projection scrubbing, you can cite our pre-print at
<https://arxiv.org/abs/2108.00319>.
