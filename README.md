
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpatialGEV

<!-- badges: start -->
<!-- badges: end -->

*Meixi Chen, Martin Lysy*

------------------------------------------------------------------------

## Description

A fast Bayesian inference method for spatial random effects modelling of
weather extremes. The latent spatial variables are efficiently
marginalized by a Laplace approximation using the
[***TMB***](https://github.com/kaskr/adcomp) library, which leverages
efficient automatic differentiation in C++. The models are compiled in
C++, whereas the optimization step is carried out in R. With this
package, users can fit spatial GEV models with different complexities to
their dataset without having to formulate the model using C++. This
package also offers a method to sample from the approximate posterior
distributions of both fixed and random effects, which will be useful for
downstream analysis.

## Installation

Before installing ***SpatialGEV***, make sure you have ***TMB***
installed following the instructions
[here](https://github.com/kaskr/adcomp/wiki/Download).

***SpatialGEV*** uses several functions from the ***INLA*** package for
SPDE approximation to the Matérn covariance as well as mesh creation on
the spatial domain. If the user would like to use the SPDE method
(i.e. `kernel="spde"` in `spatialGEV_fit()`), please first install
package ***INLA***. Since ***INLA*** is not on CRAN, it needs to be
downloaded following their instruction
[here](https://www.r-inla.org/download-install).

To download the stable version of this package, run

``` r
install.packages("SpatialGEV")
```

To download the development version of this package, run

``` r
devtools::install_github("meixichen/SpatialGEV")
```

## Example

Using the simulated data set `simulatedData2` provided in the package,
we demonstrate how to use this package. Spatial variation of the GEV
parameters are plotted below.

``` r
library(SpatialGEV)
# GEV parameters simulated from Gaussian random fields
a <- simulatedData2$a # location
logb <- simulatedData2$logb # log scale
logs <- simulatedData2$logs # log shape
locs <- simulatedData2$locs # coordinate matrix
n_loc <- nrow(locs) # number of locations
y <- Map(evd::rgev, n=sample(50:70, n_loc, replace=TRUE),
         loc=a, scale=exp(logb), shape=exp(logs)) # observations

filled.contour(unique(locs$x), unique(locs$y), matrix(a, ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", 
               main="Spatial variation of a",
               cex.lab=1,cex.axis=1)
```

<img src="man/figures/README-show-data-1.png" width="50%" />

``` r
filled.contour(unique(locs$x), unique(locs$y), matrix(exp(logb), ncol=sqrt(n_loc)), 
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude", 
               main="Spatial variation of b",
               cex.lab=1,cex.axis=1)
```

<img src="man/figures/README-show-data-2.png" width="50%" />

``` r
filled.contour(unique(locs$x), unique(locs$y), matrix(exp(logs), ncol=sqrt(n_loc)),
               color.palette = terrain.colors, xlab="Longitude", ylab="Latitude",
               main="Spatial variation of s",
               cex.lab=1,cex.axis=1)
```

<img src="man/figures/README-show-data-3.png" width="50%" />

To fit a GEV-GP model to the simulated data, use the `spatialGEV_fit()`
function. We use `random="abs"` to indicate that all three GEV
parameters are treated as random effects. The shape parameter `s` is
constrained to be positive (log transformed) by specifying
`reparam_s="positive"`. The covariance kernel function used here is the
exponential kernel `kernel="exp"`. Initial parameter values are passed
to `init_param` using a list.

``` r
fit <- spatialGEV_fit(y = y, locs = locs, random = "abs",
                      init_param = list(a = rep(60, n_loc),
                                        log_b = rep(2,n_loc),
                                        s = rep(-3,n_loc),
                                        beta_a = 60, beta_b = 2, beta_s = -2,
                                        log_sigma_a = 1.5, log_ell_a = 5,
                                        log_sigma_b = 1.5, log_ell_b = 5,
                                        log_sigma_s = -1, log_ell_s = 5),
                      reparam_s = "positive", kernel="exp", silent = T)

class(fit)
#> [1] "spatialGEVfit"
print(fit)
#> Model fitting took 165.758660793304 seconds 
#> The model has reached relative convergence 
#> The model uses a exp kernel 
#> Number of fixed effects in the model is 9 
#> Number of random effects in the model is 1200 
#> Hessian matrix is positive definite. Use spatialGEV_sample to obtain posterior samples
```

Posterior samples of the random and fixed effects are drawn using
`spatialGEV_sample()`. Specify `observation=TRUE` if we would also like
to draw from the posterior predictive distribution.

``` r
sam <- spatialGEV_sample(model = fit, n_draw = 2000, observation = T)
print(sam)
#> The samples contains 2000 draws of 1209 parameters 
#> The samples contains 2000 draws of response at 400 locations 
#> Use summary() to obtain summary statistics of the samples
```

To get summary statistics of the posterior samples, use `summary()` on
the sample object.

``` r
pos_summary <- summary(sam)
pos_summary$param_summary[1:5,]
#>        2.5%      25%      50%      75%    97.5%     mean
#> a1 59.83426 61.20700 61.92395 62.65845 64.08872 61.93718
#> a2 59.84455 61.14915 61.84187 62.49733 63.78515 61.84007
#> a3 60.00949 61.24158 61.90882 62.54333 63.83870 61.89577
#> a4 59.53783 60.77994 61.41770 62.06576 63.18262 61.40233
#> a5 58.35922 59.80460 60.49968 61.20136 62.62408 60.49384
pos_summary$y_summary[1:5,]
#>        2.5%      25%      50%      75%    97.5%     mean
#> y1 38.07496 57.15007 70.49057 87.70135 139.8428 75.25941
#> y2 36.46683 55.46055 69.76261 87.76045 142.6813 74.86213
#> y3 37.60051 55.37956 68.99792 86.42350 143.9566 73.99941
#> y4 37.21712 55.71365 69.11869 88.02201 147.2459 75.22608
#> y5 36.21457 54.40182 68.31368 85.91682 143.4205 73.37834
```
