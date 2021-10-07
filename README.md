# SpatialGEV

*Meixi Chen, Martin Lysy*

---

### Description

A fast Bayesian inference method for spatial random effects modelling of weather extremes. The latent spatial variables are efficiently marginalized by a Laplace approximation using the [***TMB***](https://github.com/kaskr/adcomp) library, which leverages efficient automatic differentiation in C++. The models are compiled in C++, whereas the optimization step is carried out in R. With this package, users can fit spatial GEV models with different complexities to their dataset without having to formulate the model using C++.  This package also offers a method to sample from the approximate posterior distributions of both fixed and random effects, which will be useful for downstream analysis. 

### Installation

Before installing ***SpatialGEV***, make sure you have ***TMB*** installed following the instructions [here](https://github.com/kaskr/adcomp/wiki/Download).

To download this package in R with a quick tutorial, run 
```r
devtools::install_github("meixichen/SpatialGEV", build_vignettes=TRUE)
```
The tutorial can be viewed by running
```r
vignette("SpatialGEV-tut")
```


### Test

To make sure the package works properly on your machine, run the built-in unit tests in ***SpatialGEV***:
```r
require(SpatialGEV)
testthat::test_package("SpatialGEV")
```
