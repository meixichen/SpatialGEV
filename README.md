# SpatialGEV
An R package that builds spatial extreme value models and does Bayesian inference

The models in this package are built using the template model builder [(TMB)](https://github.com/kaskr/adcomp) in R. TMB has the fast ability to integrate out the random effects in a hierarchical model using the Laplace approximation, which is done by leveraging C++ automatic differentiation. The model building step is performed in TMB using C++, whereas the optimization step is carried out in R. However, using TMB requires some familiarity with C++. This package allows the users to choose in the fit function which parameter(s) is considered as following the Gaussian Process, so the users can fit spatial GEV models with different complexities to their dataset without having to formulate the TMB model using C++. This package also offers a method to sample from the approximate posterior distributions of both fixed and random effects, which will be useful for downstream analysis. 

To download this package in R, run `devtools::install_github("meixichen/SpatialGEV")`.
