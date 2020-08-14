# SpatialGEV
An R package that builds spatial extreme value model and does Bayesian inference

The models in this package are built using the template model builder [(TMB)](https://github.com/kaskr/adcomp) in R, which has the fast ability to integrate out random effects using Laplace approximation. This package allows the users to choose in the fit function which parameter(s) is considered as following the GP, so the users can fit spatial GEV models with different complexities to their dataset without having to write the models in TMB by themselves. This package also offers a method to sample from both fixed and random effects posterior.
