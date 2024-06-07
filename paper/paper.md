---
title: 'SpatialGEV: Fast Bayesian inference for spatial extreme value models in \textsf{R}'
tags:
- R
- Bayesian inference
- spatial model
- extreme value
date: "\\today"
output: pdf_document
header-includes:
  - \usepackage{subcaption}
authors:
- name: Meixi Chen^[Corresponding author]
  orcid: 0000-0003-1012-5352
  affiliation: 1
- name: Martin Lysy
  orcid: 0000-0001-9974-1121
  affiliation: 1
- name: Reza Ramezan
  orcid: 0000-0003-0450-3249
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Department of Statistics and Actuarial Science, University of Waterloo
  index: 1
---

<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: { equationNumbers: { autoNumber: "AMS" } }
});
</script>

\newcommand{\bm}[1]{\boldsymbol{#1}}
\newcommand{\tx}[1]{\mathrm{#1}}
\newcommand{\xx}{{\bm{x}}}
\newcommand{\yy}{{\bm{y}}}
\newcommand{\XX}{{\bm{X}}}
\newcommand{\YY}{{\bm{Y}}}
\newcommand{\ZZ}{{\bm{Z}}}
\newcommand{\tth}{{\bm{\theta}}}
\newcommand{\pps}{{\bm{\psi}}}
\newcommand{\uu}{{\bm{u}}}
\newcommand{\SSi}{{\bm{\Sigma}}}
\newcommand{\VV}{{\bm{V}}}
\newcommand{\iid}{{\overset{\mathrm{iid}}{\sim}}}

# Summary

Extreme weather phenomena such as floods and hurricanes are of great concern due to their potential to cause extensive damage. To develop more reliable damage prevention protocols, statistical models are often used to infer the chance of observing an extreme weather event at a given location [@coles98; @cooley07; @sang-gelfand10]. Here we present **SpatialGEV**, an \textsf{R} package providing a fast and convenient toolset for analyzing spatial extreme values using a hierarchical Bayesian modeling framework. In this framework, the marginal behavior of the extremes is given by a generalized extreme value (GEV) distribution, whereas the spatial dependence between locations is captured by modeling the GEV parameters as spatially varying random effects following a Gaussian process (GP). Model inference is carried out using an efficient implementation of the Laplace approximation, which produces highly accurate posterior estimates several orders of magnitude faster than Markov Chain Monte Carlo (MCMC) methods. Users are provided with a streamlined way to build and fit various GEV-GP models in \textsf{R}, which are compiled in \textsf{C++} under the hood. For downstream analyses, the package offers methods for Bayesian parameter estimation and forecasting of extreme events. 

# Background

Let $y_{ij}$ denote observation $j$ of an extreme weather event at spatial location $i$, of which the two-dimensional spatial coordinates are $\xx_i$. The general form of the GEV-GP models that can be fit using **SpatialGEV** is
\begin{equation}
\begin{aligned}
  y_{ij} &\iid \operatorname{GEV}(a(\xx_i), \exp(b(\xx_i)), s(\xx_i)),\\
  u(\xx) &\sim \operatorname{GP}(\bm{c}_u(\xx)^T\bm{\beta}_u(\xx), \ \mathcal{K}(\xx,\xx' \mid \bm{\eta}_u)),
\end{aligned}
\end{equation}
where $\operatorname{GEV}(a, b, s)$ denotes a GEV distribution whose cumulative density function (CDF) is given by
\begin{equation}
F(y \mid a, b, s) = 
\begin{cases}
\exp\left\{-\left(1+s\cdot \frac{y-a}{b}\right)^{-1/s}\right\}, \ &s\neq 0,\\
\exp\left\{-\exp\left(-\frac{y-a}{b}\right)\right\}, \ &s=0, 
\end{cases}
\end{equation}
with location parameter $a$, positive scale parameter $b$ and shape parameter $s$,
$\operatorname{GP}(\mu(\xx), \mathcal{K}(\xx,\xx'))$ is a Gaussian process with mean function $\mu(\xx)$ (possibly depending on covariates $\bm{c}(\xx)$ via coefficients $\bm{\beta}(\xx)$) and covariance kernel $\mathcal{K}(\xx,\xx')$, and $u\in\{a,b,s\}$ is any subset of the GEV parameters which we would like to model as spatially varying.

The GEV-GP model has important applications in meteorological studies. For example, let $y=y(\xx)$ denote the amount of rainfall at a spatial location $\xx$. To forecast extreme rainfalls, it is often of interest for meteorologists to estimate the $1/p$-year rainfall return level $z_p(\xx)$, which is the value above which precipitation levels at location $\xx$ occur with probability $p$, i.e.,
\begin{equation}
\Pr\big(y(\xx) > z_p(\xx)\big) = 1-F\big(z_p(\xx)\mid a(\xx), b(\xx), s(\xx)\big) = p,
\label{eq:return-level}
\end{equation}
where $F\big(z_p(\xx)\mid a(\xx), b(\xx), s(\xx)\big)$ is the CDF of the GEV distribution specific to location $\xx$. When $p$ is chosen to be a small value, $z_p(\xx)$ indicates how extreme the precipitation level might be at location $\xx$.

# Statement of need
In a Bayesian context, the posterior distribution $p(z_p(\xx)\mid \YY)$ conditional on all data $\YY$ is very useful for forecasting extreme weather events. Traditionally, MCMC methods are used to sample from the posterior distribution of the GEV model [e.g., @cooley07; @schliep10; @dyrrdal15]. However, this can be extremely computationally intensive when the number of locations is large. **SpatialGEV** implements an approximate Bayesian inference approach as an alternative to MCMC, making large-scale spatial analyses orders of magnitude faster while achieving roughly the same accuracy as MCMC. We construct a Normal approximation to the joint posterior distribution of both GEV parameters and GP hyperparameters $p(u(\xx), \bm{\eta}_u, \bm{\beta}_u\mid \YY)$, which is then used to estimate the return level posterior. This is done via an efficient Laplace approximation to the marginal hyperparameter posterior $p(\bm{\eta}_u, \bm{\beta}_u \mid \YY)$ transforming a high-dimensional MCMC into a nested optimization problem that is faster to solve [@tierney-kadane86; @kristensen16; @chen-etal24]. The Laplace approximation is carried out using the \textsf{R}/\textsf{C++} package **TMB** [@kristensen16]. Details of the inference method can be found in @chen-etal24. 

The \textsf{R} package **SpatialExtremes** [@spatialextremes] is a popular software for fitting spatial extreme value models including GEV-GP. Although it supports a wider range of model classes, its inference for GEV-GP models relies on a basic Gibbs sampler updating each of the hyperparameters and random effects one at a time, which tends to converge very slowly since these variables are highly correlated with each other. Furthermore, GP computation in **SpatialExtremes** scales as $\mathcal{O}(n^3)$ with the number of locations $n$, whereas **SpatialGEV** offers an option for approximate GP computation scaling as $O(n^{3/2})$ [@lindgren-etal11]. Coupled with the Laplace approximation, this allows **SpatialGEV** to fit GEV-GP models to several hundreds spatial locations on a personal computer in minutes [@chen-etal24]. A more efficient MCMC algorithm for hierarchical spatial models is Hamiltonian Monte Carlo and its variants [@neal11; @hoffman-gelman14], for which a highly efficient and self-tuning implementation is provided by the \textsf{R}/\textsf{C++} package **RStan** [@rstan].
@chen-etal24 compares the speed and accuracy of the Laplace method implemented in **SpatialGEV** to **RStan**. 
It is found that, while **SpatialGEV** tends to underestimate the posterior variance of the hyperparameters, it accurately estimates the posteriors of both GEV parameters and return levels -- and does this three orders of magnitude faster than **RStan**. A well-known alternative to MCMC is the integrated Laplace approximation (INLA) method, whose \textsf{R} implementation is provided in the **R-INLA** package [@lindgren-rue15]. As an extension of the Laplace approximation, INLA is typically more accurate. However, **R-INLA** is inapplicable to GEV-GP models in which two or more GEV parameters are modeled as random effects following different Gaussian processes. In contrast, **SpatialGEV** offers more flexibility as it is straightforward for the user to choose what GEV parameters are spatial random effects. @mgcv and @youngman22 provide another means for estimating spatially varying GEV parameters via a scalable basis representation reducing the number of random effects in the model. Compared to **SpatialGEV** which keeps all random effects for inference, the basis function expansion approach is less accurate for estimating spatial processes that are not smooth or exhibit short-range correlation [@wood20; @lindgren-etal21].

# Example
## Model fitting
The main functions of the **SpatialGEV** package are `spatialGEV_fit()`, `spatialGEV_sample()`, and `spatialGEV_predict()`. This example shows how to apply these functions to analyze a simulated dataset using the GEV-GP model. The spatial domain is a $20\times 20$ regular lattice on $[0, 10] \times [0,10] \subset \mathbb{R}^2$, such that there are $n=400$ locations in total. The GEV location parameter $a(\xx)$ and the scale parameter $b(\xx)$ are generated from surfaces depicted in Figure \ref{fig:sim-par}, whereas the shape parameter $s$ is a constant $\exp(-2)$ across space. 10 to 30 observations per location are simulated from the GEV distribution conditional on the GEV parameters $(a(\xx), b(\xx), s)$. The simulated data is provided by the package as a list called `simulatedData`, whose values are calibrated to the level of total daily precipitation in mm.

\begin{figure}[h]
\centering
\begin{subfigure}{.4\textwidth}
\includegraphics[]{sim-plot-a.png}
\end{subfigure}%
\begin{subfigure}{.4\textwidth}
\includegraphics[]{sim-plot-b.png}
\end{subfigure}
\caption{The simulated GEV location parameters $a(\xx_i)$ and scale parameters $(b(\xx_i))$ plotted on regular lattices.}\label{fig:sim-par}
\end{figure}

The GEV-GP model is fit by calling `spatialGEV_fit()`. By specifying `random="ab"`, only the GEV parameters $a$ and $b$ are considered spatial random effects. Initial parameter values are passed to `init_param`, where `log_sigma_{a/b}` and `log_ell_{a/b}` are hyperparameters in the GP exponential kernel functions, and `beta_{a/b}` are regression coefficients in the GP mean function. The GP kernel function is chosen using `kernel="exp"`. Other kernel function options are the Matérn kernel and the approximate GP computation method employing an SPDE approximation to the Matérn [@lindgren-etal11]. The argument `reparam_s="positive"` means we constrain the shape parameter to be positive, i.e., its estimation is done on the log scale. Covariates to include in the mean functions can be provided in a matrix form to `X_{a/b}`. In this example, we only include the intercepts. The posterior mean estimates of the spatial random effects can be accessed from `mod_fit$report$par.random`, whereas the fixed effects can be obtained from `mod_fit$report$par.fixed`.

```r
set.seed(123)                         # set seed for reproducible results
library(SpatialGEV)                   # load package
n_loc <- 50                           # number of locations
locs <- simulatedData$locs[1:n_loc,]  # location coordinates
a <- simulatedData$a[1:n_loc]         # true GEV location parameters
logb <- simulatedData$logb[1:n_loc]   # true GEV (log) scale parameters
logs <- simulatedData$logs            # true GEV (log) shape parameter
y <- simulatedData$y[1:n_loc]         # simulated observations
# Model fitting
fit <- spatialGEV_fit(data = y, locs = locs, random = "ab",
                      init_param = list(a = rep(4, n_loc),
                                        log_b = rep(0,n_loc),
                                        s = -2,
                                        beta_a = 4, beta_b = 0,
                                        log_sigma_a = 0, log_ell_a = 1,
                                        log_sigma_b = 0, log_ell_b = 1),
                      reparam_s = "positive", kernel="exp",
                      X_a = matrix(1, nrow=n_loc, ncol=1),
                      X_b = matrix(1, nrow=n_loc, ncol=1),
                      silent=T)                
print(fit)
#> Model fitting took 8.78002285957336 seconds 
#> The model has reached relative convergence 
#> The model uses a exp kernel 
#> Number of fixed effects in the model is 7 
#> Number of random effects in the model is 100 
#> Hessian matrix is positive definite. 
#> Use spatialGEV_sample to obtain posterior samples 
```

## Sampling from the joint posterior
Now, we show how to sample 2000 times from the joint posterior distribution of the GEV parameters using the function `spatialGEV_sample()`. Only three arguments need to be passed to this function: `model` takes in the list output by `spatialGEV_fit()`, `n_draw` is the number of samples to draw from the posterior distribution, and `observation` indicates whether to draw from the posterior predictive distribution of the data at the observed locations. Call `summary()` on the sample object to obtain summary statistics of the posterior samples.

```r
sam <- spatialGEV_sample(model = fit, n_draw = 2000, observation = TRUE)
print(sam)
#> The samples contains 2000 draws of 107 parameters 
#> The samples contains 2000 draws of response at 50 locations 
#> Use summary() to obtain summary statistics of the samples
pos_summary <- summary(sam)
```

The samples are then used to calculate the posterior mean estimate of the 10-year return level $z_{10}(\xx)$ at each location, which are plotted against their true values in Figure \ref{fig:sim-return-level}.
```r
library(evd)
# True return levels
z_true <- unlist(Map(evd::qgev, p=0.1, loc=a, scale=exp(logb),
              shape=exp(logs), lower.tail=F))
# Posterior samples of return levels at all locations
return_period <- 10 
z_draws <- apply(sam$parameter_draws, 1,
                 function(all_draw){
                 mapply(evd::qgev, p=1/return_period, 
                        loc=all_draw[paste0("a", 1:n_loc)], 
                        scale=exp(all_draw[paste0("log_b", 1:n_loc)]),
                        shape=exp(all_draw["s"]), 
                        lower.tail=F)
                 })
z_mean <- apply(z_draws, 1, mean)
plot(z_true, z_mean, xlab="True", ylab="Posterior mean",
     main="10-year return levels (mm)")
abline(0, 1, lty="dashed", col="blue")
```

![Posterior mean estimates of the 10-year return level $z_{10}(\xx)$ plotted against the true values at different locations.\label{fig:sim-return-level}](sim-return-level.png){width=45%}

## Prediction at new locations
Next, we show how to predict the values of the extreme event at test locations. First, we divide the simulated dataset into training and test sets, and fit the model to the training dataset using the Matérn-SPDE kernel. We can simulate from the posterior predictive distribution of observations at the test locations using the `spatialGEV_predict()` function, which requires the fit model to the training data passed to `model`, a matrix of the coordinates of the test locations passed to `locs_new`, and the number of simulation draws passed to `n_draw`. Figure \ref{fig:sim-pred} plots the 90\% quantile values of the posterior predictive distributions against the 90\% quantile values of the observations at all test locations.

```r
set.seed(123)
n_test <- 100                               # number of test locations
test_ind <- sample(1:400, n_test)           # indices of the test locations
locs_test <- simulatedData$locs[test_ind,]  # coordinates of the test locations
y_test <- simulatedData$y[test_ind]         # observations at the test locations
locs_train <- simulatedData$locs[-test_ind,]# coordinates of the training locations
y_train <- simulatedData$y[-test_ind]       # observations at the training locations

# Fit the GEV-GP model to the training set
train_fit <- spatialGEV_fit(data = y_train, locs = locs_train, random = "ab", 
                            init_param = list(a = simulatedData$a[-test_ind], 
                                              log_b = simulatedData$logb[-test_ind], 
                                              s = -2,
                                              beta_a = 60, beta_b = 2,
                                              log_sigma_a = 0, log_kappa_a = -1,
                                              log_sigma_b = 0, log_kappa_b = -1),
                            reparam_s = "positive", kernel="spde", silent=T)
                          
# Make predictions at the test locations
pred <- spatialGEV_predict(model = train_fit, locs_new = locs_test, 
                           n_draw = 2000)
plot(sapply(y_test, quantile, probs=0.9),
     apply(pred$pred_y_draws, 2, quantile, probs=0.9), 
     xlim=c(3,10), ylim=c(3,10), 
     main="90% quantiles of responses at test locations",
     xlab="Observed",
     ylab="Predicted")
abline(0, 1, lty="dashed", col="blue")
```
![90\% quantile values of posterior predictive distributions at test locations plotted against the observed 90\% quantile values at the corresponding locations. Each circle corresponds to a test location. \label{fig:sim-pred}](sim-pred.png){width=45%}

# Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada, grant numbers RGPIN-2018-04376 (Ramezan), DGECR-2018-00349 (Ramezan) and RGPIN-2020-04364 (Lysy).

# References
