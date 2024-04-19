context("model_a_exp")

test_that("`model_a_exp` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
  n_tests <- 30 # number of test simulations
  for (ii in 1:n_tests){
    # simulate parameters and data
    n_sqrt <- sample(5:10, 1)
    n <- n_sqrt^2
    lon <- seq(0, 10, length.out = n_sqrt)
    lat <- seq(0, 10, length.out = n_sqrt)
    X <- expand.grid(x = lon, y = lat)
    dd <- as.matrix(stats::dist(X))
    log_sigma_a <- runif(1, 0, 1)
    log_ell_a <- rnorm(1, 0.5, 0.1)
    cov_a <- kernel_exp(dd, exp(log_sigma_a), exp(log_ell_a))
    mean_a <- rep(rnorm(1, 1, 1), n)
    a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
    log_b <- runif(1, -3, 0)
    beta_a <- mean(a)
    # Positive s
    s <- runif(1, 0.01, 0.1)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE),
             loc=a, scale=exp(log_b), shape=s)
    init_param <- list(a=a, log_b=log_b, s=log(s), beta_a=beta_a,
                       log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    adfun <- spatialGEV_fit(y, locs=X, random="a",
                            init_param=init_param,
                            reparam_s="positive",
                            kernel="exp",
                            sp_thres=-1,
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
                   kernel="exp", beta_a=beta_a)
    expect_equal(nll_r, nll_tmb)
    # Unconstrained s
    init_param$s <- s
    adfun <- spatialGEV_fit(y, locs=X, random="a",
                            init_param=init_param,
                            reparam_s="unconstrained",
                            kernel="exp", sp_thres=-1,
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    expect_equal(nll_r, nll_tmb)
    # Negative s
    s <- runif(1, -0.1, -0.01)
    ## y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE),
             loc=a, scale=exp(log_b), shape=s)
    init_param$s <- log(abs(s))
    adfun <- spatialGEV_fit(y, locs=X, random="a",
                            init_param=init_param,
                            reparam_s="negative",
                            kernel="exp", sp_thres=-1,
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
                   kernel="exp", beta_a=beta_a)
    expect_equal(nll_r, nll_tmb)
    # s=0
    s <- 0
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE),
             loc=a, scale=exp(log_b), shape=s)
    init_param$s <- 0
    adfun <- spatialGEV_fit(y, locs=X, random="a",
                            init_param=init_param,
                            reparam_s="zero",
                            kernel="exp", sp_thres=-1,
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    ## nll_tmb <- adfun$fn(unlist(init_param))
    nll_tmb <- adfun$fn(unlist(init_param[names(init_param) != "s"]))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
                   kernel="exp", beta_a=beta_a)
    expect_equal(nll_r, nll_tmb)
  }
})
