context("model_a")

test_that("`model_a` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    #####To be deleted before making a package ######
    mod <- "model_a" ##
    compile(paste0(mod, ".cpp")) ##
    dyn.load(dynlib(mod)) ##
    #################################
    # simulate parameters and data
    n_sqrt <- sample(5:10, 1)
    n <- n_sqrt^2 
    lon <- seq(0, 10, length.out = n_sqrt)
    lat <- seq(0, 10, length.out = n_sqrt)
    X <- expand.grid(x = lon, y = lat)
    dd <- as.matrix(stats::dist(X))
    log_sigma_a <- runif(1, 0, 1)
    log_ell_a <- rnorm(1, 0.5, 0.1)
    cov_a <- exp(log_sigma_a)*exp(-dd/exp(log_ell_a))
    mean_a <- rep(rnorm(1, 1, 1), n)
    a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
    log_b <- runif(1, -3, 0)
    
    # Positive s
    s <- runif(1, 0.01, 0.1)
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    
    init.param=list(a=a, log_b=log_b, s=log(s), log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    adfun <- spatialGEV_fit(y, X, random="a",
                            init.param=init.param,
                            reparam.s="positive",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    expect_equal(nll_r, nll_tmb)
    
    # Unconstrained s
    init.param=list(a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    adfun <- spatialGEV_fit(y, X, random="a",
                            init.param=init.param,
                            reparam.s="unconstrained",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    expect_equal(nll_r, nll_tmb)
    
    # Negative s
    s <- runif(1, -0.1, -0.01)
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    
    init.param=list(a=a, log_b=log_b, s=log(abs(s)), log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    adfun <- spatialGEV_fit(y, X, random="a",
                            init.param=init.param,
                            reparam.s="negative",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    expect_equal(nll_r, nll_tmb)
    
    # s=0
    s <- 0
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    
    init.param=list(a=a, log_b=log_b, s=0, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    adfun <- spatialGEV_fit(y, X, random="a",
                            init.param=init.param,
                            reparam.s="zero",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a)
    expect_equal(nll_r, nll_tmb)
    
  }
})
