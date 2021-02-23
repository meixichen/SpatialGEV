context("model_ab")

test_that("`model_ab` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    # simulate parameters and data
    n_sqrt <- sample(5:10, 1)
    n <- n_sqrt^2 
    lon <- seq(0, 10, length.out = n_sqrt)
    lat <- seq(0, 10, length.out = n_sqrt)
    X <- expand.grid(x = lon, y = lat)
    dd <- as.matrix(stats::dist(X))
    log_sigma_a <- runif(1, 0, 1)
    log_ell_a <- rnorm(1, 0.5, 0.1)
    log_sigma_b <- runif(1, 0,0.5)
    log_ell_b <- rnorm(1, 0.5, 0.1)
    cov_a <- exp(log_sigma_a)*exp(-dd/exp(log_ell_a))
    mean_a <- rep(rnorm(1, 1, 1), n)
    cov_b <- exp(log_sigma_b)*exp(-dd/exp(log_ell_b))
    mean_b <- rep(rnorm(1, 0.5, 0.5), n)
    a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
    log_b <- mvtnorm::rmvnorm(1, mean_b, cov_b)
    
    # Positive s
    s <- runif(1, 0.05, 0.1)
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    init.param=list(a=a, log_b=log_b, s=log(s), log_sigma_a=log_sigma_a, log_ell_a=log_ell_a,
                    log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init.param=init.param,
                            reparam.s="positive",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a, 
                   log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    expect_equal(nll_r, nll_tmb)
    
    # Unconstrained s
    init.param=list(a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a,
                    log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init.param=init.param,
                            reparam.s="unconstrained",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    expect_equal(nll_r, nll_tmb)
    
    # Negative s
    s <- runif(1, -0.1, -0.05)
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    init.param=list(a=a, log_b=log_b, s=log(abs(s)), log_sigma_a=log_sigma_a, log_ell_a=log_ell_a,
                    log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init.param=init.param,
                            reparam.s="negative",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a, 
                   log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    expect_equal(nll_r, nll_tmb)
    
    # s=0
    s <- 0
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    init.param=list(a=a, log_b=log_b, s=0, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a,
                    log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init.param=init.param,
                            reparam.s="zero",
                            sp.thres=0,
                            adfun.only=TRUE,
                            ignore.random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init.param))
    
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s, log_sigma_a=log_sigma_a, log_ell_a=log_ell_a, 
                   log_sigma_b=log_sigma_b, log_ell_b=log_ell_b)
    expect_equal(nll_r, nll_tmb)
  }
})
