context("model_ab_matern")

test_that("`model_ab_matern` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
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
    log_kappa_a <- rnorm(1, 0.5, 0.1)
    log_sigma_b <- runif(1, 0,0.5)
    log_kappa_b <- rnorm(1, 0.5, 0.1)
    cov_a <- kernel_matern(dd, exp(log_sigma_a), exp(log_kappa_a)) 
    mean_a <- rep(rnorm(1, 1, 1), n)
    cov_b <- kernel_matern(dd, exp(log_sigma_b), exp(log_kappa_b))
    mean_b <- rep(rnorm(1, 0.5, 0.5), n)
    a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
    log_b <- mvtnorm::rmvnorm(1, mean_b, cov_b)
    beta_a <- mean(a)
    beta_b <- mean(log_b)   
     
    # Positive s
    s <- runif(1, 0.05, 0.1)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE), loc=a, scale=exp(log_b), shape=s)
    init_param=list(a=a, log_b=log_b, s=log(s), beta_a=beta_a, beta_b=beta_b, 
		    log_sigma_a=log_sigma_a, log_kappa_a=log_kappa_a,
                    log_sigma_b=log_sigma_b, log_kappa_b=log_kappa_b)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init_param=init_param,
                            reparam_s="positive",
                            sp_thres=-1,
			    kernel="matern",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_kappa_a)),
		   hyperparam_b=c(exp(log_sigma_b), exp(log_kappa_b)),
                   kernel="matern", beta_a=beta_a, beta_b=beta_b) 
    expect_equal(nll_r, nll_tmb)
    
    # Unconstrained s
    init_param$s <- s
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init_param=init_param,
                            reparam_s="unconstrained",
                            sp_thres=-1,
			    kernel="matern",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    expect_equal(nll_r, nll_tmb)
    
    # Negative s
    s <- runif(1, -0.1, -0.05)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE), loc=a, scale=exp(log_b), shape=s)
    init_param$s <- log(abs(s))
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init_param=init_param,
                            reparam_s="negative",
                            sp_thres=-1,
			    kernel="matern",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_kappa_a)),
		   hyperparam_b=c(exp(log_sigma_b), exp(log_kappa_b)),
                   kernel="matern", beta_a=beta_a, beta_b=beta_b) 
    expect_equal(nll_r, nll_tmb)
    
    # s=0
    s <- 0
    y <- unlist(Map(evd::rgev, n=1, loc=a, scale=exp(log_b), shape=s))
    init_param$s <- s
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init_param=init_param,
                            reparam_s="zero",
                            sp_thres=-1,
			    kernel="matern",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_kappa_a)),
		   hyperparam_b=c(exp(log_sigma_b), exp(log_kappa_b)),
                   kernel="matern", beta_a=beta_a, beta_b=beta_b) 
    expect_equal(nll_r, nll_tmb)
    
    # Test a different value of nu
    nu <- 0.5
    s <- runif(1, 0.05, 0.1)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE), loc=a, scale=exp(log_b), shape=s)
    init_param$s <- log(s)
    adfun <- spatialGEV_fit(y, X, random="ab",
                            init_param=init_param,
                            reparam_s="positive",
                            sp_thres=-1,
                            kernel="matern",
			    nu=nu,
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_kappa_a)),
                   hyperparam_b=c(exp(log_sigma_b), exp(log_kappa_b)),
                   kernel="matern", nu=nu, beta_a=beta_a, beta_b=beta_b)
    expect_equal(nll_r, nll_tmb)

  }
})
