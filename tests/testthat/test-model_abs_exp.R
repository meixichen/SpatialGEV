context("model_abs_exp")

test_that("`model_abs_exp` gives the same likelihood as the one calculated in R under different parametrizations of shape parameter", {
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
    log_sigma_s <- runif(1, -4, -2)
    log_ell_s <- rnorm(1, 1, 2)
    cov_a <- kernel_exp(dd, exp(log_sigma_a), exp(log_ell_a)) 
    mean_a <- rep(rnorm(1, 1, 1), n)
    cov_b <- kernel_exp(dd, exp(log_sigma_b), exp(log_ell_b))
    mean_b <- rep(rnorm(1, 0.5, 0.5), n)
    cov_s <- kernel_exp(dd, exp(log_sigma_s), exp(log_ell_s))
    mean_s <- rep(rnorm(1, -3, 1), n)
    a <- mvtnorm::rmvnorm(1, mean_a, cov_a)
    log_b <- mvtnorm::rmvnorm(1, mean_b, cov_b)
    log_s <- mvtnorm::rmvnorm(1, mean_s, cov_s)
    s <- exp(log_s)
    beta_a <- mean(a)
    beta_b <- mean(log_b)   
     
    # Positive s
    beta_s <- mean(log_s)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE), loc=a, scale=exp(log_b), shape=s)
    init_param=list(
		    a=a, log_b=log_b, s=log_s, 
		    beta_a=beta_a, beta_b=beta_b, beta_s=beta_s,
		    log_sigma_a=log_sigma_a, log_ell_a=log_ell_a,
                    log_sigma_b=log_sigma_b, log_ell_b=log_ell_b,
                    log_sigma_s=log_sigma_s, log_ell_s=log_ell_s)
    adfun <- spatialGEV_fit(y, X, random="abs",
                            init_param=init_param,
                            reparam_s="positive",
                            sp_thres=-1,
			    kernel="exp",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
		   hyperparam_b=c(exp(log_sigma_b), exp(log_ell_b)),
		   hyperparam_s=c(exp(log_sigma_s), exp(log_ell_s)),
                   kernel="exp", beta_a=beta_a, beta_b=beta_b, beta_s=beta_s, 
                   f_s=function(x){log(x)}) 
    expect_equal(nll_r, nll_tmb)
    
    # Unconstrained s
    beta_s <- mean(s)
    init_param$s <- s
    init_param$beta_s <- beta_s
    adfun <- spatialGEV_fit(y, X, random="abs",
                            init_param=init_param,
                            reparam_s="unconstrained",
                            sp_thres=-1,
			    kernel="exp",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
		   hyperparam_b=c(exp(log_sigma_b), exp(log_ell_b)),
		   hyperparam_s=c(exp(log_sigma_s), exp(log_ell_s)),
                   kernel="exp", beta_a=beta_a, beta_b=beta_b, beta_s=beta_s) 
    expect_equal(nll_r, nll_tmb)
    
    # Negative s
        s <- -exp(log_s)
    y <- Map(evd::rgev, n=sample(1:20, n, replace=TRUE), loc=a, scale=exp(log_b), shape=s)
    beta_s <- mean(log_s)
    init_param$s <- log_s
    init_param$beta_s <- beta_s
    adfun <- spatialGEV_fit(y, X, random="abs",
                            init_param=init_param,
                            reparam_s="negative",
                            sp_thres=-1,
                            kernel="exp",
                            adfun_only=TRUE,
                            ignore_random=TRUE,
                            silent=TRUE)
    nll_tmb <- adfun$fn(unlist(init_param))
    nll_r <- r_nll(y, dd, a=a, log_b=log_b, s=s,
                   hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
                   hyperparam_b=c(exp(log_sigma_b), exp(log_ell_b)),
                   hyperparam_s=c(exp(log_sigma_s), exp(log_ell_s)),
                   kernel="exp", beta_a=beta_a, beta_b=beta_b, beta_s=beta_s,
                   f_s=function(x){log(abs(s))})
    expect_equal(nll_r, nll_tmb)
  }
})
