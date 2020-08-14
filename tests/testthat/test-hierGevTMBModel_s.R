context("hierGev_model_s")

test_that("The TMB model s gives same likelihood as `Likelihood calculated in R`", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    n.obs <- sample(50:100, 1)
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s
    log.sigs <- sim$log.sigs
    log.ells <- sim$log.ells
    dd <- sim$dd

    # test model s
    a <- max(exp(log.s))/exp(log.b[1]) + runif(1,0,2)
    y <- runif(n.obs, a, a+10)
    nll_r <- hierGev_nll(a=a, log.b = log.b[1], log.s = log.s, log.sigs=log.sigs,
                         log.ells = log.ells,y = y,  dd=dd)
    # f = MakeADFun(data=list(model="model_s", n=n.obs, y=y, dd=dd, sp_thres=0),
    #               parameters=list(a = a, log_b = log.b[1], log_s = log.s, 
    #                               log_sigma = log.sigs, log_ell = log.ells),
    #               DLL = "SpatialGEV_TMBExports", 
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(a, log.b[1], log.s, log.sigs, log.ells))
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, 
                               theta=list(a = a, log_b = log.b[1], log_s = log.s,
                                          log_sigma = log.sigs, log_ell = log.ells),
                               random = "s")

    if (is.na(nll_tmb)) nll_tmb <- Inf # TMB model gives NA if the loglikelihood is Inf.
    expect_equal(nll_r, nll_tmb)
    
    
  }
})
