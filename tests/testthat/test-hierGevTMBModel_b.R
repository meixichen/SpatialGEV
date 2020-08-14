context("hierGev_model_b")

test_that("The TMB model b gives same likelihood as `Likelihood calculated in R`", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    n.obs <- sample(50:100, 1)
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s
    log.sigb <- sim$log.sigb
    log.ellb <- sim$log.ellb
    dd <- sim$dd

    # test model b
    s <- exp(log.s[1])
    a <- s/(min(exp(log.b))) + runif(1,0,2)
    y <- runif(n.obs, a, a+10)
    nll_r <- hierGev_nll(a=a, log.b = log.b, log.s = log.s[1], log.sigb=log.sigb,
                         log.ellb = log.ellb,y = y,  dd=dd)
    # f = MakeADFun(data=list(model="model_b", n=n.obs, y=y, dd=dd, sp_thres=0),
    #               parameters=list(a = a, log_b = log.b, log_s = log.s[1], 
    #                               log_sigma = log.sigb, log_ell = log.ellb),
    #               DLL = "SpatialGEV_TMBExports", 
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(a, log.b, log.s[1], log.sigb, log.ellb))
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, 
                               theta=list(a = a, log_b = log.b, log_s = log.s[1],
                                          log_sigma = log.sigb, log_ell = log.ellb),
                               random = "b")
    if (is.na(nll_tmb)) nll_tmb <- Inf # TMB model gives NA if the loglikelihood is Inf.
    expect_equal(nll_r, nll_tmb)

  }
})
