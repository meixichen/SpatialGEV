context("hierGev_model_a")

test_that("The TMB model a gives same likelihood as `Likelihood calculated in R`", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    n.obs <- sample(50:100, 1)
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s
    log.siga <- sim$log.siga
    log.ella <- sim$log.ella
    dd <- sim$dd
    
    # test model a
    s <- exp(log.s[1])
    b <- exp(log.b[1])
    a <- s/b + runif(n.obs,0,2)
    y_min <- max(a-s/b)+0.5
    logit.tau <- logit(abs(a-s/b)/y_min)
    y <- runif(n.obs, max(a), max(a)+10)
    nll_r <- hierGev_nll(logit.tau = logit.tau, log.b = log.b[1], log.s = log.s[1], log.siga=log.siga,
                         log.ella = log.ella,y = y, y_min=y_min, dd=dd)
    # f = MakeADFun(data=list(model="model_a", n=n.obs, y=y, y_min=y_min, dd=dd, sp_thres=0),
    #               parameters=list(logit_tau = logit.tau, log_b = log.b[1], log_s = log.s[1], 
    #                               log_sigma = log.siga, log_ell = log.ella),
    #               DLL = "SpatialGEV_TMBExports", 
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(logit.tau, log.b[1], log.s[1], log.siga, log.ella))
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, y_min=y_min, 
                               theta=list(logit_tau = logit.tau, 
                                          log_b = log.b[1], log_s = log.s[1], 
                                          log_sigma = log.siga, log_ell = log.ella),
                               random = "a")
    if (is.na(nll_tmb)) nll_tmb <- Inf # TMB model gives NA if the loglikelihood is Inf.
    expect_equal(nll_r, nll_tmb)
  }
})
