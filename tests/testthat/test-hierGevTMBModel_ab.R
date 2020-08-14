context("hierGev_model_ab")

test_that("The TMB model ab gives same likelihood as `Likelihood calculated in R`", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    n.obs <- sample(50:100, 1)
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s
    log.siga <- sim$log.siga
    log.ella <- sim$log.ella
    log.sigb <- sim$log.sigb
    log.ellb <- sim$log.ellb
    dd <- sim$dd
    
    # test model ab
    s <- exp(log.s[1])
    a <- s/(exp(log.b)) + runif(n.obs,0,2)
    logit.tau <- logit(abs(a-s/exp(log.b)) / (max(abs(a-s/exp(log.b))) + 1))
    y <- runif(n.obs, max(a), max(a)+10)
    y_min <- max(abs(a-s/exp(log.b))) + 1
    nll_r <- hierGev_nll(logit.tau = logit.tau, log.b = log.b, log.s = log.s[1], log.siga=log.siga,
                         log.ella = log.ella, log.sigb=log.sigb, log.ellb=log.ellb, y = y, y_min=y_min, dd=dd)
    
    # f = MakeADFun(data=list(model="model_ab", n=n.obs, y=y, y_min=y_min, dd=dd, sp_thres=0),
    #               parameters=list(logit_tau = logit.tau, log_b = log.b, log_s = log.s[1],
    #                               log_sigma_a = log.siga, log_ell_a = log.ella,
    #                               log_sigma_b = log.sigb, log_ell_b = log.ellb),
    #               DLL = "SpatialGEV_TMBExports",
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(logit.tau, log.b, log.s[1], log.siga, log.ella, log.sigb, log.ellb))
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, y_min=y_min, 
                               theta=list(logit_tau = logit.tau, log_b = log.b, log_s = log.s[1],
                                          log_sigma_a = log.siga, log_ell_a = log.ella,
                                          log_sigma_b = log.sigb, log_ell_b = log.ellb),
                               random = c("a","b"))
    if (is.na(nll_tmb)) nll_tmb <- Inf # TMB model gives NA if the loglikelihood is Inf.
    expect_equal(nll_r, nll_tmb)
    
  }
})

