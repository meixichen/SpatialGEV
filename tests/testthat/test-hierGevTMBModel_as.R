context("hierGev_model_as")

test_that("The TMB model as gives same likelihood as `Likelihood calculated in R`", {
  n.tests <- 30 # number of test simulations
  
  for (ii in 1:n.tests){
    n.obs <- sample(50:100, 1)
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s
    log.siga <- sim$log.siga
    log.ella <- sim$log.ella
    log.sigs <- sim$log.sigs
    log.ells <- sim$log.ells
    dd <- sim$dd
    
    # test model ab
    b <- exp(log.b[1])
    a <- exp(log.s)/b + runif(n.obs,0,2)
    logit.tau <- logit(abs(a-exp(log.s)/b) / (max(abs(a-exp(log.s)/b)) + 1))
    y <- runif(n.obs, max(a), max(a)+10)
    y_min <- max(abs(a-exp(log.s)/b)) + 1
    nll_r <- hierGev_nll(logit.tau = logit.tau, log.b = log.b[1], log.s = log.s, log.siga=log.siga,
                         log.ella = log.ella, log.sigs=log.sigs, log.ells=log.ells, y = y, y_min=y_min, dd=dd)
    # f = MakeADFun(data=list(model="model_as", n=n.obs, y=y, y_min=y_min, dd=dd, sp_thres=0),
    #               parameters=list(logit_tau = logit.tau, log_b = log.b[1], log_s = log.s, 
    #                               log_sigma_a = log.siga, log_ell_a = log.ella,
    #                               log_sigma_s = log.sigs, log_ell_s = log.ells),
    #               DLL = "SpatialGEV_TMBExports", 
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(logit.tau, log.b[1], log.s, log.siga, log.ella, log.sigs, log.ells))
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, y_min=y_min, 
                               theta=list(logit_tau = logit.tau, log_b = log.b[1], log_s = log.s, 
                                          log_sigma_a = log.siga, log_ell_a = log.ella,
                                          log_sigma_s = log.sigs, log_ell_s = log.ells),
                               random = c("a","s"))
    if (is.na(nll_tmb)) nll_tmb <- Inf # TMB model gives NA if the loglikelihood is Inf.
    expect_equal(nll_r, nll_tmb)
    
  }
})

