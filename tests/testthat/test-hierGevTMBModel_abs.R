context("hierGev_model_abs")

test_that("The TMB model abs gives same likelihood as `Likelihood calculated in R`", {
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
    log.sigs <- sim$log.sigs
    log.ells <- sim$log.ells
    y <- sim$y
    dd <- sim$dd
    
    # test model ab
    b <- exp(log.b)
    s <- exp(log.s)
    a <- s/b + runif(n.obs,0,2)
    logit.tau <- logit(abs(a-s/b) / (max(abs(a-s/b)) + 1))
    y <- runif(n.obs, max(a), max(a)+10)
    y_min <- max(abs(a-s/b)) + 1
    
    nll_r <- hierGev_nll(logit.tau = logit.tau, log.b = log.b, log.s = log.s,
                           log.siga = log.siga, log.ella = log.ella,
                           log.sigb = log.sigb, log.ellb = log.ellb,
                           log.sigs = log.sigs, log.ells = log.ells, 
                           y = y, y_min=y_min, dd = dd)

    # f = MakeADFun(data=list(model="model_abs", n=n.obs,y=y,y_min=y_min,dd=dd,sp_thres=0),
    #               parameters=list(logit_tau = logit.tau, log_b = log.b, log_s = log.s, 
    #                               log_sigma_a = log.siga, log_ell_a = log.ella,
    #                               log_sigma_b = log.sigb, log_ell_b = log.ellb,
    #                               log_sigma_s = log.sigs, log_ell_s = log.ells),
    #               DLL = "SpatialGEV_TMBExports", 
    #               silent = TRUE)
    # nll_tmb <- f$fn(c(logit.tau, log.b, log.s, log.siga, log.ella, log.sigb, log.ellb, log.sigs, log.ells))
    
    nll_tmb <- hierGev_tmb_nll(y=y, dd=dd, y_min=y_min, 
                               theta=list(logit_tau = logit.tau, log_b = log.b, log_s = log.s, 
                                          log_sigma_a = log.siga, log_ell_a = log.ella,
                                          log_sigma_b = log.sigb, log_ell_b = log.ellb,
                                          log_sigma_s = log.sigs, log_ell_s = log.ells),
                               random = c("a", "b", "s"))
    if (is.na(nll_tmb)) nll_tmb <- Inf
    nll_r-nll_tmb
    expect_equal(nll_r, nll_tmb)
  }
})

