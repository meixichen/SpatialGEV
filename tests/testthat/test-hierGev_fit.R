context("hierGev_fit")

test_that("MLE produced by `hierGev_fit()` is able to find the at least local optimal estimates", {
  set.seed(125) 
  ntest <- 2
  for(ii in 1:ntest) {
    # Data simulation 
    n.obs <- 50
    sim <- hierGev_sim(n.obs)
    log.b <- sim$log.b 
    log.s <- sim$log.s[1]
    log.sigb <- sim$log.sigb
    log.ellb <- sim$log.ellb
    x <- sim$x
    dd <- as.matrix(sim$dd)
    s <- exp(log.s[1])
    a <- s/(min(exp(log.b))) + runif(1,0,2)
    y <- runif(n.obs, a, a+10)
    
    # Other function arguments input:
    theta = list(a = a, log_b = log.b, log_s = log.s[1], 
                 log_sigma = log.sigb, log_ell = log.ellb)
    random = "b"
    sp.thres=0
    
    # The following commands are directly copy-pasted from the `hierGev_fit` function
    candidates1 <- c("a","b","s")
    inx <- which(candidates1 %in% random)
    chosen <- candidates1[inx]
    mod <- paste("model",paste(chosen,collapse=""), sep="_")
    n <- length(y)
    dd <- as.matrix(dist(x))
    data <- list(model=mod, n=n, y=y, dd=dd, sp_thres=sp.thres)
    
    if ("a" %in% random){
      candidates2 <- c("logit_tau", "log_b", "log_s")
      random <- candidates2[inx]
      y_min <- runif(1,0, min(y))
      data[["y_min"]] <- y_min
      model <- TMB::MakeADFun(data=data,
                              parameters=theta,
                              random=random,
                              DLL = "SpatialGEV_TMBExports", 
                              silent = TRUE)
      fit <- nlminb(model$par, model$fn, model$gr)
      rep <- TMB::sdreport(model)
      #list(estimates = summary(rep), convergence = fit$convergence, y_min = y_min)
    }
    else {
      candidates2 <- c("a", "log_b", "log_s")
      random <- candidates2[inx]
      model <- TMB::MakeADFun(data=data,
                              parameters=theta,
                              random=random,
                              DLL = "SpatialGEV_TMBExports", 
                              silent = TRUE)
      fit <- nlminb(model$par, model$fn, model$gr)
      rep <- TMB::sdreport(model)
      #list(estimates = summary(rep), convergence = fit$convergence)
    }
    
    # Check numerical projection plots

    theta.names <- c('a', 'log_s','log_sigma', 'log_ell')
    oproj <- optimCheck::optim_proj(fun = model$fn,
               xsol = fit$par,
               maximize = FALSE,
               xnames = theta.names,
               plot = FALSE)
    max_diff <- abs(diff(oproj))
    expect_lt(max(pmin(max_diff[,1], max_diff[,2])), .11)
  }
})
