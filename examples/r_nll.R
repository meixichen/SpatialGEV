library(SpatialGEV)
a <- simulatedData$a
logb <- simulatedData$logb
logs <- simulatedData$logs
s <- exp(logs)
y <- simulatedData$y
locs <- simulatedData$locs
dd <- as.matrix(stats::dist(locs))
log_sigma_a <- -1; log_ell_a <- 5
log_sigma_b <- -2; log_ell_b <- 10
beta_a <- mean(a); beta_b <- mean(logb)
# Negative marginal log-likelihood produced in R using the exponential kernel
nll_r <- r_nll(y, dd, a=a, log_b=logb, s=s,
               hyperparam_a=c(exp(log_sigma_a), exp(log_ell_a)),
               hyperparam_b=c(exp(log_sigma_b), exp(log_ell_b)),
               kernel="exp", beta_a=beta_a, beta_b=beta_b)
# Negative marg loglik produced by TMB template
init_param <- list(beta_a=beta_a, beta_b=beta_b,
                   a=a, log_b=logb, s=log(s),
                   log_sigma_a=log_sigma_a,
                   log_ell_a=log_ell_a,
                   log_sigma_b=log_sigma_b,
                   log_ell_b=log_ell_b)
adfun <- spatialGEV_fit(y, locs, random="ab",
                        init_param=init_param,
                        reparam_s="positive",
                        kernel="exp",
                        adfun_only=TRUE,
                        ignore_random=TRUE,
                        silent=TRUE)
nll_tmb <- adfun$fn(unlist(init_param))
nll_r - nll_tmb
