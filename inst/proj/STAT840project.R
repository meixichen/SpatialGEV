load(system.file("extdata", "mws.RData", package = "SpatialGEV"))
start.time <- Sys.time()
# Test model with only location~GP
set.seed(322) 
n.obs <- 300
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_a <- hierGev_fit(y, real_test[,c(1,2)], 
                     list(logit_tau = rep(1,n.obs), log_b = 1, log_s = -1, 
                          log_sigma = 0, log_ell = 0),
                     random = "a")
fit_a$convergence
fit_a_rep <- fit_a$estimates
tau_est <- fit_a_rep[rownames(fit_a_rep)=="logit_tau",]
map_plot(tau_est[,1], x=real_test[,c(1,2)], file="model_a_location_estimates.png", 
         title="Spatial variation of the transformed location parameter estimates")
map_plot(tau_est[,2], x=real_test[,c(1,2)], file="model_a_location_errors.png", 
         title="Spatial variation of the transformed location parameter estimate errors")

# Test model with only scale~GP
set.seed(234) 
n.obs <- 200
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_b <- hierGev_fit(y, real_test[,c(1,2)], 
                     list(a = 1, log_b = rep(0,n.obs), log_s = -1, 
                          log_sigma = 0, log_ell = 0),
                     random = "b")
fit_a$convergence
fit_b_rep <- fit_b$estimates
logb_est <- fit_b_rep[rownames(fit_b_rep)=="log_b",]
map_plot(logb_est[,1], x=real_test[,c(1,2)], file="model_b_scale_estimates.png", 
         title="Spatial variation of the transformed scale parameter estimates")
map_plot(logb_est[,2], x=real_test[,c(1,2)], file="model_b_scale_errors.png", 
         title="Spatial variation of the transformed scale parameter estimate errors")

# Test model with shape~GP
set.seed(211) # worked on seed 211 n 150
n.obs <- 150
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_s <- hierGev_fit(y, real_test[,c(1,2)], 
                     list(a = 1, log_b = 1, log_s = rep(0,n.obs), 
                          log_sigma = 0, log_ell = 0),
                     random = "s")
fit_s$convergence
fit_s_rep <- fit_s$estimates
logs_est <- fit_s_rep[rownames(fit_s_rep)=="log_s",]
map_plot(logs_est[,1], x=real_test[,c(1,2)], file="model_s_shape_estimates.png", 
         title="Spatial variation of the transformed shape parameter estimates")
map_plot(logs_est[,2], x=real_test[,c(1,2)], file="model_s_shape_errors.png", 
         title="Spatial variation of the transformed shape parameter estimate errors")

# Test model with location,scale~GP
set.seed(222)
n.obs <- 300
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_ab <- hierGev_fit(y, real_test[,c(1,2)], 
                      list(logit_tau = rep(5,n.obs), log_b = rep(2,n.obs), log_s = 0, 
                           log_sigma_a = 2, log_ell_a = 10,
                           log_sigma_b = 2, log_ell_b = 10),
                      random = c("a","b"))
fit_ab$convergence
fit_ab_rep <- fit_ab$estimates
tau_est <- fit_ab_rep[rownames(fit_ab_rep)=="logit_tau",]
logb_est <- fit_ab_rep[rownames(fit_ab_rep)=="log_b",]
map_plot(tau_est[,1], x=real_test[,c(1,2)], file="model_ab_location_estimates.png", 
         title="Spatial variation of the transformed location parameter estimates")
map_plot(tau_est[,2], x=real_test[,c(1,2)], file="model_ab_location_errors.png", 
         title="Spatial variation of the transformed location parameter estimate errors")
map_plot(logb_est[,1], x=real_test[,c(1,2)], file="model_ab_scale_estimates.png", 
         title="Spatial variation of the transformed scale parameter estimates")
map_plot(logb_est[,2], x=real_test[,c(1,2)], file="model_ab_scale_errors.png", 
         title="Spatial variation of the transformed scale parameter estimate errors")

# Test model with location,shape~GP
set.seed(222)
n.obs <- 300
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_as <- hierGev_fit(y, real_test[,c(1,2)], 
                      list(logit_tau = rep(5, n.obs), log_b = 2, log_s = rep(0,n.obs), 
                           log_sigma_a = 2, log_ell_a = 10,
                           log_sigma_s = 2, log_ell_s = 10),
                      random = c("a", "s"))
fit_as$convergence
fit_as_rep <- fit_as$estimates
tau_est <- fit_as_rep[rownames(fit_as_rep)=="logit_tau",]
logs_est <- fit_as_rep[rownames(fit_as_rep)=="log_s",]
map_plot(tau_est[,1], x=real_test[,c(1,2)], file="model_as_location_estimates.png", 
         title="Spatial variation of the transformed location parameter estimates")
map_plot(tau_est[,2], x=real_test[,c(1,2)], file="model_as_location_errors.png", 
         title="Spatial variation of the transformed location parameter estimate errors")
map_plot(logs_est[,1], x=real_test[,c(1,2)], file="model_as_shape_estimates.png", 
         title="Spatial variation of the transformed shape parameter estimates")
map_plot(logs_est[,2], x=real_test[,c(1,2)], file="model_as_shape_errors.png", 
         title="Spatial variation of the transformed shape parameter estimate errors")


# Test model with scale,shape~GP
set.seed(222)
n.obs <- 300
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_bs <- hierGev_fit(y, real_test[,c(1,2)], 
                      list(a = 5, log_b = rep(2,n.obs), log_s = rep(0,n.obs), 
                           log_sigma_b = 2, log_ell_b = 10,
                           log_sigma_s = 2, log_ell_s = 10),
                      random = c("b", "s"))
fit_bs$convergence
fit_bs_rep <- fit_bs$estimates
logb_est <- fit_bs_rep[rownames(fit_bs_rep)=="log_b",]
logs_est <- fit_bs_rep[rownames(fit_bs_rep)=="log_s",]
map_plot(logb_est[,1], x=real_test[,c(1,2)], file="model_bs_scale_estimates.png", 
         title="Spatial variation of the transformed scale parameter estimates")
map_plot(logb_est[,2], x=real_test[,c(1,2)], file="model_bs_scale_errors.png", 
         title="Spatial variation of the transformed scale parameter estimate errors")
map_plot(logs_est[,1], x=real_test[,c(1,2)], file="model_bs_shape_estimates.png", 
         title="Spatial variation of the transformed shape parameter estimates")
map_plot(logs_est[,2], x=real_test[,c(1,2)], file="model_bs_shape_errors.png", 
         title="Spatial variation of the transformed shape parameter estimate errors")

# Test model with all parameters~GP
set.seed(222) # gave good results with seed=222 and n.obs=300
n.obs <- 300
real_test <- mws.df[sample(1:nrow(mws.df), n.obs),c(1,2,4)]
y <- real_test[,3]
fit_abs <- hierGev_fit(y, real_test[,c(1,2)], 
                       list(logit_tau = rep(5, n.obs), log_b = rep(2,n.obs), log_s = rep(0,n.obs), 
                            log_sigma_a = 5, log_ell_a = 10,
                            log_sigma_b = 2, log_ell_b = 10,
                            log_sigma_s = 2, log_ell_s = 10),
                       random = c("a", "b", "s"))
fit_abs$convergence
fit_abs_rep <- fit_abs$estimates
tau_est <- fit_abs_rep[rownames(fit_abs_rep)=="logit_tau",]
logb_est <- fit_abs_rep[rownames(fit_abs_rep)=="log_b",]
logs_est <- fit_abs_rep[rownames(fit_abs_rep)=="log_s",]
map_plot(tau_est[,1], x=real_test[,c(1,2)], file="model_abs_location_estimates.png", 
         title="Spatial variation of the transformed location parameter estimates")
map_plot(tau_est[,2], x=real_test[,c(1,2)], file="model_abs_location_errors.png", 
         title="Spatial variation of the transformed location parameter estimate errors")
map_plot(logb_est[,1], x=real_test[,c(1,2)], file="model_abs_scale_estimates.png", 
         title="Spatial variation of the transformed scale parameter estimates")
map_plot(logb_est[,2], x=real_test[,c(1,2)], file="model_abs_scale_errors.png", 
         title="Spatial variation of the transformed scale parameter estimate errors")
map_plot(logs_est[,1], x=real_test[,c(1,2)], file="model_abs_shape_estimates.png", 
         title="Spatial variation of the transformed shape parameter estimates")
map_plot(logs_est[,2], x=real_test[,c(1,2)], file="model_abs_shape_errors.png", 
         title="Spatial variation of the transformed shape parameter estimate errors")

end.time <- Sys.time()
paste("Time taken:", end.time-start.time, sep = " ")