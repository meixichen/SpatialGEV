#ifndef model_a_matern_hpp
#define model_a_matern_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a_matern(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Sigma_a(phi_a, kappa_a))
  */ 
  using namespace density;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(n_obs); // number of observations per location
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained to positive: 1 , constrained to be negative: 2, or unconstrained: 3 
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. This has default value of 1 in INLA.
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  // parameter list
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER(log_b); // log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER(log_sigma_a); // hyperparameter of the Matern (This is in fact sigma^2) 
  PARAMETER(log_kappa_a); // hyperparameter of the Matern

  int n = n_obs.size();
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  
  // construct the covariance matrix
  matrix<Type> cova(n,n);
  cov_matern<Type>(cova, dd, sigma_a, kappa_a, nu, sp_thres);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  nll_accumulator_a<Type>(nll, y, n_obs, a, log_b, s, n, reparam_s, s_mean, s_sd);  
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += MVNORM(cova)(mu_a);

  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
