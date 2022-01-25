#ifndef model_ab_exp_hpp
#define model_ab_exp_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_ab_exp(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: 
  a ~ GP(0, Sigma_a(sigma_a, ell_a))
  logb ~ GP(0, Sigma_b(sigma_b, ell_b))
  */ 
  using namespace density;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(n_obs); // number of observations per location
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=0, no thresholding is made.
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is "zero", "unconstrained", constrained to be "negative", or constrained to be "positve"
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If tail = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_sigma_b); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_b
  PARAMETER(log_ell_b); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_b

  int n = n_obs.size();
  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);
  Type sigma_b = exp(log_sigma_b);
  Type ell_b = exp(log_ell_b);
  
  // construct the covariance matrices
  matrix<Type> cova(n,n);
  matrix<Type> covb(n,n);
  cov_expo<Type>(cova, dd, sigma_a, ell_a, sp_thres);
  cov_expo<Type>(covb, dd, sigma_b, ell_b, sp_thres);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  nll_accumulator_ab<Type>(nll, y, n_obs, a, log_b, s, n, reparam_s, s_mean, s_sd);

  nll += MVNORM(cova)(a);
  nll += MVNORM(covb)(log_b);
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
