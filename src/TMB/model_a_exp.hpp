#ifndef model_a_exp_hpp
#define model_a_exp_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a_exp(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Sigma_a(sigma_a, ell_a))
  */ 
  using namespace density;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(loc_ind); // location index to which each observation in y is associated
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained to positive: 1 , constrained to be negative: 2, or unconstrained: 3 
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER(log_b); // log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a

  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);

  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  // data layer
  nll += nll_accumulator_a<Type>(y, loc_ind, a, log_b, s, reparam_s);  
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += gp_exp_nlpdf<Type>(mu_a, dd, sigma_a, ell_a, sp_thres);
  // prior
  nll += nll_accumulator_s_prior<Type>(s, s_mean, s_sd);
  nll += nll_accumulator_beta<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);

  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
