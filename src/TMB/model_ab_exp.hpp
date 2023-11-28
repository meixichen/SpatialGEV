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
  DATA_IVECTOR(loc_ind); // location index to which each observation in y is associated
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is "zero", "unconstrained", constrained to be "negative", or constrained to be "positve"
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_b_prior); // length 2 vector containing mean and sd of normal prior on beta
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If tail = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_sigma_b); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_b
  PARAMETER(log_ell_b); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_b

  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);
  Type sigma_b = exp(log_sigma_b);
  Type ell_b = exp(log_ell_b);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0);
  // data layer
  nll += nll_accumulator_ab<Type>(y, loc_ind, a, log_b, s, reparam_s); 
  // GP latent layer 
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += gp_exp_nlpdf<Type>(mu_a, dd, sigma_a, ell_a, sp_thres);
  nll += gp_exp_nlpdf<Type>(mu_b, dd, sigma_b, ell_b, sp_thres);
  // prior
  nll += nll_accumulator_s_prior<Type>(s, s_mean, s_sd);
  nll += nll_accumulator_beta<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll += nll_accumulator_beta<Type>(beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
