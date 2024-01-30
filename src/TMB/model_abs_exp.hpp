#ifndef model_abs_exp_hpp
#define model_abs_exp_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_exp(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: 
  a ~ GP(0, Sigma_a(sigma_a, ell_a))
  logb ~ GP(0, Sigma_b(sigma_b, ell_b))
  g(s) ~ GP(0, Sigma_s(sigma_s, ell_s)) where s is a transformation function of s
  */ 
  using namespace density;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(loc_ind); // location index to which each observation in y is associated
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_MATRIX(design_mat_s); // n x r design matrix for s
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  // a flag indicating whether the shape parameter is "zero", 
  // "unconstrained", constrained to be "negative", or constrained to be "positve"
  DATA_INTEGER(reparam_s);
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_b_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_s_prior); // length 2 vector containing mean and sd of normal prior on beta
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER_VECTOR(s); // if reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER_VECTOR(beta_s); // r x 1 mean vector coefficients for s
  PARAMETER(log_sigma_a); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_ell_a); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_a
  PARAMETER(log_sigma_b); // hyperparameter: log-transformed squared amplitude parameter (scalar) of the exponential covariance function in Sigma_b
  PARAMETER(log_ell_b); // hyperparameter: log-transformed smoothness parameter (scalar) of the exponential covariance function in Sigma_b
  PARAMETER(log_sigma_s); // hyperparameter for s
  PARAMETER(log_ell_s); // as above

  Type sigma_a = exp(log_sigma_a);
  Type ell_a = exp(log_ell_a);
  Type sigma_b = exp(log_sigma_b);
  Type ell_b = exp(log_ell_b);
  Type sigma_s = exp(log_sigma_s); 
  Type ell_s = exp(log_ell_s); 
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  // data layer
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y[i], a[loc_ind[i]], log_b[loc_ind[i]],
                                  s[loc_ind[i]], reparam_s);
  }
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll += nlpdf_gp_exp<Type>(mu_a, dd, sigma_a, ell_a, sp_thres);
  nll += nlpdf_gp_exp<Type>(mu_b, dd, sigma_b, ell_b, sp_thres);
  nll += nlpdf_gp_exp<Type>(mu_s, dd, sigma_s, ell_s, sp_thres);
  // prior
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior[0], beta_s_prior[1]);
  
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
