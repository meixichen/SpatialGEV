#ifndef model_abs_spde_hpp
#define model_abs_spde_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_spde(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: 
  a ~ GP(0, Matern_SPDE)
  logb ~ GP(0, Matern_SPDE)
  g(s) ~ GP(0, Matern_SPDE) where s is a transformation function of s
  */ 
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(loc_ind); // location index to which each observation in y is associated
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_MATRIX(design_mat_s); // n x r design matrix for s
  // a flag indicating whether the shape parameter is zero: 0, 
  // constrained to positive: 1 , 
  // constrained to be negative: 2, 
  // or unconstrained: 3  
  DATA_INTEGER(reparam_s); 
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. 
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
 // Type of prior on beta. 1 is weakly informative normal prior and 
 // any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_INTEGER(beta_prior); 
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_b_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_s_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_INTEGER(a_pc_prior); // 1 for using PC prior on a, 0 for using flat unif prior 
  DATA_VECTOR(range_a_prior); // length 2 vector (rho_0, p_rho) s.t. P(rho < rho_0) = p_rho
  DATA_VECTOR(sigma_a_prior); // length 2 vector (sig_0, p_sig) s.t. P(sig > sig_0) = p_sig
  DATA_INTEGER(b_pc_prior); 
  DATA_VECTOR(range_b_prior);
  DATA_VECTOR(sigma_b_prior);
  DATA_INTEGER(s_pc_prior);  
  DATA_VECTOR(range_s_prior);
  DATA_VECTOR(sigma_s_prior);
  // parameter list
  PARAMETER_VECTOR(a);  
  PARAMETER_VECTOR(log_b);
  // If reparam_s = "negative" or "positive", the initial input should be log(|s|)
  PARAMETER_VECTOR(s); 
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER_VECTOR(beta_s); // r x 1 mean vector coefficients for s
  PARAMETER(log_sigma_a); // hyperparameter for the Matern SPDE for a
  PARAMETER(log_kappa_a); // as above
  PARAMETER(log_sigma_b); // hyperparameter for the Matern SPDE for log(b)
  PARAMETER(log_kappa_b); // as above
  PARAMETER(log_sigma_s); // hyperparameter for the Matern SPDE for s
  PARAMETER(log_kappa_s); // as above

  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b); 
  Type kappa_b = exp(log_kappa_b); 
  Type sigma_s = exp(log_sigma_s); 
  Type kappa_s = exp(log_kappa_s); 
 
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
  nll += nlpdf_gp_spde<Type>(mu_a, spde, sigma_a, kappa_a, nu);
  nll += nlpdf_gp_spde<Type>(mu_b, spde, sigma_b, kappa_b, nu);
  nll += nlpdf_gp_spde<Type>(mu_s, spde, sigma_s, kappa_s, nu);
  // prior
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior[0], beta_s_prior[1]);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_a, log_sigma_a, a_pc_prior,
                                              nu, range_a_prior, sigma_a_prior);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_b, log_sigma_b, b_pc_prior,
                                               nu, range_b_prior, sigma_b_prior);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_s, log_sigma_s, s_pc_prior,
                                               nu, range_s_prior, sigma_s_prior);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

