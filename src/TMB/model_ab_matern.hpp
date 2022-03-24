#ifndef model_ab_matern_hpp
#define model_ab_matern_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_ab_matern(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Sigma_a(phi_a, kappa_a))
                 logb ~ GP(0, Sigma_b(phi_b, kappa_b))
  */ 
  using namespace density;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(n_obs); // number of observations per location
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_MATRIX(dd); // distance matrix
  DATA_SCALAR(sp_thres); // a number used to make the covariance matrix sparse by thresholding. If sp_thres=-1, no thresholding is made.
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained to positive: 1 , constrained to be negative: 2, or unconstrained: 3 
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. This has default value of 1 in INLA.
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_b_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_INTEGER(a_pc_prior); // 1 for using PC prior on a, 0 for using flat unif prior 
  DATA_VECTOR(range_a_prior); // length 2 vector (rho_0, p_rho) s.t. P(rho < rho_0) = p_rho
  DATA_VECTOR(sigma_a_prior); // length 2 vector (sig_0, p_sig) s.t. P(sig > sig_0) = p_sig
  DATA_INTEGER(b_pc_prior); 
  DATA_VECTOR(range_b_prior);
  DATA_VECTOR(sigma_b_prior);
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER_VECTOR(log_b); // random effect to be integrated out: log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER(log_sigma_a); // hyperparameter for the Matern for a (This is in fact sigma^2)
  PARAMETER(log_kappa_a); // hyperparameter 
  PARAMETER(log_sigma_b); // hyperparameter for the Matern for log(b)
  PARAMETER(log_kappa_b); // hyperparameter 

  int n = n_obs.size();
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b);
  Type kappa_b = exp(log_kappa_b);
  
  // construct the covariance matrix
  matrix<Type> cova(n,n);
  matrix<Type> covb(n,n);
  cov_matern<Type>(cova, dd, sigma_a, kappa_a, nu, sp_thres);
  cov_matern<Type>(covb, dd, sigma_b, kappa_b, nu, sp_thres);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  nll_accumulator_ab<Type>(nll, y, n_obs, a, log_b, s, n, reparam_s, s_mean, s_sd);  
 
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += MVNORM(cova)(mu_a);
  nll += MVNORM(covb)(mu_b);
  
  // prior
  nll_accumulator_beta<Type>(nll, beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll_accumulator_beta<Type>(nll, beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll_accumulator_matern_hyperpar<Type>(nll, log_kappa_a, log_sigma_a, a_pc_prior,
                                        nu, range_a_prior, sigma_a_prior);
  nll_accumulator_matern_hyperpar<Type>(nll, log_kappa_b, log_sigma_b, b_pc_prior,
                                        nu, range_b_prior, sigma_b_prior);

  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
