#ifndef model_ab_spde_hpp
#define model_ab_spde_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_ab_spde(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: 
  a ~ GP(0, Matern)
  logb ~ GP(0, Matern)
  */ 
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
  
  // data inputs
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(loc_ind); // meshlocidx, location index to which each observation in y is associated 
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for logb
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained to positive: 1 , constrained to be negative: 2, or unconstrained: 3  
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. 
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
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
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If tail = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for logb
  PARAMETER(log_sigma_a); // hyperparameter for the Matern SPDE for a
  PARAMETER(log_kappa_a); // as above
  PARAMETER(log_sigma_b); // hyperparameter for the Matern SPDE for b
  PARAMETER(log_kappa_b); // as above

  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b); 
  Type kappa_b = exp(log_kappa_b); 
 
  // calculate the negative log likelihood
  Type nll = Type(0.0);
  // data layer
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y[i], a[loc_ind[i]], log_b[loc_ind[i]], s, reparam_s);
  }
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  nll += nlpdf_gp_spde<Type>(mu_a, spde, sigma_a, kappa_a, nu);
  nll += nlpdf_gp_spde<Type>(mu_b, spde, sigma_b, kappa_b, nu);
  // prior
  nll += nlpdf_s_prior<Type>(s, s_mean, s_sd);
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_a, log_sigma_a, a_pc_prior,
                                               nu, range_a_prior, sigma_a_prior);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_b, log_sigma_b, b_pc_prior,
                                               nu, range_b_prior, sigma_b_prior);
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
