#ifndef model_a_spde_hpp
#define model_a_spde_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_a_spde(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
  Model layer 2: a ~ GP(0, Matern)
  */ 
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
 
  // data inputs
  DATA_VECTOR(y); // response vector
  DATA_IVECTOR(loc_ind); // meshlocidx, location index to which each observation in y is associated 
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_INTEGER(reparam_s); // a flag indicating whether the shape parameter is zero: 0, constrained   to positive: 1 , constrained to be negative: 2, or unconstrained: 3
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov.
  DATA_SCALAR(s_mean); // The mean of the normal prior on s or log(|s|), depending on what reparametrization is used for s. 
  DATA_SCALAR(s_sd); // The standard deviation of the normal prior on s or log(|s|). If s_sd>9999, a flat prior is imposed.
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
  DATA_VECTOR(beta_a_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_INTEGER(a_pc_prior); // 1 for using PC prior on a, 0 for using flat unif prior 
  DATA_VECTOR(range_a_prior); // length 2 vector (rho_0, p_rho) s.t. P(rho < rho_0) = p_rho
  DATA_VECTOR(sigma_a_prior); // length 2 vector (sig_0, p_sig) s.t. P(sig > sig_0) = p_sig
  // parameter list
  PARAMETER_VECTOR(a); // random effect to be integrated out. 
  PARAMETER(log_b); // log-transformed scale parameters of the GEV model  
  PARAMETER(s); // initial shape parameter of the GEV model. IMPORTANT: If reparam_s = "negative" or "postive", the initial input should be log(|s|)
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER(log_sigma_a); // hyperparameter for the Matern SPDE for a
  PARAMETER(log_kappa_a); // as above

  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  
  // calculate the negative log likelihood
  Type nll = Type(0.0); 
  // data layer
  nll += nll_accumulator_a<Type>(y, loc_ind, a, log_b, s, reparam_s);
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  nll += gp_spde_nlpdf<Type>(mu_a, spde, sigma_a, kappa_a, nu); 
  // prior
  nll += nll_accumulator_s_prior<Type>(s, s_mean, s_sd);
  nll += nll_accumulator_beta<Type>(beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll += nll_accumulator_matern_hyperpar<Type>(log_kappa_a, log_sigma_a, a_pc_prior,
                                               nu, range_a_prior, sigma_a_prior);
  return nll;  
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
