#ifndef model_abs_spde_maxsmooth_hpp
#define model_abs_spde_maxsmooth_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_spde_maxsmooth(objective_function<Type>* obj){
  /*
  Model layer 1: theta_hat ~ MVN(theta), theta_cov),
                 theta = (a, logb, g(s))
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
  DATA_MATRIX(obs); // 3 x n_loc matrix of parameter estimates.
  DATA_MATRIX(cov_obs); // 3 x (3*n_loc) matrix of corresponding variance estiates.
  DATA_MATRIX(design_mat_a); // n_mesh x r design matrix for a
  DATA_MATRIX(design_mat_b); // n_mesh x r design matrix for b
  DATA_MATRIX(design_mat_s); // n_mesh x r design matrix for s
  DATA_INTEGER(reparam_s); // currently unused
  DATA_IVECTOR(loc_ind); // n_loc vector of locations indices in the mesh matrix.
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov kernel
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  // Type of prior on beta. 1 is weakly informative normal prior and any other numbers
  // mean noninformative uniform prior U(-inf, inf).
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
  PARAMETER_VECTOR(a); // random effect to be integrated out. Length is dim(Q)[1] instead of n.
  PARAMETER_VECTOR(log_b); // random effect to be integrated out
  PARAMETER_VECTOR(s); // random effect to be integrated out
  PARAMETER_VECTOR(beta_a); // r x 1 mean vector coefficients for a
  PARAMETER_VECTOR(beta_b); // r x 1 mean vector coefficients for b
  PARAMETER_VECTOR(beta_s); // r x 1 mean vector coefficients for s
  PARAMETER(log_sigma_a); // hyperparameter for Sigma_a
  PARAMETER(log_kappa_a); // hyperparameter
  PARAMETER(log_sigma_b); // hyperparameter for Sigma_b
  PARAMETER(log_kappa_b); // hyperparameter
  PARAMETER(log_sigma_s); // hyperparameter for Sigma_s
  PARAMETER(log_kappa_s); // as above

  int n_param = 3;
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b);
  Type kappa_b = exp(log_kappa_b);
  Type sigma_s = exp(log_sigma_s);
  Type kappa_s = exp(log_kappa_s);

  // calculate the negative log likelihood
  Type nll = Type(0.0);
  // data layer: Normal distribution
  vector<Type> offset(n_param);
  matrix<Type> cov_block(n_param,n_param);
  vector<Type> mu_obs(n_param);
  for(int i=0; i<loc_ind.size(); i++) {
    offset(0) = a(loc_ind(i));
    offset(1) = log_b(loc_ind(i));
    offset(2) = s(loc_ind(i));
    cov_block = cov_obs.block(0,i*n_param,n_param,n_param);
    mu_obs = obs.col(i).array() - offset;
    nll += MVNORM(cov_block)(mu_obs);
  }
  // GP latent layer
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll += gp_spde_nlpdf<Type>(mu_a, spde, sigma_a, kappa_a, nu);
  nll += gp_spde_nlpdf<Type>(mu_b, spde, sigma_b, kappa_b, nu);
  nll += gp_spde_nlpdf<Type>(mu_s, spde, sigma_s, kappa_s, nu);
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

  // ------------- Output z -----------------------
  DATA_VECTOR(return_levels);
  int n_loc = spde.M0.rows();
  int has_returns = return_levels(0) > Type(0.0);
  matrix<Type> z(return_levels.size(), n_loc);
  if (has_returns){
    for(int i=0; i<n_loc; i++) {
      gev_reparam_quantile<Type>(z.col(i), return_levels,
                                 a(i), log_b(i), s(i), reparam_s);
    }
  }
  ADREPORT(z);

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
