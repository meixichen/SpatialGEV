#ifndef model_abs_spde_maxsmooth_hpp
#define model_abs_spde_maxsmooth_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_spde_maxsmooth(objective_function<Type>* obj){
  /* Model:
  
  obs = g(a, b, s) ~ MVN(a', b', s'), cov_obs)
  a' ~ MVN(beta_a, cov_a)
  b' ~ MVN(beta_b, cov_b)
  s' ~ MVN(beta_s, cov_s)
  theta ~ prior

  Summary:
  - Data input: obs, cov_obs
  - Random effects: a', b', s' 
  - Fixed effects: beta_* and hyperparameters theta in cov_*
  */
  using namespace density;
  using namespace SpatialGEV;
  // data inputs
  DATA_VECTOR(obs); // noisy obs of transformed c(loc_1,...,loc_n, scale_1,..., shape_n)
  DATA_VECTOR(cov_obs); // covariance matrix of the noisy observations
  DATA_MATRIX(design_mat_a); // n x r design matrix for a
  DATA_MATRIX(design_mat_b); // n x r design matrix for b
  DATA_MATRIX(design_mat_s); // n x r design matrix for s
  DATA_IVECTOR(meshidxloc); // indices of the locations in the mesh matrix
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov kernel
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  DATA_INTEGER(beta_prior); // Type of prior on beta. 1 is weakly informative normal prior and any other numbers mean noninformative uniform prior U(-inf, inf).
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
  PARAMETER_VECTOR(a); // random effect to be integrated out
  PARAMETER_VECTOR(b); // random effect to be integrated out
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

  int n = a_obs.size();
  Type sigma_a = exp(log_sigma_a);
  Type kappa_a = exp(log_kappa_a);
  Type sigma_b = exp(log_sigma_b);
  Type kappa_b = exp(log_kappa_b);
  Type sigma_s = exp(log_sigma_s);
  Type kappa_s = exp(log_kappa_s);

  // construct the covariance matrices
  matrix<Type> cova(n,n);
  matrix<Type> covb(n,n);
  matrix<Type> covs(n,n);
  cov_matern<Type>(cova, dd, sigma_a, kappa_a, nu, sp_thres);
  cov_matern<Type>(covb, dd, sigma_b, kappa_b, nu, sp_thres);
  cov_matern<Type>(covs, dd, sigma_s, kappa_s, nu, sp_thres);

  Type nll = Type(0.0);

  // spde approx
  SparseMatrix<Type> Q_a = Q_spde(spde, kappa_a);
  SparseMatrix<Type> Q_b = Q_spde(spde, kappa_b);
  SparseMatrix<Type> Q_s = Q_spde(spde, kappa_s);
  Type sigma_marg_a = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_a, 2*nu)); // marginal variance for a
  Type sigma_marg_b = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_b, 2*nu)); // marginal variance for log(b)
  Type sigma_marg_s = exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa_s, 2*nu)); // marginal variance for s
  vector<Type> mu_a = a - design_mat_a * beta_a;
  vector<Type> mu_b = log_b - design_mat_b * beta_b;
  vector<Type> mu_s = s - design_mat_s * beta_s;
  nll = SCALE(GMRF(Q_a), sigma_a/sigma_marg_a)(mu_a);
  nll += SCALE(GMRF(Q_b), sigma_b/sigma_marg_b)(mu_b);
  nll += SCALE(GMRF(Q_s), sigma_s/sigma_marg_s)(mu_s);

  // calculate the negative log likelihood
  int start_ind = 0; // index of the first observation of location i in n_obs
  int end_ind = 0; // index of the last observation of location i in n_obs
  for(int i=0;i<n;i++) {
    end_ind += n_obs[i];
    for(int j=start_ind;j<end_ind;j++){
      nll -= gumbel_lpdf<Type>(y[j], a[meshidxloc[i]], log_b[meshidxloc[i]]);
    }
    start_ind += n_obs[i];
  }

  // prior
  nll_accumulator_beta<Type>(nll, beta_a, beta_prior, beta_a_prior[0], beta_a_prior[1]);
  nll_accumulator_beta<Type>(nll, beta_b, beta_prior, beta_b_prior[0], beta_b_prior[1]);
  nll_accumulator_beta<Type>(nll, beta_s, beta_prior, beta_s_prior[0], beta_s_prior[1]);
  nll_accumulator_matern_hyperpar<Type>(nll, log_kappa_a, log_sigma_a, a_pc_prior,
                                        nu, range_a_prior, sigma_a_prior);
  nll_accumulator_matern_hyperpar<Type>(nll, log_kappa_b, log_sigma_b, b_pc_prior,
                                        nu, range_b_prior, sigma_b_prior);
  nll_accumulator_matern_hyperpar<Type>(nll, log_kappa_s, log_sigma_s, s_pc_prior,
                                        nu, range_s_prior, sigma_s_prior);
  return nll;  
   
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
