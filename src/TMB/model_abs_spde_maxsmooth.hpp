#ifndef model_abs_spde_maxsmooth_hpp
#define model_abs_spde_maxsmooth_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_spde_maxsmooth(objective_function<Type>* obj){
/// TMB specification of the Max-and-Smooth GEV-GP model with the SPDE kernel
/// Model layer 1: theta_hat ~ MVN(theta), theta_cov),
///               theta = (a, logb, g(s))
/// Model layer 2:
/// a ~ GP(0, Matern_SPDE)
/// logb ~ GP(0, Matern_SPDE)
/// g(s) ~ GP(0, Matern_SPDE) where s is a transformation function of s
/// --------- Data provided from R ---------------
/// @param[in] obs 3 x n_loc matrix of parameter estimates.
/// @param[in] cov_obs 3 x (3*n_loc) matrix of corresponding variance estiates.
/// @param[in] loc_ind n_loc vector of locations indices in the mesh matrix.
/// @param[in] reparam_s Currently unused.
/// @param[in] beta_prior Integer specifying the type of prior on the design
/// matrix coefficients. 1 is weakly informative normal prior and any other
/// numbers means Lebesgue prior `pi(beta) \propto 1`.
/// @param[in] return_periods Vector of return periods to ADREPORT. If the first
/// element of this vector is 0, then no return level calculations are performed
/// .
/// @param[in] spde Object of type `spde_t` as constructed in R by a call to
/// [INLA::inla.spde2.matern()] consisting of `n_loc` mesh locations.
/// @param[in] design_mat_a Design matrix of size
/// `n_loc x n_covariate` for parameter a.
/// @param[in] beta_a_prior Vector of length 2 containing the mean
/// and sd of the normal prior on `beta_a`.
/// @param[in] design_mat_b Design matrix of size
/// `n_loc x n_covariate` for parameter log_b.
/// @param[in] beta_b_prior Vector of length 2 containing the mean
/// and sd of the normal prior on `beta_b`.
/// @param[in] design_mat_s Design matrix of size
/// `n_loc x n_covariate` for parameter s.
/// @param[in] beta_s_prior Vector of length 2 containing the mean
/// and sd of the normal prior on `beta_s`.
/// @param[in] nu Presepecified smoothness parameter for the Matérn covariance
/// kernel applicable to all random effects.
/// @param[in] a_pc_prior Integer specifying the type of prior to
/// use on the Matérn GP on a. 1 for using PC prior on
/// a, 0 for using Lebesgue prior.
/// @param[in] range_a_prior PC prior on the range parameter for
/// the Matérn GP on
/// a. Vector of length 2 `(rho_0, p_rho)` s.t.
/// `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_a_prior PC prior on the variance parameter for
/// the Matérn GP on
/// a. Vector of length 2 `(sig_0, p_sig)` s.t.
/// `Pr(sig > sig_0) = p_sig`.
/// @param[in] b_pc_prior Integer specifying the type of prior to
/// use on the Matérn GP on log_b. 1 for using PC prior on
/// log_b, 0 for using Lebesgue prior.
/// @param[in] range_b_prior PC prior on the range parameter for
/// the Matérn GP on
/// log_b. Vector of length 2 `(rho_0, p_rho)` s.t.
/// `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_b_prior PC prior on the variance parameter for
/// the Matérn GP on
/// log_b. Vector of length 2 `(sig_0, p_sig)` s.t.
/// `Pr(sig > sig_0) = p_sig`.
/// @param[in] s_pc_prior Integer specifying the type of prior to
/// use on the Matérn GP on s. 1 for using PC prior on
/// s, 0 for using Lebesgue prior.
/// @param[in] range_s_prior PC prior on the range parameter for
/// the Matérn GP on
/// s. Vector of length 2 `(rho_0, p_rho)` s.t.
/// `Pr(rho < rho_0) = p_rho`.
/// @param[in] sigma_s_prior PC prior on the variance parameter for
/// the Matérn GP on
/// s. Vector of length 2 `(sig_0, p_sig)` s.t.
/// `Pr(sig > sig_0) = p_sig`.
///
/// --------- Parameters to estimate ------------
/// @param[in] a GEV location parameter.
/// Vector of length `n_loc`.
/// @param[in] log_b GEV scale parameter on the log scale.
/// Vector of length `n_loc`.
/// @param[in] s GEV shape parameter on the scale specified by `reparam_s`.
/// Vector of length `n_loc`.
/// @param[in] beta_a GP mean covariate coefficient vector of
/// length `n_covariate` for a.
/// @param[in] log_sigma_a GP covariance kernel variance
/// hyperparameter for a.
/// @param[in] log_kappa_a GP covariance kernel range
/// hyperparameter for a.
/// @param[in] beta_b GP mean covariate coefficient vector of
/// length `n_covariate` for log_b.
/// @param[in] log_sigma_b GP covariance kernel variance
/// hyperparameter for log_b.
/// @param[in] log_kappa_b GP covariance kernel range
/// hyperparameter for log_b.
/// @param[in] beta_s GP mean covariate coefficient vector of
/// length `n_covariate` for s.
/// @param[in] log_sigma_s GP covariance kernel variance
/// hyperparameter for s.
/// @param[in] log_kappa_s GP covariance kernel range
/// hyperparameter for s.

  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;
  // data inputs
  DATA_MATRIX(obs);
  DATA_MATRIX(cov_obs);
  DATA_MATRIX(design_mat_a);
  DATA_MATRIX(design_mat_b);
  DATA_MATRIX(design_mat_s);
  DATA_INTEGER(reparam_s);
  DATA_IVECTOR(loc_ind);
  DATA_SCALAR(nu);
  DATA_STRUCT(spde, spde_t);
  DATA_INTEGER(beta_prior);
  DATA_VECTOR(beta_a_prior);
  DATA_VECTOR(beta_b_prior);
  DATA_VECTOR(beta_s_prior);
  DATA_INTEGER(a_pc_prior);
  DATA_VECTOR(range_a_prior);
  DATA_VECTOR(sigma_a_prior);
  DATA_INTEGER(b_pc_prior);
  DATA_VECTOR(range_b_prior);
  DATA_VECTOR(sigma_b_prior);
  DATA_INTEGER(s_pc_prior);
  DATA_VECTOR(range_s_prior);
  DATA_VECTOR(sigma_s_prior);
  // parameter list
  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);
  PARAMETER_VECTOR(beta_a);
  PARAMETER_VECTOR(beta_b);
  PARAMETER_VECTOR(beta_s);
  PARAMETER(log_sigma_a);
  PARAMETER(log_kappa_a);
  PARAMETER(log_sigma_b);
  PARAMETER(log_kappa_b);
  PARAMETER(log_sigma_s);
  PARAMETER(log_kappa_s);

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
  nll += nlpdf_gp_spde<Type>(mu_a, spde, sigma_a, kappa_a, nu);
  nll += nlpdf_gp_spde<Type>(mu_b, spde, sigma_b, kappa_b, nu);
  nll += nlpdf_gp_spde<Type>(mu_s, spde, sigma_s, kappa_s, nu);
  // prior
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior, beta_a_prior[0],
      beta_a_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_b, beta_prior, beta_b_prior[0],
      beta_b_prior[1]);
  nll += nlpdf_beta_prior<Type>(beta_s, beta_prior, beta_s_prior[0],
      beta_s_prior[1]);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_a, log_sigma_a, a_pc_prior,
					   nu, range_a_prior, sigma_a_prior);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_b, log_sigma_b, b_pc_prior,
					   nu, range_b_prior, sigma_b_prior);
  nll += nlpdf_matern_hyperpar_prior<Type>(log_kappa_s, log_sigma_s, s_pc_prior,
					   nu, range_s_prior, sigma_s_prior);

  // ------------- Output return levels -----------------------
  DATA_VECTOR(return_periods);
  int n_loc = spde.M0.rows();
  int has_returns = return_periods(0) > Type(0.0);
  matrix<Type> return_levels(return_periods.size(), n_loc);
  if (has_returns){
    for(int i=0; i<n_loc; i++) {
      gev_reparam_quantile<Type>(return_levels.col(i), return_periods,
                                 a(i), log_b(i), s(i), reparam_s);
    }
  }
  ADREPORT(return_levels);

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
