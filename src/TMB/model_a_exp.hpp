#ifndef model_a_exp_hpp
#define model_a_exp_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/// TMB specification of GEV-GP models with a chosen covariance kernel.
///
/// The model is defined as follows:
///
/// y ~ GEV(a, b, s),
/// a ~ GP(log_sigma_a, log_ell_a)
/// where the GP is parameterized using the exp covariance kernel.
///
/// --------- Data provided from R ---------------
/// @param[in] y Response vector of length `n_obs`.  Assumed to be > 0.
/// @param[in] loc_ind Location vector of length `n_obs` of integers
/// `0 <= i_loc < n_loc` indicating to which locations each element of `y` is
/// associated.
/// @param[in] reparam_s Integer indicating the type of shape parameter. 0:
/// `s = 0`, i.e., use Gumbel instead of GEV distribution.  1: `s > 0`, in which
/// case we operate on `log(s)`.  2: `s < 0`, in which case we operate on
/// `log(-s)`.  3: unconstrained.
/// @param[in] beta_prior Integer specifying the type of prior on the design
/// matrix coefficients. 1 is weakly informative normal prior and any other
/// numbers means Lebesgue prior `pi(beta) \propto 1`.
/// @param[in] return_periods Vector of return periods to ADREPORT. If the first
/// element of this vector is 0, then no return level calculations are performed
/// .
/// @param[in] dist_mat `n_loc x n_loc` distance matrix typically constructed
/// via `stats::dist(coordinates)`.
/// @param[in] sp_thres Scalar number used to make the covariance matrix sparse
/// by thresholding. If sp_thres=-1, no thresholding is made.
/// @param[in] design_mat_a Design matrix of size
/// `n_loc x n_covariate` for parameter a.
/// @param[in] beta_a_prior Vector of length 2 containing the mean
/// and sd of the normal prior on `beta_a`.
/// @param[in] s_mean Scalar for Normal prior mean on s.
/// @param[in] s_sd Scalar for Normal prior sd on s.
///
/// --------- Parameters to estimate ------------
/// @param[in] a GEV location parameter.
/// Vector of length `n_loc`.
/// @param[in] log_b GEV scale parameter on the log scale.
/// Vector of length 1.
/// @param[in] s GEV shape parameter on the scale specified by `reparam_s`.
/// Vector of length 1.
/// @param[in] beta_a GP mean covariate coefficient vector of
/// length `n_covariate` for a.
/// @param[in] log_sigma_a GP covariance kernel variance
/// hyperparameter for a.
/// @param[in] log_ell_a GP covariance kernel range
/// hyperparameter for a.
template<class Type>
Type model_a_exp(objective_function<Type>* obj){
  using namespace density;
  using namespace R_inla;
  using namespace Eigen;
  using namespace SpatialGEV;

  // ------ Data inputs ------------
  DATA_VECTOR(y);
  DATA_IVECTOR(loc_ind);
  DATA_INTEGER(reparam_s);
  DATA_INTEGER(beta_prior);
  DATA_VECTOR(return_periods);
  int has_returns = return_periods(0) > Type(0.0);
  DATA_MATRIX(dist_mat);
  DATA_SCALAR(sp_thres);
  int n_loc = dist_mat.rows(); // number of spatial dimensions

  // Inputs for a
  DATA_MATRIX(design_mat_a);
  DATA_VECTOR(beta_a_prior);
  DATA_SCALAR(s_mean);
  DATA_SCALAR(s_sd);

  // ------------ Parameters ----------------------

  PARAMETER_VECTOR(a);
  PARAMETER_VECTOR(log_b);
  PARAMETER_VECTOR(s);

  PARAMETER_VECTOR(beta_a);
  PARAMETER(log_sigma_a);
  PARAMETER(log_ell_a);

  // Initialize the negative log likelihood
  Type nll = Type(0.0);

  // ---------- Likelihood contribution from a ------------------
  // GP latent layer
  vector<Type> mu_a = a -
    design_mat_a * beta_a;
  nll += nlpdf_gp_exp<Type>(mu_a, dist_mat,
				   exp(log_sigma_a),
				   exp(log_ell_a),
                                   sp_thres);
  // Priors
  nll += nlpdf_beta_prior<Type>(beta_a, beta_prior,
      beta_a_prior(0), beta_a_prior(1));
  // FIXME: rename this to not depend on `s`
  nll += nlpdf_s_prior<Type>(s(0), s_mean, s_sd);

  // ------------- Data layer -----------------
  for(int i=0;i<y.size();i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), a(loc_ind(i)), log_b(0),
	s(0), reparam_s);
  }

  // ------------- Output return levels -----------------------
  if(has_returns) {
    matrix<Type> return_levels(return_periods.size(), n_loc);
    for(int i=0; i<n_loc; i++) {
      gev_reparam_quantile<Type>(return_levels.col(i), return_periods,
                                 a(i), log_b(0), s(0), reparam_s);
    }
    ADREPORT(return_levels);
  }

  return nll;
}
#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif


