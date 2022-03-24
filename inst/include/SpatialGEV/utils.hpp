/// @file utils.hpp
///
/// @brief Utilities for `SpatialGEV`.

#ifndef SPATIALGEV_UTILS_HPP
#define SPATIALGEV_UTILS_HPP

namespace SpatialGEV {

  /// @typedef
  /// @brief Standard typedefs for arguments to Eigen functions.
  template <class Type>
  using RefMatrix_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  template <class Type>
  using cRefMatrix_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> >;
  template <class Type>
  using RefVector_t = Eigen::Ref <Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using cRefVector_t = const Eigen::Ref <const Eigen::Matrix<Type, Eigen::Dynamic, 1> >;
  template <class Type>
  using RefRowVector_t = Eigen::Ref <Eigen::Matrix<Type, 1, Eigen::Dynamic> >;
  template <class Type>
  using cRefRowVector_t = const Eigen::Ref <const Eigen::Matrix<Type, 1, Eigen::Dynamic> >;
  
  /// Calculates the log-density of the Gumbel distribution.
  ///
  /// @param[in] x Argument to the density.
  /// @param[in] a Location parameter.
  /// @param[in] log_b Log of scale parameter.
  ///
  /// @return Log-density of Gumbel distribution evaluated at its inputs.
  template <class Type>
  Type gumbel_lpdf(Type x, Type a, Type log_b) {
    Type t = (x - a) / exp(log_b);
    return -exp(-Type(1.0) * t) - t - log_b;
  }
  
  /// Calculates the log-density of the GEV distribution.
  ///
  /// @param[in] x Argument to the density.
  /// @param[in] a Location parameter.
  /// @param[in] log_b Log of scale parameter.
  /// @param[in] s Shape parameter.
  ///
  /// @return Log-density of GEV distribution evaluated at its inputs.
  template <class Type>
  Type gev_lpdf(Type x, Type a, Type log_b, Type s) {
    Type log_t = log(Type(1.0) + s * (x - a) / exp(log_b));
    return -exp(-Type(1.0) * log_t/s) - (s + Type(1.0))/s * log_t - log_b;
    // return pow(t - Type(1.0)/s) + (s + Type(1.0))/s + log(t);
  }
 
  /// Compute the exponential kernel function
  ///
  /// @param[in] x Value to evaluate at
  /// @param[in] sigma Amplitude parameter
  /// @param[in] ell Smoothness parameter
  ///
  /// @return A scalar of squared exponential kernel function value
  template <class Type>
  Type kernel_exp(Type x, Type sigma, Type ell) {
    return sigma * exp(- x / ell);
  }

  /// Compute the variance matrix for the exponential kernel.
  ///
  /// @param[out] cov Matrix into which to store the output.
  /// @param[in] dd Distance matrix.
  /// @param[in] sigma Scale parameter.
  /// @param[in] ell Length parameter.
  /// @param[in] sp_thres Threshold parameter.
  template <class Type>
  void cov_expo(RefMatrix_t<Type> cov, cRefMatrix_t<Type>& dd,
  	       Type sigma, Type ell, Type sp_thres) {
    int i,j;
    int n = dd.rows();
    if (sp_thres == -1){
      cov = - dd / ell;
      cov = cov.array().exp();
      cov *= sigma;
    }else {
      for (i = 0; i < n; i++){
	cov(i,i) = sigma;
	for (j = 0; j < i; j++){
	  if (dd(i,j) >= sp_thres) {
	    cov(i,j) = 0;
	    cov(j,i) = 0;
	  } else {
	    cov(i,j) = kernel_exp(dd(i,j), sigma, ell);  
	    cov(j,i) = cov(i,j);
	  }
	}
      }
    }
    return;
  }

  /// Compute the variance matrix for the matern kernel.
  ///
  /// @param[out] cov Matrix into which to store the output.
  /// @param[in] dd Distance matrix.
  /// @param[in] sigma Hyperparameter of the Matern. It is in fact sigma^2.
  /// @param[in] kappa Hyperparameter of the Matern. Positive.
  /// @param[in] nu Smoothness parameter of the Matern.
  /// @param[in] sp_thres Threshold parameter.
  template <class Type>
  void cov_matern(RefMatrix_t<Type> cov, cRefMatrix_t<Type>& dd,
  	       Type sigma, Type kappa, Type nu, Type sp_thres) {
    int i,j;
    int n = dd.rows();
    if (sp_thres == -1){
      for (i = 0; i < n; i++){
	cov(i,i) = sigma;
	for (j = 0; j < i; j++){
	    cov(i,j) = sigma*matern(dd(i,j), 1/kappa, nu);  
	    cov(j,i) = cov(i,j);
	}
      }
    }else {
      for (i = 0; i < n; i++){
	cov(i,i) = sigma;
	for (j = 0; j < i; j++){
	  if (dd(i,j) >= sp_thres) {
	    cov(i,j) = 0;
	    cov(j,i) = 0;
	  } else {
	    cov(i,j) = sigma*matern(dd(i,j), 1/kappa, nu);  
	    cov(j,i) = cov(i,j);
	  }
	}
      }
    } //end else
    return;
  }

  /// Add negative log-likelihood contributed by prior on beta
  ///
  /// @param[out] nll Negative log-likelihood accumulator.
  /// @param[in] beta Beta parameter in the model.
  /// @param[in] prior Type of prior. 1 is weakly informative normal prior and any other numbers 
  /// mean noninformative prior.
  /// @param[in] mean Mean of the normal prior. Only relevant if prior=1.
  /// @param[in] sd Standard deviation of the normal prior. Only relevant if prior=1.
  template <class Type>
  void nll_accumulator_beta(Type &nll, RefVector_t<Type> beta, Type prior,
                            Type mean, Type sd) {
    if (prior == 1) { 
      for(int i=0; i< (beta.size()); i++){
        nll -= dnorm(beta[i], mean, sd, 1);
      }
    }
    return;
  }

  
  /// Add negative log-likelihood contributed by prior on Matern hyperparameters
  ///
  /// @param[out] nll Negative log-likelihood accumulator.
  /// @param[in] log_kappa Log of inverse range parameter of Matern
  /// @param[in] log_sigma Log of marginal variance parameter of Matern
  /// @param[in] prior Type of prior. 1 is weakly penalized complexity (PC) prior and any other
  /// number means noninformative prior.
  /// @param[in] nu Matern smoothness hyperparameter.
  /// @param[in] range_prior. Length 2 vector (rho_0, p_rho) s.t. P(rho < rho_0) = p_rho. 
  /// Only relevant if prior=1.
  /// @param[in] sigma_prior. Length 2 vector (sig_0, p_sig) s.t. P(sig > sig_0) = p_sig.  
  /// Only relevant if prior=1.
  template <class Type>
  void nll_accumulator_matern_hyperpar(Type &nll, Type log_kappa, Type log_sigma,
                                       Type prior, Type nu,                
                                       cRefVector_t<Type> range_prior, 
				       cRefVector_t<Type> sigma_prior) {
    if (prior == 1) {
       // See Theorem 6 of Fuglstad et al. (2017) https://arxiv.org/pdf/1503.00256.pdf 
       Type log_rho = 0.5*log(8.0*nu) - log_kappa; // get range parameter
       Type rho = exp(log_rho);
       Type log_sig = 0.5*log_sigma; // get sd parameter
       Type sig = exp(log_sig);
       Type rho_0 = range_prior[0];
       Type p_rho = range_prior[1];
       Type sig_0 = sigma_prior[0];
       Type p_sig = sigma_prior[1];
       Type lam1 = -1.0 * log(p_rho) * rho_0;
       Type lam2 = -1.0 * log(p_sig) / sig_0; 
       // PC prior log density 
       Type logpi = log(lam1) + log(lam2) - 2.0 * log_rho - lam1 / rho - lam2 * sig; 
       // Jacobian adjustment = log(1 / (|dlogkappa/drho| * |dlogsigma/dsig|))
       logpi += log_kappa + 0.5 * log_sigma - log(2.0) - 0.5*log(8.0*nu); 
       nll -= logpi;
    }
    return;
  }

  /// Add negative log-likelihood contributed by the data layer for model_a.
  ///
  /// @param[out] nll negative log-likelihood accumulator.
  /// @param[in] y Data.
  /// @param[in] n_obs Vector of number of observations per location 
  /// @param[in] a GEV Location parameter vector.
  /// @param[in] log_b GEV (log) scale parameter.
  /// @param[in] s GEV Shape parameter (possibly transformed).
  /// @param[in] n Number of locations.
  /// @param[in] reparam_s Flag indicating reparametrization of s
  /// @param[in] s_mean Mean of s prior distn.
  /// @param[in] s_sd SD of s prior distn.
  template <class Type>
  void nll_accumulator_a(Type &nll, cRefVector_t<Type>& y, vector<int> n_obs, 
		      RefVector_t<Type> a, Type log_b, Type s,
		      Type n, Type reparam_s, Type s_mean, Type s_sd) {
    int start_ind = 0; // index of the first observation of location i in n_obs
    int end_ind = 0; // index of the last observation of location i in n_obs
    if (reparam_s == 0){ // this is the case we are using Gumbel distribution
      for(int i=0;i<n;i++) {
	end_ind += n_obs[i]; 
	for (int j=start_ind;j<end_ind;j++){
          nll -= gumbel_lpdf<Type>(y[j], a[i], log_b);
	}
        start_ind += n_obs[i];	
      }
    } else{ // the case where we are using GEV distribution with nonzerio shape parameter
      if (s_sd<9999){ // put a prior on s, or log(s), or log(|s|)
	nll -= dnorm(s, s_mean, s_sd, true);
      }
      if (reparam_s == 1){ // if we have stated that s is constrained to be positive, this implies that we are optimizing log(s)
	s = exp(s);
      } else if (reparam_s == 2){ // if we have stated that s is constrained to be negative, this implies that we are optimizing log(-s)
	s = -exp(s);
      } // if we don't use any reparametrization, then s is unconstrained
      for(int i=0;i<n;i++) {
	end_ind += n_obs[i];
	for (int j=start_ind;j<end_ind;j++){
          nll -= gev_lpdf<Type>(y[j], a[i], log_b, s);
	} 
	start_ind += n_obs[i];
      }
    } // end else
    return;
  }

  /// Add negative log-likelihood contributed by the data layer for model_ab.
  ///
  /// @param[out] nll negative log-likelihood accumulator.
  /// @param[in] y Data.
  /// @param[in] n_obs Vector of number of observations per location 
  /// @param[in] a GEV location parameter vector.
  /// @param[in] log_b GEV (log) scale parameter vector.
  /// @param[in] s GEV shape parameter (possibly transformed).
  /// @param[in] n Number of locations.
  /// @param[in] reparam_s Flag indicating reparametrization of s
  /// @param[in] s_mean Mean of s prior distn.
  /// @param[in] s_sd SD of s prior distn.
  template <class Type>
  void nll_accumulator_ab(Type &nll, cRefVector_t<Type>& y, vector<int> n_obs, 
		      RefVector_t<Type> a, RefVector_t<Type> log_b, Type s,
		      Type n, Type reparam_s, Type s_mean, Type s_sd) {
    int start_ind = 0; // index of the first observation of location i in n_obs
    int end_ind = 0; // index of the last observation of location i in n_obs
    if (reparam_s == 0){ // this is the case we are using Gumbel distribution
      for(int i=0;i<n;i++) {
	end_ind += n_obs[i]; 
	for (int j=start_ind;j<end_ind;j++){
          nll -= gumbel_lpdf<Type>(y[j], a[i], log_b[i]);
	}
        start_ind += n_obs[i];	
      }
    } else{ // the case where we are using GEV distribution with nonzerio shape parameter
      if (s_sd<9999){ // put a prior on s, or log(s), or log(|s|)
	nll -= dnorm(s, s_mean, s_sd, true);
      }
      if (reparam_s == 1){ // if we have stated that s is constrained to be positive, this implies that we are optimizing log(s)
	s = exp(s);
      } else if (reparam_s == 2){ // if we have stated that s is constrained to be negative, this implies that we are optimizing log(-s)
	s = -exp(s);
      } // if we don't use any reparametrization, then s is unconstrained
      for(int i=0;i<n;i++) {
	end_ind += n_obs[i];
	for (int j=start_ind;j<end_ind;j++){
          nll -= gev_lpdf<Type>(y[j], a[i], log_b[i], s);
	} 
	start_ind += n_obs[i];
      }
    } // end else
    return;
  }

  /// Add negative log-likelihood contributed by the data layer for model_abs.
  ///
  /// @param[out] nll negative log-likelihood accumulator.
  /// @param[in] y Data.
  /// @param[in] n_obs Vector of number of observations per location 
  /// @param[in] a GEV location parameter vector.
  /// @param[in] log_b GEV (log) scale parameter vector.
  /// @param[in] s GEV shape parameter (possibly transformed).
  /// @param[in] n Number of locations.
  /// @param[in] reparam_s Flag indicating reparametrization of s
  template <class Type>
  void nll_accumulator_abs(Type &nll, cRefVector_t<Type>& y, vector<int> n_obs, 
		      RefVector_t<Type> a, RefVector_t<Type> log_b, RefVector_t<Type> s,
		      Type n, Type reparam_s) {
    int start_ind = 0; // index of the first observation of location i in n_obs
    int end_ind = 0; // index of the last observation of location i in n_obs
    if (reparam_s == 1){ // if we have stated that s is constrained to be positive, this implies that we are optimizing log(s)
      s = s.array().exp();
    } else if (reparam_s == 2){ // if we have stated that s is constrained to be negative, this implies that we are optimizing log(-s)
      s = -s.array().exp();
    } // if we don't use any reparametrization, then s is unconstrained
    for(int i=0;i<n;i++) {
      end_ind += n_obs[i];
      for (int j=start_ind;j<end_ind;j++){
	nll -= gev_lpdf<Type>(y[j], a[i], log_b[i], s[i]);
      } 
      start_ind += n_obs[i];
    }
    return;
  } 
} // end namespace SpatialGEV

#endif
