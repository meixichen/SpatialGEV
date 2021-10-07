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
    if (sp_thres == 0){
      cov = -dd/ell;
      cov = cov.array().exp();
      cov *= sigma;
    } else {
      // construct the covariance matrix
      for (i = 0; i < n; i++){
        cov(i,i) = sigma;
        for (j = 0; j < i; j++){
          if (dd(i,j) >= sp_thres) {
            cov(i,j) = 0;
            cov(j,i) = 0;
          } else {
            cov(i,j) = sigma*exp(-dd(i,j)/ell);  
            cov(j,i) = cov(i,j);
          }
        }
      }
    }
    return;
  }

} // end namespace SpatialGEV

#endif
