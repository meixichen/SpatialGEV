#ifndef SpatialGEV_model_gev_hpp
#define SpatialGEV_model_gev_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// @param[in] prior_s Vector of length two specifying the mean and sd of a normal prior on the shape parameter `s.`  For `sd > 9999` a flat prior is used.
template<class Type>
Type model_gev(objective_function<Type>* obj) {
  using namespace SpatialGEV;

  // data inputs
  DATA_VECTOR(y);
  DATA_INTEGER(reparam_s);
  DATA_VECTOR(s_prior);

  // parameter list
  PARAMETER(a);  
  PARAMETER(log_b);
  PARAMETER(s);

  Type nll = Type(0.0);

  for(int i=0; i<y.size(); i++) {
    nll -= gev_reparam_lpdf<Type>(y(i), a, log_b, s,
				  reparam_s);
  }
  nll += nlpdf_s_prior<Type>(s, s_prior(0), s_prior(1));


  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

