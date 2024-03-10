// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_SpatialGEV_TMBExports
#include <TMB.hpp>
#include "model_a_exp.hpp"
#include "model_a_matern.hpp"
#include "model_a_spde.hpp"
#include "model_ab_exp.hpp"
#include "model_ab_matern.hpp"
#include "model_ab_spde.hpp"
#include "model_abs_exp.hpp"
#include "model_abs_matern.hpp"
#include "model_abs_spde_maxsmooth.hpp"
#include "model_abs_spde.hpp"
#include "model_gev_spde.hpp"
#include "model_gev.hpp"
#include "model_ptp_spde.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "model_a_exp") {
    return model_a_exp(this);
  } else if(model == "model_a_matern") {
    return model_a_matern(this);
  } else if(model == "model_a_spde") {
    return model_a_spde(this);
  } else if(model == "model_ab_exp") {
    return model_ab_exp(this);
  } else if(model == "model_ab_matern") {
    return model_ab_matern(this);
  } else if(model == "model_ab_spde") {
    return model_ab_spde(this);
  } else if(model == "model_abs_exp") {
    return model_abs_exp(this);
  } else if(model == "model_abs_matern") {
    return model_abs_matern(this);
  } else if(model == "model_abs_spde_maxsmooth") {
    return model_abs_spde_maxsmooth(this);
  } else if(model == "model_abs_spde") {
    return model_abs_spde(this);
  } else if(model == "model_gev_spde") {
    return model_gev_spde(this);
  } else if(model == "model_gev") {
    return model_gev(this);
  } else if(model == "model_ptp_spde") {
    return model_ptp_spde(this);
  } else {
    Rf_error("Unknown model.");
  }
  return 0;
}
