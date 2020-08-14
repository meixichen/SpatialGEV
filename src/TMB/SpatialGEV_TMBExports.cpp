// Generated by TMBtools: do not edit by hand

#define TMB_LIB_INIT R_init_SpatialGEV_TMBExports
#include <TMB.hpp>
#include "model_a.hpp"
#include "model_ab.hpp"
#include "model_abs.hpp"
#include "model_as.hpp"
#include "model_b.hpp"
#include "model_bs.hpp"
#include "model_s.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "model_a") {
    return model_a(this);
  } else if(model == "model_ab") {
    return model_ab(this);
  } else if(model == "model_abs") {
    return model_abs(this);
  } else if(model == "model_as") {
    return model_as(this);
  } else if(model == "model_b") {
    return model_b(this);
  } else if(model == "model_bs") {
    return model_bs(this);
  } else if(model == "model_s") {
    return model_s(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
