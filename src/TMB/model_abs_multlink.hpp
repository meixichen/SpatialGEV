#ifndef model_abs_multlink_hpp
#define model_abs_multlink_hpp

#include "SpatialGEV/utils.hpp"

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type model_abs_multlink(objective_function<Type>* obj){
  /*
  Model layer 1: y ~ GEV(a, b, s)
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
  DATA_VECTOR(y); // response vector: mws. Assumed to be > 0
  DATA_IVECTOR(loc_ind); // location index to which each observation in y is associated
  DATA_MATRIX(design_mat_psi); // n x r design matrix
  DATA_MATRIX(design_mat_tau); 
  DATA_MATRIX(design_mat_phi); 
  DATA_IVECTOR(meshidxloc); // indices of the locations in the mesh matrix
  DATA_SCALAR(nu); // Smoothness parameter for the Matern cov. 
  DATA_STRUCT(spde, spde_t); // take the returned object by INLA::inla.spde2.matern in R
  // Type of prior on beta. 1 is weakly informative normal prior and any other numbers 
  // mean noninformative uniform prior U(-inf, inf).
  DATA_INTEGER(beta_prior);   
  DATA_VECTOR(beta_psi_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_tau_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_VECTOR(beta_phi_prior); // length 2 vector containing mean and sd of normal prior on beta
  DATA_INTEGER(psi_pc_prior); // 1 for using PC prior on a, 0 for using flat unif prior 
  DATA_VECTOR(range_psi_prior); // length 2 vector (rho_0, p_rho) s.t. P(rho < rho_0) = p_rho
  DATA_VECTOR(sigma_psi_prior); // length 2 vector (sig_0, p_sig) s.t. P(sig > sig_0) = p_sig
  DATA_INTEGER(tau_pc_prior); 
  DATA_VECTOR(range_tau_prior);
  DATA_VECTOR(sigma_tau_prior);
  DATA_INTEGER(phi_pc_prior);  
  DATA_VECTOR(range_phi_prior);
  DATA_VECTOR(sigma_phi_prior);
  // parameter list
  PARAMETER_VECTOR(psi); // random effect to be integrated out. 
  PARAMETER_VECTOR(tau);  
  PARAMETER_VECTOR(phi); 
  PARAMETER_VECTOR(beta_psi); 
  PARAMETER_VECTOR(beta_tau);
  PARAMETER_VECTOR(beta_phi);
  PARAMETER(log_sigma_psi); 
  PARAMETER(log_kappa_psi);
  PARAMETER(log_sigma_tau); 
  PARAMETER(log_kappa_tau); 
  PARAMETER(log_sigma_phi);
  PARAMETER(log_kappa_phi);

  Type sigma_psi = exp(log_sigma_psi);
  Type kappa_psi = exp(log_kappa_psi);
  Type sigma_tau = exp(log_sigma_tau); 
  Type kappa_tau = exp(log_kappa_tau); 
  Type sigma_phi = exp(log_sigma_phi); 
  Type kappa_phi = exp(log_kappa_phi); 
 
  // Transform to GEV parameters
  Type xi00 = 0.;
  Type alp3 = 0.8;
  Type sig3 = (log(1 - pow(xi00+0.5, alp3)))*
    (1 - pow(xi00+0.5, alp3))*
    (-(1./alp3)*pow(xi00+0.5, (-alp3+1)));
  Type b3 = -sig3*(log(-log(1-pow(0.5, alp3))));
  vector<Type> a = exp(psi);
  vector<Type> log_b = tau + psi;
  vector<Type> s = pow(1. - exp(-exp( (phi - b3) / sig3)), 1./alp3) - .5;
  
  // calculate the negative log likelihood
  Type nll = Type(0.0);
  // data layer 
  int reparam_s = 3; // unconstrained s
  nll += nll_accumulator_abs<Type>(y, loc_ind, a, log_b, s, reparam_s);
  // GP latent layer
  vector<Type> mu_psi = psi - design_mat_psi * beta_psi;
  vector<Type> mu_tau = tau - design_mat_tau * beta_tau;
  vector<Type> mu_phi = phi - design_mat_phi * beta_phi;
  nll += gp_spde_nlpdf<Type>(mu_psi, spde, sigma_psi, kappa_psi, nu);
  nll += gp_spde_nlpdf<Type>(mu_tau, spde, sigma_tau, kappa_tau, nu);
  nll += gp_spde_nlpdf<Type>(mu_phi, spde, sigma_phi, kappa_phi, nu);
  // prior
  nll += nll_accumulator_beta<Type>(beta_psi, beta_prior, beta_psi_prior[0], beta_psi_prior[1]);
  nll += nll_accumulator_beta<Type>(beta_tau, beta_prior, beta_tau_prior[0], beta_tau_prior[1]);
  nll += nll_accumulator_beta<Type>(beta_phi, beta_prior, beta_phi_prior[0], beta_phi_prior[1]);
  nll += nll_accumulator_matern_hyperpar<Type>(log_kappa_psi, log_sigma_psi, psi_pc_prior,
                                        nu, range_psi_prior, sigma_psi_prior);
  nll += nll_accumulator_matern_hyperpar<Type>(log_kappa_tau, log_sigma_tau, tau_pc_prior,
                                        nu, range_tau_prior, sigma_tau_prior);
  nll += nll_accumulator_matern_hyperpar<Type>(log_kappa_phi, log_sigma_phi, phi_pc_prior,
                                        nu, range_phi_prior, sigma_phi_prior);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

