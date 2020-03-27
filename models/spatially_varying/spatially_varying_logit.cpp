//
// Author: Anita Nandi
// Date: 2019-07-24
//
// Simple example of a TMB model - linear regression with field
//
// Data: prevalence survey data and covariate data
//
// The model: prev = invlogit(intercept + x * slope + logit_prevalence_field.array())
//

#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()()
{
  
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;
  
  // ------------------------------------------------------------------------ //
  // Spatial field data
  // ------------------------------------------------------------------------ //
  
  // The A matrices are for projecting the mesh to a point for the pixel and point data respectively.
  DATA_SPARSE_MATRIX(Apixel);
  DATA_STRUCT(spde, spde_t);
  
  // ------------------------------------------------------------------------ //
  // Input data
  // ------------------------------------------------------------------------ //
  
  // Covariates: Environmental data
  DATA_MATRIX(x);
  
  // Response: Prevalence point data
  DATA_VECTOR(positive_cases);
  DATA_VECTOR(examined_cases);
  
  // ------------------------------------------------------------------------ //
  // Parameters
  // ------------------------------------------------------------------------ //
  
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  Type priormean_intercept = -4.0;
  Type priorsd_intercept = 2.0; //priormean_intercept from data entry
  Type priormean_slope = 0.0;
  Type priorsd_slope = 0.5;
  
  // spde hyperparameters
  PARAMETER(log_kappa);
  PARAMETER(log_tau);
  
  // Priors on spde hyperparameters
  Type priormean_log_kappa = -2;
  Type priorsd_log_kappa = 0.1;
  Type priormean_log_tau = -5;
  Type priorsd_log_tau = 0.5;
  
  // Convert hyperparameters to natural scale
  Type tau = exp(log_tau);
  Type kappa = exp(log_kappa);
  
  PARAMETER_VECTOR(nodemean);
  
  // Number of prev points
  int n_points = positive_cases.size();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  // Likelihood of hyperparameters for field
  nll -= dnorm(log_kappa, priormean_log_kappa, priorsd_log_kappa, true);
  nll -= dnorm(log_tau, priormean_log_tau, priorsd_log_tau, true);
  
  // Build spde matrix
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  // Likelihood of the random field.
  nll += SCALE(GMRF(Q), 1.0 / tau)(nodemean);
  
  Type nllpriors = nll;
  
  // ------------------------------------------------------------------------ //
  // Calculate random field effects
  // ------------------------------------------------------------------------ //
  
  // Calculate field for pixel data
  vector<Type> logit_prevalence_field;
  logit_prevalence_field = Apixel * nodemean;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from data
  // ------------------------------------------------------------------------ //
  
  vector<Type> pixel_linear_pred(n_points);
  vector<Type> pixel_pred(n_points);
  vector<Type> reportnll(n_points);
  
  pixel_linear_pred = intercept + x * slope  + logit_prevalence_field.array();
  pixel_pred = invlogit(pixel_linear_pred);
  
  for (int point = 0; point < n_points; point++) {
    nll -= dbinom(positive_cases[point], examined_cases[point], pixel_pred[point], true);
    reportnll[point] = -dbinom(positive_cases[point], examined_cases[point], pixel_pred[point], true);
  }
  
  REPORT(pixel_pred);
  REPORT(reportnll);
  REPORT(positive_cases);
  REPORT(examined_cases);
  REPORT(nllpriors);
  REPORT(nll);
  
  return nll;
}
