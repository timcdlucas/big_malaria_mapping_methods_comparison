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
  
  // rgression priors
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  
  // spde hyperparameters
  PARAMETER(log_sigma);
  PARAMETER(log_rho);
  Type sigma = exp(log_sigma);
  Type rho = exp(log_rho);
  
  // Priors on spde hyperparameters
  DATA_SCALAR(prior_rho_min);
  DATA_SCALAR(prior_rho_prob);
  DATA_SCALAR(prior_sigma_max);
  DATA_SCALAR(prior_sigma_prob);
  
  // Convert hyperparameters to natural scale
  DATA_SCALAR(nu);
  Type kappa = sqrt(8.0) / rho;
  
  // Random effect parameters
  PARAMETER_VECTOR(nodemean);
  
  // Number of pixels
  int n_points = x.rows();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
  // Likelihood from priors
  // ------------------------------------------------------------------------ //
  
  nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }

  // Likelihood of hyperparameters for field. 
  // From https://www.tandfonline.com/doi/full/10.1080/01621459.2017.1415907 (Theorem 2.6)
  Type lambdatilde1 = -log(prior_rho_prob) * prior_rho_min;
  Type lambdatilde2 = -log(prior_sigma_prob) / prior_sigma_max;
  Type log_pcdensity = log(lambdatilde1) + log(lambdatilde2) - 2*log_rho - lambdatilde1 * pow(rho, -1) - lambdatilde2 * sigma;
  // log_rho and log_sigma from the Jacobian
  nll -= log_pcdensity + log_rho + log_sigma;
  
  // Build spde matrix
  SparseMatrix<Type> Q = Q_spde(spde, kappa);
  
  // From Lindgren (2011) https://doi.org/10.1111/j.1467-9868.2011.00777.x, see equation for the marginal variance
  Type scaling_factor = sqrt(exp(lgamma(nu)) / (exp(lgamma(nu + 1)) * 4 * M_PI * pow(kappa, 2*nu)));
  
  // Likelihood of the random field.
  nll += SCALE(GMRF(Q), sigma / scaling_factor)(nodemean);
  
  Type nll_priors = nll;
  
  
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
  REPORT(nll_priors);
  REPORT(nll);
  
  return nll;
}
