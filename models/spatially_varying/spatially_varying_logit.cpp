//
  // Author: Anita Nandi
// Date: 2019-02-14

// Data: Spatial field mesh and matrices, polygon data, covariate pixel data


#define TMB_LIB_INIT R_init_disaggregation
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
    // Polygon level data
  // ------------------------------------------------------------------------ //
    
    // Covariate pixel data
  DATA_MATRIX(x);
  
  // Shape data. Cases and region id.
  DATA_VECTOR(response_data);
  DATA_VECTOR(response_sample_size);
  

  // ------------------------------------------------------------------------ //
    // Parameters
  // ------------------------------------------------------------------------ //
    
  PARAMETER(intercept);
  PARAMETER_VECTOR(slope);
  
  DATA_SCALAR(priormean_intercept);
  DATA_SCALAR(priorsd_intercept);
  DATA_SCALAR(priormean_slope);
  DATA_SCALAR(priorsd_slope);
  
  // Priors for likelihood
  PARAMETER(log_tau_gaussian);
  Type tau_gaussian = exp(log_tau_gaussian);
  Type gaussian_sd = 1 / sqrt(tau_gaussian);
  

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
  
  // Model component flags
  DATA_INTEGER(field);

  // Number of pixels
  int n_pixels = x.rows();
  
  Type nll = 0.0;
  
  // ------------------------------------------------------------------------ //
    // Likelihood from priors
  // ------------------------------------------------------------------------ //
    
    nll -= dnorm(intercept, priormean_intercept, priorsd_intercept, true);
  for (int s = 0; s < slope.size(); s++) {
    nll -= dnorm(slope[s], priormean_slope, priorsd_slope, true);
  }
  
  if(field) {
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
  }
  
  Type nll_priors = nll;
  
  // ------------------------------------------------------------------------ //
    // Likelihood from data
  // ------------------------------------------------------------------------ //
    
  vector<Type> pixel_linear_pred(n_pixels);
  pixel_linear_pred = intercept + x * slope;
  pixel_linear_pred = invlogit(pixel_linear_pred);
  
  
  if(field) {
    // Calculate field for pixel data
    vector<Type> linear_pred_field(n_pixels);
    linear_pred_field = Apixel * nodemean;
    pixel_linear_pred += linear_pred_field.array();
  }
  

  Type response;
  vector<Type> pixel_pred;
  vector<Type> reportprediction_rate(n_pixels);
  vector<Type> reportnll(n_pixels);

  // For each shape get pixel predictions within and aggregate to polygon level
  for (int polygon = 0; polygon < n_pixels; polygon++) {
    
    nll -= dbinom(response_data[polygon], response_sample_size[polygon], pixel_pred[polygon], true);
    reportnll[polygon] = -dbinom(response_data[polygon], response_sample_size[polygon], pixel_pred[polygon], true);
    
    
  }
  
  REPORT(reportprediction_rate);
  REPORT(reportnll);
  REPORT(response_data);
  REPORT(nll_priors);
  REPORT(nll);

  return nll;
}

