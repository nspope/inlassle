#include "Weighted.hpp"

namespace Covariance {

mat Weighted::covariance (void)
{
  mat d = arma::eye<mat>(dim, dim) * std::pow(nugget,2);
  for (uword i=0; i<npars; ++i)
    d += std::pow(stddev,2) * weight(i) * D.slice(i);
  return d;
}

mat Weighted::differential (void) 
{
  mat d (dim * dim, npars+1);
  mat I = arma::eye<mat>(dim, dim) * std::pow(nugget,2);

  /* dC/ds */
  d.col(0) = 2./stddev * arma::vectorise(C - I); 

  /* dC/dp */ 
  for (uword i=0; i<npars; ++i)
    d.col(i+1) = std::pow(stddev,2) * arma::vectorise(D.slice(i));

  return d;
}

cube Weighted::differential_dist (void)
{
  cube d (dim, dim, npars, arma::fill::zeros);

  /* dC/dD */  // TODO: consolidate with Weighted::differential
  for (uword i=0; i<npars; ++i)
    d.slice(i).fill(weight(i) * std::pow(stddev, 2));

  return d;
}

vec Weighted::parameters (void)
{
  arma::vec out = {stddev};
  return arma::join_vert(out, weight);
}

} /* namespace Covariance */

//------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Weighted_C")]]
arma::mat test_Weighted_C (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Weighted cov (D, nu, delta, pars);
  return cov.C;
}

// [[Rcpp::export("inlassle_test_Weighted_dC_dt")]]
arma::mat test_Weighted_dC_dt (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Weighted cov (D, nu, delta, pars);
  return cov.dC_dt;
}

// [[Rcpp::export("inlassle_test_Weighted_dC_dD")]]
arma::cube test_Weighted_dC_dD (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Weighted cov (D, nu, delta, pars);
  return cov.dC_dD;
}

