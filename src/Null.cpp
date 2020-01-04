#include "Null.hpp"

namespace Covariance {

mat Null::covariance (void)
{
  mat d = arma::eye<mat>(dim, dim) * std::pow(nugget,2);
  return d;
}

mat Null::differential (void) 
{
  mat d (dim * dim, npars+1, arma::fill::zeros);
  return d;
}

cube Null::differential_dist (void)
{
  cube d (dim, dim, npars, arma::fill::zeros);
  return d;
}

vec Null::parameters (void)
{
  arma::vec out = {stddev};
  return arma::join_vert(out, weight);
}

} /* namespace Covariance */

//------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Null_C")]]
arma::mat test_Null_C (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Null cov (D, nu, delta, pars);
  return cov.C;
}

// [[Rcpp::export("inlassle_test_Null_dC_dt")]]
arma::mat test_Null_dC_dt (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Null cov (D, nu, delta, pars);
  return cov.dC_dt;
}

// [[Rcpp::export("inlassle_test_Null_dC_dD")]]
arma::cube test_Null_dC_dD (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Null cov (D, nu, delta, pars);
  return cov.dC_dD;
}

