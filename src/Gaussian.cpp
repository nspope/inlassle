#include "Gaussian.hpp"

//------------------------------------------- implementation

namespace Covariance {

mat Gaussian::distances (void)
{
  mat d = arma::zeros<mat>(dim, dim);
  for (uword i=0; i<npars; ++i)
    d += D.slice(i) / range(i);
  return d;
}

mat Gaussian::covariance (void)
{
  return pow(stddev, 2) * arma::exp(-arma::pow(Dk, rate)) + arma::eye<mat>(dim,dim)*pow(nugget,2);
}

mat Gaussian::differential (void)
{
  arma::mat  I = arma::eye<mat>(dim, dim) * pow(nugget, 2);
  arma::cube d (dim, dim, 1 + npars);
  d.slice(0) = 2./stddev * (C - I); // dC/dstddev
  for (uword i=0; i<npars; ++i)
  {
    d.slice(1+i) = pow(range(i), -2) * rate * D.slice(i) % C % arma::pow(Dk, rate - 1); // dC/range(i) 
    d.slice(1+i).diag().zeros();
  }

  // make Jacobian
  mat out (dim * dim, npars+1);
  for (uword k=0; k<d.n_slices; ++k)
    out.col(k) = arma::vectorise(d.slice(k));

  return out;

//  d.slice(0) = -C % arma::pow(Dk, rate) % arma::log(Dk + arma::eye<mat>(dim,dim)); // dC/drate: needing if rate is ever treated as not fixed
}

cube Gaussian::differential_dist (void) 
{
  arma::cube d (dim, dim, npars);
  for (uword i=0; i<npars; ++i)
  {
    d.slice(i) = -rate/range(i) * C % arma::pow(Dk, rate - 1); // dC/D(i)
    d.slice(i).diag().zeros();
  }
  return d;
}

vec Gaussian::parameters (void)
{
  arma::vec out = {stddev};
  return arma::join_vert(out, range);
}

} // namespace Covariance

//------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Gaussian_C")]]
arma::mat test_Gaussian_C (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Gaussian cov (D, nu, delta, pars);
  return cov.C;
}

// [[Rcpp::export("inlassle_test_Gaussian_dC_dt")]]
arma::mat test_Gaussian_dC_dt (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Gaussian cov (D, nu, delta, pars);
  return cov.dC_dt;
}

// [[Rcpp::export("inlassle_test_Gaussian_dC_dD")]]
arma::cube test_Gaussian_dC_dD (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Gaussian cov (D, nu, delta, pars);
  return cov.dC_dD;
}
