#include "Matern.hpp"

namespace Covariance {

mat Matern::distances (void)
  /* sum_k \sqrt{2v}/p_k d_{i,j,k} */
{
  mat d = arma::zeros<mat>(dim, dim);
  for (uword i=0; i<npars; ++i)
    d += D.slice(i) * sqrt(2*rate) / range(i);
  return d;
}

mat Matern::besselfn (void)
  /* K_v (sum_k \sqrt{2v}/p_k d_{i,j,k}) */
{
  mat d = arma::zeros<mat>(dim, dim);
  for (uword i=0; i<dim-1; ++i)
    for (uword j=i+1; j<dim; ++j)
      d.at(i,j) = cyl_bessel_k(rate, Dk.at(i,j));
  return d;
}

mat Matern::covariance (void)
{
  double constant = stddev*stddev * exp(-lgamma(rate) + (1-rate)*log(2.));
  mat    d        = arma::eye<mat>(dim, dim) * stddev*stddev + nugget*nugget;
  for (uword i=0; i<dim-1; ++i)
    for (uword j=i+1; j<dim; ++j)
      d.at(i,j)  = d.at(j,i) 
                 = pow(Dk.at(i,j), rate) * Bk.at(i,j) * constant;
  return d;
}

mat Matern::differential (void) 
{
  double constant;
  arma::mat  I = arma::eye<mat>(dim, dim) * pow(nugget, 2);
  arma::cube d (dim, dim, 1 + npars, arma::fill::zeros);

  /* dC/ds */
  d.slice(0) = 2./stddev * (C - I); 

  /* dC/dp */ 
  for (uword i=0; i<dim-1; ++i)
    for (uword j=i+1; j<dim; ++j)
    {
      // Simplifying things:
      // stuff * x^v B_k(x) * (2*v/x - B_{k+1}(x)/B_{k}(x)) =
      //  stuff * x^v * B_k(x) * 2*v/x - stuff * x^v * B_{k+1}(x)

      constant = pow(stddev, 2) * exp(-lgamma(rate) + (1-rate)*log(2.)) * pow(Dk.at(i,j), rate);
      constant = 2 * C.at(i,j) * rate / Dk.at(i,j) - constant * cyl_bessel_k(rate + 1., Dk.at(i,j));

      for (uword k=0; k<npars; ++k)
        d.at(i,j,1+k) = d.at(j,i,1+k) 
                      = -pow(range(k),-2) * sqrt(2*rate) * constant * D.at(i,j,k);
    }

  /* as Jacobian matrix */
  mat out (dim * dim, npars+1);
  for (uword k=0; k<d.n_slices; ++k)
    out.col(k) = arma::vectorise(d.slice(k));

  return out;
}

cube Matern::differential_dist (void)
{
  double constant;
  arma::cube d (dim, dim, npars, arma::fill::zeros);

  /* dC/dD */  // TODO: consolidate with Matern::differential
  for (uword i=0; i<dim-1; ++i)
    for (uword j=i+1; j<dim; ++j)
    {
      // We can have underflow in the Bessel functions here. The most stable way
      // to defend against this is to directly approximate the ratio K_{v+1}(x)/K_{v}(x)
      // for small x. Asymptotics give:
      //    K_m(x)/K_n(x) = 1 + (m-n)(m+n)/(2x) + (m-n)(m+n)(m^2 - n^2 - 2)/(8x^2) + O(1/x^3)
      // this should be fine for, say, x > 100.

      constant = pow(stddev, 2) * exp(-lgamma(rate) + (1-rate)*log(2.)) * pow(Dk.at(i,j), rate);
      constant = 2 * C.at(i,j) * rate / Dk.at(i,j) - constant * cyl_bessel_k(rate + 1., Dk.at(i,j));

      for (uword k=0; k<npars; ++k)
        d.at(i,j,k) = d.at(j,i,k) 
                    = constant * sqrt(2*rate)/range(k);
    }

  return d;
}

vec Matern::parameters (void)
{
  arma::vec out = {stddev};
  return arma::join_vert(out, range);
}

} /* namespace Covariance */

//------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Matern_C")]]
arma::mat test_Matern_C (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Matern cov (D, nu, delta, pars);
  return cov.C;
}

// [[Rcpp::export("inlassle_test_Matern_dC_dt")]]
arma::mat test_Matern_dC_dt (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Matern cov (D, nu, delta, pars);
  return cov.dC_dt;
}

// [[Rcpp::export("inlassle_test_Matern_dC_dD")]]
arma::cube test_Matern_dC_dD (arma::cube D, double nu, double delta, arma::vec pars)
{
  Covariance::Matern cov (D, nu, delta, pars);
  return cov.dC_dD;
}

