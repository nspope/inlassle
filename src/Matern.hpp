#ifndef MATERN_H
#define MATERN_H

#include "Inlassle.hpp"

namespace Covariance {

struct Matern
 /* Matern covariance structure,
  *   \delta_{i,j} n^2 + s^2 2^{1-v}/\Gamma(v) (\sum_k \sqrt{2v}/p_k d_{i,j,k})^v K_v(\sum_k \sqrt{2v}/p_k d_{i,j,k})
  * where "stddev" equals s, "range" equals p, "rate" equals v, and "nugget" equals n */
{
  const cube   D;      // distance matrices
  const uword  dim,    // number of populations
               npars;  // number of distance matrices
  const double nugget, // (square root of) additional spatial variance at distance 0
               stddev, // marginal standard deviation
               rate;   // rate of decay in correlation with distance
  const vec    range;  // scale of correlation
  const mat    Dk,     // summed (weighted) distance matrices
               Bk,     // Bessel function part of kernel
               C;      // covariance
  const mat    dC_dt;  // differential of covariance with pars
  const cube   dC_dD;  // differential of covariance with distances 

  Matern (const cube& D, const double nu, const double delta, const vec& pars) 
    : D (D)
    , dim (D.n_rows)
    , npars (D.n_slices)
    , nugget (delta)
    , stddev (pars(0))
    , range (pars.tail(npars))
    , rate (nu)
    , Dk (distances())
    , Bk (besselfn())
    , C (covariance())
    , dC_dt (differential())
    , dC_dD  (differential_dist())
  {
    if (D.n_cols != dim             || pars.n_elem != 1+npars      || 
        (stddev < 0 && nugget <= 0) || (stddev <= 0 && nugget < 0) || 
        arma::any(range <= 0)       || rate <= 0                   )
      Rcpp::stop("Matern: invalid parameters");

    for (uword i=0; i<npars; ++i)
      if (!arma::is_finite(D.slice(i)))
        Rcpp::stop("Matern: non-finite distances");
      else if (arma::any(arma::any(D.slice(i) < 0)))
        Rcpp::warning("Matern: negative distances");
  }

  mat distances (void);
  mat besselfn (void);
  mat covariance (void);
  mat differential (void);
  cube differential_dist (void); 
  vec parameters (void);
};

} // namespace Covariance

#endif /* MATERN_H */
