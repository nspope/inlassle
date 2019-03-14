#ifndef WEIGHTED_H
#define WEIGHTED_H

#include "Inlassle.hpp"

namespace Covariance {

struct Weighted
 /* Covariance structure where the "distance matrices" are actually covariance matrices; 
  *   \delta_{i,j}^2 + s^2 (\sum_k p_k d_{i,j,k})
  * e.g. global covariance is a weighted sum.
  */ 
{
  const cube   D;      // covariance matrices
  const uword  dim,    // number of populations
               npars;  // number of distance matrices
  const double nugget, // (square root of) additional spatial variance at distance 0
               stddev; // marginal standard deviation (not used)
  const vec    weight; // weights for covariance matrices
  const mat    C;      // covariance
  const mat    dC_dt;  // differential of global covariance with pars
  const cube   dC_dD;  // differential of global covariance with individual covariance matrices

  Weighted (const cube& D, const double nu, const double delta, const vec& pars) // nu is ignored
    : D (D)
    , dim (D.n_rows)
    , npars (D.n_slices)
    , nugget (delta)
    , stddev (pars(0))
    , weight (pars.tail(npars))
    , C (covariance())
    , dC_dt (differential())
    , dC_dD  (differential_dist())
  {
    if (D.n_cols != dim             || pars.n_elem != 1+npars      || 
        (stddev < 0 && nugget <= 0) || (stddev <= 0 && nugget < 0) ||
        arma::any(weight <= 0)       )
      Rcpp::stop("Weighted: invalid parameters");

    for (uword i=0; i<npars; ++i)
    {
      vec d_check = D.slice(i).diag();
      if (!arma::is_finite(D.slice(i)))
        Rcpp::stop("Weighted: non-finite covariances");
      else if (arma::any(D.slice(i).diag() < 0.0))
        Rcpp::warning("Weighted: negative variances");
    }
  }

  mat covariance (void);
  mat differential (void);
  cube differential_dist (void); 
  vec parameters (void);
};

} // namespace Covariance

#endif 
