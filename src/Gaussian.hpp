#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include "Inlassle.hpp"

namespace Covariance {

struct Gaussian
  /* Gaussian covariance structure,
   *   \delta_{i,j} n^2 + s^2 \exp \{ -(\sum_k 1/p_k d_{i,j,k})^v \}
   * where the "stddev" is s, the "range" is p, the "rate" is v, and the "nugget" is n */
{
  const cube   D;      // distance matrices
  const uword  dim,    // number of populations
               npars;  // number of distance matrices
  const double nugget, // (square root of) additional variance distance 0
               stddev, // marginal standard deviation
               rate;   // rate of decay in correlation with distance
  const vec    range;  // scale of correlation
  const mat    Dk,     // summed (weighted) distance matrices
               C;      // covariance
  const mat    dC_dt;  // differential of covariance with pars
  const cube   dC_dD;  // differential of covariance with distance

  Gaussian (const cube& D, const double nu, const double delta, const vec& pars) :
    D      (D),
    dim    (D.n_rows),
    npars  (D.n_slices),
    nugget (delta),
    stddev (pars(0)),
    rate   (nu),
    range  (pars.tail(npars)),
    Dk     (distances()),
    C      (covariance()),
    dC_dt  (differential()),
    dC_dD  (differential_dist())
  {
    if (D.n_cols != dim             || pars.n_elem != 1+npars      || 
        (stddev < 0 && nugget <= 0) || (stddev <= 0 && nugget < 0) || 
        arma::any(range <= 0)       || rate <= 0                   )
      Rcpp::stop("Gaussian: invalid parameters");

    for (uword i=0; i<npars; ++i)
      if (!arma::is_finite(D.slice(i)))
        Rcpp::stop("Gaussian: non-finite distances");
      else if (arma::any(arma::any(D.slice(i) < 0)))
        Rcpp::warning("Gaussian: negative distances");
  }

  mat distances (void);
  mat covariance (void);
  mat differential (void);
  cube differential_dist (void);
  vec parameters (void);
};

} // namespace Covariance

#endif /* GAUSSIAN_H */
