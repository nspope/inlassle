#ifndef PROBLEM_H
#define PROBLEM_H

#include "Inlassle.hpp"

template <class Prior> struct Parameters;

// PURPOSE
//   "Problem" represents the primary optimization problem of this package: fitting a spatially-explicit model
//   to SNP data across loci, where the mean frequency and overdispersion of variant alleles are functions of location-specific
//   covariates, and the covariance in variant counts across locations is the sum of a (Matern) spatial covariance 
//   structure and locus-specific random effects.
//
// RESPONSIBILITIES
//   ...
//
// CONSTRUCTOR
//   N : matrix containing chromosome counts for each locus (column) and location (row)
//   Y : matrix containing variant counts for each locus (column) and location (row)
//   X : design matrix for fixed effects, e.g. covariates (column) that may relate to allele frequencies across locations (row)
//   Z : design matrix for dispersion, e.g. covariates (column) that may relate to overdispersion in allele frequencies across locations (row)
//   D : array of distance matrices for spatial covariance among locations (row x column), across different distance types (slice/third index)
//   tol : double controlling convergence tolerance
//   parallel : boolean controlling whether multiple threads are used

struct Problem
{
  /* data */
  const mat  N, // chromosomes
             Y, // allele counts
             X, // design matrix (fixed effects)
             Z; // design matrix (dispersions)
        cube D; // distances
  const vec  hyperparameters; // vector of prior parameters FIXME this is just a placeholder, hyperparameters are currently hardcoded

  /* dimensions */
  const uword n_popul,  // populations
              n_fixef,  // fixed effects
              n_dispr,  // dispersion parameters
              n_vcomp,  // variance components
              n_sppar,  // spatial parameters
              n_loci,   // loci
              dim,
              npars;

  /* derivatives */
  vec gradient;           // gradient of loglikelihood w.r.t problem parameters
  mat hessian;            // (approximate) hessian of loglikelihood w.r.t problem parameters
  cube gradient_distance; // gradient of loglikelihood w.r.t distance matrix

  /* settings */
  const bool   parallel = false;  // multithreaded?
  const double nu = 2;            // "order" of spatial covariance
  const uword  maxiter = 100;     // maximum number of steps in optimization line search
        uword  iter = 0;          // iteration count

  Problem (const mat&, const mat&, const mat&, const mat&, const cube&, const double, const bool); // constructor; see above
  template <class Spatial, class Prior> double likelihood (const Parameters<Prior>&, const uvec&); // calculate likelihood/gradient given parameter values, boolean for fixed parameters
  template <class Spatial, class Prior> vec optimize (const vec&, const uvec&, const double); // optimize given starting values, boolean for fixed parameters, convergence tolerance
  template <class Spatial, class Prior> vec optimize (const vec&, const double); // optimize given starting values, convergence tolerance
//  template <class Spatial, class Prior> vec optimize_fixed (const vec&, const double, const uvec&);
};

#endif /* PROBLEM_H */
