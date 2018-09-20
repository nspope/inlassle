#ifndef PROBLEM_H
#define PROBLEM_H

#include "Inlassle.hpp"

template <class Prior> struct Parameters;

struct Problem
{
  /* data */
  const mat  N, // chromosomes
             Y, // allele counts
             X; // covariates
        cube D; // distances
  const vec  hyperparameters; // vector of prior parameters TODO this is just a placeholder

  /* dimensions */
  const uword n_popul,  // populations
              n_fixef,  // fixed effects
              n_vcomp,  // variance components
              n_sppar,  // spatial parameters
              n_loci,   // loci
              dim,
              npars;

  /* derivatives */
  vec gradient;
  mat hessian;
  cube gradient_distance;

  /* settings */
  const bool   parallel = false;
  const double nu = 2;            // "order" of spatial covariance
  const uword  maxiter = 100;
        uword  iter = 0;

  Problem (const mat&, const mat&, const mat&, const cube&, const double, const bool);
  template <class Spatial, class Prior> double likelihood (const Parameters<Prior>&);
  template <class Spatial, class Prior> vec optimize (const vec&, const double);
  template <class Spatial, class Prior> vec optimize_fixed (const vec&, const double, const uvec&);
};

#endif /* PROBLEM_H */
