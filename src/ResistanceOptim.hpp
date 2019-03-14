#ifndef RESISTANCEOPTIM_H
#define RESISTANCEOPTIM_H

#include "Inlassle.hpp"
#include "Problem.hpp"
#include "ResistanceSolver.hpp"

// PURPOSE
//   "ResistanceOptim" represents the meta-problem in this package: optimizing a resistance surface
//   using the likelihood provided by Problem::optimize. Thus it incorporates the solve for
//   resistance distances ("ResistanceSolver") and the random field likelihood ("Problem")

struct ResistanceOptim
{
  const double inner_tol = 1e-10; // convergence tolerance for subproblem
  const double outer_tol = 1e-5;  // convergence tolerance for outer problem
  const uword  maxiter = 100,
               verbose = 1;       // 0 == no info, 1 == meta-problem, 2 == meta and sub-problem
        uword  iter = 0;          // iteration counter

  Problem subproblem;
  ResistanceSolver resistance;

  vec start;    // starting values for subproblem, retained across iterations

  ResistanceOptim (const Problem&, const MatrixXd&, const UiVector&, const MatrixXi&, const double, const double, const uword, const uword);
//  template <class Spatial, class Prior> double fixed_likelihood (const vec&, const vec&); 
//  template <class Spatial, class Prior> vec fixed_optimize (const vec&, const vec&, const double); 
  double BFGS (vec&, vec&, mat&, const uvec&, double&);
  double linesearch(vec&, vec&, vec&, const uvec&, double&);
  double optimize_cov (vec&, vec&, mat&, const uvec&);
  double likelihood_cov (const vec&, vec&, const uvec&);
  template <class Spatial, class Prior> double likelihood (const vec&, vec&, const uvec&); 
  template <class Spatial, class Prior> vec optimize (const vec&, const uvec&); 
  template <class Spatial, class Prior> vec optimize_global (const vec&, const vec&, const uvec&); 
  template <class Spatial, class Prior> vec grid (const mat&, const uvec&, cube&); 
  void print (const vec&, const vec&, const vec&) const;
};

#endif /* RESISTANCEOPTIM_H */
