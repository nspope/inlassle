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
  const uword  verbose = 1,       // 0 == no info, 1 == meta-problem, 2 == meta and sub-problem
               maxiter = 100;     // iterations for meta-problem
        uword  iter = 0;          // iteration counter

  Problem subproblem;
  ResistanceSolver resistance;

  vec start,    // starting values for subproblem, retained across iterations
      gradient; // gradient for meta-problem

  ResistanceOptim (const Problem&, const MatrixXd&, const UiVector&, const MatrixXi&, const double, const uword, const uword);
//  template <class Spatial, class Prior> double fixed_likelihood (const vec&, const vec&); 
//  template <class Spatial, class Prior> vec fixed_optimize (const vec&, const vec&, const double); 
  template <class Spatial, class Prior> double likelihood (const vec&, const uvec&); 
  template <class Spatial, class Prior> vec optimize (const vec&, const uvec&, const double); 
  template <class Spatial, class Prior> vec optimize_global (const vec&, const vec&, const uvec&); 
  template <class Spatial, class Prior> vec grid (const mat&, const uvec&, cube&); 
  void print (const vec&, const vec&) const;
};

#endif /* RESISTANCEOPTIM_H */
