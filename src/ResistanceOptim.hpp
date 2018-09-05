#ifndef RESISTANCEOPTIM_H
#define RESISTANCEOPTIM_H

#include "Inlassle.hpp"
#include "Problem.hpp"
#include "ResistanceSolver.hpp"

struct ResistanceOptim
{
  const double inner_tol = 1e-10; // convergence tolerance for subproblem
  const uword  verbose = 1,
               maxiter = 100;
        uword  iter = 0;

  Problem subproblem;
  ResistanceSolver resistance;

  vec start,
      gradient;

  ResistanceOptim (const Problem&, const MatrixXd&, const UiVector&, const MatrixXi&, const double, const uword, const uword);
  template <class Spatial, class Prior> double likelihood (const vec&); 
  template <class Spatial, class Prior> vec optimize (const vec&, const double); 
  template <class Spatial, class Prior> vec optimize_global (const vec&, const vec&); 
  template <class Spatial, class Prior> vec grid (const mat&, cube&); 
  void print (const vec&, const vec&) const;
};

#endif /* RESISTANCEOPTIM_H */
