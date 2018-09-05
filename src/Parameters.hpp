#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Inlassle.hpp"

struct Problem;

template <class Prior>
struct Parameters
{

  /* container that holds parameter values and is responsible for transforming between
   * constrained and unconstrained parameterizations, as well as incorporating prior information */

  vec t, v, s, b;
  mat LLt, dC_dv;
  vec tuc, vuc, suc, buc;
  Prior prior;

  Parameters (const Problem&, const vec&);
  vec t_constrained (vec);
  vec v_constrained (vec);
  vec s_constrained (vec);
  vec b_constrained (vec);
  vec gradient_unconstrained (void) const;
  mat hessian_unconstrained (void) const;
  vec lower_bounds (void) const;
  vec upper_bounds (void) const;
  vec get_constrained (void) const;
  vec get_unconstrained (void) const;
  double logprob_prior (void) const;
  vec gradient_prior (void) const;
  mat hessian_prior (void) const;
};

#endif /* PARAMETERS_H */
