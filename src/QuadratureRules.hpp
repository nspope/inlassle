#ifndef QUADRULES_H
#define QUADRULES_H

#include <cmath>
#include "Inlassle.hpp"

//--------------------------------------- wrapper for a multidimensional rule

struct CompositeRule
{
  uword dim = 0, order = 0;
  mat abscissa = {};
  vec weights = {};

  CompositeRule (void) {}

  template<class Lhs, class... Rhs>
    static CompositeRule make (const Lhs& lhs, const Rhs&... rhs)
    {
      if (lhs.dim != 1)
        Rcpp::stop("CompositeRule: dimension mismatch");
      CompositeRule inp = CompositeRule::make<Rhs...>(rhs...);
      CompositeRule out;
      out.dim = inp.dim + lhs.dim;
      out.order = inp.order * lhs.order;
      out.abscissa = arma::join_horiz(
          arma::kron(inp.abscissa, arma::ones<vec>(lhs.order)), 
          arma::kron(arma::ones<vec>(inp.order), lhs.abscissa));
      out.weights = arma::kron(inp.weights, lhs.weights);
      return out;
    }

  template<class Rhs>
    static CompositeRule make (const Rhs& rhs)
    {
      if (rhs.dim != 1) 
        Rcpp::stop("CompositeRule: dimension mismatch");
      CompositeRule out;
      out.dim = rhs.dim;
      out.order = rhs.order;
      out.abscissa = rhs.abscissa;
      out.weights = rhs.weights;
      return out;
    }
};

// we are integrating a function that is
//   int_0^inf gamma/pi * 1/(x^2 + gamma^2) * h(x)
// change variables, u = x^2
//   int_0^inf gamma/pi * u^{-0.5} / (u + gamma^2) * h(u^{0.5})
// apply Mobius transform u = gamma^2 (1+v)/(1-v)
//   gamma/pi 1/( (gamma^2)^0.5 ) int_-1^1 (1 + v)^-0.5 (1 - v)^-0.5 * h( (gamma^2 (1 + v)/(1 - v))^0.5 )
//   1/pi int_-1^1 (1 + v)^-0.5 (1 - v)^-0.5 * h( gamma (1 - v)^-0.5 / (1 + v)^-0.5 )
// this can be integrated with Jacobi rule, [-1, 1] range, alpha = beta = -0.5
// (according to Homeier & Steinborn 1996 Computer Physics Communications)
// (this ref just basically compares a stupid approach against Gaussian quadrature with change of varibles)

#endif
