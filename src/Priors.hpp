#ifndef PRIORS_H
#define PRIORS_H

#include "Inlassle.hpp"
#include "Distributions.hpp"

//---------------------------------------------------- implementation

namespace Prior {

struct MLE
{
  // No priors or penalization
  Flat t;
  Flat v;
  Flat s;
  Flat b;

  MLE (const vec& pars, const vec& t, const vec& v, const vec& s, const vec& b) : t(t), v(v), s(s), b(b) {}
};

struct Penalized
{
  // Penalization of variance parameters so as to get non-zero estimates
  ImproperGamma     t; // only the first parameter in vector (marginal std. deviation) is penalized, see constructor
  ImproperWishart   v;
  Flat              s; 
  Flat              b;

  Penalized (const vec& pars, const vec& t, const vec& v, const vec& s, const vec& b) : t(t, uvec({0})), v(v), s(s), b(b) {}
};

struct Inlassle
{
  // Nate's preferred priors
  private:
//    const double t_mu = 0.,      
    const double t_mu = -0.5,      
//                 t_lambda = 0.1,
                 t_lambda = 1,
                 v_eta = 2.,
                 v_gamma = 1.,
                 s_mu = 0.5,
                 s_delta = 1,
                 b_mu = 0.,
                 b_lambda = 0.1;

  public:
    Gaussian     t;
    LKJ          v;
    LogBetaPrime s;
    Gaussian     b;

  Inlassle (const vec& pars, const vec& t, const vec& v, const vec& s, const vec& b) 
    : t (t, t_mu*arma::ones<vec>(t.n_elem), t_lambda*arma::eye<mat>(t.n_elem, t.n_elem))
    , v (v, v_eta, v_gamma*arma::ones<vec>(b.n_elem))
    , s (s, s_mu*arma::ones<vec>(s.n_elem), s_delta*arma::ones<vec>(s.n_elem))
    , b (b, b_mu*arma::ones<vec>(b.n_elem), b_lambda*arma::eye<mat>(b.n_elem, b.n_elem)) 
  {}
};

struct Bedassle
{
  // Priors used in bedassle
  // TODO
};

} /* namespace Prior */

#endif /* PRIORS_H */

