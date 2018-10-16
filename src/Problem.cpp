#include "Problem.hpp"
#include "Matern.hpp"
#include "Parameters.hpp"
#include "Likelihood.hpp"
#include "Priors.hpp"
#include <dlib/optimization.h> 
#include "dlib_modified_newton.hpp"

// ----------------------------------------------- implementation

Problem::Problem (const mat& N, const mat& Y, const mat& X, const mat& Z, const cube& D, const double nu, const bool parallel = false) 
  : N (arma::mat(N))
  , Y (arma::mat(Y))
  , X (arma::mat(X))
  , Z (arma::mat(Z))
  , D (arma::cube(D))
  , n_popul (Y.n_rows)
  , n_fixef (X.n_cols)
  , n_dispr (Z.n_cols)
  , n_vcomp (X.n_cols*(X.n_cols+1)/2)
  , n_sppar (D.n_slices+1)
  , n_loci  (Y.n_cols)
  , dim (n_loci)
  , npars (n_popul + n_fixef + n_vcomp + n_sppar)
  , parallel (parallel)
  , nu (nu)
{ 
  if (N.n_rows != n_popul || N.n_cols != n_loci  || X.n_rows != n_popul ||
      D.n_rows != n_popul || D.n_cols != n_popul || Z.n_rows != n_popul )
    Rcpp::stop ("Problem: inconsistent dimensions");
}

template <class Spatial, class Prior>
double Problem::likelihood (const Parameters<Prior>& param, const uvec& variable)
{
  ++iter;

  Likelihood<Spatial, Prior> loglik (*this, param);

  gradient = loglik.gradient(param);
  //hessian  = loglik.fisher(param);
  gradient_distance = loglik.gradient_distance (param);

  gradient %= arma::conv_to<vec>::from(variable); // variable == 1 if allowed to be optimized, == 0 else

  return loglik.loglikelihood;
}

//template <class Spatial, class Prior>
//vec Problem::optimize_fixed (const vec& start, const double tol, const uvec& fix)
//{
//  if (fix.n_elem != start.n_elem)
//    Rcpp::stop ("optimize_fixed: inconsistent dimensions");
//
//  iter = 0;
//
//  uvec var_pars = arma::find (fix == 0);
//  vec  arma_inp = start;
//
//  /* objective, gradient */
//  auto loglik = [&](const dlib_mat& inp) 
//  { 
//    arma_inp.elem (var_pars) = dlib_to_arma(inp);
//    Parameters<Prior> p(*this, arma_inp); 
//    return likelihood<Spatial, Prior>(p); 
//  };             
//  auto grad = [&](const dlib_mat& inp) 
//  { 
//    return arma_to_dlib(gradient.elem (var_pars)); 
//  }; 
//
//  /* Either Newton or BFGS is used to minimize objective. */
//  auto pars = arma_to_dlib (arma_inp.elem(var_pars));
//  auto result = dlib::find_min
//    (dlib::bfgs_search_strategy (),
//     dlib::objective_delta_stop_strategy (tol/double(n_loci), maxiter), // NOTE: rescaled by dimension
//     loglik, grad, pars,
//     0.); // minimum function value for dlib::find_min()
//
//  /* return vector of length start.n_elem + 1, where first element is loglikelihood */
//  vec out = { result };
//  arma_inp.elem(var_pars) = dlib_to_arma(pars);
//  return arma::join_cols(out, arma_inp); 
//}

template <class Spatial, class Prior>
vec Problem::optimize (const vec& start, const double tol)
{
  return optimize<Spatial, Prior>(start, arma::ones<uvec>(start.n_elem), tol);
}

template <class Spatial, class Prior>
vec Problem::optimize (const vec& start, const uvec& variable, const double tol)
{
  if (variable.n_elem != start.n_elem || variable.max() > 1)
    Rcpp::stop ("Problem::optimize -- fixed indicator has inconsistent dimensions and/or is nonbinary");

  iter = 0;

  /* objective, gradient, approximated Hessian */
  auto loglik    = [&](const dlib_mat& inp) { Parameters<Prior> p(*this, dlib_to_arma(inp)); return likelihood<Spatial, Prior>(p, variable); };             
  auto grad      = [&](const dlib_mat& inp) { return arma_to_dlib(gradient); }; 
  //auto inv_hess  = [&](const dlib_mat& inp) { return arma_to_dlib(hessian); }; // low rank approx ... already inverted via likelihood ...

  /* Either Newton or BFGS is used to minimize objective. */
  auto pars = arma_to_dlib (start);
  auto result = dlib::find_min
    (dlib::bfgs_search_strategy (),
     dlib::objective_delta_stop_strategy (tol/double(n_loci), maxiter), // NOTE: rescaled by dimension
     loglik, grad, pars,
     0.); // minimum function value for dlib::find_min()

  /* return vector of length start.n_elem + 1, where first element is loglikelihood */
  vec out = { result };
  return arma::join_cols(out, dlib_to_arma(pars)); // FIXME: currently returning unconstrained estimates ...
}

// ----------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Problem_likelihood")]]
Rcpp::List test_Problem_likelihood (arma::mat N, arma::mat Y, arma::mat X, arma::mat Z, arma::cube D, arma::vec t, arma::vec v, arma::vec s, arma::vec b, bool parallel) 
{
  Problem prob (N, Y, X, Z, D, 2, parallel);
  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
  uvec f = arma::ones<uvec>(p.n_elem);
  Parameters<Prior::MLE> parm (prob, p);
  double ll = prob.likelihood<Covariance::Matern, Prior::MLE>(parm, f);
  return Rcpp::List::create(
      Rcpp::_["loglik"] = ll,
      Rcpp::_["gradient"] = prob.gradient,
      Rcpp::_["hessian"] = prob.hessian
      );
}

//Rcpp::List test_Problem_plikelihood (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, arma::vec t, arma::vec v, arma::vec s, arma::vec b, bool parallel) 
//{
//  Problem prob (N, Y, X, D, 2, parallel);
//  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
//  Parameters<Prior::Penalized> parm (prob, p);
//  double ll = prob.likelihood<Covariance::Matern, Prior::Penalized>(parm);
//  return Rcpp::List::create(
//      Rcpp::_["loglik"] = ll,
//      Rcpp::_["gradient"] = prob.gradient,
//      Rcpp::_["hessian"] = prob.hessian
//      );
//}

// [[Rcpp::export("inlassle_test_Problem_optimize")]]
arma::vec test_Problem_optimize (arma::mat N, arma::mat Y, arma::mat X, arma::vec Z, arma::cube D, arma::uvec f, bool parallel, double tol) 
{
  Problem prob (N, Y, X, Z, D, 2, parallel);
  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
  vec out = prob.optimize<Covariance::Matern, Prior::MLE>(start, f, tol);
  return out;
}

//arma::vec test_Problem_optimize_fixed (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol, double nu, arma::vec start, arma::uvec fix) 
//{
//  Problem prob (N, Y, X, D, nu, parallel);
//  vec out = prob.optimize_fixed<Covariance::Matern, Prior::MLE>(start, tol, fix);
//  return out;
//}

//arma::vec test_Problem_penalize (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol) 
//{
//  Problem prob (N, Y, X, D, 2, parallel);
//  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
//  vec out = prob.optimize<Covariance::Matern, Prior::Penalized>(start, tol);
//  return out;
//}

//arma::vec test_Problem_priors (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol) 
//{
//  Problem prob (N, Y, X, D, 2, parallel);
//  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
//  vec out = prob.optimize<Covariance::Matern, Prior::Inlassle>(start, tol);
//  return out;
//}
