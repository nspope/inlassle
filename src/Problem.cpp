#include "Problem.hpp"
#include "Matern.hpp"
#include "Parameters.hpp"
#include "Likelihood.hpp"
#include "Priors.hpp"
#include <dlib/optimization.h> 
#include "dlib_modified_newton.hpp"

// ----------------------------------------------- implementation

Problem::Problem (const mat& N, const mat& Y, const mat& X, const cube& D, const double nu, const bool parallel = false) 
  : N (arma::mat(N))
  , Y (arma::mat(Y))
  , X (arma::mat(X))
  , D (arma::cube(D))
  , n_popul (Y.n_rows)
  , n_fixef (X.n_cols)
  , n_vcomp (X.n_cols*(X.n_cols+1)/2)
  , n_sppar (D.n_slices+1)
  , n_loci  (Y.n_cols)
  , dim (n_loci)
  , npars (n_popul + n_fixef + n_vcomp + n_sppar)
  , parallel (parallel)
  , nu (nu)
{ 
  if (N.n_rows != n_popul || N.n_cols != n_loci  || X.n_rows != n_popul ||
      D.n_rows != n_popul || D.n_cols != n_popul )
    Rcpp::stop ("Problem: inconsistent dimensions");
}

template <class Spatial, class Prior>
double Problem::likelihood (const Parameters<Prior>& param)
{
  ++iter;

  Likelihood<Spatial, Prior> loglik (*this, param);

  gradient = loglik.gradient(param);
  //hessian  = loglik.fisher(param);
  gradient_distance = loglik.gradient_distance (param);

  return loglik.loglikelihood;
}

template <class Spatial, class Prior>
vec Problem::optimize (const vec& start, const double tol)
{
  iter = 0;

  /* objective, gradient, approximated Hessian */
  auto loglik    = [&](const dlib_mat& inp) { Parameters<Prior> p(*this, dlib_to_arma(inp)); return likelihood<Spatial, Prior>(p); };             
  auto grad      = [&](const dlib_mat& inp) { return arma_to_dlib(gradient); }; 
  auto inv_hess  = [&](const dlib_mat& inp) { return arma_to_dlib(hessian); }; // low rank approx ... already inverted via likelihood ...

  /* Either Newton or BFGS is used to minimize objective. */
  auto pars = arma_to_dlib (start);
  auto result = dlib::find_min
    (//dlib::modified_newton_search_strategy (inv_hess), 
     //dlib::lbfgs_search_strategy (10),
     dlib::bfgs_search_strategy (),
     dlib::objective_delta_stop_strategy (tol/double(n_loci), maxiter), // NOTE: rescaled by dimension
     //dlib::gradient_norm_stop_strategy (tol/double(n_loci), maxiter), // NOTE: rescaled by dimension
     loglik, grad, pars,
     0.); // minimum function value for dlib::find_min()

  /* return vector of length start.n_elem + 1, where first element is loglikelihood */
  vec out = { result };
  return arma::join_cols(out, dlib_to_arma(pars)); // TODO: currently returning unconstrained estimates ...
}

// ----------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Problem_likelihood")]]
Rcpp::List test_Problem_likelihood (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, arma::vec t, arma::vec v, arma::vec s, arma::vec b, bool parallel) 
{
  Problem prob (N, Y, X, D, 2, parallel);
  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
  Parameters<Prior::MLE> parm (prob, p);
  double ll = prob.likelihood<Covariance::Matern, Prior::MLE>(parm);
  return Rcpp::List::create(
      Rcpp::_["loglik"] = ll,
      Rcpp::_["gradient"] = prob.gradient,
      Rcpp::_["hessian"] = prob.hessian
      );
}

// [[Rcpp::export("inlassle_test_Problem_plikelihood")]]
Rcpp::List test_Problem_plikelihood (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, arma::vec t, arma::vec v, arma::vec s, arma::vec b, bool parallel) 
{
  Problem prob (N, Y, X, D, 2, parallel);
  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
  Parameters<Prior::Penalized> parm (prob, p);
  double ll = prob.likelihood<Covariance::Matern, Prior::Penalized>(parm);
  return Rcpp::List::create(
      Rcpp::_["loglik"] = ll,
      Rcpp::_["gradient"] = prob.gradient,
      Rcpp::_["hessian"] = prob.hessian
      );
}

// [[Rcpp::export("inlassle_test_Problem_optimize")]]
arma::vec test_Problem_optimize (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol) 
{
  Problem prob (N, Y, X, D, 2, parallel);
  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
  vec out = prob.optimize<Covariance::Matern, Prior::MLE>(start, tol);
  return out;
}

// [[Rcpp::export("inlassle_test_Problem_penalize")]]
arma::vec test_Problem_penalize (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol) 
{
  Problem prob (N, Y, X, D, 2, parallel);
  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
  vec out = prob.optimize<Covariance::Matern, Prior::Penalized>(start, tol);
  return out;
}

// [[Rcpp::export("inlassle_test_Problem_priors")]]
arma::vec test_Problem_priors (arma::mat N, arma::mat Y, arma::mat X, arma::cube D, bool parallel, double tol) 
{
  Problem prob (N, Y, X, D, 2, parallel);
  vec start = arma::zeros<vec>(prob.n_sppar + prob.n_vcomp + prob.n_popul + prob.n_fixef);
  vec out = prob.optimize<Covariance::Matern, Prior::Inlassle>(start, tol);
  return out;
}
