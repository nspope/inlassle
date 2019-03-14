#include "ResistanceOptim.hpp"
#include "Priors.hpp"
#include "Matern.hpp"
#include "Weighted.hpp"
#include "Parameters.hpp"
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

//--------------------------------------- utils

mat dlib2arma (const dlib_mat& x)
{
  mat out (x.nr(), x.nc());
  for (uword i=0; i<x.nr(); ++i)
    for (uword j=0; j<x.nc(); ++j)
      out(i,j) = x(i,j);
  return out;
}

dlib::matrix<double, 0L, 1L> arma2dlib_vec (const vec& x)
{
  // this is annoying
  dlib::matrix<double, 0L, 1L> out;
  out.set_size(x.n_elem);
  for (uword i=0; i<x.n_elem; ++i)
      out(i) = x(i);
  return out; 
}

dlib_mat arma2dlib (const mat& x)
{
  dlib_mat out (x.n_rows, x.n_cols);
  for (uword i=0; i<x.n_rows; ++i)
    for (uword j=0; j<x.n_cols; ++j)
      out(i,j) = x(i,j);
  return out; 
}

MatrixXd arma2eigen (const mat& x)
{
  MatrixXd out (x.n_rows, x.n_cols);
  for (uword i=0; i<x.n_rows; ++i)
    for (uword j=0; j<x.n_cols; ++j)
      out(i,j) = x(i,j);
  return out;
}

mat eigen2arma (const MatrixXd& x)
{
  mat out (x.rows(), x.cols());
  for (uword i=0; i<x.rows(); ++i)
    for (uword j=0; j<x.cols(); ++j)
      out(i,j) = x(i,j);
  return out;
}

MatrixXd dlib2eigen (const dlib_mat& x)
{
  MatrixXd out (x.nr(), x.nc());
  for (uword i=0; i<x.nr(); ++i)
    for (uword j=0; j<x.nc(); ++j)
      out(i,j) = x(i,j);
  return out;
}

dlib_mat eigen2dlib (const MatrixXd& x)
{
  dlib_mat out (x.rows(), x.cols());
  for (uword i=0; i<x.rows(); ++i)
    for (uword j=0; j<x.cols(); ++j)
      out(i,j) = x(i,j);
  return out;
}

//--------------------------------------- implementation

//TODO need to refactor this stuff, possibly?
ResistanceOptim::ResistanceOptim (const Problem& subproblem, const MatrixXd& spatial, const UiVector& targets, const MatrixXi& adjacency,
                                  const double inner_tol, const double outer_tol, const uword maxiter, const uword verbose)
  : subproblem (subproblem)
  , resistance (spatial, targets, adjacency, subproblem.parallel)
  , start (arma::zeros<vec>(subproblem.npars))
  , inner_tol (inner_tol)
  , outer_tol (outer_tol)
  , maxiter (maxiter)
  , verbose (verbose)
{
  if (targets.size() != subproblem.N.n_rows)
    Rcpp::stop ("ResistanceOptim: inconsistent dimensions");
}

double ResistanceOptim::BFGS (vec& par, vec& grad, mat& hess, const uvec& variable, double& loglik)
{
  vec descent = -hess * grad,
      grad_diff = -grad;
  double alpha = linesearch(par, grad, descent, variable, loglik);
  vec step = descent * alpha;
  grad_diff += grad;

  double p = arma::dot(grad_diff, step);
  if (p < 0.1)
    hess = arma::eye<mat>(arma::size(hess));
  else
    hess = (arma::eye<mat>(arma::size(hess)) - step*grad_diff.t()/p) * hess * (arma::eye<mat>(arma::size(hess)) - grad_diff*step.t()/p) + step*step.t()/p;

  double change = alpha > 0.0 ? sqrt(arma::dot(step, step)) : arma::datum::inf;

  if (verbose)
    Rcpp::Rcout << "Iter " << iter << ", loglik = " << loglik << ", ||x'-x|| = " << change << std::endl;

  return change;
}

double ResistanceOptim::linesearch (vec& par, vec& grad, vec& descent, const uvec& variable, double& loglik)
{
  const double c1 = 1e-4, c2 = 0.9; // Wolfe constants
  const double l = loglik;
  const vec p = par, g = grad;
  double alpha = 5., b = 0., a = arma::dot(grad, descent);
  if (a >= 0.)
  {
    Rcpp::warning("\tResistanceOptim: not a descent direction, resetting Hessian");
    return 0.; //restarts Hessian
  }

  for (uword i=0; i<maxiter; ++i)
  {
    alpha *= 0.2;
    par = p + alpha * descent;
    loglik = likelihood_cov(par, grad, variable);
    b = arma::dot(grad, descent);
    if (arma::is_finite(b) && (std::fabs(b) <= c2*std::fabs(a) || b <= 0.))
      return alpha;
  }
  Rcpp::warning("\tResistanceOptim: failed to find sufficient decrease, resetting Hessian");
  par = p; grad = g; loglik = l;
  return 0.; //restarts Hessian
}

double ResistanceOptim::optimize_cov (vec& par, vec& grad, mat& hess, const uvec& variable)
{
  grad = arma::zeros<vec>(par.n_elem);
  hess = arma::eye<mat>(par.n_elem, par.n_elem);
  double loglik = likelihood_cov(par, grad, variable); 

  for(iter=0; iter<maxiter; ++iter)
    if (BFGS(par, grad, hess, variable, loglik) < outer_tol)
      return loglik;

  Rcpp::warning("ResistanceOptim: BFGS did not converge");
  return loglik;
}

// using covariance rather than distance
double ResistanceOptim::likelihood_cov (const vec& par, vec& grad, const uvec& variable)
{
  // the first <> parameters are for the resistance surface
  // the remaining parameters are for the binomial random field
  vec resistance_par = par.head(resistance.npars);
  vec field_par = par.tail(par.n_elem - resistance.npars);

  // calculate resistance distance
  subproblem.D.slice(0) = eigen2arma(resistance.resistance_covariance<Link::Log>(arma2eigen(resistance_par)));

  vec tmp = {1., 1.};
  field_par = arma::join_horiz(tmp, field_par); 
  Parameters<Prior::MLE> p(subproblem, field_par);
  double result = subproblem.likelihood<Covariance::Weighted, Prior::MLE>(p, variable);

  // of resistance stuff
  vec resistance_gradient = eigen2arma(resistance.rd_resistance_covariance<Link::Log>(arma2eigen(subproblem.gradient_distance.slice(0))));

  // of other stuff
  vec field_gradient = subproblem.gradient;
  field_gradient.shed_rows(0,1);//this drops what normally we be the stddev and range

  //
  grad = arma::join_vert(resistance_gradient, field_gradient);

  // verbose
  Rcpp::Rcout << iter << " " << double(subproblem.dim)*result << " " << par << std::endl;

  return result;
}

template <class Spatial, class Prior>
double ResistanceOptim::likelihood (const vec& par, vec& grad, const uvec& variable)
{
  ++iter;

  // calculate resistance distance
  subproblem.D.slice(0) = eigen2arma(resistance.resistance_distances<Link::ReciprocalLogit>(arma2eigen(par)));

  // solve subproblem, e.g. find (possibly penalized) maximum likelihood estimate of spatial parameters
  auto result = subproblem.optimize<Spatial, Prior>(start, variable, inner_tol);

  // use parameter estimates as starting values in next step: this may improve efficieny,
  // especially if taking small steps. However, if estimates are too close to a boundary
  // the optimizer can get "stuck" (this is especially a problem with the rates).
  start = arma::clamp(result.tail(result.n_elem - 1), -3, 3); // restricted to interval [-3, 3]

  grad = eigen2arma(resistance.rd_resistance_distances<Link::ReciprocalLogit>(arma2eigen(subproblem.gradient_distance.slice(0))));

  // verbose
  print (par, grad, result);

  return result(0);
}

template <class Spatial, class Prior>
vec ResistanceOptim::optimize (const vec& start, const uvec& variable)
{
  if (start.n_elem != resistance.npars)
    Rcpp::stop("ResistanceOptim: inconsistent dimensions");

  iter = 0;

  vec gradient = arma::zeros<vec>(size(start));

  auto loglik = [&] (const dlib_mat& par) { return likelihood<Spatial, Prior>(dlib2arma(par), gradient, variable); };
  auto grad   = [&] (const dlib_mat& par) { return arma2dlib(gradient); };

  auto pars = arma2dlib (start);

  auto result = dlib::find_min
    (dlib::bfgs_search_strategy (),
     dlib::objective_delta_stop_strategy (outer_tol/double(subproblem.dim), maxiter), // NOTE: rescaled by dimension
     loglik, grad, pars,
     0.); // minimum function value for dlib::find_min()

  vec out = { result };
  return arma::join_cols(out, dlib2arma(pars)); 
}

template <class Spatial, class Prior>
vec ResistanceOptim::optimize_global (const vec& lower, const vec& upper, const uvec& variable) 
{
  if (lower.n_elem != resistance.npars || upper.n_elem != resistance.npars)
    Rcpp::stop("ResistanceOptim: inconsistent dimensions");

  iter = 0;
  vec gradient (lower.n_elem);

  auto loglik = [&] (const dlib::matrix<double,0L,1L>& par) { return likelihood<Spatial, Prior>(dlib2arma(par), gradient, variable); };

  auto result = dlib::find_min_global
    (loglik, arma2dlib_vec(lower), arma2dlib_vec(upper),
     dlib::max_function_calls(maxiter)); 

  vec out = { result.y };
  return arma::join_cols(out, dlib2arma(result.x)); 
}

template <class Spatial, class Prior>
vec ResistanceOptim::grid (const mat& par, const uvec& variable, cube& rd)
{
  if (par.n_cols != resistance.npars)
    Rcpp::stop("ResistanceOptim: inconsistent dimensions");

  vec out (par.n_rows);
  vec gradient (par.n_cols);
  rd = arma::zeros<cube>(subproblem.n_popul, subproblem.n_popul, par.n_rows);
  for (uword i=0; i<par.n_rows; ++i)
  {
    out(i) = likelihood<Spatial, Prior>(par.row(i).t(), gradient, variable);
    rd.slice(i) = subproblem.D.slice(0);
  }

  return out;
}

void ResistanceOptim::print (const vec& par, const vec& grad, const vec& result) const
{
  if (verbose > 0)
  {
    Rcpp::Rcout << 
      "Iter " << iter << 
      ", |f(x)| = " << result(0) << 
      ", |f'(x)| = " << arma::norm(grad) <<
      ", subproblem solved in " << subproblem.iter << 
      " iterations with |g'(x)| = " << arma::norm(subproblem.gradient) <<
      std::endl;
  }
  if (verbose > 1)
  {
    par.t().print("Parameters:");
    grad.t().print("Gradient:");
    result.tail(result.n_elem-1).t().print("Subproblem parameters:");
    subproblem.gradient.t().print("Subproblem gradient:");
  }
}

//--------------------------------------- tests

// [[Rcpp::export("inlassle_test_ResistanceOptim_likelihood")]]
Rcpp::List test_ResistanceOptim_likelihood (const arma::vec& pars, const arma::uvec& variable, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::mat& Z, 
                                            const arma::cube& D, const Eigen::MatrixXd& S, const std::vector<unsigned>& T, const Eigen::MatrixXi& A)
{
  Problem subproblem (N, Y, X, Z, D, 2, true);
  ResistanceOptim model (subproblem, S, T, A, 1e-10, 1e-4, 100, 1);
  vec grad = arma::zeros<vec>(pars.n_elem);
  double ll = model.likelihood<Covariance::Matern, Prior::MLE> (pars, grad, variable);
  return Rcpp::List::create (Rcpp::_["loglik"] = ll,
                             Rcpp::_["start"] = model.start,
                             Rcpp::_["gradient"] = grad,
                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
}

//Rcpp::List test_ResistanceOptim_fixed_optimize (const arma::vec& pars, const arma::vec& fix, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::cube& D, 
//                                                const Eigen::MatrixXd& Z, const std::vector<unsigned>& T, const Eigen::MatrixXi& A, const arma::uword verbose, const double tol)
//{
//  Problem subproblem (N, Y, X, D, 2, true);
//  ResistanceOptim model (subproblem, Z, T, A, 1e-10, 100, verbose);
//  auto result = model.fixed_optimize<Covariance::Matern, Prior::MLE> (pars, fix, tol); //tol was fixed to 1e-7
//  return Rcpp::List::create (Rcpp::_["loglik"] = result(0),
//                             Rcpp::_["par"] = result.tail(result.n_elem - 1),
//                             Rcpp::_["start"] = fix,
//                             Rcpp::_["gradient"] = model.gradient,
//                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
//}

// [[Rcpp::export("inlassle_test_ResistanceOptim_optimize")]]
Rcpp::List test_ResistanceOptim_optimize (const arma::vec& pars, const arma::uvec& variable, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::mat& Z,
                                          const arma::cube& D, const Eigen::MatrixXd& S, const std::vector<unsigned>& T, const Eigen::MatrixXi& A, const arma::uword verbose)
{
  Problem subproblem (N, Y, X, Z, D, 2, true);
  ResistanceOptim model (subproblem, S, T, A, 1e-10, 1e-4, 100, verbose);
  auto result = model.optimize<Covariance::Matern, Prior::MLE> (pars, variable);
  return Rcpp::List::create (Rcpp::_["loglik"] = result(0),
                             Rcpp::_["par"] = result.tail(result.n_elem - 1),
                             Rcpp::_["start"] = model.start,
                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
}

// [[Rcpp::export("inlassle_test_ResistanceOptim_optimize_global")]]
Rcpp::List test_ResistanceOptim_optimize_global (const arma::vec& lb, const arma::vec& ub, const arma::uvec& variable, const arma::mat& N, const arma::mat& Y, const arma::mat& X, 
                                                 const arma::mat& Z,  const arma::cube& D, const Eigen::MatrixXd& S, const std::vector<unsigned>& T, const Eigen::MatrixXi& A, 
                                                 const arma::uword verbose)
{
  Problem subproblem (N, Y, X, Z, D, 2, true);
  ResistanceOptim model (subproblem, S, T, A, 1e-10, 1e-4, 100, verbose);
  auto result = model.optimize_global<Covariance::Matern, Prior::MLE> (lb, ub, variable);
  return Rcpp::List::create (Rcpp::_["loglik"] = result(0),
                             Rcpp::_["par"] = result.tail(result.n_elem - 1),
                             Rcpp::_["start"] = model.start,
                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
}

// [[Rcpp::export("inlassle_test_ResistanceOptim_grid")]]
Rcpp::List test_ResistanceOptim_grid (const arma::mat& pars, const arma::uvec& variable, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::mat& Z, const arma::cube& D, 
                                      const Eigen::MatrixXd& S, const std::vector<unsigned>& T, const Eigen::MatrixXi& A)
{
  Problem subproblem (N, Y, X, Z, D, 2, true);
  ResistanceOptim model (subproblem, S, T, A, 1e-10, 1e-3, 100, 1);
  cube rd;
  auto result = model.grid<Covariance::Matern, Prior::MLE> (pars, variable, rd);
  return Rcpp::List::create (Rcpp::_["loglik"] = result,
                             Rcpp::_["pars"] = pars,
                             Rcpp::_["rd"] = rd);
}

//Rcpp::List test_ResistanceOptim_priors_likelihood (const arma::vec& pars, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::cube& D, 
//                                            const Eigen::MatrixXd& Z, const std::vector<unsigned>& T, const Eigen::MatrixXi& A)
//{
//  Problem subproblem (N, Y, X, D, 2, false);
//  ResistanceOptim model (subproblem, Z, T, A, 1e-10, 100, 1);
//  double ll = model.likelihood<Covariance::Matern, Prior::Inlassle> (pars);
//  return Rcpp::List::create (Rcpp::_["loglik"] = ll,
//                             Rcpp::_["start"] = model.start,
//                             Rcpp::_["gradient"] = model.gradient,
//                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
//}

//Rcpp::List test_ResistanceOptim_priors_optimize (const arma::vec& pars, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::cube& D, 
//                                          const Eigen::MatrixXd& Z, const std::vector<unsigned>& T, const Eigen::MatrixXi& A, const arma::uword verbose)
//{
//  Problem subproblem (N, Y, X, D, 2, true);
//  ResistanceOptim model (subproblem, Z, T, A, 1e-10, 100, verbose);
//  auto result = model.optimize<Covariance::Matern, Prior::Inlassle> (pars, 1e-7);
//  return Rcpp::List::create (Rcpp::_["loglik"] = result(0),
//                             Rcpp::_["par"] = result.tail(result.n_elem - 1),
//                             Rcpp::_["start"] = model.start,
//                             Rcpp::_["gradient"] = model.gradient,
//                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
//}

//Rcpp::List test_ResistanceOptim_priors_optimize_global (const arma::vec& lb, const arma::vec& ub, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::cube& D, 
//                                          const Eigen::MatrixXd& Z, const std::vector<unsigned>& T, const Eigen::MatrixXi& A, const arma::uword verbose)
//{
//  Problem subproblem (N, Y, X, D, 2, true);
//  ResistanceOptim model (subproblem, Z, T, A, 1e-10, 100, verbose);
//  auto result = model.optimize_global<Covariance::Matern, Prior::Inlassle> (lb, ub);
//  return Rcpp::List::create (Rcpp::_["loglik"] = result(0),
//                             Rcpp::_["par"] = result.tail(result.n_elem - 1),
//                             Rcpp::_["start"] = model.start,
//                             Rcpp::_["gradient"] = model.gradient,
//                             Rcpp::_["subproblem_gradient"] = model.subproblem.gradient);
//}

//Rcpp::List test_ResistanceOptim_priors_grid (const arma::mat& pars, const arma::mat& N, const arma::mat& Y, const arma::mat& X, const arma::cube& D, 
//                                          const Eigen::MatrixXd& Z, const std::vector<unsigned>& T, const Eigen::MatrixXi& A)
//{
//  Problem subproblem (N, Y, X, D, 2, false);
//  ResistanceOptim model (subproblem, Z, T, A, 1e-10, 100, 1);
//  cube rd;
//  auto result = model.grid<Covariance::Matern, Prior::Inlassle> (pars, rd);
//  return Rcpp::List::create (Rcpp::_["loglik"] = result,
//                             Rcpp::_["pars"] = pars,
//                             Rcpp::_["rd"] = rd);
//}
