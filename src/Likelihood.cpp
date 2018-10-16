#include "Likelihood.hpp"
#include "Parameters.hpp"
#include "Problem.hpp"
#include "Field.hpp"
#include "Matern.hpp"
#include "Priors.hpp"

// ----------------------------------------------- implementation

template <class Spatial, class Prior>
Likelihood<Spatial, Prior>::Likelihood (const Problem& data, const Parameters<Prior>& par) 
  : data (data)
  , cov (data.D, data.nu, 0.1, par.t) // FIXME hardcoded delta = 0.1, make set-able
//  , s (par.s)
  , s (arma::exp(-data.Z * par.s))
  , mu (data.X * par.b)
//  , Q (arma::inv_sympd(cov.C + data.X * par.LLt * data.X.t()))
  , Q (precision(data, par))
  , dC_dv (par.dC_dv)
  , loglik (data.n_loci, arma::fill::zeros)
  , ll (loglik)
  , dl_dC (pow(data.n_popul,2), data.n_loci, arma::fill::zeros)
  , dl_ds (data.n_popul, data.n_loci, arma::fill::zeros)
  , dl_dmu (data.n_popul, data.n_loci, arma::fill::zeros)
  , dlp_dC (dl_dC)
  , dlp_ds (dl_ds)
  , dlp_dmu (dl_dmu)
{
  if (!singular)
    loglikelihood = likelihood(par);
  else
  {
    Rcpp::stop ("Likelihood: field covariance is singular");
    loglikelihood = arma::datum::inf;
  }
}

template <class Spatial, class Prior>
Likelihood<Spatial, Prior>::Likelihood (const Likelihood<Spatial, Prior>& rhs, RcppParallel::Split) 
  : data (rhs.data)
  , cov (rhs.cov)
  , s (rhs.s)
  , mu (rhs.mu)
  , Q (rhs.Q)
  , dC_dv (rhs.dC_dv)
  , ll (rhs.ll)
  , dlp_dC (rhs.dlp_dC)
  , dlp_ds (rhs.dlp_ds)
  , dlp_dmu (rhs.dlp_dmu)
{}

template <class Spatial, class Prior>
mat Likelihood<Spatial, Prior>::precision (const Problem& data, const Parameters<Prior>& par)
{
  mat out;
  singular = !arma::inv_sympd(out, cov.C + data.X * par.LLt * data.X.t());
  return out;
}

template <class Spatial, class Prior>
void Likelihood<Spatial, Prior>::operator () (std::size_t begin, std::size_t end)
{
  for (auto i=begin; i!=end; ++i)
    fit_field (i);
}

template <class Spatial, class Prior>
double Likelihood<Spatial, Prior>::likelihood (const Parameters<Prior>& par)
{
  if (data.parallel)
    RcppParallel::parallelFor (0, data.n_loci, *this);
  else
    (*this)(0, data.n_loci);
  return arma::mean(loglik) + par.logprob_prior()/double(data.n_loci); 
}

template <class Spatial, class Prior>
void Likelihood<Spatial, Prior>::fit_field (const uword i)
{
  /* fit latent field and calculate Laplace approximation at locus */
  Field field (data.Y.col(i), data.N.col(i), mu, s, Q);

  ll(i) = field.lapapp; 

  dlp_dC.col(i) = arma::vectorise(field.dlp_dC);
  dlp_ds.col(i) = field.dlp_ds;
  dlp_dmu.col(i) = field.dlp_dmu;
}

template <class Spatial, class Prior>
vec Likelihood<Spatial, Prior>::gradient (const Parameters<Prior>& par) const
{
  // This function calculates the average gradient for parameters across i.i.d loci.

  vec grad_C = arma::mean(dlp_dC, 1),
      grad_s = arma::mean(dlp_ds, 1),
      grad_mu = arma::mean(dlp_dmu, 1);

  // let s = exp(-S), S = Z * z
  // by the chain rule, we have
  // dlp/dS = dlp/ds % d(exp(-S))/dS = dlp/ds % -exp(-S) = dlp/ds % -s
  // dlp/dz = dlp/dS % d(Z * z)/dz = Z.t() * dlp/dS
  grad_s %= -s; 

  return arma::join_vert(cov.dC_dt.t() * grad_C,
         arma::join_vert(par.dC_dv.t() * grad_C,
         arma::join_vert(data.Z.t() * grad_s, data.X.t() * grad_mu))) %
         par.gradient_unconstrained() + 
         par.gradient_prior()/double(data.n_loci);
}

template <class Spatial, class Prior>
cube Likelihood<Spatial, Prior>::gradient_distance (const Parameters<Prior>& par) const
{
  // This function calculates the average gradient for distances across i.i.d loci.

  mat grad_C = arma::reshape(arma::mean(dlp_dC, 1), data.n_popul, data.n_popul);

  return cov.dC_dD.each_slice() % grad_C; 
}

template <class Spatial, class Prior>
mat Likelihood<Spatial, Prior>::fisher (const Parameters<Prior>& par) const
{
  // This function provides an estimate of the Fisher information from the empirical
  // covariance of the gradients across i.i.d loci. The estimate is potentially low 
  // rank, as very small singular values are eliminated, to avoid numerical issues.
  
  mat dlp_dt = cov.dC_dt.t() * dlp_dC,
      dlp_dv = par.dC_dv.t() * dlp_dC,
      dlp_dS = dlp_ds.each_col() % -s, // see ::gradient above
      dlp_db = data.X.t() * dlp_dmu;

  dlp_dS = data.Z.t() * dlp_dS; // see ::gradient above

  mat scatter = arma::join_vert(dlp_dt,
                arma::join_vert(dlp_dv,
                arma::join_vert(dlp_dS, dlp_db)));

  scatter = (scatter * scatter.t()) % par.hessian_unconstrained() + par.hessian_prior(); 
  scatter /= double(data.n_loci);

  mat eigvec;
  vec eigval;
  if (!arma::eig_sym(eigval, eigvec, scatter))
    Rcpp::stop ("Likelihood: eigendecomposition unsuccessful");

  vec eigval_mod = arma::clamp(eigval, arma::max(eigval) * sqrt(arma::datum::eps), arma::datum::inf); // tolerance from MASS::ginv

  return eigvec * arma::diagmat(1./eigval_mod) * eigvec.t();
}

// ----------------------------------------------- explicit instantiations

template class Likelihood<Covariance::Matern, Prior::MLE>;
template class Likelihood<Covariance::Matern, Prior::Penalized>;
template class Likelihood<Covariance::Matern, Prior::Inlassle>;

// ----------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Likelihood")]]
Rcpp::List test_Likelihood (arma::mat N, arma::mat Y, arma::mat X, arma::mat Z, arma::cube D, arma::vec t, arma::vec v, arma::vec s, arma::vec b, bool parallel) 
{
  Problem prob (N, Y, X, Z, D, 2, parallel);
  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
  Parameters<Prior::MLE> parm (prob, p);
  Likelihood<Covariance::Matern, Prior::MLE> lik (prob, parm);
  return Rcpp::List::create (
           Rcpp::_["Q"] = lik.precision(prob, parm),
           Rcpp::_["loglik"] = lik.loglikelihood,
           Rcpp::_["gradient"] = lik.gradient(parm),
           Rcpp::_["gradient_distance"] = lik.gradient_distance(parm),
           Rcpp::_["hessian"] = lik.fisher(parm));
}
