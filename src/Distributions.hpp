#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include "Inlassle.hpp"

struct Flat
{
  private:
    const vec x;
    const uword n;

  public:
    const double logprob = 0.;

  Flat (const vec& x) : x(x), n(x.n_elem)
  {
    if (!arma::is_finite(x)) Rcpp::stop("Flat: nonfinite input");
  }

  vec gradient (void) const
  {
    return arma::zeros<vec>(n);
  }

  mat hessian (void) const
  {
    return arma::zeros<mat>(n,n);
  }

  double normalizing_constant (void) const
  {
    return -arma::datum::inf; // improper prior
  }

  const vec& unconstrained (void) const
  {
    return x;
  }

  const vec& constrained (void) const
  {
    return x;
  }
};

struct ImproperGamma
{
  private:
    const vec x;
    const uword n;
    const uvec which;

  public:
    double logprob;

  ImproperGamma (const vec& x, const uvec& which) : x(x), n(x.n_elem), which(which)
  {
    // DESCRIPTION
    //    Implementation of the penalization strategy discussed in Chung et al, 2013, "A nondegenerate penalized likelihood
    //    estimator for variance parameters in multilevel models", to get non-zero estimates of the mode of scale parameters.
    //
    //    For a scale parameter \sigma, this corresponds to the improper gamma prior,
    //      p(\sigma) \propto \sigma \approx \sigma \sim gamma(2, \lambda), \lambda \rightarrow 0
    //
    //    Although this functor operates on an unconstrained (log) parameterization of the standard deviation, we do NOT include
    //    the Jacobian in this density.
    //
    //    This should be used for penalized MLE, and NOT for full Bayesian inference.
    //
    //    The optional argument "which" determines which parameters in the vector recieve the penalty. The default is all.

    if (!arma::is_finite(x)) 
      Rcpp::stop("ImproperGamma: non-finite input");

    if (which.n_elem > n || !arma::all(which < n)) 
      Rcpp::stop("ImproperGamma: inconsistent dimensions");

    logprob = -arma::accu(x.elem(which));
  }

  ImproperGamma (const vec& x) : ImproperGamma (x, arma::regspace<uvec>(0, x.n_elem-1)) 
  {}

  vec gradient (void) const
  {
    vec gr = arma::zeros<vec>(n);
    gr.elem(which).ones();
    return -gr; // the gradient is linear, thus pulling estimates away from 0
  }

  mat hessian (void) const
  {
    // the Hessian is zero. 
    return arma::zeros<mat>(n, n);
  }

  double normalizing_constant (void) const
  {
    return -arma::datum::inf; // it's a improper prior
  }

  const vec& unconstrained (void) const
  {
    return x;
  }

  vec constrained (void) const
  {
    return arma::exp(x);
  }
};

struct ImproperWishart
{
  private:
    const vec x;
    const uword n;

  public:
    double logprob;

  ImproperWishart (const vec& x) : x(x), n((sqrt(1. + 8.*double(x.n_elem))-1.)/2.)
  {
    // DESCRIPTION
    //    Implementation of the penalization strategy discussed in Chung et al., 2015, ... to get non-zero
    //    estimates of the mode of a covariance matrix.
    //
    //    For an n-dimensional covariance matrix \Sigma, this corresponds to an improper Wishart prior,
    //      p(\Sigma) \propto |\Sigma|^{1/2}
    //
    //    Although this functor operates on an unconstrained log-LDL' parameterization, we do NOT include
    //    the Jacobian in this density.
    //
    //    This should be used for penalized MLE, and NOT for full Bayesian inference.

    if (!arma::is_finite(x)) Rcpp::stop("ImproperWishart: non-finite input");

    logprob = 0.;
    for (uword i=0; i<n; ++i)
      logprob += -x(i*n-i*(i-1)/2); // the input is the log'd sqrt of the D in the LDL'
  }

  vec gradient (void) const
  {
    vec grad = arma::zeros<vec>(n*(n+1)/2);
    for (uword i=0; i<n; ++i)
      grad(i*n-i*(i-1)/2) = -1; // the gradient is linear (on the diagonal elements), thus pulling estimates away from 0
    return grad;
  }

  mat hessian (void) const
  {
    // The Hessian is actually 0. However, this is bad for optimization. So we include a regularizer here.
    //return arma::zeros<mat>(n*(n+1)/2, n*(n+1)/2); // correct but blows up
    mat hess = arma::zeros<mat>(n*(n+1)/2, n*(n+1)/2);
    for (uword i=0; i<n; ++i)
      hess(i*n-i*(i-1)/2,i*n-i*(i-1)/2) = 1; 
    return hess;
  }

  double normalizing_constant (void) const
  {
    return -arma::datum::inf; // it's a improper prior
  }

  const vec& unconstrained (void) const
  {
    return x;
  }

  const vec& constrained (void) const
  {
    return x; // TODO
  }
};

struct Gaussian
{
  private:
    const vec x,
              mu;
    const mat Lambda; // Precision matrix
    const uword n;

  public:
    double logprob;

  Gaussian (const vec& x, const vec& mu, const mat& Lambda) : x(x), mu(mu), Lambda(Lambda), n(x.n_elem)
  {
    // DESCRIPTION
    //   Your friendly, locally-curated multivariate Gaussian density,
    //         p(x; \mu, \sigma) \propto \exp \{ -0.5 (x - \mu)' \Lambda (x - \mu) \}
    //
    // ARGUMENTS
    //   x  :: vector of length n, Gaussian random variables
    //   mu :: vector of length n, location parameters
    //   Lambda :: n by n matrix, precisions, where det(Lambda) > 0 (this is NOT checked)
    //
    // RETURNS
    //   The negative log density without the normalizing constant (that depends on pars)

    if (!arma::is_finite(x)) Rcpp::stop("Gaussian: non-finite input");
    if (mu.n_elem != n || Lambda.n_cols != n || Lambda.n_rows != n) Rcpp::stop ("Gaussian: inconsistent dimensions");

    logprob = 0.5 * arma::dot(x - mu, Lambda * (x - mu));
  }

  vec gradient (void) const
  {
    return Lambda * (x - mu);
  }

  mat hessian (void) const
  {
    return Lambda;
  }

  double normalizing_constant (void) const
  {
    return 0.; // TODO
  }

  const vec& constrained (void) const
  {
    return x;
  }

  const vec& unconstrained (void) const
  {
    return x;
  }
};

struct LKJ
{
  private:
    const vec x;
    const double eta;
    const vec gamma;
    const uword n;
    vec S, D;

  public:
    double logprob;

  LKJ (const vec& x, const double eta, const vec& gamma) : x(x), eta(eta), gamma(gamma), n(gamma.n_elem) 
  { 
    // DESCRIPTION
    //   A combination of a LKJ prior on n*(n-1)/2 correlation parameters and a half-Cauchy prior on n standard deviations is,
    //         p(C, S; \eta) \propto |C|^{\eta} \prod_i (S_i^2 + \gamma_i^2)^{-1}
    //
    //   As we work in an unconstrained space (e.g. SCS = L \exp{D}^2 L'), we use a transform that maps the n(n-1)/2 lower-triangular elements
    //   of L and the n diagonal elements of D onto n(n-1)/2 correlations and n standard deviations. The log Jacobian determinant for this transform is:
    //         \log |J(L, D \rightarrow C, S)| = 2 \sum_i (n - i + 1) D_i - n \sum_i \log S_i
    //    
    // ARGUMENTS
    //    x :: vector of length n*(n+1)/2, elements of log-Cholesky factor, such that we recover L, D, S by 
    //          "L.elem(lower_tri) = v; D = diagmat(L); L.diag().ones(); D.diag().exp(); S = L * D * D * L.t();"
    //    eta :: double. eta > 0. Penalization on determinant of correlation matrix, larger ==> diagonal.
    //    gamma :: vector. all(gamma > 0). Scale parameters for Cauchy prior on standard deviations.
    //
    // RETURNS
    //    The negative log density in log-LDL' space, sin a normalizing constant that depends on pars, saved as member "logprob"

    if (!arma::is_finite(x)) Rcpp::stop("LKJ: non-finite input");
    if(x.n_elem != n*(n+1)/2) Rcpp::stop("LKJ: inconsistent dimensions");
    if(!arma::all(gamma > 0) || !(eta > 0)) Rcpp::stop("LKJ: invalid parameter values"); 

    // Reparameterize LDL' 
    D = arma::zeros<vec>(n);
    for (uword i=0; i<n; ++i)
      D(i) = x[i*n-i*(i-1)/2];

    S = arma::exp(2 * D);
    for (uword j=n; j>0; --j) // rows
      for (uword i=0; i<j-1; ++i) // cols
        S(j-1) += S(i) * pow(x[i*(n-1)-i*(i-1)/2+j-1],2);
    S = arma::sqrt(S);

    // Negative prior density
    logprob  = -2 * (eta - 1) * arma::accu(D - arma::log(S)); // LKJ
    logprob += arma::accu(arma::log(arma::pow(S,2) + arma::pow(gamma,2))); // Cauchy
    logprob += -2*arma::accu(arma::regspace<vec>(n,1) % D) + double(n)*arma::accu(arma::log(S)); // Jacobian
  }

  vec gradient (void) const
  {
    // Gradient
    mat dS_dv = _dS_dv (S, x);
    vec dl_dS = _dl_dS (S);
    return dS_dv.t() * dl_dS + _dD_dv().t() * _dl_dD();
  }

  mat hessian (void) const
  {
    // Hessian
    mat dS_dv = _dS_dv (S, x);
    vec dl_dS = _dl_dS (S);
    mat hess = dS_dv.t() * arma::diagmat(_d2l_dS2(S)) * dS_dv;
    for (uword i=0; i<n; ++i) // cols
      for (uword j=i; j<n; ++j) // rows
        hess.col(i*(n-1)-i*(i-1)/2+j) += _d2S_dv2(S, x, dS_dv, j, i).t() * dl_dS;
    return hess;
  }

  double normalizing_constant (void) const
  {
    return 0.; //TODO
  }

  const vec& unconstrained (void) const
  {
    return x;
  }

  const vec& constrained (void) const
  {
    return x; //TODO
  }

  mat test_dS_dv (void)
  {
    // Reparameterize LDL' 
    D = arma::zeros<vec>(n);
    for (uword i=0; i<n; ++i)
      D(i) = x[i*n-i*(i-1)/2];

    S = arma::exp(2 * D);
    for (uword j=n; j>0; --j) // rows
      for (uword i=0; i<j-1; ++i) // cols
        S(j-1) += S(i) * pow(x[i*(n-1)-i*(i-1)/2+j-1],2);
    S = arma::sqrt(S);

    return _dS_dv (S, x);
  }

  private:

  // Jacobians
  mat _dD_dv (void) const
  {
    mat dD_dv = arma::zeros<mat>(n, n*(n+1)/2);
    for (uword i=0; i<n; ++i)
      dD_dv(i, i*n-i*(i-1)/2) = 1.;
    return dD_dv;
  }

  vec _dl_dD (void) const
  {
    return -2 * (eta - 1) - 2 * arma::regspace<vec>(n,1);
  }

  vec _dl_dS (const vec& _S) const
  {
    return 1./_S * (double(n) + 2*(eta-1)) + 2 * _S/(arma::pow(_S, 2) + arma::pow(gamma,2)); // latter term is Cauchy .. could replace?
  }

  mat _d2S_dv2 (const vec& _S, const vec& _x, const mat& dS_dv, const uword r, const uword c) const
  {
    mat d2S_dv2 = arma::zeros<mat>(n, n*(n+1)/2);
    uword tr = c*(n-1)-c*(c-1)/2+r;

    if (r == c) // element of D
    {
      for (uword j=r; j<n; ++j)
      {
        d2S_dv2.row(j) = -dS_dv.row(j) * dS_dv(j,tr) / _S(j);
        for (uword i=r; i<n; ++i)
        {
          uword li = c*(n-1)-c*(c-1)/2+i;
          d2S_dv2 (j, li) += 2 * dS_dv(j, li);
        }
      }
    }
    else // element of L
    {
      uword di = c*n-c*(c-1)/2;
      if (_x(tr) == 0.)
        d2S_dv2(r,tr) = exp(2 * _x(di)) / _S(r);
      else
      {
        d2S_dv2.row(r) = -dS_dv.row(r) * dS_dv(r,tr) / _S(r);
        d2S_dv2(r,tr) += dS_dv(r,tr) / _x(tr); // problems
        d2S_dv2(r,di) += 2 * dS_dv(r,di) / _x(tr); // problems when _x = 0
      }
    }

    return d2S_dv2;
  }

  mat _dS_dv (const vec& _S, const vec& _x) const
  {
    mat dS_dv = arma::zeros<mat>(n, n*(n+1)/2);
    for (uword i=0; i<n; ++i)
    {
      uword le = i*n-i*(i-1)/2;
      double eD = exp(2 * _x(le));
      dS_dv (i,le) = eD/_S(i);
      for (uword j=i+1; j<n; ++j)
      {
        uword li = i*(n-1)-i*(i-1)/2+j;
        dS_dv (j,li) = eD*_x(li)/_S(j);
        dS_dv (j,le) = _x(li)*dS_dv(j,li);
      }
    }
    return dS_dv;
  }

  mat _d2l_dS2 (const vec& _S) const
  {
    // latter term is Hessian of Cauchy--could replace
    mat d2l_dS2 = arma::zeros<mat>(n,n);
    d2l_dS2.diag() = -(double(n) + 2*(eta-1))/arma::pow(_S,2) + 2/(arma::pow(_S,2) + arma::pow(gamma,2)) % (1 - 2 * arma::pow(_S,2)/(arma::pow(_S,2) + arma::pow(gamma,2)));
    return d2l_dS2;
  }

//  double everything (const vec& v) 
//  {
//    // DESCRIPTION
//    //   A combination of a LKJ prior on n*(n-1)/2 correlation parameters and a half-Cauchy prior on n standard deviations is,
//    //         p(C, S; \eta) \propto |C|^{\eta} \prod_i (S_i^2 + \gamma_i^2)^{-1}
//    //
//    //   As we work in an unconstrained space (e.g. SCS = L \exp{D}^2 L'), we use a transform that maps the n(n-1)/2 lower-triangular elements
//    //   of L and the n diagonal elements of D onto n(n-1)/2 correlations and n standard deviations. The log Jacobian determinant for this transform is:
//    //         \log |J(L, D \rightarrow C, S)| = 2 \sum_i (n - i + 1) D_i - n \sum_i \log S_i
//    //    
//    // ARGUMENTS
//    //    v    :: vector of length n*(n+1)/2, elements of log-Cholesky factor, such that we recover L, D, S by 
//    //               "L.elem(lower_tri) = v; D = diagmat(L); L.diag().ones(); D.diag().exp(); S = L * D * D * L.t();"
//    //
//    // RETURNS
//    //    The negative log density in log-LDL' space, sin a normalizing constant that depends on pars.
//    //    The gradient and hessian are saved as members.
//
//    if (n*(n+1)/2 != v.n_elem) Rcpp::stop ("LKJ: inconsistent dimension");
//
//    // Reparameterize LDL' 
//    vec D = arma::zeros<vec>(n);
//    for (uword i=0; i<n; ++i)
//      D(i) = v[i*n-i*(i-1)/2];
//
//    vec S = arma::exp(2 * D);
//    for (uword j=n; j>0; --j) // rows
//      for (uword i=0; i<j-1; ++i) // cols
//        S(j-1) += S(i) * pow(v[i*(n-1)-i*(i-1)/2+j-1],2);
//    S = arma::sqrt(S);
//
//    // Negative prior density
//    double lp = 2 * (eta - 1) * arma::accu(D - arma::log(S)); // LKJ
//    lp += -arma::accu(arma::log(arma::pow(S,2) + arma::pow(gamma,2))); // Cauchy
//    lp += 2*arma::accu(arma::regspace<vec>(n,1) % D) - double(n)*arma::sum(arma::log(S)); // Jacobian
//
//    // Gradient
//    mat dS_dv = _dS_dv (S, v);
//    vec dl_dS = _dl_dS (S);
//    grad = dS_dv.t() * dl_dS + _dD_dv().t() * _dl_dD();
//
//    // Hessian
//    hess = dS_dv.t() * arma::diagmat(_d2l_dS2(S)) * dS_dv;
//    for (uword i=0; i<n; ++i) // cols
//      for (uword j=i; j<n; ++j) // rows
//        hess.col(i*(n-1)-i*(i-1)/2+j) += _d2S_dv2(S, v, dS_dv, j, i).t() * dl_dS;
//
//    return -lp;
//  }
};

struct LogBetaPrime 
{
  private:
    const vec   mu,
                delta,
                x,
                s,
                f;
    const uword n;

  public:
    double logprob;

  LogBetaPrime (const vec& x, const vec& mu, const vec& delta) : x(x), s(exp(x)), f(1/(s + 1)), mu(mu), delta(delta), n(x.n_elem)
  {
    // DESCRIPTION
    //    A beta random variable f \in (0,1) has the probability density (in mean-dispersion parameterization),
    //         p(f; \mu, \delta) \propto f^{\mu \delta - 1} (1-f)^{\delta (1-\mu) - 1}
    //   
    //    An unconstrained parameterization for f is s = \log (1/f - 1) \in (-\infty, \infty). The transformed density is,
    //         p(s; \mu, \delta) \propto p(f; \mu, \delta) |\frac{\D f}{\D s}|
    //
    //    Where |\frac{\D f}{\D s}| = f (1 - f). Thus,
    //         p(s; \mu, \delta) \propto f^{\mu \delta} (1-f)^{\delta (1-\mu)}
    //    
    //    This is a "sort-of" reparameterization of the classic beta-prime distribution, hence the name in this module.
    //
    // ARGUMENTS
    //    x     :: vector of unconstrained random variables of dimension n
    //    mu    :: vector of mean parameters, of length n, where all(mu > 0 && mu < 1)
    //    delta :: vector of concentration parameters, of length n, where all(delta > 0)
    //
    // RETURNS
    //    The negative log density without the normalizing constant (that depends on pars), saved in member "logprob"

    if (!arma::is_finite(x)) Rcpp::stop("LogBetaPrime: non-finite input");
    
    if (mu.n_elem != n || delta.n_elem != n) 
      Rcpp::stop("LogBetaPrime: inconsistent dimensions");

    if (arma::any(mu <= 0) || arma::any(delta <= 0) || arma::any(mu >= 1)) 
      Rcpp::stop("LogBetaPrime: invalid parameter values");

    logprob  = -arma::accu(mu % delta % arma::log(f) + (1.-mu) % delta % arma::log(1.-f)); 
  }

  vec gradient (void) const
  {
    return mu % delta % s / (s + 1.) - (1.-mu) % delta / (s + 1.);
  }

  mat hessian (void) const
  {
    return arma::diagmat(mu % delta % s / (s + 1.) % (1. - s/(s + 1.)) + (1.-mu) % delta % s / arma::pow(s + 1., 2));
  }

  double normalizing_constant (void) const
  {
    double nc = 0;
    for (uword i=0; i<n; ++i)
      nc += R::lbeta(mu(i)*delta(i), (1.-mu(i))*delta(i));
    return nc;
  }

  const vec& unconstrained (void) const
  {
    return x;
  }

  const vec& constrained (void) const
  {
    return s;
  }
};

#endif /* DISTRIBUTIONS_H */
