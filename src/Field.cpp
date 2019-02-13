#include "Field.hpp"

// ----------------------------------------------- implementation

Field::Field (const vec& y, const vec& n, const vec& mu, const mat& Q)  
  : y (y)
  , n (n)
  , dim (y.n_elem)
  , Q (adjust_missing(Q)) 
  , mu (mu)
  , mode (dim, arma::fill::zeros)
  , freq (dim, arma::fill::zeros)
  , dlp_dx (dim, arma::fill::zeros)
  , dlp_dmu (dim, arma::fill::zeros)
  , dlp_dC (dim, dim, arma::fill::zeros)
  , dl_dx (dim, arma::fill::zeros)
  , d2l_dx2 (dim, arma::fill::zeros)
  , d3l_dx3 (dim, arma::fill::zeros)
  , grad (dim, arma::fill::zeros)
  , hess (dim, dim, arma::fill::zeros)
  , tol (1e-8)
{
  if (n.n_elem != dim || mu.n_elem != dim || Q.n_rows != dim || Q.n_cols != dim)
    Rcpp::stop("Field: inconsistent dimensions from input");
  lapapp = newton_raphson();
}

mat Field::adjust_missing (mat Qi) const
{
  /* adjust precision matrix for missing data */
  uvec miss = arma::find(n <= 0);
  if (miss.n_elem)
  {
    mat I = arma::eye<mat>(arma::size(Qi)),
        E = I.cols(miss);
    Qi   -= Qi * E * arma::solve(E.t()*Qi*E, E.t()) * Qi;
    Qi.diag() += arma::sum(E, 1);
  }
  return Qi;
}

void Field::guess_mode (void)
{
  /* starting values for mode */
  for (uword i=0; i<dim; ++i)
  {
    if (n[i] > 0)
      mode[i] = log(double(y[i]) + 0.5) - log(double(n[i] - y[i]) + 1.);
    else
      mode[i] = mu[i]; // ensures that missing values stay fixed at mode
  }
  freq = 1./(1.+ arma::exp(-mode));

  // gradient, Hessian
  score();
  grad = Q * (mode - mu) - dl_dx;
  hess = Q - arma::diagmat(d2l_dx2); 
}

//-------------------------------------- Newton-Raphson with line search

double Field::dbinom (void) const
{
  double lp = 0;
  for (uword i=0; i<dim; ++i)
    if (n[i] > 0) 
      lp += y[i] * log(freq[i]) + (n[i] - y[i]) * log(1. - freq[i]);
  return lp;
}

void Field::score (void)
{
  for (uword i=0; i<dim; ++i)
    if (n[i] > 0)
    {
      dl_dx[i] = y[i] * (1.- freq[i]) - (n[i] - y[i]) * freq[i];
      d2l_dx2[i] = -n[i] * (1.- freq[i]) * freq[i];
      d3l_dx3[i] = (1.-2.* freq[i]) * d2l_dx2[i];
    }
}

void Field::linesearch (vec&& descent)
{
  /* a gradient-only variant of line search */
  const unsigned maxit = 100;//get rid after debugging
  const double c1 = 1e-4, c2 = 0.9;
  double alpha = 5., a = arma::dot(descent, grad), b = 0.;
  vec origin = mode;
  if (a > 0.) // not a descent direction, U-turn
    descent *= -1.;
  for (unsigned i=0; i<maxit; ++i)
  {
    alpha *= 0.2;
    mode = origin + alpha * descent;
    freq = 1./(1.+ arma::exp(-mode));
    score ();
    grad = Q * (mode - mu) - dl_dx;
    b = arma::dot(descent, grad);
    if (arma::is_finite(b) && (std::fabs(b) <= c2*std::fabs(a) || b < 0.))
      return;
  } 
  Rcpp::warning ("Field: linesearch failed to find sufficient decrease");
}

double Field::newton_raphson (void)
{
  /* optimize (negative, unnormalized) log density of latent field and return log-likelihood */
  converge = false;
  grad.fill(arma::datum::inf);
  guess_mode (); //also initializes gradient and hessian

  for (iter=0; iter<maxiter; ++iter)
  {
    converge = arma::all(arma::abs(grad) < tol);
    if (converge)
      break;
    linesearch(-arma::solve(hess, grad)); //updates mode/freq/grad/score
    hess = Q - arma::diagmat(d2l_dx2); 
  }

  if (!converge)
    Rcpp::warning ("Field: Newton-Raphson failed to converge after 100 iterations");

  /* Laplace approximation */
  residual = mode - mu;
  loglik = 0.5 * arma::dot(residual, Q * residual) - dbinom();
  logdet = log(arma::det(Q)); //what if no missing data? Use default? TODO
  loghes = log(arma::det(hess));

  /* gradient with regard to parameters */
  fisher = arma::inv_sympd (hess);
  dlp_dx = fisher * (fisher.diag() % d3l_dx3);
  dlp_dmu = -Q * (residual + 0.5 * dlp_dx);
  dlp_dC = 0.5 * Q - 0.5 * dlp_dmu * dlp_dmu.t() + 0.5 * Q * (0.25 * dlp_dx * dlp_dx.t() - fisher) * Q; 

  return loglik - 0.5 * (logdet - loghes);
}

//-------------------------------------- tests

// [[Rcpp::export("inlassle_test_Field")]]
Rcpp::List test_Field (arma::vec y, arma::vec n, arma::vec mu, arma::mat Q)
{
  Field field (y, n, mu, Q);
  return Rcpp::List::create(
      Rcpp::_["mode"] = field.mode,
      Rcpp::_["loglik"] = field.loglik,
      Rcpp::_["logdet"] = field.logdet,
      Rcpp::_["loghes"] = field.loghes,
      Rcpp::_["lapapp"] = field.lapapp,
      Rcpp::_["dlp_dmu"] = field.dlp_dmu,
      Rcpp::_["dlp_dC"] = field.dlp_dC,
      Rcpp::_["Q"] = field.Q,
      Rcpp::_["iter"] = field.iter);
}
