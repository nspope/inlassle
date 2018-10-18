#include "Field.hpp"
#include <dlib/optimization.h> 
#include "dlib_modified_newton.hpp"

// ----------------------------------------------- utils

mat psd_inv (const mat& m, double c = 1)
{
  mat eigvec;
  vec eigval;
  arma::eig_sym(eigval, eigvec, m);
  for (uword i=0; i<eigval.n_elem; ++i)
    if (eigval(i) <= 0.) 
      eigval(i) = c; 
  return eigvec * arma::diagmat(eigval) * eigvec.t();
}

// ----------------------------------------------- implementation

Field::Field (const vec& y, const vec& n, const vec& mu, const vec& s, const mat& Q)  
  : y (y)
  , n (n)
  , dim (y.n_elem)
  , Q (adjust_missing(Q)) 
  , s (s)
  , mu (mu)
  , mode (dim, arma::fill::zeros)
  , dlp_dmu (dim, arma::fill::zeros)
  , dlp_ds (dim, arma::fill::zeros)
  , dlp_dC (dim, dim, arma::fill::zeros)
  , dl_dx (dim, arma::fill::zeros)
  , d2l_dx2 (dim, arma::fill::zeros)
  , dl_ds (dim, arma::fill::zeros)
  , d2l_dsdx (dim, arma::fill::zeros)
  , d3l_dsdx2 (dim, arma::fill::zeros)
  , d3l_dx3 (dim, arma::fill::zeros)
  , grad (dim, arma::fill::zeros)
  , hess (dim, dim, arma::fill::zeros)
  , tol (1e-8)
{
  if (n.n_elem != dim || mu.n_elem != dim || Q.n_rows != dim || Q.n_cols != dim || s.n_elem != dim)
    Rcpp::stop("Field: inconsistent dimensions");
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

//-------------------------------------- hand-rolled Newton-Raphson 

double Field::dbetabinom (void) const
{
  /* vectorized (log) beta-binomial density that incorporates logit-link */
  double z, m, a, b, lp = 0;
  for (uword i=0; i<dim; ++i)
    if (n[i] > 0) 
    {
      z = n[i] - y[i]; m = plogis(mode[i]);
      a = s[i]*m;      b = s[i]*(1-m);
      lp += lgamma(n[i]+1) + lgamma(y[i]+a) + lgamma(z+b) + lgamma(a+b) -
            lgamma(y[i]+1) - lgamma(z+1) - lgamma(n[i]+a+b) - lgamma(a) - lgamma(b);
    }
  return lp;
}

void Field::score (void) 
{
  /* score and Fisher information of log beta-binomial density with respect to latent field x */
  double z, m, a, b, dl_dm, dm_dx, d2l_dm2, d2m_dx2;
  for (uword i=0; i<dim; ++i)
    if (n[i] > 0)
    {
      m = plogis(mode[i]); z = n[i] - y[i];
      a = s[i] * m;        b = s[i] * (1-m);

      dg0[0] = polygamma(0,y[i]+a); dg1[0] = polygamma(1,y[i]+a);
      dg0[1] = polygamma(0,a);      dg1[1] = polygamma(1,a);           
      dg0[2] = polygamma(0,z+b);    dg1[2] = polygamma(1,z+b);   
      dg0[3] = polygamma(0,b);      dg1[3] = polygamma(1,b);     
      
      dl_dm      = s[i]*(dg0[0] - dg0[1] - dg0[2] + dg0[3]);
      dm_dx      = m * (1. - m);
      dl_dx[i]   = dl_dm*dm_dx;
      d2l_dm2    = s[i]*s[i]*(dg1[0] - dg1[1] + dg1[2] - dg1[3]); 
      d2m_dx2    = dm_dx * (1. - 2 * m);
      d2l_dx2[i] = d2m_dx2*dl_dm + d2l_dm2*pow(dm_dx,2);
    }
}

void Field::derivatives (void)
{
  double z, 
         m, 
         a, 
         b, 
         dl_dm, 
         dm_dx, 
         d2l_dm2, 
         d2m_dx2, 
         d2l_dsdm, 
         d3l_dsdm2, 
         d3l_dm3, 
         d3m_dx3;

  for (uword i=0; i<dim; ++i)
    if (n[i] > 0)
    {
      m = plogis(mode[i]);
      z = n[i] - y[i];
      a = s[i] * m;
      b = s[i] * (1-m);

      dg0[0] = polygamma(0,y[i]+a); dg1[0] = polygamma(1,y[i]+a); dg2[0] = polygamma(2,y[i]+a);
      dg0[1] = polygamma(0,a);      dg1[1] = polygamma(1,a);      dg2[1] = polygamma(2,a);
      dg0[2] = polygamma(0,z+b);    dg1[2] = polygamma(1,z+b);    dg2[2] = polygamma(2,z+b);
      dg0[3] = polygamma(0,b);      dg1[3] = polygamma(1,b);      dg2[3] = polygamma(2,b);

      dl_dm        = s[i]*(dg0[0] - dg0[1] - dg0[2] + dg0[3]);
      dm_dx        = exp(-mode[i])*pow(1 + exp(-mode[i]), -2);
      dl_dx        = dl_dm*dm_dx;
      d2l_dm2      = s[i]*s[i]*(dg1[0] - dg1[1] + dg1[2] - dg1[3]);
      d2m_dx2      = dm_dx*(1. - 2*m);
      d2l_dx2[i]   = d2m_dx2*dl_dm + d2l_dm2*pow(dm_dx,2);
      dl_ds[i]     = polygamma(0,s[i]) - polygamma(0,n[i]+s[i]) + m*(dg0[0] - dg0[1]) + (1-m)*(dg0[2] - dg0[3]);
      d2l_dsdm     = a*(dg1[0] - dg1[1]) - b*(dg1[2] - dg1[3]) + dg0[0] - dg0[1] - dg0[2] + dg0[3];
      d2l_dsdx[i]  = d2l_dsdm*dm_dx;
      d3l_dsdm2    = s[i]*a*(dg2[0] - dg2[1]) + s[i]*b*(dg2[2] - dg2[3]) + 2*s[i]*(dg1[0] - dg1[1] + dg1[2] - dg1[3]);
      d3l_dsdx2[i] = d2l_dsdm*d2m_dx2 + d3l_dsdm2*pow(dm_dx,2);
      d3l_dm3      = pow(s[i],3)*(dg2[0] - dg2[1] - dg2[2] + dg2[3]);
      d3m_dx3      = dm_dx*(1. - 6*dm_dx);
      d3l_dx3[i]   = 3*d2l_dm2*d2m_dx2*dm_dx + dl_dm*d3m_dx3 + d3l_dm3*pow(dm_dx,3);
    }
}

double Field::newton_raphson (void)
{
  /* optimize (negative, unnormalized) log density of latent field and return log-likelihood */
  converge = false;
  grad.fill(arma::datum::inf);

  /* initial values */
  // this could be improved? By somehow incorporating the mean/precision matrix.
  // In general this will work well with more observations.
  for (uword i=0; i<dim; ++i)
  {
    if (n[i] > 0)
      mode[i] = log(double(y[i]) + 0.5) - log(double(n[i] - y[i]) + 1.);
    else
      mode[i] = mu[i]; // ensures that missing values stay fixed at mode
  }
  mode.print("starting mode");//DEBUG

  // start dlib
  auto ll = [&](const dlib_mat& inp) 
  { 
    double z, m, a, b, lp = 0;
    mode = dlib_to_arma(inp);
    for (uword i=0; i<dim; ++i)
      if (n[i] > 0) 
      {
        z = n[i] - y[i]; m = plogis(mode[i]);
        a = s[i]*m;      b = s[i]*(1-m);
        lp += lgamma(n[i]+1) + lgamma(y[i]+a) + lgamma(z+b) + lgamma(a+b) -
              lgamma(y[i]+1) - lgamma(z+1) - lgamma(n[i]+a+b) - lgamma(a) - lgamma(b);
        Rcpp::Rcout << "Obj" << std::endl << lp << " " << mode[i] << " " << s[i] << " " << y[i] << " " << n[i] << std::endl; //DEBUG
      }
    return 0.5 * arma::dot(mode-mu, Q * (mode-mu)) - lp;
  };

  // key point: the gradient and Hessian vanish when m or (1-m) = 0.
  // But: digamma fnc will throw errors. So to make safe we should
  // escape before then, setting dl_dx to 0.
  auto gr = [&](const dlib_mat& inp) 
  { 
    double z, m, a, b;
    mode = dlib_to_arma(inp);
    for (uword i=0; i<dim; ++i)
      if (n[i] > 0)
      {
        m = plogis(mode[i]); 
        z = n[i] - y[i];
        a = s[i] * m;        
        b = s[i] * (1-m);
        dl_dx[i] = m * (1. - m) * s[i] * (polygamma(0,y[i]+a) - polygamma(0,a) - polygamma(0,z+b) + polygamma(0,b));      
      }
    grad = Q * (mode - mu) - dl_dx;

    grad.t().print("Grad");//DEBUG

    return arma_to_dlib(grad);
  }; 

  auto he = [&](const dlib_mat& inp) 
  { 
    double z, m, a, b;
    mode = dlib_to_arma(inp);
    for (uword i=0; i<dim; ++i)
      if (n[i] > 0)
      {
        m = plogis(mode[i]); 
        z = n[i] - y[i];
        a = s[i] * m;        
        b = s[i] * (1-m);
        d2l_dx2[i] = dl_dx[i] * (1. - 2*m) + 
          pow(m * (1. - m), 2) * s[i] * s[i] * (polygamma(1,y[i]+a) - polygamma(1,a) + polygamma(1,z+b) - polygamma(1,b)); 
      }
    hess = Q - arma::diagmat(d2l_dx2);
    return arma_to_dlib(hess);
  }; 

  vec lb = mode; lb.fill(-arma::datum::inf);
  vec ub = mode; ub.fill(arma::datum::inf);

  auto pars = arma_to_dlib(mode);

  auto result = dlib::find_min_box_constrained // use b/c only does backtracking, no gradient calculation (that could blow up) in line search
    (dlib::bfgs_search_strategy (),
     //dlib::newton_search_strategy (he),
     dlib::objective_delta_stop_strategy (1e-9, maxiter), 
     ll, gr, pars, arma_to_dlib(lb), arma_to_dlib(ub));

  mode = dlib_to_arma(pars);
  // end dlib

//  // The big problem here is that if mu is outside of the radius of convergence, then 
//  // Newton-Raphson will fail miserably. This could perhaps be solved via line-search?
//  // Another problem is that it's very possible for the Hessian to be indefinite at mode = mu.
//  //
//  // I've got an implementation of a modified Cholesky that should work for this, in that it always returns a descent direction is is well conditioned.
//
//  double min_lambda;
//  
//  /* iterations */
//  for (iter=0; iter<maxiter; ++iter)
//  {
//    if (arma::all(arma::abs(grad) < tol))
//    {
//      converge = true;
//      break;
//    }
//    if (iter)
//      mode -= arma::solve(hess, grad);
//    score (); 
//    hess = Q - arma::diagmat(d2l_dx2); // orig
//    grad = Q * (mode - mu) - dl_dx;
//  }
//
//  if (!converge)
//    Rcpp::warning ("Newton-Raphson failed to converge after 100 iterations");

  derivatives ();

  /* Laplace approximation */
  vec residual = mode - mu;
  hess = Q - arma::diagmat(d2l_dx2); 
  loglik = 0.5 * arma::dot(residual, Q * residual) - dbetabinom();
  logdet = log(arma::det(Q)); //what if no missing data? Use default? TODO
  loghes = log(arma::det(hess));

  /* gradient with regard to parameters */
  fisher = arma::inv_sympd (hess);
  dlp_dx = fisher * (fisher.diag() % d3l_dx3);
  dlp_dmu = -Q * (residual + 0.5 * dlp_dx);
  dlp_ds = -dl_ds - 0.5 * fisher.diag() % d3l_dsdx2 - 0.5 * dlp_dx % d2l_dsdx;
  dlp_dC = 0.5 * Q - 0.5 * dlp_dmu * dlp_dmu.t() + 0.5 * Q * (0.25 * dlp_dx * dlp_dx.t() - fisher) * Q; 

  return loglik - 0.5 * (logdet - loghes);
}

//-------------------------------------- tests

// [[Rcpp::export("inlassle_test_Field_mode")]]
arma::vec test_Field_mode (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.mode;
}

// [[Rcpp::export("inlassle_test_Field_loglik")]]
double test_Field_loglik (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.loglik;
}

// [[Rcpp::export("inlassle_test_Field_logdet")]]
double test_Field_logdet (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.logdet;
}

// [[Rcpp::export("inlassle_test_Field_loghes")]]
double test_Field_loghes (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.loghes;
}

// [[Rcpp::export("inlassle_test_Field_lapapp")]]
double test_Field_lapapp (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.lapapp;
}

// [[Rcpp::export("inlassle_test_Field_dlp_dmu")]]
arma::vec test_Field_dlp_dmu (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.dlp_dmu;
}

// [[Rcpp::export("inlassle_test_Field_dlp_ds")]]
arma::vec test_Field_dlp_ds (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.dlp_ds;
}

// [[Rcpp::export("inlassle_test_Field_dlp_dC")]]
arma::mat test_Field_dlp_dC (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.dlp_dC;
}

// [[Rcpp::export("inlassle_test_Field_Q")]]
arma::mat test_Field_Q (arma::vec y, arma::vec n, arma::vec mu, arma::vec s, arma::mat Q)
{
  Field field (y, n, mu, s, Q);
  return field.Q;
}
