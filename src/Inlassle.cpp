#include "Inlassle.hpp"

bool uvecComp::operator() (uvec const& lhs, uvec const& rhs) const
{
  uword i;
  for (i=0; i<lhs.n_elem; ++i)
  {
    if (i+1 > rhs.n_elem) 
    { 
      return false; 
    } 
    else if (lhs[i] > rhs[i]) 
    {
      return false;
    }
    else if (lhs[i] < rhs[i]) 
    {
      return true;
    }
  }
  if (i < rhs.n_elem) 
  {
    return true;
  } 
  else 
  {
    return false;
  }
}

vec logit (vec p)
{
  return arma::log(p) - arma::log(1-p);
}

double logit (double p)
{
  return log(p) - log(1-p);
}

vec plogis (vec x)
{
  x = arma::exp(x);
  return x / (1. + x);
}

double plogis (double x)
{
  x = exp(x);
  return x / (1. + x);
}

uword tri_diag_ind (uword i, uword n)
{
  return (n*(n+1) - (n-i)*(n-i+1))/2;
}

dlib_mat arma_to_dlib (const mat& x)
{
  dlib_mat y;
  y.set_size (x.n_rows, x.n_cols);
  for (uword i = 0; i < x.n_rows; ++i)
    for (uword j = 0; j < x.n_cols; ++j)
      y(i,j) = x(i,j);
  return y;
}

mat dlib_to_arma (const dlib_mat& x)
{
  mat y (x.nr(), x.nc()); 
  for (uword i = 0; i < x.nr(); ++i)
    for (uword j = 0; j < x.nc(); ++j)
      y(i,j) = x(i,j);
  return y;
}

void cholesky_rev_smith (mat& F, const mat& L)
  /* reverse-mode autodiff for (lower) Cholesky decomposition ala Smith 1995,
   *   e.g. for L = arma::chol(..., "lower") 
   * done in place, with input matrix F of differentials per Cholesky element */
{
  for (int k=F.n_rows; k-- > 0; )
  {
    for (uword j=k+1; j<F.n_rows; ++j)
    {
      for (uword i=j; i<F.n_rows; ++i)
      {
        F.at(i,k) -= F.at(i,j)*L.at(j,k);
        F.at(j,k) -= F.at(i,j)*L.at(i,k);
      }
      F.at(j,k) /= L.at(k,k);
      F.at(k,k) -= L.at(j,k)*F.at(j,k);
    }
    F.at(k,k) /= 2.*L.at(k,k);
  }
}

