#ifndef SHRINKAGE_H
#define SHRINKAGE_H

#include "Inlassle.hpp"

struct Problem;

struct Gradient
{
  /* convenience wrapper for gradient of model parameters */

  template <typename T>
  struct Block
  {
    uword n_loci;
    T     x_bar;

    Block (const T& v)
      : n_loci (0)
      , x_bar  (arma::zeros<T>(arma::size(v)))
    {}

    Block<T>& operator+= (const T& input)
      /* stable online calculation of mean, Welford */
    {
      n_loci += 1;
      x_bar  += (input - x_bar)/double(n_loci);
      return *this;
    }

    Block<T>& operator+= (const Block<T>& rhs)
      /* stable parallel calculation of mean and variance, Chan et al. */
    {
      x_bar   = (double(n_loci) * x_bar + double(rhs.n_loci) * rhs.x_bar) / double(n_loci + rhs.n_loci);
      n_loci += rhs.n_loci;
      return *this;
    }

    Block<T> operator<< (const Block<T>& rhs)
      /* vertical join */
    {
      Block<T> lhs (*this);
      if (lhs.n_loci != rhs.n_loci)
        Rcpp::warning ("GradientBlock: vertical join with different numbers of loci");
      lhs.x_bar = arma::join_cols (lhs.x_bar, rhs.x_bar);
      return lhs;
    }

    T operator() (void)
    {
      return x_bar;
    }
  };

  Block<mat>  dC;
  Block<cube> dD;
  Block<vec>  dt, dv, ds, db;

  Gradient (const Problem&);
  vec operator() (void);
  Gradient& operator+= (const Gradient&);
};

struct Hessian
{
  /* estimated Fisher information of parameters */

  struct Block 
  {
    uword n_loci;
    mat   m_bar,  // running sum of squares
          x_bar;  // running mean

    Block (const mat& m)
      : n_loci (0)
      , m_bar  (arma::zeros<mat>(arma::size(m)))
      , x_bar  (arma::zeros<mat>(arma::size(m)))
    {}

    Block& operator+= (const mat& input)
      /* stable online calculation of mean and variance, Welford */
    {
      n_loci += 1;
      mat r_n = input - x_bar,
          x_n = x_bar + r_n/double(n_loci);
      m_bar += r_n % (input - x_n);
      x_bar  = x_n;
      return *this;
    }

    Block& operator+= (const Block& rhs)
      /* stable parallel calculation of mean and variance, Chan et al. */
    {
      m_bar  += rhs.m_bar + arma::pow(x_bar - rhs.x_bar, 2) * double(n_loci * rhs.n_loci) / double(n_loci + rhs.n_loci); // is this really stable?
      x_bar   = (double(n_loci) * x_bar + double(rhs.n_loci) * rhs.x_bar) / double(n_loci + rhs.n_loci);
      n_loci += rhs.n_loci;
      return *this;
    }

    Block operator<< (const Block& rhs)
      /* vertical join */
    {
      Block lhs (*this);
      if (lhs.n_loci != rhs.n_loci)
        Rcpp::warning ("Block: vertical join with different numbers of loci");
      lhs.m_bar   = arma::join_cols (lhs.m_bar, rhs.m_bar);
      lhs.x_bar   = arma::join_cols (lhs.x_bar, rhs.x_bar);
      return lhs;
    }

    Block operator>> (const Block& rhs)
      /* horizontal join */
    {
      Block lhs (*this);
      if (lhs.n_loci != rhs.n_loci)
        Rcpp::warning ("Block: horizontal join with different numbers of loci");
      lhs.m_bar   = arma::join_rows (lhs.m_bar, rhs.m_bar);
      lhs.x_bar   = arma::join_rows (lhs.x_bar, rhs.x_bar);
      return lhs;
    }

    Block t (void)
    {
      Block lhs (*this);
      lhs.m_bar = arma::trans(lhs.m_bar);
      lhs.x_bar = arma::trans(lhs.x_bar);
      return lhs;
    }

    mat shrinkage_estimator (bool type = false)
    {
      mat v_bar = m_bar * double(n_loci) / pow(double(n_loci - 1), 3);
      double lambda;
      if (type)
        lambda = arma::accu(arma::trimatu(v_bar,1))/arma::accu(arma::pow(arma::trimatu(x_bar,1),2)); // shrink toward unaltered diagonal
      else
        lambda = arma::accu(v_bar)/arma::accu(arma::pow(x_bar - arma::eye<mat>(arma::size(x_bar)),2)); // shrink toward unit diagonal
      lambda = std::max(0., std::min(1., lambda));
      return x_bar * (1 - lambda) + arma::diagmat(x_bar) * lambda;
    }

    mat operator() (void)
    {
      return x_bar;
    }
  };

  Block dt2, dv2, ds2, db2, dvdt, dsdt, dbdt, dsdv, dbdv, dbds;
  
  Hessian (const Problem&);
  mat operator() (void);
  mat raw (void);
  Hessian& operator+= (const Hessian&);
};

#endif /* SHRINKAGE_H */
