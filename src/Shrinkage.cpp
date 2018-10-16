#include "Shrinkage.hpp"
#include "Problem.hpp"

// FIXME: I'm pretty sure this is NOT being used anymore, as there's no shrinkage estimator of the Hessian being used

Gradient::Gradient (const Problem& data) 
  : dD (arma::zeros<cube>(data.n_popul, data.n_popul, data.n_sppar-1))
  , dC (arma::zeros<mat>(data.n_popul, data.n_popul))
  , dt (arma::zeros<vec>(data.n_sppar))
  , dv (arma::zeros<vec>(data.n_vcomp))
  , ds (arma::zeros<vec>(data.n_popul))
  , db (arma::zeros<vec>(data.n_fixef))
{}   

vec Gradient::operator() (void)
{
  auto gradient = (dt << dv << ds << db);
  return gradient.x_bar;
}

Gradient& Gradient::operator+= (const Gradient& rhs)
{
  dt += rhs.dt; dv += rhs.dv; ds += rhs.ds;
  db += rhs.db; dC += rhs.dC; dD += rhs.dD;
  return *this;
}

Hessian::Hessian (const Problem& data) 
  : dt2 (arma::zeros<mat>(data.n_sppar, data.n_sppar))
  , dv2 (arma::zeros<mat>(data.n_vcomp, data.n_vcomp))
  , ds2 (arma::zeros<mat>(data.n_popul, data.n_popul))
  , db2 (arma::zeros<mat>(data.n_fixef, data.n_fixef))
  , dvdt (arma::zeros<mat>(data.n_vcomp, data.n_sppar))
  , dsdt (arma::zeros<mat>(data.n_popul, data.n_sppar))
  , dbdt (arma::zeros<mat>(data.n_fixef, data.n_sppar))
  , dsdv (arma::zeros<mat>(data.n_popul, data.n_vcomp))
  , dbdv (arma::zeros<mat>(data.n_fixef, data.n_vcomp))
  , dbds (arma::zeros<mat>(data.n_fixef, data.n_popul))
{}   

mat Hessian::operator() (void)
{
  auto hessian = (dt2      << dvdt     << dsdt     << dbdt) >>
                 (dvdt.t() << dv2      << dsdv     << dbdv) >>
                 (dsdt.t() << dsdv.t() << ds2      << dbds) >>
                 (dbdt.t() << dbdv.t() << dbds.t() << db2 );
  return hessian.shrinkage_estimator();
}

mat Hessian::raw (void)
{
  auto hessian = (dt2      << dvdt     << dsdt     << dbdt) >>
                 (dvdt.t() << dv2      << dsdv     << dbdv) >>
                 (dsdt.t() << dsdv.t() << ds2      << dbds) >>
                 (dbdt.t() << dbdv.t() << dbds.t() << db2 );
  return hessian.x_bar;
}

Hessian& Hessian::operator+= (const Hessian& rhs)
{
  dt2  += rhs.dt2;  dv2 += rhs.dv2;   ds2 += rhs.ds2;
  db2  += rhs.db2;  dvdt += rhs.dvdt; dsdt += rhs.dsdt;
  dbdt += rhs.dbdt; dsdv += rhs.dsdv; dbdv += rhs.dbdv;
  dbds += rhs.dbds;
  return *this;
}

// [[Rcpp::export]]
arma::cube inlassle_test_HessianBlock_inc1 (arma::mat inp)
{
  /* test online incrementing with arma::mat and calculation of shrinkage estimator */
  Hessian::Block block (inp.col(0) * inp.col(0).t());
  for (uword i=0; i<inp.n_cols; ++i)
    block += inp.col(i) * inp.col(i).t();

  arma::cube out (inp.n_rows, inp.n_rows, 3);
  out.slice(0) = block.x_bar;
  out.slice(1) = block.m_bar;
  out.slice(2) = block.shrinkage_estimator();

  return out;
}

// [[Rcpp::export]]
arma::cube inlassle_test_HessianBlock_inc2 (arma::mat inp1, arma::mat inp2)
{
  /* test parallel incrementing with Block and calculation of shrinkage estimator */
  Hessian::Block block1 (inp1.col(0) * inp1.col(0).t());
  Hessian::Block block2 (inp2.col(0) * inp2.col(0).t());
  for (uword i=0; i<inp1.n_cols; ++i)
    block1 += inp1.col(i) * inp1.col(i).t();
  for (uword i=0; i<inp2.n_cols; ++i)
    block2 += inp2.col(i) * inp2.col(i).t();

  block1 += block2;

  arma::cube out (inp1.n_rows, inp1.n_rows, 3);
  out.slice(0) = block1.x_bar;
  out.slice(1) = block1.m_bar;
  out.slice(2) = block1.shrinkage_estimator();

  return out;
}

// [[Rcpp::export]]
arma::cube inlassle_test_HessianBlock_join (arma::mat inp1, arma::mat inp2)
{
  /* test transpose, horizontal and vertical joins */
  Hessian::Block block1 (inp1.col(0) * inp1.col(0).t());
  Hessian::Block block2 (inp2.col(0) * inp2.col(0).t());
  Hessian::Block block3 (inp1.col(0) * inp2.col(0).t());
  for (uword i=0; i<inp1.n_cols; ++i)
    block1 += inp1.col(i) * inp1.col(i).t();
  for (uword i=0; i<inp2.n_cols; ++i)
    block2 += inp2.col(i) * inp2.col(i).t();
  for (uword i=0; i<inp2.n_cols; ++i)
    block3 += inp1.col(i) * inp2.col(i).t();

  auto block = (block1 << block3.t()) >> (block3 << block2);

  arma::cube out (inp1.n_rows+inp2.n_rows, inp1.n_rows+inp2.n_rows, 3);
  out.slice(0) = block.x_bar;
  out.slice(1) = block.m_bar;
  out.slice(2) = block.shrinkage_estimator();

  return out;
}

// [[Rcpp::export]]
arma::vec inlassle_test_GradientBlock_inc1 (arma::mat inp)
{
  /* test online incrementing with arma::vec */
  Gradient::Block<vec> block (inp.col(0));
  for (uword i=0; i<inp.n_cols; ++i)
    block += inp.col(i);
  return block.x_bar;
}

// [[Rcpp::export]]
arma::vec inlassle_test_GradientBlock_inc2 (arma::mat inp1, arma::mat inp2)
{
  /* test parallel incrementing with Block */
  Gradient::Block<vec> block1 (inp1.col(0));
  Gradient::Block<vec> block2 (inp2.col(0));
  for (uword i=0; i<inp1.n_cols; ++i)
    block1 += inp1.col(i);
  for (uword i=0; i<inp2.n_cols; ++i)
    block2 += inp2.col(i);
  block1 += block2;
  return block1.x_bar;
}

// [[Rcpp::export]]
arma::vec inlassle_test_GradientBlock_join (arma::mat inp1, arma::mat inp2)
{
  /* test vertical join */
  Gradient::Block<vec> block1 (inp1.col(0));
  Gradient::Block<vec> block2 (inp2.col(0));
  for (uword i=0; i<inp1.n_cols; ++i)
    block1 += inp1.col(i);
  for (uword i=0; i<inp2.n_cols; ++i)
    block2 += inp2.col(i);
  auto block = (block1 << block2);
  return block.x_bar;
}

// [[Rcpp::export]]
arma::vec inlassle_test_Gradient (arma::mat dt, arma::mat dv, arma::mat ds, arma::mat db)
{
  mat Y = arma::ones<mat>(ds.n_rows, ds.n_cols),
      N = arma::ones<mat>(ds.n_rows, ds.n_cols),
      Z = arma::eye<mat>(ds.n_rows, ds.n_rows),
      X = arma::ones<mat>(ds.n_rows, db.n_rows);
  cube D = arma::zeros<cube>(ds.n_rows, ds.n_rows, dt.n_rows - 1);
  Problem prob (N, Y, X, Z, D, 2, false);
  Gradient grad (prob);

  for (uword i = 0; i < ds.n_cols; ++i)
  {
    grad.dt += dt.col(i);
    grad.dv += dv.col(i);
    grad.ds += ds.col(i);
    grad.db += db.col(i);
  }

  return grad();
}

// [[Rcpp::export]]
arma::mat inlassle_test_Gradient_dC (arma::cube dC)
{
  mat Y = arma::ones<mat>(dC.n_rows, dC.n_slices),
      N = arma::ones<mat>(dC.n_rows, dC.n_slices),
      Z = arma::eye<mat>(dC.n_rows, dC.n_rows),
      X = arma::ones<mat>(dC.n_rows, 1);
  cube D = arma::zeros<cube>(dC.n_rows, dC.n_rows, 1);
  Problem prob (N, Y, X, Z, D, 2, false);
  Gradient grad (prob);

  for (uword i = 0; i < dC.n_slices; ++i)
    grad.dC += dC.slice(i);

  return grad.dC();
}

// [[Rcpp::export]]
arma::mat inlassle_test_Hessian (arma::mat dt, arma::mat dv, arma::mat ds, arma::mat db)
{
  mat Y = arma::ones<mat>(ds.n_rows, ds.n_cols),
      N = arma::ones<mat>(ds.n_rows, ds.n_cols),
      Z = arma::eye<mat>(ds.n_rows, ds.n_rows),
      X = arma::ones<mat>(ds.n_rows, db.n_rows);
  cube D = arma::zeros<cube>(ds.n_rows, ds.n_rows, dt.n_rows - 1);
  Problem prob (N, Y, X, Z, D, 2, false);
  Hessian hess (prob);

  for (uword i = 0; i < ds.n_cols; ++i)
  {
    hess.dt2  += dt.col(i) * dt.col(i).t();
    hess.dvdt += dv.col(i) * dt.col(i).t();
    hess.dsdt += ds.col(i) * dt.col(i).t();
    hess.dbdt += db.col(i) * dt.col(i).t();
    hess.dv2  += dv.col(i) * dv.col(i).t();
    hess.dsdv += ds.col(i) * dv.col(i).t();
    hess.dbdv += ds.col(i) * dv.col(i).t();
    hess.ds2  += ds.col(i) * ds.col(i).t();
    hess.dbds += db.col(i) * ds.col(i).t();
    hess.db2  += db.col(i) * db.col(i).t();
  }

  return hess();
}
