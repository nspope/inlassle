#ifndef INLASSLE_H
#define INLASSLE_H

#include <RcppArmadillo.h> 
#include <RcppParallel.h> 
#include <boost/math/special_functions/polygamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <dlib/matrix.h> 

// [[Rcpp::depends(RcppArmadillo,RcppParallel,BH)]]
// [[Rcpp::plugins(cpp11)]]

using vec   = arma::vec;
using uvec  = arma::uvec;
using mat   = arma::mat;
using uword = arma::uword;
using cube  = arma::cube;

using dlib_mat = dlib::matrix<double>;

using boost::math::polygamma;
using boost::math::lgamma;
using boost::math::cyl_bessel_k;

using Worker = RcppParallel::Worker;
using Split  = RcppParallel::Split;

struct uvecComp
/* 'less-than' comparison functor for uvec, where each vector
 * is the positions of "set bits" in a possibly infinite but 
 * sparse binary vector */
{
  bool operator() (uvec const&, uvec const&) const;
};

using uvec_map_uword = std::map<uvec, uword, uvecComp>;
using uvec_map_mat   = std::map<uvec,   mat, uvecComp>;
using uvec_map_vec   = std::map<uvec,   vec, uvecComp>;
using uvec_list      = std::vector<uvec>;

/* transforms */
vec logit (vec);
double logit (double);
vec plogis (vec);
double plogis (double);
uword tri_diag_ind (uword, uword);
dlib_mat arma_to_dlib (const mat&);
mat dlib_to_arma (const dlib_mat&);

/* utilities */
void cholesky_rev_smith (mat&, const mat&); // reverse autodiff of Cholesky decomposition

#endif /* INLASSLE_H */
