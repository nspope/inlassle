#ifndef GENOTYPELIKELIHOOD_H
#define GENOTYPELIKELIHOOD_H

#include <cmath>
#include "Inlassle.hpp"

// algorithm for computing site allele frequency likelihood
// from Nielsen et al. 2012 PLoS One
// the input are genotype likelihoods for (AA, Aa, aa); not log-likelihoods; scaling is arbitrary
vec saf (const mat& gl)
{
  uword n = gl.n_rows*2;
  vec h = arma::zeros<vec>(n + 1);
  h.at(0) = gl.at(0,0); 
  h.at(1) = 2*gl.at(0,1); 
  h.at(2) = gl.at(0,2); 
  for (uword i=2; i<=gl.n_rows; ++i) // samples
    for (uword j=2*i; j>=2; --j)//check indices
    {
      h.at(j) = gl.at(i-1,2)*h.at(j-2) + 2*gl.at(i-1,1)*h.at(j-1) + gl.at(i-1,0)*h.at(j);
      h.at(1) = gl.at(i-1,0)*h.at(1) + gl.at(i-1,1)*h.at(0);
      h.at(0) = gl.at(i-1,0)*h.at(0);
    }
  for (uword j=0; j<h.n_elem; ++j)
    h.at(j) *= double(n + 1) * std::beta(n-j+1, j+1); //choose(n,k) == 1/((n+1)*beta(n-k+1,k+1)
  return h;
}
//needs to be normalized?

struct safReader
{
  //TODO
  //GL input is 3 matrices, each nsites x nind, that contain Major/Major Major/Minor Minor/Minor liks
  //parse GL input
  //print number of individuals
  //determine missing data
  //output std::vector with vec of saf for each site. n for site is vec.n_elem - 1.

  //parallelize
  
};

#endif
