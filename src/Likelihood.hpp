#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "Inlassle.hpp"
#include "Shrinkage.hpp"

struct Problem;
template <class Prior> struct Parameters;

template <class Spatial, class Prior>
struct Likelihood : public RcppParallel::Worker
{
  // WHAT THIS CLASS REPRESENTS
  //   TODO

  bool singular = true;
  double loglikelihood = arma::datum::inf;

  private:
  /* inputs */
  const Problem &data;   
  const Spatial cov;

  /* parameters */
  const vec sigma, mu;
  const mat Q, dC_dv; //is dC_dv used for anything?? I think it has been move into struct Parameters

  /* output */
  vec loglik, &ll;
  mat dl_dC, &dlp_dC,
      dl_dsigma, &dlp_dsigma,
      dl_dmu, &dlp_dmu;

  public:

  Likelihood (const Problem&, const Parameters<Prior>&);
  Likelihood (const Likelihood&, RcppParallel::Split);

  void operator () (std::size_t, std::size_t);
  void fit_field (const uword);
  double likelihood (const Parameters<Prior>&);
  vec gradient (const Parameters<Prior>&) const;
  cube gradient_distance (const Parameters<Prior>&) const;
  mat fisher (const Parameters<Prior>&) const;
  mat precision (const Problem&, const Parameters<Prior>&);
};

#endif /* LIKELIHOOD_H */
