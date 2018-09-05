#ifndef FIELD_H
#define FIELD_H

#include "Inlassle.hpp"

struct Field
{
  /* functor for the objective function of the latent field */

  /* settings */
  const double tol      = 1e-4;  // convergence tolerence in optimization
  const uword  maxiter  = 100;   // maximum iterations allowed for Newton-Raphson
        bool   converge = false;
        uword  iter     = 0;

  /* data */
  const vec&   y, n;    // successes, trials

  /* params */
  const uword  dim;     // number of points in field
  const mat    Q;       // precision matrix of field, corrected for missing data 
  const vec&   s,       // beta-binomial dispersion parameters
               mu;      // mean (logit) frequency in field

  /* output */
  double lapapp = 0,    // Laplace approximation 
         loglik = 0,    // loglikelihood
         logdet = 0,    // log determinant of precision TODO this *could* be const and only recalculated if necessary!
         loghes = 0;    // log determinant of Hessian 
  vec    mode,          // values of field
         dlp_dmu,       // derivative of Laplace approx w.r.t mu 
         dlp_ds;        // derivative of Laplace approx w.r.t s 
  mat    dlp_dC;        // derivative of Laplace approx w.r.t inv(Q)

  private:
    /* temporaries */
    vec::fixed<4> dg0, dg1, dg2;
    mat hess, fisher; 
    vec grad, dl_dx, d2l_dx2, dl_ds, d2l_dsdx, d3l_dsdx2, d3l_dx3, dlp_dx;

  public:
    Field (const vec&, const vec&, const vec&, const vec&, const mat&); 
    mat adjust_missing (mat) const; 
    double dbetabinom (void) const;
    void score (void);
    double newton_raphson (void);
    void derivatives (void);
};

#endif /* FIELD_H */
