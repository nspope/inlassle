#ifndef FIELD_H
#define FIELD_H

#include "Inlassle.hpp"

struct Field
{
  // WHAT THIS OBJECT REPRESENTS
  //   TODO  

  /* settings */
  const double tol      = 1e-8;  // convergence tolerence in optimization
  const uword  maxiter  = 100;   // maximum iterations allowed for Newton-Raphson
        bool   converge = false;
        uword  iter     = 0;

  /* data */
  const vec&   y, n;    // successes, trials

  /* params */
  const uword  dim;     // number of points in field
  const mat    Q;       // precision matrix of field, corrected for missing data 
  const vec&   mu;      // mean (logit) frequency in field

  /* output */
  double lapapp = 0,    // Laplace approximation 
         loglik = 0,    // loglikelihood
         logdet = 0,    // log determinant of precision TODO this *could* be const and only recalculated if necessary!
         loghes = 0;    // log determinant of Hessian 
  vec    mode,          // values of field
         freq,          // probability parameter of binomial
         dlp_dmu;       // derivative of Laplace approx w.r.t mu 
  mat    dlp_dC;        // derivative of Laplace approx w.r.t inv(Q)

  private:
    /* temporaries */
    mat hess, fisher; 
    vec grad, dl_dx, d2l_dx2, d3l_dx3, dlp_dx, residual;

  public:
    Field (const vec&, const vec&, const vec&, const mat&); 
    mat adjust_missing (mat) const; 
    void guess_mode (void);
    double dbinom (void) const;
    void score (void);
    void linesearch (const vec&);
    double newton_raphson (void);
};

#endif /* FIELD_H */
