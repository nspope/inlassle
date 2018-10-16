#include "Parameters.hpp"
#include "Problem.hpp"
#include "Priors.hpp"

//---------------------------------------------------- utils

uword linear_index (const uword r, const uword c, const uword n)
{
  return arma::sub2ind (arma::size(n, n), r, c);
}

uword lower_triangular_index (const uword r, const uword c, const uword n)
{
  // r >= c
  return n*(n+1)/2 - (n-c)*(n-c+1)/2 + r - c;
}

//---------------------------------------------------- implementation

template <class Prior>
Parameters<Prior>::Parameters (const Problem& data, const vec& p) 
  : t (data.n_sppar, arma::fill::zeros)
  , v (data.n_vcomp, arma::fill::zeros)
  , s (data.n_dispr, arma::fill::zeros)
  , b (data.n_fixef, arma::fill::zeros)
  , LLt (data.n_fixef, data.n_fixef, arma::fill::zeros)
  , dC_dv (data.n_popul*data.n_popul, data.n_vcomp, arma::fill::zeros)
  , tuc (p.subvec(0, arma::size(t)))
  , vuc (p.subvec(t.n_elem, arma::size(v)))
  , suc (p.subvec(t.n_elem + v.n_elem, arma::size(s)))
  , buc (p.subvec(t.n_elem + v.n_elem + s.n_elem, arma::size(b)))
  , prior (data.hyperparameters, tuc, vuc, suc, buc) // prior is ALWAYS in unconstrained parameterization
{
  if (p.n_elem != t.n_elem + v.n_elem + s.n_elem + b.n_elem)
    Rcpp::stop ("Parameters: wrong length");

  // map parameters from unconstrained to constrained
  t = t_constrained(tuc);
  v = v_constrained(vuc);
  s = s_constrained(suc);
  b = b_constrained(buc);

  // form log-LDL' decomposition of the variance components in matrix form
  uvec nonz = arma::find(arma::trimatl(arma::ones(arma::size(LLt))));
  LLt.elem(nonz) = v; 
  mat P = arma::diagmat(LLt); 
  LLt.diag().ones(); 

  // form Jacobian matrix of the isomorphism f(L,D) = LDL' = C
  for (uword i=0; i<LLt.n_rows; ++i) //cols
    for (uword j=i; j<LLt.n_rows; ++j) //rows
    {
      uword li = linear_index (j, i, LLt.n_rows), // linear index of lower element of symmetric pair
            ui = linear_index (i, j, LLt.n_rows), // linear index of upper element of symmetric pair
            lt = lower_triangular_index (j, i, LLt.n_rows); // lower triangular index
      if (i==j)
        dC_dv.col(lt) = arma::kron(2 * data.X * LLt.col(i), data.X * LLt * P.col(j));
      else
        dC_dv.col(lt) = arma::kron(data.X.col(j), data.X * LLt * P * P.col(i)) +
            arma::kron(data.X * LLt * P * P.col(i), data.X.col(j));
    }

  // form covariance matrix for variance components, via isomorphism above
  LLt *= arma::diagmat(P);
  LLt *= LLt.t();
}

template <class Prior>
vec Parameters<Prior>::t_constrained (vec tu)
  /* transform unconstrained spatial parameters to constrained spatial parameters:
   *   exp ()
   */
{
  if (!arma::is_finite(tu))
    Rcpp::stop ("Parameters: invalid value for spatial parameters");
  return arma::exp (tu);
}

template <class Prior>
vec Parameters<Prior>::v_constrained (vec vu)
  /* transform unconstrained variance components to constrained variance components
   *   log-Cholesky, e.g. exp() if diagonal
   */
{
  if (!arma::is_finite(vu))
    Rcpp::stop ("Parameters: invalid value for variance components");
  for (uword i=0; i<b.n_elem; ++i)
    vu(tri_diag_ind(i, b.n_elem)) = exp(vu(tri_diag_ind(i, b.n_elem))); 
//    if (vu(tri_diag_ind(i, b.n_elem)) < 0)
//      Rcpp::stop ("Parameters: invalid value for variance components");
  return vu;
}

template <class Prior>
vec Parameters<Prior>::s_constrained (vec su)
  /* transform unconstrained concentrations to constrained concentrations
   *   .. currently none ...
   */
{
  if (!arma::is_finite(su)) 
    Rcpp::stop ("Parameters: invalid value for concentration parameters");
  return su;
}

template <class Prior>
vec Parameters<Prior>::b_constrained (vec bu)
  /* transform unconstrained fixed effects to constrained fixed effects
   *   .. currently none ...
   */
{
  if (!arma::is_finite(bu)) 
    Rcpp::stop ("Parameters: invalid value for fixed effects");
  return bu;
}

template <class Prior>
vec Parameters<Prior>::gradient_unconstrained (void) const
{
  // DESCRIPTION
  //   If f(z) = x is a transformation from unconstrained parameters z to constrained parameters x,
  // and the transformations involve only single variables, then the Jacobian is diagonal with the
  // diagonal elements J_{i,i} = f'_i(z_i). This is need to, e.g. map a gradient from constrained to
  // unconstrained space.
  //
  // RETURNS
  //   A vector with the diagonal of the Jacobian transformation from unconstrained to constrained.
  // Ordered as: spatial parameters, variance components, concentration parameters, fixed effects.
  //
  // FIXME make these individual functions, e.g. gradients of the mappings t_constrained etc.

  vec dt_dtu = t,                              // exp(tu)
      dv_dvu = arma::ones<vec>(arma::size(v)), // 1 if not diagonal element
      ds_dsu = arma::ones<vec>(arma::size(s)), // 1
      db_dbu = arma::ones<vec>(arma::size(b)); // 1

  for (uword i=0; i<b.n_elem; ++i)
    dv_dvu(tri_diag_ind(i, b.n_elem)) = v(tri_diag_ind(i, b.n_elem)); // exp(vu) if diagonal element

  return arma::join_vert(dt_dtu, arma::join_vert(dv_dvu, arma::join_vert(ds_dsu, db_dbu)));
}

template <class Prior>
mat Parameters<Prior>::hessian_unconstrained (void) const
  /* chain rule for mapping estimated Fisher information in constrained space to
   * unconstrained space. This is, e.g.:
   *  sum (dl_dt % dt_dtu) * (dl_dt % dt_dtu).t() =
   *  (dt_dtu * dt_dtu.t()) % sum (dl_dt * dl_dt.t())
   */
{
  vec gr = gradient_unconstrained ();
  return gr * gr.t();
}

template <class Prior>
vec Parameters<Prior>::lower_bounds (void) const
{
  // FIXME deprecated, delete
  vec lb_t = arma::ones<vec>(arma::size(t)) * -arma::datum::inf, // -Inf
      lb_v = arma::ones<vec>(arma::size(v)) * -arma::datum::inf, // -Inf
      lb_s = arma::ones<vec>(arma::size(s)) * -arma::datum::inf, // -Inf
      lb_b = arma::ones<vec>(arma::size(b)) * -arma::datum::inf; // -Inf
  return arma::join_vert(lb_t, arma::join_vert(lb_v, arma::join_vert(lb_s, lb_b)));
}

template <class Prior>
vec Parameters<Prior>::upper_bounds (void) const
{
  // FIXME deprecated, delete
  vec ub_t = arma::ones<vec>(arma::size(t)) * arma::datum::inf, // Inf
      ub_v = arma::ones<vec>(arma::size(v)) * arma::datum::inf, // Inf
      ub_s = arma::ones<vec>(arma::size(s)) * arma::datum::inf, // Inf
      ub_b = arma::ones<vec>(arma::size(b)) * arma::datum::inf; // Inf
  return arma::join_vert(ub_t, arma::join_vert(ub_v, arma::join_vert(ub_s, ub_b)));
}

template <class Prior>
vec Parameters<Prior>::get_constrained (void) const
{
  // DESCRIPTION
  //   Gives the constrained parameterization -- defined on subspaces of (-Inf, Inf).
  //
  // RETURNS
  //   A vector with the constrained parameters, ordered as: spatial parameters, variance components,
  // concentration parameters, fixed effects.

  return arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
}

template <class Prior>
vec Parameters<Prior>::get_unconstrained (void) const
{
  // DESCRIPTION
  //   Gives the unconstrained parameterization -- defined on (-Inf, Inf) -- for all parameters.
  //
  // RETURNS
  //   A vector with the unconstrained parameters, ordered as: spatial parameters, variance components,
  // concentration parameters, fixed effects.

  return arma::join_vert(tuc, arma::join_vert(vuc, arma::join_vert(suc, buc)));
}

template <class Prior>
double Parameters<Prior>::logprob_prior (void) const
{
  // DESCRIPTION
  //   Gives the value of the prior distribution in unconstrained space. The prior(s) reflect the
  // transformation between constrained and unconstrained space: that is, if p(x) is the prior in
  // constrained space (e.g. the positive reals, or PSD matrices), and f(z) = x is a transformation
  // between unconstrained variable z and constrained variable x, then g(z) = p(x) |J_{f(z)}| is the
  // transformed distribution, where |J_.| is the determinant of the Jacobian of the transformation.
  //
  // RETURNS
  //   The _negative log prior probability_ of the unconstrained parameters.

  return prior.t.logprob + 
         prior.v.logprob + 
         prior.s.logprob + 
         prior.b.logprob;
}

template <class Prior>
vec Parameters<Prior>::gradient_prior (void) const
{
  // DESCRIPTION
  //   Gives the gradient of the prior distribution in unconstrained space. For definition, see
  // Parameters::logprob_prior.
  //   
  // RETURNS
  //   The gradient vector for the _negative log prior probability_ w.r.t unconstrained variables. 
  // The vector is ordered: spatial parameters, variance components, concentration parameters, then 
  // fixed effects.

  return arma::join_vert(prior.t.gradient(),
         arma::join_vert(prior.v.gradient(),
         arma::join_vert(prior.s.gradient(),
                         prior.b.gradient())));
}

template <class Prior>
mat Parameters<Prior>::hessian_prior (void) const
{
  // DESCRIPTION
  //   Gives the Hessian of the prior distribution in unconstrained space. For definition, see
  // Parameters::logprob_prior.
  //
  // RETURNS
  //   The Hessian matrix for the _negative log prior probability_ w.r.t the unconstrained variables. 
  // The matrix is ordered: spatial parameters, variance components, concentration parameters, then 
  // fixed effects.

  mat out = arma::zeros<mat>(t.n_elem+v.n_elem+s.n_elem+b.n_elem, t.n_elem+v.n_elem+s.n_elem+b.n_elem);
  out.submat (0, 0, arma::size(t.n_elem, t.n_elem)) = prior.t.hessian();
  out.submat (t.n_elem, t.n_elem, arma::size(v.n_elem, v.n_elem)) = prior.v.hessian();
  out.submat (t.n_elem+v.n_elem, t.n_elem+v.n_elem, arma::size(s.n_elem, s.n_elem)) = prior.s.hessian();
  out.submat (t.n_elem+v.n_elem+s.n_elem, t.n_elem+v.n_elem+s.n_elem, arma::size(b.n_elem, b.n_elem)) = prior.b.hessian();
  return out;
}

//---------------------------------------------------- explicit instantiations

template class Parameters<Prior::MLE>;
template class Parameters<Prior::Penalized>;
template class Parameters<Prior::Inlassle>;

//---------------------------------------------------- tests

// [[Rcpp::export("inlassle_test_Parameters")]]
Rcpp::List test_Parameters (arma::vec t, arma::vec v, arma::vec s, arma::vec b)
{
  /* helper that creates Parameters object used for tests */
  mat X = arma::ones<mat> (s.n_elem, b.n_elem),
      Z = arma::eye<mat> (s.n_elem, s.n_elem),
      Y = arma::ones<mat> (s.n_elem, 1),
      N = arma::ones<mat> (s.n_elem, 1);
  cube D = arma::ones<cube> (s.n_elem, s.n_elem, t.n_elem-1);

  Problem prob (N, Y, X, Z, D, 2, false);
  vec p = arma::join_vert(t, arma::join_vert(v, arma::join_vert(s, b)));
  Parameters<Prior::MLE> parm (prob, p);

  return Rcpp::List::create (
           Rcpp::_["get_constrained"] = parm.get_constrained(),
           Rcpp::_["get_unconstrained"] = parm.get_unconstrained (),
           Rcpp::_["lower_bounds"] = parm.lower_bounds (),
           Rcpp::_["upper_bounds"] = parm.upper_bounds (),
           Rcpp::_["gradient_constrained"] = parm.gradient_unconstrained (),
           Rcpp::_["hessian_unconstrained"] = parm.hessian_unconstrained (),
           Rcpp::_["LLt"] = parm.LLt,
           Rcpp::_["dC_dv"] = parm.dC_dv,
           Rcpp::_["prior"] = parm.logprob_prior(),
           Rcpp::_["gradient_prior"] = parm.gradient_prior(),
           Rcpp::_["hessian_prior"] = parm.hessian_prior());
}
