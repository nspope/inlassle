# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

inlassle_test_Field <- function(y, n, mu, Q) {
    .Call('_inlassle_test_Field', PACKAGE = 'inlassle', y, n, mu, Q)
}

inlassle_test_Gaussian_C <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Gaussian_C', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Gaussian_dC_dt <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Gaussian_dC_dt', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Gaussian_dC_dD <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Gaussian_dC_dD', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Likelihood <- function(N, Y, X, Z, D, t, v, s, b, parallel) {
    .Call('_inlassle_test_Likelihood', PACKAGE = 'inlassle', N, Y, X, Z, D, t, v, s, b, parallel)
}

inlassle_test_Likelihood_cov <- function(N, Y, X, Z, D, v, s, b, parallel) {
    .Call('_inlassle_test_Likelihood_cov', PACKAGE = 'inlassle', N, Y, X, Z, D, v, s, b, parallel)
}

inlassle_test_Matern_C <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Matern_C', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Matern_dC_dt <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Matern_dC_dt', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Matern_dC_dD <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Matern_dC_dD', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Null_C <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Null_C', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Null_dC_dt <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Null_dC_dt', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Null_dC_dD <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Null_dC_dD', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Parameters <- function(t, v, s, b) {
    .Call('_inlassle_test_Parameters', PACKAGE = 'inlassle', t, v, s, b)
}

inlassle_test_ImproperGamma <- function(x) {
    .Call('_inlassle_test_ImproperGamma', PACKAGE = 'inlassle', x)
}

inlassle_test_ImproperWishart <- function(x) {
    .Call('_inlassle_test_ImproperWishart', PACKAGE = 'inlassle', x)
}

inlassle_test_Gaussian <- function(x, mu, Lambda) {
    .Call('_inlassle_test_Gaussian', PACKAGE = 'inlassle', x, mu, Lambda)
}

inlassle_test_LKJ <- function(x, eta, gamma) {
    .Call('_inlassle_test_LKJ', PACKAGE = 'inlassle', x, eta, gamma)
}

inlassle_test_LKJ_jac <- function(x, eta, gamma) {
    .Call('_inlassle_test_LKJ_jac', PACKAGE = 'inlassle', x, eta, gamma)
}

inlassle_test_LogBetaPrime <- function(x, mu, delta) {
    .Call('_inlassle_test_LogBetaPrime', PACKAGE = 'inlassle', x, mu, delta)
}

inlassle_test_Problem_likelihood <- function(N, Y, X, Z, D, t, v, s, b, parallel) {
    .Call('_inlassle_test_Problem_likelihood', PACKAGE = 'inlassle', N, Y, X, Z, D, t, v, s, b, parallel)
}

inlassle_test_Problem_likelihood_cov <- function(N, Y, X, Z, D, v, s, b, parallel) {
    .Call('_inlassle_test_Problem_likelihood_cov', PACKAGE = 'inlassle', N, Y, X, Z, D, v, s, b, parallel)
}

inlassle_test_Problem_optimize <- function(N, Y, X, Z, D, f, parallel, tol) {
    .Call('_inlassle_test_Problem_optimize', PACKAGE = 'inlassle', N, Y, X, Z, D, f, parallel, tol)
}

inlassle_test_HomeierRule <- function(order, gamma) {
    .Call('_inlassle_test_HomeierRule', PACKAGE = 'inlassle', order, gamma)
}

inlassle_test_ResistanceOptim_likelihood <- function(pars, variable, N, Y, X, Z, D, S, T, A) {
    .Call('_inlassle_test_ResistanceOptim_likelihood', PACKAGE = 'inlassle', pars, variable, N, Y, X, Z, D, S, T, A)
}

inlassle_test_ResistanceOptim_optimize <- function(pars, variable, N, Y, X, Z, D, S, T, A, verbose) {
    .Call('_inlassle_test_ResistanceOptim_optimize', PACKAGE = 'inlassle', pars, variable, N, Y, X, Z, D, S, T, A, verbose)
}

inlassle_test_ResistanceOptim_optimize_global <- function(lb, ub, variable, N, Y, X, Z, D, S, T, A, verbose) {
    .Call('_inlassle_test_ResistanceOptim_optimize_global', PACKAGE = 'inlassle', lb, ub, variable, N, Y, X, Z, D, S, T, A, verbose)
}

inlassle_test_ResistanceOptim_grid <- function(pars, variable, N, Y, X, Z, D, S, T, A) {
    .Call('_inlassle_test_ResistanceOptim_grid', PACKAGE = 'inlassle', pars, variable, N, Y, X, Z, D, S, T, A)
}

testlink <- function() {
    invisible(.Call('_inlassle_testlink', PACKAGE = 'inlassle'))
}

testrd <- function(spdat, adj, targ, pars, diff, parallel) {
    .Call('_inlassle_testrd', PACKAGE = 'inlassle', spdat, adj, targ, pars, diff, parallel)
}

#' @export ResistanceSolver
NULL

inlassle_test_HessianBlock_inc1 <- function(inp) {
    .Call('_inlassle_inlassle_test_HessianBlock_inc1', PACKAGE = 'inlassle', inp)
}

inlassle_test_HessianBlock_inc2 <- function(inp1, inp2) {
    .Call('_inlassle_inlassle_test_HessianBlock_inc2', PACKAGE = 'inlassle', inp1, inp2)
}

inlassle_test_HessianBlock_join <- function(inp1, inp2) {
    .Call('_inlassle_inlassle_test_HessianBlock_join', PACKAGE = 'inlassle', inp1, inp2)
}

inlassle_test_GradientBlock_inc1 <- function(inp) {
    .Call('_inlassle_inlassle_test_GradientBlock_inc1', PACKAGE = 'inlassle', inp)
}

inlassle_test_GradientBlock_inc2 <- function(inp1, inp2) {
    .Call('_inlassle_inlassle_test_GradientBlock_inc2', PACKAGE = 'inlassle', inp1, inp2)
}

inlassle_test_GradientBlock_join <- function(inp1, inp2) {
    .Call('_inlassle_inlassle_test_GradientBlock_join', PACKAGE = 'inlassle', inp1, inp2)
}

inlassle_test_Gradient <- function(dt, dv, ds, db) {
    .Call('_inlassle_inlassle_test_Gradient', PACKAGE = 'inlassle', dt, dv, ds, db)
}

inlassle_test_Gradient_dC <- function(dC) {
    .Call('_inlassle_inlassle_test_Gradient_dC', PACKAGE = 'inlassle', dC)
}

inlassle_test_Hessian <- function(dt, dv, ds, db) {
    .Call('_inlassle_inlassle_test_Hessian', PACKAGE = 'inlassle', dt, dv, ds, db)
}

inlassle_test_Weighted_C <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Weighted_C', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Weighted_dC_dt <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Weighted_dC_dt', PACKAGE = 'inlassle', D, nu, delta, pars)
}

inlassle_test_Weighted_dC_dD <- function(D, nu, delta, pars) {
    .Call('_inlassle_test_Weighted_dC_dD', PACKAGE = 'inlassle', D, nu, delta, pars)
}

