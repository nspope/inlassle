#include "Priors.hpp"

//---------------------------------------------------- tests

// [[Rcpp::export("inlassle_test_ImproperGamma")]]
Rcpp::List test_ImproperGamma (arma::vec x)
{
  ImproperGamma iga (x);
  return Rcpp::List::create (Rcpp::_["logprob"] = iga.logprob,
                             Rcpp::_["gradient"] = iga.gradient(),
                             Rcpp::_["hessian"] = iga.hessian(),
                             Rcpp::_["normalizing_constant"] = iga.normalizing_constant());
}

// [[Rcpp::export("inlassle_test_ImproperWishart")]]
Rcpp::List test_ImproperWishart (arma::vec x)
{
  ImproperWishart iwi (x);
  return Rcpp::List::create (Rcpp::_["logprob"] = iwi.logprob,
                             Rcpp::_["gradient"] = iwi.gradient(),
                             Rcpp::_["hessian"] = iwi.hessian(),
                             Rcpp::_["normalizing_constant"] = iwi.normalizing_constant());
}

// [[Rcpp::export("inlassle_test_Gaussian")]]
Rcpp::List test_Gaussian (arma::vec x, arma::vec mu, arma::mat Lambda)
{
  Gaussian gauss (x, mu, Lambda);
  return Rcpp::List::create (Rcpp::_["logprob"] = gauss.logprob,
                             Rcpp::_["gradient"] = gauss.gradient(),
                             Rcpp::_["hessian"] = gauss.hessian(),
                             Rcpp::_["normalizing_constant"] = gauss.normalizing_constant());
}

// [[Rcpp::export("inlassle_test_LKJ")]]
Rcpp::List test_LKJ (arma::vec x, double eta, arma::vec gamma)
{
  LKJ lkj (x, eta, gamma);
  return Rcpp::List::create (Rcpp::_["logprob"] = lkj.logprob,
                             Rcpp::_["gradient"] = lkj.gradient(),
                             Rcpp::_["hessian"] = lkj.hessian(),
                             Rcpp::_["normalizing_constant"] = lkj.normalizing_constant());
}

// [[Rcpp::export("inlassle_test_LKJ_jac")]]
Rcpp::List test_LKJ_jac (arma::vec x, double eta, arma::vec gamma)
{
  LKJ lkj (x, eta, gamma);
  return Rcpp::List::create (Rcpp::_["logprob"] = lkj.test_dS_dv(),
                             Rcpp::_["gradient"] = lkj.gradient(),
                             Rcpp::_["hessian"] = lkj.hessian(),
                             Rcpp::_["normalizing_constant"] = lkj.normalizing_constant());
}

// [[Rcpp::export("inlassle_test_LogBetaPrime")]]
Rcpp::List test_LogBetaPrime (arma::vec x, arma::vec mu, arma::vec delta)
{
  LogBetaPrime lbp (x, mu, delta);
  return Rcpp::List::create (Rcpp::_["logprob"] = lbp.logprob,
                             Rcpp::_["gradient"] = lbp.gradient(),
                             Rcpp::_["hessian"] = lbp.hessian(),
                             Rcpp::_["normalizing_constant"] = lbp.normalizing_constant());
}
