#ifndef RESISTANCESOLVER_H
#define RESISTANCESOLVER_H

#include "Inlassle.hpp"

// [[Rcpp::depends(RcppEigen,RcppParallel)]]
// [[Rcpp::plugins(cpp11)]]

using MatrixXd = Eigen::MatrixXd;
using MatrixXi = Eigen::MatrixXi;
using VectorXd = Eigen::VectorXd;
using VectorXi = Eigen::VectorXi;
using SpMatrix = Eigen::SparseMatrix<double>;
using SpSolver = Eigen::SimplicialLDLT<SpMatrix>;
using LuDecomp = Eigen::PartialPivLU<MatrixXd>;
using UiVector = std::vector<unsigned>;
using UiMap    = std::map<unsigned,unsigned>;

namespace Link
{

  /* link functions that map linear predictor to conductance */

  struct Log
  {
    double operator() (const double x) const
    {
      return exp(x);
    }

    double deriv (const double p) const
    {
      return p;
    }
  };

  struct Logit
  {
    double operator() (const double x) const
    {
      return 1. / (1. + exp(-x));
    }

    double deriv (const double p) const
    {
      return p * (1. - p);
    }
  };

  struct Identity
  {
    double operator() (const double x) const
    {
      return x;
    }

    double deriv (const double p) const
    {
      return 1;
    }
  };

  struct ReciprocalLogit
  {
    //equivalent to logit on resistances; seems to work better than just conductance
    double operator() (const double x) const
    {
      return 1. + exp(-x);
    }

    double deriv (const double p) const
    {
// d(1 + exp(-x))/dx = -exp(-x)
// p = 1 + e^-x
// p - 1 = e^-x
// -(p - 1) = -e^-x
      return -(p - 1);
    }
  };

  struct Softplus
  {
    double operator() (const double x) const
    {
      return 1. + log(1. + exp(x));
    }

    double deriv (const double p) const
    {
      return 1. - 1./exp(p - 1.);
    }
  };

} // namespace Link

//' @export ResistanceSolver
struct ResistanceSolver
{
  const MatrixXd    spatial_data;
  const unsigned    dim,
                    npars;
  const long double constant;
  const SpMatrix    targets;

  UiMap    ground;    // cells that are "grounded" to the missing node
  SpMatrix adjacency; // adjacency matrix for reduced graph (e.g. with one missing node)
  VectorXd rates,     // diagonal of Laplacian for reduced graph
           conductances,   // conductances per node
           rd_conductances,// differential per node (used in reverse-autodiff)
           rd_parameters,  // differential per parameter (used in reverse-autodiff)
           LinvCSum;  // currently not used, see deprecated code
  MatrixXd Linv,      // solve(Laplacian, ei)
           LinvD,     // differential %*% t(solve(Laplacian, ei))
           Rd,        // scaled resistance distances
           Rdinv;
  SpSolver solver;    // sparse direct solver
  double   logdet = 0.;    // log determinant (scaled by dimension) of rd matrix

  // settings
  bool parallel = false;
  bool use_iterative_solver = false;

  ResistanceSolver (const MatrixXd data, const UiVector& targ, const MatrixXi& adj, const bool parallel = false) :
    spatial_data (MatrixXd(data)),
    dim (spatial_data.rows()),
    npars (spatial_data.cols()),
    constant (1e-300 * double(dim)),
    targets (make_targets(targ)),
    adjacency (make_adjacency(adj)),
    rates (VectorXd::Ones(dim-1)),
    conductances (VectorXd::Ones(dim)),
    rd_conductances (VectorXd::Zero(dim)),
    rd_parameters (VectorXd::Zero(spatial_data.cols())),
    LinvCSum (VectorXd::Zero(targets.cols())),
    Linv (MatrixXd::Zero(dim-1, targets.cols())),
    LinvD (MatrixXd::Zero(targets.cols(), dim-1)),
    Rd (MatrixXd::Zero(targets.cols(), targets.cols())),
    Rdinv (MatrixXd::Zero(targets.cols(), targets.cols())),
    parallel (parallel)
  {
    if (adjacency.rows() != adjacency.cols() ||
        adjacency.rows() != dim-1            )
      Rcpp::stop ("Dimension mismatch");
    
    // symbolic Cholesky decomposition
    solver.analyzePattern(adjacency);
  }

  SpMatrix make_adjacency (const MatrixXi&);
  SpMatrix make_targets (const UiVector&);
  template <class LinkFn> MatrixXd resistance_distances (const VectorXd);
  template <class LinkFn> VectorXd rd_resistance_distances (MatrixXd);
  template <class LinkFn> MatrixXd resistance_covariance (const VectorXd);
  template <class LinkFn> VectorXd rd_resistance_covariance (MatrixXd);
  MatrixXd getAdjacency (void);
  VectorXd getLaplacianDiagonal (void);
  VectorXd getConductance (void);
  VectorXd getGradConductance (void);

  // wrappers for R API
  MatrixXd resistance_distances_log (const VectorXd);
  VectorXd rd_resistance_distances_log (MatrixXd);
  MatrixXd resistance_distances_logit (const VectorXd);
  VectorXd rd_resistance_distances_logit (MatrixXd);
  MatrixXd resistance_distances_rlogit (const VectorXd);
  VectorXd rd_resistance_distances_rlogit (MatrixXd);
  MatrixXd resistance_distances_softplus (const VectorXd);
  VectorXd rd_resistance_distances_softplus (MatrixXd);
  MatrixXd resistance_distances_identity (const VectorXd);
  VectorXd rd_resistance_distances_identity (MatrixXd);

  MatrixXd resistance_covariance_log (const VectorXd);
  VectorXd rd_resistance_covariance_log (MatrixXd);
  MatrixXd resistance_covariance_logit (const VectorXd);
  VectorXd rd_resistance_covariance_logit (MatrixXd);
  MatrixXd resistance_covariance_rlogit (const VectorXd);
  VectorXd rd_resistance_covariance_rlogit (MatrixXd);
  MatrixXd resistance_covariance_softplus (const VectorXd);
  VectorXd rd_resistance_covariance_softplus (MatrixXd);
  MatrixXd resistance_covariance_identity (const VectorXd);
  VectorXd rd_resistance_covariance_identity (MatrixXd);
};

template <class LinkFn>
struct map_data_to_conductance : public RcppParallel::Worker
{
  const MatrixXd &spatial_data;
  const VectorXd &pars;
        VectorXd &conductances;
        LinkFn   link;

  map_data_to_conductance (ResistanceSolver& rhs, const VectorXd& pars) :
    spatial_data (rhs.spatial_data),
    pars (pars),
    conductances (rhs.conductances)
  {
    conductances.setZero();
  }

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto j = begin; j != end; ++j)
      conductances(j) = link(spatial_data.row(j) * pars);
  }
};

template <class LinkFn>
struct map_conductance_to_parameters
{
  const MatrixXd &spatial_data;
  const VectorXd &rd_conductances,
                 &conductances;
        VectorXd _parameters;
        LinkFn   link;

  map_conductance_to_parameters (ResistanceSolver& rhs) :
    spatial_data (rhs.spatial_data),
    rd_conductances (rhs.rd_conductances),
    conductances (rhs.conductances),
    _parameters (VectorXd::Zero(spatial_data.cols()))
  {}

  map_conductance_to_parameters (map_conductance_to_parameters<LinkFn>& rhs, RcppParallel::Split) :
    spatial_data (rhs.spatial_data),
    rd_conductances (rhs.rd_conductances),
    conductances (rhs.conductances),
    _parameters (VectorXd::Zero(spatial_data.cols()))
  {}

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto j = begin; j != end; ++j)
    {
      double d = rd_conductances(j) * link.deriv(conductances(j));
      _parameters += d * spatial_data.row(j).transpose();
    }
  }

  void join (const map_conductance_to_parameters<LinkFn>& rhs)
  {
    _parameters += rhs._parameters;
  }
};

#endif /* RESISTANCESOLVER_H */
