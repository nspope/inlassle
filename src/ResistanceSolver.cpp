#include "ResistanceSolver.hpp"

//double link_fn(const double x)
//{
//  // use a logistic link as this helps prevent the Laplacian from getting ill-conditioned
//  return 1. / (1. + exp(-x));
//}
//
//double link_fn_deriv(const double p)
//{
//  return p * (1. - p);
//}
              
/* parallel workhorses */

//struct map_data_to_conductance : public RcppParallel::Worker
//{
//  const MatrixXd &spatial_data;
//  const VectorXd &pars;
//        VectorXd &conductances;
//
//  map_data_to_conductance (ResistanceSolver& rhs, const VectorXd& pars) :
//    spatial_data (rhs.spatial_data),
//    pars (pars),
//    conductances (rhs.conductances)
//  {
//    conductances.setZero();
//  }
//
//  void operator() (std::size_t begin, std::size_t end)
//  {
//    for (auto j = begin; j != end; ++j)
//      conductances(j) = link_fn(spatial_data.row(j) * pars);
//  }
//};

// explicit instantiation
//template class map_data_to_conductance<Link::Logit>;
//template class map_data_to_conductance<Link::Identity>;

struct map_conductance_to_adjacency : public RcppParallel::Worker
{ 
  const VectorXd &conductances;
        VectorXd _rates;
        SpMatrix &adjacency;

  map_conductance_to_adjacency (ResistanceSolver& rhs) :
    conductances (rhs.conductances),
    _rates (VectorXd::Zero(rhs.rates.rows())),
    adjacency (rhs.adjacency)
  {
    // ground the Laplacian
    for (auto& i : rhs.ground)
      _rates(i.first) += (conductances(rhs.dim-1) + conductances(i.first))/2.;
  }

  map_conductance_to_adjacency (map_conductance_to_adjacency& rhs, RcppParallel::Split) :
    conductances (rhs.conductances),
    _rates (VectorXd::Zero(rhs._rates.rows())),
    adjacency (rhs.adjacency)
  {}

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto j = begin; j != end; ++j)
    {
      SpMatrix::InnerIterator i(adjacency, j); ++i; //skip diagonal
      for (i; i; ++i)
      {
        // map spatial data onto conductance
        i.valueRef() = -(conductances(i.col()) + conductances(i.row()))/2.;

        // subtract conductance for edge from Laplacian diagonal
        _rates(i.col()) -= i.value();
        _rates(i.row()) -= i.value();
      }
    }
  }

  void join (const map_conductance_to_adjacency& rhs)
  {
    _rates += rhs._rates;
  }
};

struct map_gradient_to_cells : public RcppParallel::Worker
{
  const MatrixXd &Linv;
  const VectorXd &LinvCSum;
        MatrixXd &LinvD;
        VectorXd &rd_conductances;

  map_gradient_to_cells (ResistanceSolver& rhs, const MatrixXd& Diff) :
    Linv (rhs.Linv),
    LinvCSum (rhs.LinvCSum),
    LinvD (rhs.LinvD),
    rd_conductances (rhs.rd_conductances)
  {
    LinvD = -Diff * Linv.transpose(); // possibly could be parallelized ...
    rd_conductances.setZero();
  }

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto j = begin; j != end; ++j)
      rd_conductances (j) = 0.5 * Linv.row(j) * LinvD.col(j);
  }
};

struct map_cells_to_conductance : public RcppParallel::Worker
{

  const SpMatrix &adjacency;
  const MatrixXd &Linv,
                 &LinvD;
  const VectorXd &rd_conductances;
        VectorXd _conductance;

  map_cells_to_conductance (ResistanceSolver& rhs) :
    adjacency (rhs.adjacency),
    Linv (rhs.Linv),
    LinvD (rhs.LinvD),
    rd_conductances (rhs.rd_conductances),
    _conductance (VectorXd::Zero(rhs.dim))
  {}

  map_cells_to_conductance (map_cells_to_conductance& rhs, RcppParallel::Split) :
    adjacency (rhs.adjacency),
    Linv (rhs.Linv),
    LinvD (rhs.LinvD),
    rd_conductances (rhs.rd_conductances),
    _conductance (VectorXd::Zero(rhs._conductance.rows()))
  {}

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto j = begin; j != end; ++j)
    {
      SpMatrix::InnerIterator i(adjacency, j); ++i; //skip diagonal
      for (; i; ++i)
      {
        // visually:
        //   -nonzeros on row/col k, +diagonal elements corresponding to nonzeros, +diagonal k,k (precomputed)
        double tmp = Linv.row(i.col()) * LinvD.col(i.row()); // this *requires* D to be symmetric so that inner project is symmetric in ordering of (i,j)
        _conductance(i.col()) += -tmp + rd_conductances(i.row()) + rd_conductances(i.col()); 
        _conductance(i.row()) += -tmp + rd_conductances(i.col()) + rd_conductances(i.row());
      }
    }
  }

  void join (const map_cells_to_conductance& rhs)
  {
    _conductance += rhs._conductance;
  }
};

//struct map_conductance_to_parameters
//{
//  const MatrixXd &spatial_data;
//  const VectorXd &rd_conductances,
//                 &conductances;
//        VectorXd _parameters;
//
//  map_conductance_to_parameters (ResistanceSolver& rhs) :
//    spatial_data (rhs.spatial_data),
//    rd_conductances (rhs.rd_conductances),
//    conductances (rhs.conductances),
//    _parameters (VectorXd::Zero(spatial_data.cols()))
//  {}
//
//  map_conductance_to_parameters (map_conductance_to_parameters& rhs, RcppParallel::Split) :
//    spatial_data (rhs.spatial_data),
//    rd_conductances (rhs.rd_conductances),
//    conductances (rhs.conductances),
//    _parameters (VectorXd::Zero(spatial_data.cols()))
//  {}
//
//  void operator() (std::size_t begin, std::size_t end)
//  {
//    for (auto j = begin; j != end; ++j)
//    {
//      double d = rd_conductances(j) * link_fn_deriv(conductances(j));
//      _parameters += d * spatial_data.row(j).transpose();
//    }
//  }
//
//  void join (const map_conductance_to_parameters& rhs)
//  {
//    _parameters += rhs._parameters;
//  }
//};

// explicit instantiation
//template class map_conductance_to_parameters<Link::Identity>;
//template class map_conductance_to_parameters<Link::Logit>;

struct solve_laplacian : public RcppParallel::Worker
{
  const unsigned dim;
  const long double constant;
  const SpSolver &solver;
  const SpMatrix &targets;
        MatrixXd &Linv,
                 &Lp;

  solve_laplacian (ResistanceSolver& rhs, MatrixXd& Lp) :
    dim (rhs.dim),
    constant (rhs.constant),
    solver (rhs.solver),
    targets (rhs.targets),
    Linv (rhs.Linv),
    Lp (Lp)
  {}

  void operator() (std::size_t begin, std::size_t end)
  {
    for (auto i = begin; i != end; ++i)
    {
      VectorXd ei  = VectorXd::Constant (dim-1, -1./double(dim));//*constant
      ei          += targets.col(i);//*constant
      Linv.col(i)  = solver.solve(ei);
      double tmp   = Linv.col(i).sum() / double(dim);
      Lp.col(i)    = targets.transpose() * Linv.col(i) - VectorXd::Ones(targets.cols()) * tmp;
    }
  }
};

/* member functions */
SpMatrix ResistanceSolver::make_adjacency (const MatrixXi& adj)
  /* construct the adjacency matrix for the graph */
{
  VectorXi a = adj.rowwise().maxCoeff(),
           b = adj.rowwise().minCoeff();
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve (adj.rows() + dim - 1);
  for (unsigned i=0; i<adj.rows(); ++i)
  {
    if (a(i) > dim-1 || b(i) < 0)
      Rcpp::stop ("ResistanceSolver: adjacency list out of bounds");

    if (a(i)==dim-1) // ground
    { 
      ground[b(i)] = dim-1; 
    }
    else
    {
      triplets.push_back (
          Eigen::Triplet<double>(
            a(i), // lower triangle is filled, so that the cell with the
            b(i), // lesser index is used as a row.
            1.)
          );
    }
  }
  for (unsigned i=0; i<dim-1; ++i) // diagonal
    triplets.push_back (Eigen::Triplet<double>(i, i, 1.));
  SpMatrix out (dim-1, dim-1);
  out.setFromTriplets(triplets.begin(), triplets.end(), [] (const double&, const double &b) { return b; }); // ensures that only a single duplicate survives
  out.makeCompressed();
  return out;
}

SpMatrix ResistanceSolver::make_targets (const UiVector& targ)
{
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(targ.size());
  for (unsigned i=0; i<targ.size(); ++i)
  {
    if (targ.at(i) > dim - 1 || targ.at(i) < 0)
      Rcpp::stop("ResistanceSolver: targets out of bounds");
    if (targ.at(i) == dim - 1) 
      Rcpp::stop("ResistanceSolver: cannot use last vertex of graph as a target");

    triplets.push_back (
        Eigen::Triplet<double>(targ.at(i), i, 1.));
  }
  SpMatrix out (dim-1, targ.size());
  out.setFromTriplets(triplets.begin(), triplets.end(), [] (const double&, const double &b) { return b; }); // ensures that only a single duplicate survives
  out.makeCompressed();
  return out;
}

template <class LinkFn>
MatrixXd ResistanceSolver::resistance_distances (const VectorXd pars)
{
  // calculate conductances
  map_data_to_conductance<LinkFn> map1 (*this, pars);
  if (parallel)
    RcppParallel::parallelFor (0, dim, map1);
  else
    map1 (0, dim);

  // compute cell values in the Laplacian
  map_conductance_to_adjacency map2 (*this);
  if (parallel)
    RcppParallel::parallelReduce (0, dim-1, map2);
  else
    map2 (0, dim-1);
  rates = map2._rates;
  adjacency.diagonal() = rates;

  // factorize Laplacian
  solver.factorize(adjacency); 

  // solve multiple right hand sides: 
  // this should be parallelized, but the question is whether the solver can be copied.
  MatrixXd Lp = MatrixXd::Zero(targets.cols(), targets.cols());
  solve_laplacian solver (*this, Lp);
  if (parallel)
    RcppParallel::parallelFor (0, targets.cols(), solver);
  else
    solver (0, targets.cols());
//  for (unsigned i=0; i<targets.cols(); ++i)
//  {
//    // I'm leaving the "constant" used by gdistance package out ... perhaps there's a better way to control overflow
//    VectorXd ei  = VectorXd::Constant (dim-1, -1./double(dim));//*constant
//    ei          += targets.col(i);//*constant
//    Linv.col(i)  = solver.solve(ei);
//    double tmp   = Linv.col(i).sum() / double(dim);
//    Lp.col(i)    = targets.transpose() * Linv.col(i) - VectorXd::Ones(targets.cols()) * tmp;
//  }
//  Lp          /= constant;
  Rd            = -2.*Lp;
  Rd.rowwise() += Lp.diagonal().transpose();
  Rd.colwise() += Lp.diagonal();

  // rescale by absolute value of determinant;
  // for distance matrices the sign depends on whether there is an odd or even number of entries.
  // we do calculations in a log-space via the LU decomposition for stability (as entries in Rd can be very large)
  logdet = 0;
  LuDecomp lu(Rd);
  auto& lu_mat = lu.matrixLU();
  for (unsigned i=0; i < lu_mat.rows(); ++i) 
    logdet += log(std::fabs(lu_mat(i,i)));
  logdet *= -1./double(lu_mat.rows());
  for (unsigned i=0; i < lu_mat.rows(); ++i) 
    for (unsigned j=i+1; j < lu_mat.rows(); ++j)
      Rd(i,j) = Rd(j,i) = exp(log(Rd(j,i)) + logdet);
  Rdinv = lu.inverse();

  return Rd;
}

template <class LinkFn>
VectorXd ResistanceSolver::rd_resistance_distances (MatrixXd D)
{
  if (D.cols() != Rd.cols() || D.rows() != Rd.rows())
    Rcpp::stop("Dimension mismatch");

  // make sure D is symmetric ...
  D += D.transpose().eval();
  D /= 2.;

  // d = abs(det(R))^(-1/n), thus
  //   d(d)/dR = -1/n * abs(det(R))^(-1/n) * solve(R)
  // if passing back a matrix of sensitivities D for f(d) = A = X * d, then
  //   df/R_ij = d(d)/d(R_ij) * sum(X * D), e.g. we have D' = -1/n * solve(R) * abs(det(R))^(-1/n) * sum(X * D)
  // now, if X = R, then 
  //   by product rule we have update D' = -1/n * solve(R) * abs(det(R))^(-1/n) * sum(X * D) + D * abs(det(R))^(-1/n)
  //                                     = abs(det(R))^(-1/n) * D - solve(R) * sum(A * D)/n)
  MatrixXd D1 = exp(logdet) * D - Rdinv * (Rd.cwiseProduct(D)).sum()/double(Rd.rows());

  // R = -2*L + ones%*%diag(L)' + diag(L)%*%ones'
  MatrixXd D2 = -2*D1;
  D2.diagonal() += 2.*D1.rowwise().sum(); //b/c symmetric

  // Ei = c*(E - ones ones'/n)
  // E.t() * (inv(L) * Ei - ones ones'/n inv(L) * Ei)
  // (E.t() - ones ones'/n) * (inv(L) * Ei)
  // 1/c * Ei.t() inv(L) Ei
  // meaning the differential is, 
  // -1/c * inv(L) * Ei * D * Ei.t() * inv(L)
  // -1/c * solve(L, Ei) %*% D %*% t(solve(L, Ei)) // the problem is this will underflow, probably
  //
  // Note that we haven't stored Xi=solve(L,Ei)/c directly b/c of overflow/underflow. Instead we have
  //   Linv.col(i) = (Xi - sum(Xi/dim))/c = (Xi - LinvCSum(i))/c
  // Problem is ... Xi/c ==> overflow, 1/(c*dim) ===> underflow, Xi.sum() ===> overflow
  //   Xi/c = Linv.col(i) + LinvCSum(i)/c
  // Is this really true? Rather, it seems to me that Xi/c *should* be OK.
  //   solve(L, n/v*(e - 1/n)) ==> n/v * solve(L, e - 1/n)
  //   and I think n/v is chosen so that n << v. Thus n/v will be *small*.
  //   then we have
  //     n/v * solve(L, e - 1/n) - sum(n/v * solve(L, e - 1/n))/n =
  //     n/v * solve(L, e - 1/n) - 1/v * sum(solve(L, e - 1/n))   =
  //     n/v * solve(L, e - 1/n) - n/v * mean(solve(L, e- 1/n)) 
  //   then we divide by n/v, so as to get
  //     solve(L, e-1/n) - mean(L, e-1/n)
  //   reading the gdistance notes, it seems to me that the issue are the entries of Lp being too large.
  //   recall that Lp = Ei.t() * Linv * Ei. So the issue is that (solve(L, e-1/n) - mean(L, e-1/n)) may be too big.
  //
  // clearly this is a large, dense matrix. The differential at a cell, say (i,j), is
  //  solve(L, Ei).row(i) %*% D %*% solve(L, Ei).col(i);
  // thus we can set, 
  //  Q = solve(L, Ei)/c %*% D (parallelize this operation; say over blocks?)
  // and set
  //  d(cell_ij) = Q.row(i) * solve(L, Ei).col(i);
  // then loop over cells of sparse matrix, and the contribution to the differential for the
  // conductance from cell (k,l) is
  //  d(conductance_k) = d(conductance_l) = -0.5*d(cell_kl)
  // and from cell kk is
  //  d(conductance_k/neighbours of k) += 0.5*d(cell_kk)
  // then we can map back to the parameters.

  // compute LinvD = sensitivities * Linv * Ei.t(), then calculate ...
  map_gradient_to_cells map1 (*this, D2);
  if (parallel)
    RcppParallel::parallelFor (0, dim-1, map1);
  else
    map1 (0, dim-1);

  // loop over non-zeros of adjacency matrix, calculate differential for each element;
  // and then pass back to vector of differential for conductances
  map_cells_to_conductance map2 (*this); 
  if (parallel)
    RcppParallel::parallelReduce (0, dim-1, map2);
  else
    map2 (0, dim-1);

  // the conductance of edges leading to the "dropped" node "j" must also be included.
  // for neighbour i of j,
  //   L(i,i) += 0.5*c(i) + 0.5*c(j)
  // so we have,
  //   d(c(i)) += 0.5 * d(L(i,i))
  //   d(c(j)) += 0.5 * d(L(i,i))
  // where d(L(i,i)) = Linv.row(i) * LinvD.col(i); this is just "rd_conductance" before adding contribution from map2.
  for (auto& i : ground)
  {
    map2._conductance (dim-1)   += rd_conductances(i.first);
    map2._conductance (i.first) += rd_conductances(i.first);
  }

  // copy over from Worker
  rd_conductances = map2._conductance;

  // finally map back onto the parameters; given
  //  f = exp(x' %*% p)
  // then
  //  d(f)/d(x'p) d(x'p)/d(p) = f*x
  map_conductance_to_parameters<LinkFn> map3 (*this);
  if (parallel)
    RcppParallel::parallelReduce (0, dim, map3);
  else
    map3 (0, dim);

  return map3._parameters;
}

MatrixXd ResistanceSolver::getAdjacency (void)
{
  return MatrixXd(adjacency); //at some point can we return sparse matrix?
}

VectorXd ResistanceSolver::getLaplacianDiagonal (void)
{
  return rates;
}

VectorXd ResistanceSolver::getConductance (void)
{
  return conductances;
}

VectorXd ResistanceSolver::getGradConductance (void)
{
  return rd_conductances;
}

// wrappers around link functions for R API
// TODO I would love to figure out a cleaner way to do this
MatrixXd ResistanceSolver::resistance_distances_logit (const VectorXd input)
{
  return resistance_distances<Link::Logit> (input);
}

VectorXd ResistanceSolver::rd_resistance_distances_logit (MatrixXd input)
{
  return rd_resistance_distances<Link::Logit> (input);
}

MatrixXd ResistanceSolver::resistance_distances_logit_sill (const VectorXd input)
{
  return resistance_distances<Link::LogitSill> (input);
}

VectorXd ResistanceSolver::rd_resistance_distances_logit_sill (MatrixXd input)
{
  return rd_resistance_distances<Link::LogitSill> (input);
}

MatrixXd ResistanceSolver::resistance_distances_softplus_sill (const VectorXd input)
{
  return resistance_distances<Link::SoftplusSill> (input);
}

VectorXd ResistanceSolver::rd_resistance_distances_softplus_sill (MatrixXd input)
{
  return rd_resistance_distances<Link::SoftplusSill> (input);
}

MatrixXd ResistanceSolver::resistance_distances_identity (const VectorXd input)
{
  return resistance_distances<Link::Identity> (input);
}

VectorXd ResistanceSolver::rd_resistance_distances_identity (MatrixXd input)
{
  return rd_resistance_distances<Link::Identity> (input);
}

// [[Rcpp::export]]
void testlink ()
{
  Link::Identity id_link;
  Link::Logit logit_link;
  Rcpp::Rcout << id_link(2.) << " " << id_link.deriv(2.) << std::endl;
  Rcpp::Rcout << logit_link(2.) << " " << logit_link.deriv(2.) << std::endl;
}

// [[Rcpp::export]]
Eigen::MatrixXd testrd (Eigen::MatrixXd spdat, Eigen::MatrixXi adj, std::vector<unsigned> targ, Eigen::VectorXd pars, Eigen::MatrixXd diff, bool parallel)
{
  //deep copy ...
  auto spdat2 = MatrixXd(spdat);
  auto diff2 = MatrixXd(diff);
  auto pars2 = VectorXd(pars);
  ResistanceSolver hi (spdat2, targ, adj, parallel);
  auto out = hi.resistance_distances<Link::Logit>(pars2);
//  std::cout << out << std::endl << std::endl;
//  std::cout << MatrixXd(hi.adjacency) << std::endl << std::endl;
//  std::cout << hi.rates << std::endl << std::endl;
//  for (auto i : hi.ground) std::cout << i << " "; std::cout << std::endl << std::endl;
  auto out2 = hi.rd_resistance_distances<Link::Logit>(diff2);
//  std::cout << out2 << std::endl;
  return out;
}

// next steps:
//  -does it work for general graphs? Yes
//  -does it "overflow" if graph is too large? ???
//  -does it work in parallel? Yes

RCPP_EXPOSED_CLASS(ResistanceSolver)
RCPP_MODULE(inlassle) {
  using namespace Rcpp;
  class_<ResistanceSolver>("ResistanceSolver")
    .constructor<Eigen::MatrixXd, std::vector<unsigned>, Eigen::MatrixXi, bool>()
    .method("resistance_distances_logit", &ResistanceSolver::resistance_distances_logit)
    .method("resistance_distances_logit_sill", &ResistanceSolver::resistance_distances_logit_sill)
    .method("resistance_distances_softplus_sill", &ResistanceSolver::resistance_distances_softplus_sill)
    .method("resistance_distances_identity", &ResistanceSolver::resistance_distances_identity)
    .method("rd_resistance_distances_logit", &ResistanceSolver::rd_resistance_distances_logit)
    .method("rd_resistance_distances_logit_sill", &ResistanceSolver::rd_resistance_distances_logit_sill)
    .method("rd_resistance_distances_softplus_sill", &ResistanceSolver::rd_resistance_distances_softplus_sill)
    .method("rd_resistance_distances_identity", &ResistanceSolver::rd_resistance_distances_identity)
    .method("getAdjacency", &ResistanceSolver::getAdjacency)
    .method("getLaplacianDiagonal", &ResistanceSolver::getLaplacianDiagonal)
    .method("getConductance", &ResistanceSolver::getConductance)
    .method("getGradConductance", &ResistanceSolver::getGradConductance)
    ;
}
