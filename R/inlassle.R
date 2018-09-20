#' @import methods Rcpp raster
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib inlassle
NULL

Rcpp::loadModule("inlassle", TRUE)
#Rcpp::loadModule("randomfield", TRUE)

### tests
test_ResistanceSolver <- function(dim=5, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  adj <- t(apply(adj,1,sort))
  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  pars <- c(1,1)
  cond <- plogis(spd%*%pars)
  L <- matrix(0,dim*dim,dim*dim)
  for(i in 1:nrow(adj)) L[adj[i,1]+1,adj[i,2]+1] <- -0.5*cond[adj[i,1]+1] - 0.5*cond[adj[i,2]+1]
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  Linv <- MASS::ginv(L)
  E <- diag(dim*dim)[,targ+1]
  Lp <- t(E)%*%Linv%*%E
  Rd <- t(-2*Lp + diag(Lp)) + diag(Lp)
  Rd <- Rd*abs(det(Rd))^(-1/length(targ))
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  Rd2 <- hi$resistance_distances_logit(c(1,1))
#  testrd(spd, adj, targ, c(1,1), matrix(0,length(targ),length(targ))) #prints debugging information
  list(Rd, Rd2)
}

test_ResistanceSolverNL <- function(dim=5, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  adj <- t(apply(adj,1,sort))
  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  pars <- c(1,1)
  cond <- 1/plogis(spd%*%pars)
  L <- matrix(0,dim*dim,dim*dim)
  for(i in 1:nrow(adj)) L[adj[i,1]+1,adj[i,2]+1] <- -0.5*cond[adj[i,1]+1] - 0.5*cond[adj[i,2]+1]
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  Linv <- MASS::ginv(L)
  E <- diag(dim*dim)[,targ+1]
  Lp <- t(E)%*%Linv%*%E
  Rd <- t(-2*Lp + diag(Lp)) + diag(Lp)
  Rd <- Rd*abs(det(Rd))^(-1/length(targ))
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  Rd2 <- hi$resistance_distances_rlogit(c(1,1))
#  testrd(spd, adj, targ, c(1,1), matrix(0,length(targ),length(targ))) #prints debugging information
  list(Rd, Rd2)
}

test_ResistanceSolver2 <- function(dim=5, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  #adj <- cbind(pmin(adj[,1],adj[,2]), pmax(adj[,1],adj[,2]))
  #adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  pars <- c(1,1)
  Rd2 <- hi$resistance_distances_logit(pars)
  Rd2
}

test_ResistanceSolver3 <- function(dim=5, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
#  adj <- cbind(pmin(adj[,1],adj[,2]), pmax(adj[,1],adj[,2]))
#  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  dif <- matrix(rnorm(3*3), 3, 3); dif = dif + t(dif); diag(dif) = 0;
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  pars <- c(1,1)
  Rd2 <- hi$resistance_distances_logit(pars)
  dRd2 <- hi$rd_resistance_distances_logit(dif)
  emp <- as.vector(dif) %*% numDeriv::jacobian(hi$resistance_distances_logit, pars)
  list(Rd2, dRd2, emp)
}

test_ResistanceSolver3NL <- function(dim=5, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
#  adj <- cbind(pmin(adj[,1],adj[,2]), pmax(adj[,1],adj[,2]))
#  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  dif <- matrix(rnorm(3*3), 3, 3); dif = dif + t(dif); diag(dif) = 0;
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  pars <- c(1,1)
  Rd2 <- hi$resistance_distances_rlogit(pars)
  dRd2 <- hi$rd_resistance_distances_rlogit(dif)
  emp <- as.vector(dif) %*% numDeriv::jacobian(hi$resistance_distances_rlogit, pars)
  list(Rd2, dRd2, emp)
}

test_RandomField <- function(seed=1, miss=FALSE)
{
  set.seed(seed)

  pars <- c(1, 0.5, 1.2)
  dis  <- as.matrix(dist(cbind(runif(3),runif(3))))

  sp <- rep(1,3)
  sm <- rep(0,3)

  Y <- matrix(rbinom(5*3, size=10, prob=0.5),3,5)
  N <- matrix(10,3,5)
  if (miss) 
  {
    N[2,2:4] <- Y[2,2:4] <- 0
  }

  hi <- RandomField$new(Y, N, matrix(1,3,1), sp, sm, FALSE)
  hi$refit (pars, c(0.7), array(dis,dim=c(3,3,1)))
  other <- hi$get_field()

  xc <- other$covariance

  dbb <- function(y, n, a, b)
  {
    sum(ifelse(n>0,
           lbeta(y + a, n - y + b) - lbeta(a, b) + lchoose(n, y),
           0))
  }

  ll <- function(y, n, x, s, xc, xm, sp, sm)
  {
    miss <- n==0
    a = plogis(x)*exp(s)
    b = (1-plogis(x))*exp(s)
    xp <- solve(xc[!miss,!miss])
    0.5 * t(x[!miss] - xm) %*% xp %*% (x[!miss] - xm) - dbb(y, n, a, b) - 0.5*determinant(xp)$modulus
  }

  fullLL <- function(x, s, m)
  {
    x  <- matrix(x, nrow(Y), ncol(Y))
    lp <- 0
    for (i in 1:ncol(Y))
      lp <- lp + ll(Y[,i], N[,i], x[,i], s, xc=xc, xm=m, sp=sp, sm=sm)
    lp <- lp + sum(0.5*sp*(s-sm)^2)
    lp
  }

  fit   <- optim(rep(0, 6*3+1), function(par) fullLL(par[1:(5*3)], par[5*3+1:3], par[6*3+1]), method="BFGS")
  hess  <- numDeriv::hessian(function(par) fullLL(par[1:(5*3)], par[5*3+1:3], par[6*3+1]), fit$par)[!c(N==0),!c(N==0)]
  grad  <- numDeriv::grad(function(par) fullLL(par[1:(5*3)], par[5*3+1:3], par[6*3+1]), fit$par)

  loglik  <- fit$value
  logdet1 <- 0.5*determinant(hess)$modulus

  out <- list(loglik=loglik, logdet=logdet1, lapapp=loglik+logdet1, field=matrix(fit$par[1:15],3,5), dispersion=exp(fit$par[16:18]), fixef=fit$par[19], grad=grad)

  list(numderiv=out, inlassle=other)
}

test_RandomField2 <- function(seed=1, loci=5, pop=3, miss=FALSE, parallel=FALSE)
{
  set.seed(seed)

  pars <- c(1, 0.5, 1.2)
  dis  <- as.matrix(dist(cbind(runif(pop),runif(pop))))

  sp <- rep(1,pop)
  sm <- rep(0,pop)

  Y <- t(MASS::mvrnorm(loci, rep(0,pop), Sigma=exp(-dis)*0.7))
  Y <- matrix(rbeta(loci*pop, plogis(Y), 1-plogis(Y)), nrow(Y), ncol(Y))
  Y <- matrix(rbinom(loci*pop, size=10, prob=Y),nrow(Y),ncol(Y))
  N <- Y*0 + 10
  if (miss) 
  {
    gone <- sample.int(1:(loci*pop), loci*pop*miss)
    N[gone] <- 0
    Y[gone] <- 0
  }

  hi <- RandomField$new(Y, N, matrix(1,pop,1), sp, sm, parallel)
  hi$refit (pars, c(0.7), array(dis,dim=c(pop,pop,1)))
  other <- hi$get_field()

  # test
  fit <- optim(c(0,0,0), function(pars){ hi$refit(c(1,pars[1:2]), pars[3], array(dis,dim=c(pop,pop,1))) }, method="L-BFGS-B", control=list(trace=1), lower=c(0.01, 0.01, 0), upper=c(Inf, Inf, Inf))

  #
  
  return(list(init_field=other, fit=fit))
}

test_RandomField3 <- function(miss=FALSE, parallel=FALSE)
{
  set.seed(1)

  pars <- c(1, 0.5, 1.2)
  dis  <- as.matrix(dist(cbind(runif(3),runif(3))))

  sp <- rep(1,3)
  sm <- rep(0,3)

  Y <- matrix(rbinom(5*3, size=10, prob=0.5),3,5)
  N <- matrix(10,3,5)
  if (miss) 
    N[2,2:4] <- Y[2,2:4] <- 0

  hi <- RandomField$new(Y, N, matrix(1,3,1), sp, sm, parallel)
  hi$refit (pars, c(0.7), array(dis,dim=c(3,3,1)))
  emp <- numDeriv::grad(function(pars)
                        {
                 hi2 <- RandomField$new(Y, N, matrix(1,3,1), sp, sm, FALSE)
                 hi2$refit (pars, c(0.7), array(dis,dim=c(3,3,1)))
                        },
                 pars
                 )
  act <- hi$gradient_covariance_hyperparameters()

  emp2 <- numDeriv::grad(function(l)
                        {
                 hi2 <- RandomField$new(Y, N, matrix(1,3,1), sp, sm, FALSE)
                 hi2$refit (pars, l, array(dis,dim=c(3,3,1)))
                        },
                 c(0.7)
                 )

  emp3 <- matrix(numDeriv::grad(function(D)
                        {
                 hi2 <- RandomField$new(Y, N, matrix(1,3,1), sp, sm, FALSE)
                 hi2$refit (pars, c(0.7), array(D,dim=c(3,3,1)))
                        },
                 dis
                 ),3,3)
  diag(emp3) <- 0

  act3 <- (hi$gradient_distances())[,,1]
  act3 <- (act3  + t(act3))
  act3[upper.tri(act3)] <- 0
  
  return(list(list(emp,emp2), act, emp3, act3))
}

## test set
get_test_set <- function(seed=1, loci=100, pop=10)
{
  set.seed(seed)
  coord <- cbind(runif(pop),runif(pop))
  dis  <- as.matrix(dist(coord))
  rr <- raster(nrows=30, ncols=30, xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- cellFromXY(rr, coord) - 1
  if( any(duplicated(targ)) | any(targ==30*30-1) ) stop("bad")

  adj <- adjacent(rr, 1:(30*30), target=1:(30*30))-1

  sp <- rep(1,pop)
  sm <- rep(0,pop)

  spd <- scale(cbind(c(volcano[1:30,1:30]), c(volcano[31:60,31:60])))
  X <- matrix(1,pop,1)

  Y <- t(MASS::mvrnorm(loci, rep(0,pop), Sigma=exp(-dis)*0.7))
  Y <- matrix(rbeta(loci*pop, plogis(Y), 1-plogis(Y)), nrow(Y), ncol(Y))
  Y <- matrix(rbinom(loci*pop, size=10, prob=Y),nrow(Y),ncol(Y))
  N <- Y*0 + 10

  return(list(spd=spd,adj=adj,targ=targ, dis=array(dis,c(nrow(Y),nrow(Y),1)), Y=Y, N=N, X=X, sp=sp, sm=sm))

}
