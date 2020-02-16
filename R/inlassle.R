#' @import methods Rcpp raster
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib inlassle
NULL

Rcpp::loadModule("inlassle", TRUE)
#Rcpp::loadModule("randomfield", TRUE)

### tests
### TODO: document purpose ...
.test_ResistanceCov <- function(dim=5, seed=1, parallel=FALSE)
{#test covariance calculation
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  adj <- t(apply(adj,1,sort))
  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  pars <- c(1,1)
  cond <- exp(spd%*%pars)
  L <- matrix(0,dim*dim,dim*dim)
  for(i in 1:nrow(adj)) L[adj[i,1]+1,adj[i,2]+1] <- -0.5*cond[adj[i,1]+1] - 0.5*cond[adj[i,2]+1]
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  Linv <- MASS::ginv(L)
  E <- diag(dim*dim)[,targ+1]
  Lp <- t(E)%*%Linv%*%E 
  hi <- ResistanceSolver$new(spd, targ, adj, parallel)
  Lp2 <- hi$resistance_covariance_log(pars)
  list(Lp, Lp2)
}

.test_ResistanceCov2 <- function(dim=5, seed=1, parallel=FALSE)
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
  Rd2 <- hi$resistance_covariance_log(pars)
  Rd2
}

.test_ResistanceCov3 <- function(dim=5, seed=1, parallel=FALSE)
{#gradient calculation
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
  Rd2 <- hi$resistance_covariance_log(pars)
  dRd2 <- hi$rd_resistance_covariance_log(dif)
  emp <- as.vector(dif) %*% numDeriv::jacobian(hi$resistance_covariance_log, pars)
  list(Rd2, dRd2, emp)
}

.test_ResistanceSolver <- function(dim=5, seed=1, parallel=FALSE)
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

.get_ResistanceNL <- function(dim=5, pts=3, seed=1, do_grid=FALSE, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  if(do_grid)
    coords <- as.matrix(expand.grid(seq(0.2,0.8,length.out=pts),seq(0.2,0.8,length.out=pts)))
  else
    coords <- cbind(runif(pts,0.2,0.8),runif(pts,0.2,0.8))
  targ <- unique(cellFromXY(rr, coords)-1)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  adj <- t(apply(adj,1,sort))
  adj <- adj[!duplicated(adj),]
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))
  pars <- c(0,0)
  cond <- 1/plogis(spd%*%pars)
  L <- matrix(0,dim*dim,dim*dim)
  for(i in 1:nrow(adj)) L[adj[i,1]+1,adj[i,2]+1] <- -0.5*cond[adj[i,1]+1] - 0.5*cond[adj[i,2]+1]
  L <- L + t(L)
  diag(L) <- -rowSums(L)
  Linv <- MASS::ginv(L)
  E <- diag(dim*dim)[,targ+1]
  Lp <- t(E)%*%Linv%*%E
  Rd <- t(-2*Lp + diag(Lp)) + diag(Lp)
  list(Rd=Rd, E=E, Lp=Lp, Linv=Linv, r=rr, coords=coords)
}

.test_ResistanceSolverNL <- function(dim=5, pts=3, seed=1, parallel=FALSE)
{
  set.seed(seed)
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
#  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
  targ <- sample(0:(dim*dim-2), pts)
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

.test_ResistanceSolver2 <- function(dim=5, seed=1, parallel=FALSE)
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

.test_ResistanceSolver3 <- function(dim=5, seed=1, parallel=FALSE)
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

.test_ResistanceSolver3NL <- function(dim=5, seed=1, parallel=FALSE)
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

.test_RandomField <- function(seed=1, miss=FALSE)
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

.test_RandomField2 <- function(seed=1, loci=5, pop=3, miss=FALSE, parallel=FALSE)
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

.test_RandomField3 <- function(miss=FALSE, parallel=FALSE)
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
.get_test_set <- function(seed=1, loci=100, pop=10)
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

## rudimentary API
## TODO: clean, add documentation, proper tests

## TODO: handle missing data in FstFromCounts
FstFromCounts <- function(Y, N)
{
  if (!all(dim(Y)==dim(N)))
    stop("Dimension mismatch")
  if (any(Y<0) || any(N<0))
    stop("Cannot have negative counts")
  if (any(Y>N))
    stop("Allele counts cannot exceed number of haploids sampled")

  Fr <- Y/N
  f2 <- apply(apply(Fr, 2, function(x) outer(x,x,"-")^2), 1, mean, na.rm=TRUE)
  cornum <- apply(Fr*(1-Fr)/(N-1), 1, mean, na.rm=TRUE)
  corden <- apply(Fr*(1-Fr)*N/(N-1), 1, mean, na.rm=TRUE)
  fst <- (f2 - c(outer(cornum, cornum, "+")))/(f2 - c(outer(cornum, cornum, "+")) + c(outer(corden, corden, "+")))
  fst <- matrix(fst, nrow(Y), nrow(Y))
  diag(fst) <- 0
  fst
}

ResistanceSurface <- function(covariates, coords, directions=4, saveStack=TRUE)
{
  if (class(covariates) != "RasterStack")
    stop("Covariate input must be class 'RasterStack'")
  if (class(coords) != "SpatialPoints")
    stop("Coordinate input must be class 'SpatialPoints'")

  # share missing cells across layers
  spdat <- raster::getValues(covariates)
  spdat[is.na(rowSums(spdat)),] <- NA

  cat("Pruning disconnected components ...\n")
  for (i in 1:dim(covariates)[3])
  {
    raster::values(covariates[[i]]) <- spdat[,i]
    cr <- raster::clump(covariates[[i]], directions = directions)

    #remove secondary clumps
    connected_component <- names(which.max(table(raster::getValues(cr))))
    bad <- raster::getValues(cr) != as.integer(connected_component)
    spdat[bad,] <- NA
    raster::values(covariates[[i]]) <- spdat[,i]
  }

  # final check that graph is connected after pruning
  cr <- raster::clump(covariates[[1]], directions = directions)
  if (length(na.omit(unique(raster::getValues(cr)))) > 1)
    stop("Pruning failed, disconnected components remain")
  if (length(na.omit(unique(raster::getValues(cr)))) == 0)
    stop("Pruning failed, no non-missing cells remain")

  # extract raster "data"
  cat("Extracting adjacency list ...\n")
  adj <- raster::adjacent(covariates[[1]], which(!is.na(spdat[,1])), target=which(!is.na(spdat[,1])), directions=directions)

  # find cells of demes
  cells <- unmapped_cells <- raster::cellFromXY(covariates[[1]], coords)

  # remove NAs and remap indices of adjacency list to be contiguous
  bad     <- is.na(spdat[,1])
  map     <- cbind(1:length(which(!bad)), (1:raster::ncell(cr))[which(!bad)])
  spdat   <- spdat[!bad,,drop=FALSE]
  adj[,1] <- map[match(adj[,1], map[,2]),1] - 1
  adj[,2] <- map[match(adj[,2], map[,2]),1] - 1
  cells   <- map[match(cells, map[,2]),1] - 1

  cat("Finding duplicated coordinates and lumping into demes ...\n")
  uniq_cells     <- unique(cells)
  uniq_cells_map <- match(cells, uniq_cells)

  # check that cells lie on connected portion of raster
  if(any(is.na(getValues(cr)[unmapped_cells])))
    stop("At least one deme is located on a missing cell")

  out <- list("solver"       = ResistanceSolver$new(spdat, uniq_cells, adj, TRUE),
              "demes"        = unique(unmapped_cells),
              "coords2demes" = uniq_cells_map,
              "covariates"   = names(covariates), 
              "stack"        = if(saveStack) covariates else NULL)
  class(out) <- "inlassle_linear_system"
  out
}

simulate_inlassleBinomial <- function(linear_system, rpar, npar, nsnp, nchr, seed = 0)
{
  if (length(npar) != 4)
    stop ("Length of nuisance parameters must be 4")
  if (length(rpar) != length(linear_system$covariates))
    stop ("Length of resistance parameters must equal number of covariates in linear system")

  covariance <- linear_system$solver$resistance_covariance_log(rpar)

  random_intercepts_per_snp <- rnorm(nsnp, npar[4], exp(npar[2]))
  true_frequencies <- random_intercepts_per_snp + 
    MASS::mvrnorm(nsnp, rep(0,length(linear_system$coords2demes)), 
                  exp(2*npar[1])*covariance + diag(rep(exp(npar[3]),length(linear_system$coords2demes))))
  true_frequencies <- t(plogis(true_frequencies))
  sampled_chromosomes <- matrix(nchr, length(linear_system$demes), nsnp) 
  observed_frequencies <- matrix(rbinom(nsnp*length(linear_system$coords2demes), 
                                        size=c(sampled_chromosomes), 
                                        prob=c(true_frequencies)), 
                                 length(linear_system$coords2demes), nsnp) #allele counts

  list(snp_count=observed_frequencies,
       chr_count=sampled_chromosomes)
}

inlassleBinomial <- function(snp, chrom, linear_system, start=rep(0,length(linear_system$covariates)+4), maxit=100, parallel=TRUE)
{
  if (class(snp) != "matrix" || class(chrom) != "matrix" || class(linear_system) != "inlassle_linear_system")
    stop("Inputs must be: a matrix of allele counts, a matrix of sample sizes, an 'inlassle_linear_system' object")
  if (ncol(snp) < nrow(snp))
    stop("Fewer SNPs than demes")
  if (!(all(dim(snp)==dim(chrom))))
    stop("Dimension mismatch in SNP inputs")
  if (length(start) != length(linear_system$covariates) + 4)
    stop("Parameter start values are of incorrect length")
  if (nrow(snp) != length(linear_system$coords2demes))
    stop("Number of demes does not match size of linear system")
  if (any(snp > chrom))
    stop("SNP count cannot exceed number of sampled chromosomes")
  if (any(snp < 0 | chrom < 0))
    stop("SNP/chromosome counts cannot be negative")

  #aggregate SNP matrix by coord2deme
  snp2 <- matrix(0, length(linear_system$demes), ncol(snp))
  chrom2 <- matrix(0, length(linear_system$demes), ncol(snp))
  for(i in 1:nrow(snp))
  {
    snp2[linear_system$coords2demes[i],] <- snp2[linear_system$coords2demes[i],] + snp[i,]
    chrom2[linear_system$coords2demes[i],] <- chrom2[linear_system$coords2demes[i],] + chrom[i,]
  }
  snp <- snp2
  chrom <- chrom2

  ...grad <<- NA
  N <- ncol(snp)
  n <- nrow(snp)
  indr <- 1:length(linear_system$covariates)
  indn <- (max(indr)+1):(max(indr)+4)
  obj <- function(par)
  {
    msd  <- exp(2. * par[indn][1]) # scaling of covariance
    rcov <- msd * array(linear_system$solver$resistance_covariance_log(par[indr]), c(n,n,1))
    ll   <- inlassle:::inlassle_test_Likelihood_cov(chrom, snp, matrix(1,n,1), matrix(1,n,1), rcov, par[indn][2], par[indn][3], par[indn][4], FALSE)#TODO: if Field::linesearch throws warnings and parallel=TRUE this can really f*** up the stack
    rd   <- linear_system$solver$rd_resistance_covariance_log(msd * ll$gradient_distance[,,1]) #backpropagate
    ...grad <<- c(rd, sum(2. * rcov[,,1] * ll$gradient_distance[,,1]), ll$gradient) * N #save gradient for later use
    return (ll$loglik * N) #multiply by N because inlassle returns avg loglik
  }
  gra <- function(par)
  {
    if (any(is.na(...grad))) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }
  fit <- optim(start, fn=obj, gr=gra, method="L-BFGS-B", hessian=TRUE, control=list(maxit=maxit))

  converged <- fit$convergence == 0
  if (!converged)
    warning("Optimizer did not converge")

  hess <- fit$hessian

  # fitted values for resistance parameters
  rtable <- c()
  rtable <- cbind(rtable, "Est"  = fit$par[indr])
  rtable <- cbind(rtable, "SE"   = sqrt(diag(solve(hess)))[indr])
  rtable <- cbind(rtable, "Z"    = rtable[,"Est"]/rtable[,"SE"])
  rtable <- cbind(rtable, "pval" = 2.*(1.-pnorm(abs(rtable[,"Z"]))))
  rownames(rtable) <- linear_system$covariates

  # fitted values for nuisance parameters
  ntable <- cbind("Est" = fit$par[indn],
                  "SE"  = sqrt(diag(solve(hess)))[indn])
  rownames(ntable) <- c("log Std.Dev (Resistance)", "log Std.Dev (SNP)", "log Std.Dev (Deme)", "(Intercept)")

  rownames(hess) <- colnames(hess) <- c(rownames(rtable),rownames(ntable))

  # loglikelihood, etc.
  loglik <- -fit$value
  AIC    <- -2*loglik + 2*nrow(hess)

  list("Resistance parameters"=rtable, 
       "Nuisance parameters"=ntable, 
       "Hessian"=hess,
       "Loglik"=loglik, "AIC"=AIC, 
       "Converged"=converged)
}

inlassleBinomialGrid <- function(snp, chrom, linear_system, grid, maxit=100, parallel=TRUE)
{
  # calculate profile likelihood across a parameter grid

  if (class(snp) != "matrix" || class(chrom) != "matrix" || class(linear_system) != "inlassle_linear_system")
    stop("Inputs must be: a matrix of allele counts, a matrix of sample sizes, an 'inlassle_linear_system' object, a matrix of parameter values")
  if (ncol(snp) < nrow(snp))
    stop("Fewer SNPs than demes")
  if (!(all(dim(snp)==dim(chrom))))
    stop("Dimension mismatch in SNP inputs")
  if (ncol(grid) != length(linear_system$covariates))
    stop("Wrong number of parameters in grid")
  if (nrow(snp) != length(linear_system$coords2demes))
    stop("Number of demes does not match size of linear system")
  if (any(snp > chrom))
    stop("SNP count cannot exceed number of sampled chromosomes")
  if (any(snp < 0 | chrom < 0))
    stop("SNP/chromosome counts cannot be negative")

  #aggregate SNP matrix by coord2deme
  snp2 <- matrix(0, length(linear_system$demes), ncol(snp))
  chrom2 <- matrix(0, length(linear_system$demes), ncol(snp))
  for(i in 1:nrow(snp))
  {
    snp2[linear_system$coords2demes[i],] <- snp2[linear_system$coords2demes[i],] + snp[i,]
    chrom2[linear_system$coords2demes[i],] <- chrom2[linear_system$coords2demes[i],] + chrom[i,]
  }
  snp <- snp2
  chrom <- chrom2

  ...grad   <<- NA
  ...rgrad <<- NA
  N <- ncol(snp)
  n <- nrow(snp)
  indn <- 1:4
  obj <- function(par, rcov)
  {
    msd  <- exp(2. * par[indn][1]) # scaling of covariance
    rcov <- msd * rcov
    ll   <- inlassle:::inlassle_test_Likelihood_cov(chrom, snp, matrix(1,n,1), matrix(1,n,1), rcov, par[indn][2], par[indn][3], par[indn][4], FALSE)#TODO: if Field::linesearch throws warnings and parallel=TRUE this can really f*** up the stack
    ...grad <<- c(sum(2. * rcov[,,1] * ll$gradient_distance[,,1]), ll$gradient) * N #save gradient for later use
    # TODO: make "resistance gradient" optional, as it comes with some computational cost
    rd <- linear_system$solver$rd_resistance_covariance_log(msd * ll$gradient_distance[,,1]) #backpropagate
    ...rgrad <<- rd * N #save gradient with regard to resistance pars
    return (ll$loglik * N)
  }
  gra <- function(par, rcov)
  {
    if (any(is.na(...grad))) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }

  llout <- rep(NA, nrow(grid))
  grout <- matrix(NA, nrow(grid), ncol(grid))
  rdout <- array(NA, dim=c(n,n,nrow(grid)))
  cvout <- rep(NA, nrow(grid))
  start <- rep(0, length(indn))
  for (i in 1:nrow(grid))
  {
    rcov <- array(linear_system$solver$resistance_covariance_log(grid[i,]), c(n,n,1))
    fit  <- optim(start, fn=obj, gr=gra, rcov=rcov, method="L-BFGS-B", control=list(maxit=maxit))

    converged <- fit$convergence == 0
    if (!converged)
      warning("Optimizer did not converge")
    else
      start <- fit$par

    cvout[i] <- converged
    llout[i] <- -fit$value
    grout[i,] <- -...rgrad
    rdout[,,i] <- rcov
  }

  list("Parameters"=grid, "Convergence"=cvout, "Loglik"=llout, "Gradient"=grout, "ResistanceCovariance"=rdout)
}

inlassleMLPE <- function(gdist, linear_system, start=rep(0,length(linear_system$covariates)), maxit = 100)
{
  if (class(gdist) != "matrix" || class(linear_system) != "inlassle_linear_system")
    stop("Inputs must be: a matrix of genetic distances, an 'inlassle_linear_system' object")
  if (length(start) != length(linear_system$covariates))
    stop("Wrong number of parameters in start")
  if (nrow(gdist) != ncol(gdist))
    stop("Genetic distances must be a square matrix")
  if (any(gdist < 0) || any(diag(gdist)!=0))
    stop("Genetic distance matrix must be non-negative with zeros on diagonal")
  if (nrow(gdist) != length(linear_system$coords2demes))
    stop("Number of demes does not match size of linear system (perhaps some demes are in same raster cell?)")

  pairs <- which(lower.tri(gdist), arr.ind=TRUE)
  dis   <- gdist[lower.tri(gdist)]

  #gradient wrt -loglik is
  #  -coef[2] * solve(sigma^2 * corMatrix) %*% (y - fitted)
  #to get solve(..., y-fitted), consider that
  #  resid(fit, "response") returns r = y-fitted
  #so, simply need to call .recalc_cpp twice on r and multiply by sigma, coef

  jacc_prod <- function(dz, z, s) #jacobian product, mapping gradient of scaled to gradient of unscaled
    1/s * (dz - z/(length(z)-1) * c(crossprod(z, dz)) - mean(dz))

#  dd <- function(x) { n <- length(x); s <- sd(x); m <- mean(x); (x - m) / ((n-1)*s) }
#  dd2 <- function(x) { n <- length(x); s <- sd(x); m <- mean(x); -1/s^2 * (x - m) / ((n-1)*s) }
#  dd3 <- function(x) { n <- length(x); s <- sd(x); m <- mean(x);
#      diag(1/rep(s,n)) + outer(x, dd2(x)) -
#        outer(rep(1/s,n), rep(1/n, n)) - outer(rep(m,n), dd2(x)) }
#
#  jacc <- function(z, s, m) { n <- length(z);
#    diag(rep(1/s, n)) + outer(z*s, -1/((n-1)*s^2) * z) - matrix(1/(n*ss),n,n)
#  }

  ...grad  <<- NA
  ...model <<- NA
  obj <- function(par)
  {
    rdis <- linear_system$solver$resistance_distances_log(par)
    rdis <- rdis[lower.tri(rdis)]
    rsd  <- sd(rdis)
    rdis <- (rdis - mean(rdis)) / rsd
    dat  <- data.frame(y = gdist[lower.tri(gdist)], x = rdis, pop1 = pairs[,1], pop2 = pairs[,2])
    fit  <- nlme::gls(y ~ x, correlation = corMLPE(form=~pop1 + pop2), data = dat, method = "ML")

    # gradient calculation
    conLin <- list("Xy" = cbind(resid(fit, "response")), "logLik" = NA)
    conLin <- corMLPE:::recalc.corMLPE(fit$modelStruct$corStruct, conLin) #multiply by inverse sqrt once ...
    conLin <- corMLPE:::recalc.corMLPE(fit$modelStruct$corStruct, conLin) #... and twice
    grad   <- c(coef(fit)[2] * conLin[[1]] / sigma(fit)^2) #gradient wrt scaled RD
    grad   <- jacc_prod (grad, rdis, rsd) #gradient wrt unscaled RD

    tmp    <- matrix(0, nrow(gdist), ncol(gdist))
    tmp[lower.tri(tmp)] <- grad
    tmp    <- 0.5 * (tmp + t(tmp))

    ...grad  <<- -linear_system$solver$rd_resistance_distances_log(tmp)
    ...model <<- fit

    return(-logLik(fit))
  }
  gra <- function(par)
  {
    if (any(is.na(...grad))) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }
  fit <- optim(start, fn=obj, gr=gra, method="L-BFGS-B", hessian=TRUE, control=list(maxit=maxit))

  converged <- fit$convergence == 0
  if (!converged)
    warning("Optimizer did not converge")

  hess <- fit$hessian

  # fitted values for resistance parameters
  rtable <- c()
  rtable <- cbind(rtable, "Est"  = fit$par)
  rtable <- cbind(rtable, "SE"   = sqrt(diag(solve(hess))))
  rtable <- cbind(rtable, "Z"    = rtable[,"Est"]/rtable[,"SE"])
  rtable <- cbind(rtable, "pval" = 2.*(1.-pnorm(abs(rtable[,"Z"]))))
  rownames(rtable) <- linear_system$covariates

  rownames(hess) <- colnames(hess) <- rownames(rtable)

  # loglikelihood, etc.
  loglik <- -fit$value
  AIC    <- -2*loglik + 2*nrow(hess)
  
  list("Resistance parameters"=rtable, 
       "Fitted regression"=...model,
       "Hessian"=hess, 
       "Loglik"=loglik, "AIC"=AIC, 
       "Converged"=converged)
}

inlassleMLPEGrid <- function(gdist, linear_system, grid, maxit = 100)
{
  if (class(gdist) != "matrix" || class(linear_system) != "inlassle_linear_system" || class(grid) != "matrix")
    stop("Inputs must be: a matrix of genetic distances, an 'inlassle_linear_system' object, a matrix of parameter values")
  if (nrow(gdist) != ncol(gdist))
    stop("Genetic distances must be a square matrix")
  if (any(gdist < 0) || any(diag(gdist)!=0))
    stop("Genetic distance matrix must be non-negative with zeros on diagonal")
  if (ncol(grid) != length(linear_system$covariates))
    stop("Wrong number of parameters in grid")
  if (nrow(gdist) != length(linear_system$coords2demes))
    stop("Number of demes does not match size of linear system (perhaps some demes are in same raster cell?)")

  pairs <- which(lower.tri(gdist), arr.ind=TRUE)
  dis   <- gdist[lower.tri(gdist)]
  n     <- nrow(gdist)

  jacc_prod <- function(dz, z, s) #jacobian product, mapping gradient of scaled to gradient of unscaled
    1/s * (dz - z/(length(z)-1) * c(crossprod(z, dz)) - mean(dz))

  ...grad  <<- NA
  ...model <<- NA
  obj <- function(par, rdis)
  {
    dat  <- data.frame(y = gdist[lower.tri(gdist)], x = rdis, pop1 = pairs[,1], pop2 = pairs[,2])
    fit  <- nlme::gls(y ~ x, correlation = corMLPE(form=~pop1 + pop2), data = dat, method = "ML")

    # gradient calculation
    conLin <- list("Xy" = cbind(resid(fit, "response")), "logLik" = NA)
    conLin <- corMLPE:::recalc.corMLPE(fit$modelStruct$corStruct, conLin) #multiply by inverse sqrt once ...
    conLin <- corMLPE:::recalc.corMLPE(fit$modelStruct$corStruct, conLin) #... and twice
    grad   <- c(coef(fit)[2] * conLin[[1]] / sigma(fit)^2) #gradient wrt scaled RD
    grad   <- jacc_prod (grad, rdis, rsd) #gradient wrt unscaled RD

    tmp    <- matrix(0, nrow(gdist), ncol(gdist))
    tmp[lower.tri(tmp)] <- grad
    tmp    <- 0.5 * (tmp + t(tmp))

    ...grad  <<- -linear_system$solver$rd_resistance_distances_log(tmp)
    ...model <<- fit

    return (-logLik(fit))
  }
  gra <- function(par, rdis)
  {
    if (any(is.na(...grad))) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }

  llout <- rep(NA, nrow(grid))
  grout <- matrix(NA, nrow(grid), ncol(grid))
  rdout <- array(NA, dim=c(n,n,nrow(grid)))
  cvout <- rep(TRUE, nrow(grid))
  for(i in 1:nrow(grid))
  {
    rdis <- tmp <- linear_system$solver$resistance_distances_log(grid[i,])
    rdis <- rdis[lower.tri(rdis)]
    rsd  <- sd(rdis)
    rdis <- (rdis - mean(rdis)) / rsd
    llout[i] <- -obj(grid[i,], rdis)
    grout[i,] <- -...grad

    tmp[] <- 0
    tmp[lower.tri(tmp)] <- rdis
    tmp <- tmp + t(tmp)
    rdout[,,i] <- tmp
  }

  list("Parameters"=grid, "Convergence"=cvout, "Loglik"=llout, "Gradient"=grout, "ResistanceDistance"=rdout)
}

#TODO:
#inlassleWishart <- function(gdist, nsnp, linear_system)
#{
#  ...grad <- NA
#  obj <- function(par)
#}

#gradient for wishart with nu fixed
# L(V) = -0.5 tr(inv(V) X) - 0.5 nu logdet(V)
# L'(V) = 0.5 inv(V) X inv(V) - 0.5 nu inv(V)

#gradient for projected wishart with nu fixed
# L(V) = d/2 logdet(inv(V) Q) + d/4 tr(inv(V) Q D)
# Q = I - X inv(X' inv(V) X) X' inv(V)
# X is design matrix for mean ... if X = 1 then
# Q = I - 1/sum(inv(V)) X colsums(inv(V))

# if U are gradiants of f(a W 1 1' W)
# c(U) %*% (kronecker(rowSums(W) %*% t(rep(1,3)), diag(3)) + kronecker(diag(3),rowSums(W) %*% t(rep(1,3)))) * a
# ==>
# M = rowSums(W) %*% t(rep(1,3))
# U %*% M + t(M) %*% U
# ==>
# dL - a dL W 1 1' - a 1 1' W dL ... this is part holding a = 1/sum(W) fixed
# sum(dL * (a^2 W 1 1' W)) ... this is part holding W fixed (it is dxd matrix)
# ==>
# dW = dL - a dL W 1 1' - a 1 1' W dL + sum(dL * (a^2 W 1 1' W))
# but note that dL in above is the gradient of matrix wrt pseudo-determinant ... how to get this?
# conjecture: logDet(WQ) is always 1 .. .seems to be true. So only need to worry aobut grad of trace

### workspace
#cc = sum(U * (1/sum(W)^2 * W %*% matrix(1,3,1) %*% matrix(1,1,3) %*% W))
#bb = 1/sum(W) * U %*% W %*% matrix(1,3,1) %*% matrix(1,1,3) + 
#     1/sum(W) * matrix(1,3,1) %*% matrix(1,1,3) %*% W
# U - bb + cc
###


# fn <- function(V) log(det(solve(V) %*% (diag(nrow(V)) - 1/sum(solve(V)) * matrix(1,3,1) %*% matrix(1,1,3) %*% solve(V))))


