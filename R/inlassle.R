#' @import methods Rcpp raster
#' @importFrom RcppParallel RcppParallelLibs
#' @useDynLib inlassle
NULL

Rcpp::loadModule("inlassle", TRUE)
#Rcpp::loadModule("randomfield", TRUE)

### tests
test_ResistanceCov <- function(dim=5, seed=1, parallel=FALSE)
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

test_ResistanceCov2 <- function(dim=5, seed=1, parallel=FALSE)
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

test_ResistanceCov3 <- function(dim=5, seed=1, parallel=FALSE)
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

#
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

get_ResistanceNL <- function(dim=5, pts=3, seed=1, do_grid=FALSE, parallel=FALSE)
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

test_ResistanceSolverNL <- function(dim=5, pts=3, seed=1, parallel=FALSE)
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

### covariance model on 101x101 grid
#coordinates start at 0
#deme grid is thus (17 33 50 67 83) - 1
setup_splatche <- function(Y, N, spd, k, parallel=FALSE)
{
  rr <- raster::raster(matrix(spd[,1], 101, 101),
                       xmn=0, xmx=100, ymn=0, ymx=100)
  coords <- expand.grid(Y=c(17,33,50,67,83),X=c(17,33,50,67,83))
  coords <- coords[,2:1] #b/c/ splatche does Lat then Lon
  plot(rr)
  points(coords)

  targ <- cellFromXY(rr, coords[,1:2]) - 1

  adj <- raster::adjacent(rr, 1:(101*101))-1

  res <- ResistanceSolver$new(spd, targ, adj, parallel)
  list(Y=Y, N=N, targ=targ, adj=adj, spd=spd, res=res, k=k,parallel=parallel)
}

sim_splatche <- function(par, spd)
{
  rr <- raster::raster(matrix(spd[,1], 101, 101),
                       xmn=0, xmx=100, ymn=0, ymx=100)
  coords <- expand.grid(Y=c(17,33,50,67,83),X=c(17,33,50,67,83))
  coords <- coords[,2:1] #b/c/ splatche does Lat then Lon

  targ <- cellFromXY(rr, coords[,1:2]) - 1

  adj <- raster::adjacent(rr, 1:(101*101))-1

  res <- ResistanceSolver$new(spd, targ, adj, parallel)
  cv <- res$resistance_covariance_log(par)
  cv <- solve(cv + 1.0 + diag(rep(1.0,25)))
  
  list(Y=Y, N=N, targ=targ, adj=adj, spd=spd, res=res, k=k,parallel=parallel)
}

obj_splatche <- function(spat, pars)
{
  #assume that spd includes intercept!
  X <- cbind(spat$spd[spat$targ+1,spat$k,drop=F])
  Z <- cbind(spat$spd[spat$targ+1,spat$k,drop=F])

  nrp <- ncol(spat$spd)
  nvp <- choose(ncol(X), 2) + ncol(X) 
  nsp <- ncol(Z)
  nbp <- ncol(X)

  r <- pars[1:nrp]
  v <- pars[(nrp+1):(nrp+nvp)]
  s <- pars[(nrp+nvp+1):(nrp+nvp+nsp)]
  b <- pars[(nrp+nvp+nsp+1):(nrp+nvp+nsp+nbp)]

  Rd <- array(spat$res$resistance_covariance_log(r),dim=c(nrow(spat$Y),nrow(spat$Y),1))
  LL <- inlassle_test_Likelihood_cov(spat$N, spat$Y, X, Z, Rd, v, s, b, spat$parallel)

  rd <- spat$res$rd_resistance_covariance_log(LL$gradient_distance[,,1])
  der <- c(rd, LL$gradient)
  list(ll=LL$loglik, gr=der, pr=LL$Q)
}

fit_splatche <- function(Y, N, spd, start, k, method="BFGS", parallel=FALSE)
{
  stp <- setup_splatche(Y, N, spd, k, parallel)
  obj <- function(par) {print(par);fit<-obj_splatche(stp, par);..grad <<- fit$gr; return(fit$ll)}
  grd <- function(par) {return(..grad)}
  optim(start, obj, gr=grd, method=method)
}

#make_gf <- function(dim, seed, type=0)
#{
#  #0 = equal
#  #1 = random
#  #2 = barrier
#  set.seed(seed)
#  rr <- raster::raster(nrows=dim, ncols=dim, 
#                       xmn=0, xmx=1, ymn=0, ymx=1)
#  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
#  adj <- raster::adjacent(rr, 1:(dim*dim))-1
#  adj <- t(apply(adj,1,sort))
#  adj <- adj[!duplicated(adj),]
#  if(type==0)
#  {
#    cond <- rep(1,dim*dim)
#  } else if(type==1)
#  {
#    cond <- exp(rnorm(dim*dim,0,0.5))
#  }
#  else if(type==2)
#  {
#    condmat <- matrix(1, dim, dim)
#    for(i in 1:ncol(condmat))
#      if (i %in% (floor(dim/2)):(ceiling(dim/2)+1))
#        condmat[,i] <- 1/10#sample(c(1,1/10),dim,prob=c(0.1,0.9),replace=T)
#    condmat[6,] <- 1
#    cond <- c(condmat)
#  }
#  condmat <- matrix(cond,dim,dim)
#  L <- matrix(0,dim*dim,dim*dim)
#  for(i in 1:nrow(adj)) L[adj[i,1]+1,adj[i,2]+1] <- -0.5*cond[adj[i,1]+1] - 0.5*cond[adj[i,2]+1]
#  L <- L + t(L)
#  AA <- -L
#  DD <- -rowSums(L)
#  diag(L) <- -rowSums(L)
#  Linv <- MASS::ginv(L)
#  Linv2 <- Linv %*% t(Linv)
#  L2 <- diag(nrow(L)) - diag(1/DD)%*%AA
#  L2inv <- MASS::ginv(L2)
#  L2inv2 <- L2inv %*% t(L2inv)
#  library(ggplot2)
#  library(reshape2)
#  cnt = 1
#  out <- c()
#  out2 <- c()
#  for(i in 1:dim)
#  for(j in 1:dim)
#  {
#    posx = ((i-1)*dim + 1):((i-1)*dim + dim)
#    posy = ((j-1)*dim + 1):((j-1)*dim + dim)
#    pntx = ((i-1)*dim + i)
#    pnty = ((j-1)*dim + j)
#    out2 <- rbind(out2, cbind(pnty, pntx))
#    out <- rbind(out, cbind(expand.grid(y=posy, x=posx), v=L2inv2[cnt,]))
#    cnt = cnt + 1
#  }
#  rmat <- matrix(rep(1:dim,dim),dim,dim)
#  cmat <- matrix(rep(1:dim,each=dim),dim,dim)
#  Adj <- diag(diag(L)) - L
#  Adj <- melt(Adj)
#  Adj <- Adj[Adj$value != 0.0,]
#  Adj$rowstart <- rmat[Adj$Var1]
#  Adj$rowend <- rmat[Adj$Var2]
#  Adj$colstart <- cmat[Adj$Var1]
#  Adj$colend <- cmat[Adj$Var2]
#  Adj2 <- Adj
#  # if row1 > row2
#  Adj2$rowstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$rowstart + 0.3, Adj2$rowstart)
#  Adj2$rowstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$rowstart - 0.3, Adj2$rowstart)
#  Adj2$rowend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$rowend - 0.3, Adj2$rowend)
#  Adj2$rowend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$rowend + 0.3, Adj2$rowend)
#  Adj2$colstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$colstart - 0.1, Adj2$colstart)
#  Adj2$colstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$colend + 0.1, Adj2$colstart)
#  Adj2$colend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$colstart - 0.1, Adj2$colend)
#  Adj2$colend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$colend + 0.1, Adj2$colend)
##
#  Adj2$colstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$colstart + 0.3, Adj2$colstart)
#  Adj2$colstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$colstart - 0.3, Adj2$colstart)
#  Adj2$colend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$colend - 0.3, Adj2$colend)
#  Adj2$colend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$colend + 0.3, Adj2$colend)
#  Adj2$rowstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$rowstart - 0.1, Adj2$rowstart)
# Adj2$rowstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$rowend + 0.1, Adj2$rowstart)
#Adj2$rowend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$rowstart - 0.1, Adj2$rowend)
#Adj2$rowend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$rowend + 0.1, Adj2$rowend)
#  p1 <- ggplot(melt(condmat), aes(y=Var1, x=Var2)) + geom_tile(aes(fill=value), color="black") + 
#    scale_x_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) + 
#    scale_y_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) +
#    geom_segment(data=Adj2, aes(x=colstart, xend=colend, y=rowstart, yend=rowend, size=exp(2*value)), arrow=arrow(length=unit(0.016,"npc"))) +
#    scale_size(range=c(0.2,1)) + 
#    xlab("X Coordinate") + ylab("Y Coordinate") +
#    #theme(axis.title=element_blank()) +
#    theme(legend.position="none", axis.text=element_blank())
#  if(type>0)
#    p1 = p1 + scale_fill_gradient(low="gray50", high="white")
#  else 
#    p1 = p1 + scale_fill_gradient(low="white", high="white")
#  p2 <- ggplot(melt(Linv2), aes(y=Var1, x=Var2, fill=value)) + geom_tile() + 
#    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
#    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0)
#  p3 <- ggplot(out, aes(y=y, x=x)) + geom_tile(aes(fill=v)) + 
#    scale_x_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) + 
#    scale_y_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) +
#    geom_hline(yintercept=seq(0.5,dim*dim+0.5,dim)) +
#    geom_vline(xintercept=seq(0.5,dim*dim+0.5,dim)) +
#    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) + 
#    xlab("X Coordinate") + ylab("Y Coordinate") +
#    #theme(axis.title=element_blank()) +
#    theme(legend.position="none", axis.text=element_blank())
#  list(p1,p2,p3)
#  # what I want to do: take row 1, 1:100 and make these 1:10,1:10
#}
#make_gf(8,1,2)
#
#make_gf_aniso <- function(dim, seed, aniso=c(0.1,0.2,1.2,1.2))
#{
#  #0 = equal
#  #1 = random
#  #2 = barrier
#  set.seed(seed)
#  rr <- raster::raster(nrows=dim, ncols=dim, 
#                       xmn=0, xmx=1, ymn=0, ymx=1)
#  targ <- c(0, floor(dim*dim/2)-1, dim*dim-2)
#  adj <- raster::adjacent(rr, 1:(dim*dim))-1
#  adj <- t(apply(adj,1,sort))
#  adj <- adj[!duplicated(adj),]
#  rmat <- matrix(rep(1:dim,dim),dim,dim)
#  cmat <- matrix(rep(1:dim,each=dim),dim,dim)
#  L <- matrix(0,dim*dim,dim*dim)
#  for(i in 1:nrow(adj)) 
#    {
#      a1 = adj[i,1]+1
#      a2 = adj[i,2]+1
#      if(rmat[a1] == rmat[a2]) #same row
#      {#greater col
#        if(cmat[a1] > cmat[a2])
#        {
#          L[a1,a2] <- -aniso[1]
#          L[a2,a1] <- -aniso[2]
#        }
#        else
#        {
#          L[a1,a2] <- -aniso[2]
#          L[a2,a1] <- -aniso[1]
#        }
#      }
#      else
#      {#same col
#        if(rmat[a1] > rmat[a2])#greater row
#        {
#          L[a1,a2] <- -aniso[3]
#          L[a2,a1] <- -aniso[4]
#        }
#        else
#        {
#          L[a1,a2] <- -aniso[4]
#          L[a2,a1] <- -aniso[3]
#        }
#      }
#    }
#  #L <- L + t(L)
#  AA <- -L
#  DD <- -rowSums(L)
#  diag(L) <- -rowSums(L)
#  Linv <- MASS::ginv(L)
#  Linv2 <- Linv %*% t(Linv)
#  L2 <- diag(nrow(L)) - diag(1/DD)%*%AA
#  L2inv <- MASS::ginv(L2)
#  L2inv2 <- L2inv %*% t(L2inv)
#  library(ggplot2)
#  library(reshape2)
#  cnt = 1
#  out <- c()
#  out2 <- c()
#  for(i in 1:dim)
#  for(j in 1:dim)
#  {
#    posx = ((i-1)*dim + 1):((i-1)*dim + dim)
#    posy = ((j-1)*dim + 1):((j-1)*dim + dim)
#    pntx = ((i-1)*dim + i)
#    pnty = ((j-1)*dim + j)
#    out2 <- rbind(out2, cbind(pnty, pntx))
##    out <- rbind(out, cbind(expand.grid(y=posy, x=posx), v=Linv2[cnt,]))
#    out <- rbind(out, cbind(expand.grid(y=posy, x=posx), v=L2inv2[cnt,]))
#    cnt = cnt + 1
#  }
#  Adj <- diag(diag(L)) - L
#  Adj <- melt(Adj)
#  Adj <- Adj[Adj$value != 0.0,]
#  Adj$rowstart <- rmat[Adj$Var1]
#  Adj$rowend <- rmat[Adj$Var2]
#  Adj$colstart <- cmat[Adj$Var1]
#  Adj$colend <- cmat[Adj$Var2]
#  Adj2 <- Adj
#  # if row1 > row2
#  Adj2$rowstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$rowstart + 0.3, Adj2$rowstart)
#  Adj2$rowstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$rowstart - 0.3, Adj2$rowstart)
#  Adj2$rowend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$rowend - 0.3, Adj2$rowend)
#  Adj2$rowend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$rowend + 0.3, Adj2$rowend)
#  Adj2$colstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$colstart - 0.1, Adj2$colstart)
#  Adj2$colstart = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$colend + 0.1, Adj2$colstart)
#  Adj2$colend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart < Adj$rowend), Adj$colstart - 0.1, Adj2$colend)
#  Adj2$colend = ifelse( (Adj$colstart == Adj$colend) & (Adj$rowstart > Adj$rowend), Adj$colend + 0.1, Adj2$colend)
##
#  Adj2$colstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$colstart + 0.3, Adj2$colstart)
#  Adj2$colstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$colstart - 0.3, Adj2$colstart)
#  Adj2$colend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$colend - 0.3, Adj2$colend)
#  Adj2$colend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$colend + 0.3, Adj2$colend)
#  Adj2$rowstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$rowstart - 0.1, Adj2$rowstart)
# Adj2$rowstart = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$rowend + 0.1, Adj2$rowstart)
#Adj2$rowend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart < Adj$colend), Adj$rowstart - 0.1, Adj2$rowend)
#Adj2$rowend = ifelse( (Adj$rowstart == Adj$rowend) & (Adj$colstart > Adj$colend), Adj$rowend + 0.1, Adj2$rowend)
##
##what about barrier effects, say we trim
#  condmat <- matrix(0,dim,dim)
#  p1 <- ggplot(melt(condmat), aes(y=Var1, x=Var2)) + geom_tile(color="black", fill="white") + 
#    scale_x_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) + 
#    scale_y_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) +
#    geom_segment(data=Adj2, aes(x=colstart, xend=colend, y=rowstart, yend=rowend, size=exp(value)), arrow=arrow(length=unit(0.016,"npc"))) +
#    scale_size(range=c(0.2,1.0)) + 
#    scale_fill_gradient(low="gray50", high="white") +
#    xlab("X Coordinate") + ylab("Y Coordinate") +
#    #theme(axis.title=element_blank()) +
#    theme(legend.position="none", axis.text=element_blank())
#  p2 <- ggplot(melt(Linv2), aes(y=Var1, x=Var2, fill=value)) + geom_tile() + 
#    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
#    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0)
#  p3 <- ggplot(out, aes(y=y, x=x)) + geom_tile(aes(fill=v)) + 
#    scale_x_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) + 
#    scale_y_continuous(expand=c(0,0), breaks=seq(0.5,dim*dim,dim)) +
#    geom_hline(yintercept=seq(0.5,dim*dim+0.5,dim)) +
#    geom_vline(xintercept=seq(0.5,dim*dim+0.5,dim)) +
#    scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) + 
#    xlab("X Coordinate") + ylab("Y Coordinate") +
#    #theme(axis.title=element_blank()) +
#    theme(legend.position="none", axis.text=element_blank())
#    #geom_point(data=data.frame(out2),aes(x=pntx,y=pnty),pch=3, size=1)
#  list(p1,p2,p3)
#  # what I want to do: take row 1, 1:100 and make these 1:10,1:10
#}
#make_gf_aniso(8,1,c(0.6,1.2,0.6,0.9))
#make_gf_aniso(11,1,c(1,1,1,1))
#
#p11 <- make_gf(9,1,0)[[1]]
#p21 <- make_gf(9,1,0)[[3]]
#p12 <- make_gf(9,1,2)[[1]]
#p22 <- make_gf(9,1,2)[[3]]
#p13 <- make_gf_aniso(9,1,c(0.6,1.2,0.6,0.9))[[1]]
#p23 <- make_gf_aniso(9,1,c(0.6,1.2,0.6,0.9))[[3]]
#library(gridExtra)
#library(extrafont)
#pdf("Ch5figS1.pdf", height=7, width=4.81, family="CM Roman")
#grid.arrange(arrangeGrob(p11, top="Graph",left="IBD"), arrangeGrob(p21, top="Covariance/Vertex"), arrangeGrob(p12,left="Barrier"), p22, arrangeGrob(p13,left="Anisotropic"), p23, ncol=2)
#dev.off()
#embed_fonts("Ch5figS1.pdf")

make_linear_system <- function(covariates, coords, directions=4, saveStack=TRUE)
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
    bad <- raster::getValues(cr) != connected_component
    spdat[bad,] <- NA
    raster::values(covariates[[i]]) <- spdat[,i]
  }

  # final check that graph is connected after pruning
  cr <- raster::clump(covariates[[1]], directions = directions)
  if (length(unique(raster::getValues(cr))) > 1)
    stop("Pruning failed, disconnected components remain")

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
  if(any(is.na(getValues(cr)[uniq_cells])))
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
  if (length(npar) != 3)
    stop ("Length of nuisance parameters must be 3")
  if (length(rpar) != length(linear_system$covariates))
    stop ("Length of resistance parameters must equal number of covariates in linear system")

  covariance <- linear_system$solver$resistance_covariance_log(rpar)

  random_intercepts_per_snp <- rnorm(nsnp, npar[3], exp(npar[1]))
  true_frequencies <- random_intercepts_per_snp + 
    MASS::mvrnorm(nsnp, rep(0,length(linear_system$coords2demes)), 
                  covariance + diag(rep(exp(npar[2]),length(linear_system$coords2demes))))
  true_frequencies <- t(plogis(true_frequencies))
  sampled_chromosomes <- matrix(nchr, length(linear_system$demes), nsnp) 
  observed_frequencies <- matrix(rbinom(nsnp*length(linear_system$coords2demes), 
                                        size=c(sampled_chromosomes), 
                                        prob=c(true_frequencies)), 
                                 length(linear_system$coords2demes), nsnp) #allele counts

  list(snp_count=observed_frequencies,
       chr_count=sampled_chromosomes)
}

inlassleBinomial <- function(snp, chrom, linear_system, start=rep(0,length(linear_system$covariates)+3), maxit=100, parallel=TRUE)
{
  if (class(snp) != "matrix" | class(chrom) != "matrix" | class(linear_system) != "inlassle_linear_system")
    stop("Inputs are of invalid type")
  if (ncol(snp) < nrow(snp))
    stop("Fewer SNPs than demes")
  if (!(all(dim(snp)==dim(chrom))))
    stop("Dimension mismatch in SNP inputs")
  if (length(start) != length(linear_system$covariates) + 3)
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
  n <- nrow(snp)
  indr <- 1:length(linear_system$covariates)
  indn <- (max(indr)+1):(max(indr)+3)
  obj <- function(par)
  {
    rcov <- array(linear_system$solver$resistance_covariance_log(par[indr]), c(n,n,1))
    ll   <- inlassle:::inlassle_test_Likelihood_cov(chrom, snp, matrix(1,n,1), matrix(1,n,1), rcov, par[indn][1], par[indn][2], par[indn][3], FALSE)#TODO: if Field::linesearch throws warnings and parallel=TRUE this can really f*** up the stack
    rd   <- linear_system$solver$rd_resistance_covariance_log(ll$gradient_distance[,,1]) #backpropagate
    ...grad <<- c(rd, ll$gradient) #save gradient for later use
    return (ll$loglik)
  }
  gra <- function(par)
  {
    if (any(is.na(...grad))) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }
  fit <- optim(start, fn=obj, gr=gra, method="BFGS", hessian=TRUE, control=list(maxit=maxit))

  converged <- fit$convergence == 0
  if (!converged)
    warning("Optimizer did not converge")

  hess <- fit$hessian

  # fitted values for resistance parameters
  rtable <- cbind("Est"  = fit$par[indr],
                  "SE"   = sqrt(diag(solve(hess)))[indr],
                  "Z"    = fit$par[indr]/sqrt(diag(solve(hess)))[indr],
                  "pval" = 2*pnorm(abs(fit$par[indr]/sqrt(diag(solve(hess)))[indr]), 0., .1, lower=FALSE))
  rownames(rtable) <- linear_system$covariates

  # fitted values for nuisance parameters
  ntable <- cbind("Est" = fit$par[indn],
                  "SE"  = sqrt(diag(solve(hess)))[indn])
  rownames(ntable) <- c("log Std.Dev (SNP)", "log Std.Dev (Deme)", "(Intercept)")

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

inlassleMLPE <- function(gdist, linear_system, start=rep(0,length(linear_system$covariates)), maxit = 100)
{
  browser()

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
    rdis <- linear_system$solver$resistance_distance_log(par)
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
    tmp    <- tmp + t(tmp)

    ...grad  <<- -linear_system$solver$rd_resistance_distance_log(grad)
    ...model <<- fit

    return -logLik(fit)
  }
  gra <- function(par)
  {
    if (is.na(...grad)) #if objective has not been called, call it
      obj(par)
    grad <- ...grad
    ...grad <<- NA
    return (grad)
  }
  fit <- optim(start, fn=obj, gr=gra, method="BFGS", hessian=TRUE, control=list(maxit=maxit))

  converged <- fit$converged == 0
  if (!converged)
    warning("Optimizer did not converge")

  hess <- fit$hessian

  # fitted values for resistance parameters
  rtable <- cbind("Est"  = fit$par,
                  "SE"   = sqrt(diag(solve(hess))),
                  "Z"    = fit$par/sqrt(diag(solve(hess))),
                  "pval" = 2*pnorm(abs(fit$par/sqrt(diag(solve(hess)))), 0., .1, lower=FALSE))
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

#TODO:
#inlassleWishart <- function(gdist, nsnp, linear_system)
#{
#  ...grad <- NA
#  obj <- function(par)
#}
