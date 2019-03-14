.near <- function(x, y, tol = 1e-5)
  all(abs(x-y) < tol)

.symm <- function(X)
  (X + t(X))/2

.test_Field <- function(seed = 1, n = 10, p = 10, miss = 0)
{
  set.seed(seed)
  Q <- rWishart(1, p, diag(p))[,,1]
  C <- solve(Q)
  Y <- rbinom(p, prob=rbeta(p,1,1), size=n)
  N <- rep(n,p)
  S <- rexp(p,1)
  mu <- rnorm(p, 0, 1)
  miss <- sample(c(TRUE,FALSE), p, prob = c(miss, 1-miss), replace=TRUE)
  N[miss] <- Y[miss] <- 0
  Qm <- solve(C[!miss,!miss,drop=FALSE])
obj <- function(x) 
    0.5*t(x-mu[!miss])%*%Qm%*%(x-mu[!miss]) - 
      sum(-lchoose(N[!miss], Y[!miss]) + dbinom(Y[!miss], N[!miss], plogis(x), log=TRUE))
  obj_lp <- function(mu, C) inlassle_test_Field (Y, N, mu, solve(C))$lapapp

  fit <- optim(rep(0,sum(!miss)), obj, method="L-BFGS-B")
  hess <- numDeriv::hessian (obj, fit$par)

  gs <- list()
  gs[["Q"]] <- Qm
  gs[["mode"]] <- fit$par
  gs[["loglik"]] <- fit$value
  gs[["logdet"]] <- determinant(Qm)$modulus
  gs[["loghes"]] <- determinant(hess)$modulus
  gs[["dlp_dmu"]] <- numDeriv::grad (function(x) obj_lp(x, C), mu)
  gs[["dlp_dC"]] <- .symm(matrix(numDeriv::grad (function(x) obj_lp(mu, x), solve(Q)),p,p))

  out <- inlassle_test_Field(Y, N, mu, Q)
  cs <- list()
  cs[["Q"]] <- out$Q[!miss,!miss]
  cs[["mode"]] <- out$mode[!miss]
  cs[["loglik"]] <- out$loglik
  cs[["logdet"]] <- out$logdet
  cs[["loghes"]] <- out$loghes
  cs[["dlp_dmu"]] <- out$dlp_dmu
  cs[["dlp_dC"]] <- out$dlp_dC

  tst <- sapply(1:length(cs), function(i) .near(gs[[i]], cs[[i]]))
  tst 
}

.test_Likelihood_cov <- function(seed = 1, n = 10, p = 10, l = 100, miss = 0)
{
  set.seed(seed)
  Y <- matrix(rbinom(p*l, prob=rbeta(p*l,1,1), size=n), p, l)
  N <- matrix(n, p, l)
  X <- cbind(rep(1,p),rnorm(p))
  Z <- X
  D <- rWishart(1,n,diag(n))[,,1]
  D <- array(D,c(n,n,1))
  t <- c(0, 0)
  v <- c(0.3, -0.4, 0.2)
  s <- rnorm(2)
  b <- c(-0.2, 0.2)
  miss <- sample(c(TRUE,FALSE), prob = c(miss, 1-miss), p * l, replace=TRUE)
  Y[miss] <- N[miss] <- 0

  obj <- function(t, v, s, b, D)
  {
    t <- exp(t)
    C <- inlassle_test_Weighted_C(D, 2, 0.1, t)
    L <- matrix(0, length(b), length(b))
    L[lower.tri(L, diag=TRUE)] <- v
    P <- diag(exp(diag(L))) #log LDL
    diag(L) <- 1
    L <- L %*% P
    S <- c(exp(2 * Z %*% s))
    Q <- solve(C + X %*% L %*% t(L) %*% t(X) + diag(S))
    mu <- X %*% b
    ll <- sum(sapply(1:ncol(N), function(l)
           inlassle_test_Field (Y[,l], N[,l], mu, Q)$lapapp))
    list(ll=ll/ncol(N), Q=Q)
  }

  obj2 <- function(v, s, b, D)
    inlassle_test_Likelihood_cov (N, Y, X, Z, D, v, s, b, TRUE)[["loglik"]]

  obj3 <- function(i)
    inlassle_test_Likelihood_cov (N[,i,drop=F], Y[,i,drop=F], X, Z, D, v, s, b, TRUE)[["gradient"]]
  grad_by_locus <- sapply(1:ncol(N), obj3)

  gs <- list()
  gs[["Q"]] <- obj (t, v, s, b, D)$Q
  gs[["loglik"]] <- obj (t, v, s, b, D)$ll
  gs[["gradient"]] <- numDeriv::grad(function(x) obj2(
                                                      x[1:length(v)], 
                                                      x[(length(v)+1):(length(v)+length(s))], 
                                                      x[(length(v)+length(s)+1):length(x)], D), c(v,s,b))
  gs[["gradient_distance"]] <- array(.symm(matrix(numDeriv::grad(function(D) obj2(v, s, b, D), D),p,p)), dim=c(p, p, 1))

  cs <- inlassle_test_Likelihood_cov (N, Y, X, Z, D, v, s, b, TRUE)

  tst <- sapply(1:length(gs), function(i) .near(gs[[i]], cs[[i]], tol=1e-5/l))
  tst 
}

.test_Likelihood <- function(seed = 1, n = 10, p = 10, l = 100, miss = 0)
{
  set.seed(seed)
  Y <- matrix(rbinom(p*l, prob=rbeta(p*l,1,1), size=n), p, l)
  N <- matrix(n, p, l)
  X <- cbind(rep(1,p),rnorm(p))
  Z <- X
  D <- array(as.matrix(dist(cbind(runif(p),runif(p)))), c(p,p,1))
  t <- c(0.7, 0.35)
  v <- c(0.3, -0.4, 0.2)
  s <- rnorm(2)
  b <- c(-0.2, 0.2)
  miss <- sample(c(TRUE,FALSE), prob = c(miss, 1-miss), p * l, replace=TRUE)
  Y[miss] <- N[miss] <- 0

  obj <- function(t, v, s, b, D)
  {
    t <- exp(t)
    C <- inlassle_test_Matern_C (D, 2, 0.1, t) 
    L <- matrix(0, length(b), length(b))
    L[lower.tri(L, diag=TRUE)] <- v
    P <- diag(exp(diag(L))) #log LDL
    diag(L) <- 1
    L <- L %*% P
    S <- c(exp(2 * Z %*% s))
    Q <- solve(C + X %*% L %*% t(L) %*% t(X) + diag(S))
    mu <- X %*% b
    ll <- sum(sapply(1:ncol(N), function(l)
           inlassle_test_Field (Y[,l], N[,l], mu, Q)$lapapp))
    list(ll=ll/ncol(N), Q=Q)
  }

  obj2 <- function(t, v, s, b, D)
    inlassle_test_Likelihood (N, Y, X, Z, D, t, v, s, b, TRUE)[["loglik"]]

  obj3 <- function(i)
    inlassle_test_Likelihood (N[,i,drop=F], Y[,i,drop=F], X, Z, D, t, v, s, b, TRUE)[["gradient"]]
  grad_by_locus <- sapply(1:ncol(N), obj3)

  gs <- list()
  gs[["Q"]] <- obj (t, v, s, b, D)$Q
  gs[["loglik"]] <- obj (t, v, s, b, D)$ll
  gs[["gradient"]] <- numDeriv::grad(function(x) obj2(x[1:length(t)], 
                                                      x[(length(t)+1):(length(t)+length(v))], 
                                                      x[(length(t)+length(v)+1):(length(t)+length(v)+length(s))], 
                                                      x[(length(t)+length(v)+length(s)+1):length(x)], D), c(t,v,s,b))
  gs[["gradient_distance"]] <- array(.symm(matrix(numDeriv::grad(function(D) obj2(t, v, s, b, D), D),p,p)), dim=c(p, p, 1))
  gs[["hessian"]] <- MASS::ginv(grad_by_locus %*% t(grad_by_locus) / ncol(N))

  cs <- inlassle_test_Likelihood (N, Y, X, Z, D, t, v, s, b, TRUE)

  tst <- sapply(1:length(gs), function(i) .near(gs[[i]], cs[[i]], tol=1e-5/l))
  tst 
}

.test_Matern <- function(seed = 1, n = 10, p = 10, d = 1, nu = 2, delta = 0.1)
{
  set.seed(seed)
  t <- rexp(d + 1, 1)
  D <- array(as.matrix(dist(cbind(runif(p),runif(p)))), c(p,p,d))

  Dk <- function(D, nu, t)
    apply(D, c(1,2), function(x) sum(sqrt(2*nu) * x / t[-c(1)]))

  maternC <- function(Dk, nu, delta, t)
    ifelse(Dk > 0, t[1]^2 * 2^(1-nu)/gamma(nu) * Dk^nu * besselK(Dk, nu), t[1]^2 + delta^2)

  gs <- list()
  gs[["C"]] <- maternC (Dk(D, nu, t), nu, delta, t)
  gs[["dC_dt"]] <- numDeriv::jacobian (function(x) maternC(Dk(D, nu, x), nu, delta, x), t)
  gs[["dC_dD"]] <- array(apply(D, c(1,2), function(x) numDeriv::grad(function(y) maternC(sum(sqrt(2*nu)* y / t[-c(1)]), nu, delta, t), x)), c(p,p,d))
  for(i in 1:d) diag(gs[["dC_dD"]][,,d]) <- 0

  cs <- list()
  cs[["C"]] <- inlassle_test_Matern_C (D, nu, delta, t)
  cs[["dC_dt"]] <- inlassle_test_Matern_dC_dt (D, nu, delta, t)
  cs[["dC_dD"]] <- inlassle_test_Matern_dC_dD (D, nu, delta, t)

  tst <- sapply(1:length(cs), function(i) .near(gs[[i]], cs[[i]]))
  tst 
}

.test_Gaussian <- function(seed = 1, n = 10, p = 10, d = 1, nu = 2, delta = 0.1)
{
  set.seed(seed)
  t <- rexp(d + 1, 1)
  D <- array(as.matrix(dist(cbind(runif(p),runif(p)))), c(p,p,d))

  Dk <- function(D, t)
    apply(D, c(1,2), function(x) sum(x / t[-c(1)]))

  gaussC <- function(Dk, nu, delta, t)
    ifelse(Dk > 0, t[1]^2 * exp(-Dk^nu), t[1]^2 + delta^2)

  gs <- list()
  gs[["C"]] <- gaussC (Dk(D, t), nu, delta, t)
  gs[["dC_dt"]] <- numDeriv::jacobian (function(x) gaussC(Dk(D, x), nu, delta, x), t)
  gs[["dC_dD"]] <- array(apply(D, c(1,2), function(x) numDeriv::grad(function(y) gaussC(sum(y / t[-c(1)]), nu, delta, t), x)), c(p,p,d))
  for(i in 1:d) diag(gs[["dC_dD"]][,,d]) <- 0

  cs <- list()
  cs[["C"]] <- inlassle_test_Gaussian_C (D, nu, delta, t)
  cs[["dC_dt"]] <- inlassle_test_Gaussian_dC_dt (D, nu, delta, t)
  cs[["dC_dD"]] <- inlassle_test_Gaussian_dC_dD (D, nu, delta, t)

  tst <- sapply(1:length(cs), function(i) .near(gs[[i]], cs[[i]]))
  tst 
}


.test_problem <- function(seed = 1, n = 10, p = 10, l = 100, miss = 0, parallel = FALSE, newton = TRUE, tol = 1e-5)
{
  set.seed(seed)
  Y <- matrix(rbinom(p*l, prob=rbeta(p*l,1,1), size=n), p, l)
  N <- matrix(n, p, l)
  X <- cbind(rep(1,p),rnorm(p))
  D <- array(as.matrix(dist(cbind(runif(p),runif(p)))), c(p,p,1))
  t <- c(0.7, 0.35)
  v <- c(0.3, -0.4, 0.2)
  s <- rexp(p,1)
  b <- c(-0.2, 0.2)
  miss <- sample(c(TRUE,FALSE), prob = c(miss, 1-miss), p * l, replace=TRUE)
  Y[miss] <- N[miss] <- 0

  obj <- function(x)
    inlassle_test_Problem_likelihood (N, Y, X, D, x[1:length(t)], 
                                      x[(length(t)+1):length(c(t,v))], 
                                      x[(length(c(t,v))+1):length(c(t,v,s))], 
                                      x[(length(c(t,v,s))+1):length(c(t,v,s,b))], 
                                      parallel)[["loglik"]]
  gra <- function(x)
    inlassle_test_Problem_likelihood (N, Y, X, D, x[1:length(t)], 
                                      x[(length(t)+1):length(c(t,v))], 
                                      x[(length(c(t,v))+1):length(c(t,v,s))], 
                                      x[(length(c(t,v,s))+1):length(c(t,v,s,b))], 
                                      parallel)[["gradient"]]

  objp <- function(x)
    inlassle_test_Problem_plikelihood (N, Y, X, D, x[1:length(t)], 
                                      x[(length(t)+1):length(c(t,v))], 
                                      x[(length(c(t,v))+1):length(c(t,v,s))], 
                                      x[(length(c(t,v,s))+1):length(c(t,v,s,b))], 
                                      parallel)[["loglik"]]
  grap <- function(x)
    inlassle_test_Problem_plikelihood (N, Y, X, D, x[1:length(t)], 
                                      x[(length(t)+1):length(c(t,v))], 
                                      x[(length(c(t,v))+1):length(c(t,v,s))], 
                                      x[(length(c(t,v,s))+1):length(c(t,v,s,b))], 
                                      parallel)[["gradient"]]

  grah <- function(x)
    inlassle_test_Problem_plikelihood (N, Y, X, D, x[1:length(t)], 
                                      x[(length(t)+1):length(c(t,v))], 
                                      x[(length(c(t,v))+1):length(c(t,v,s))], 
                                      x[(length(c(t,v,s))+1):length(c(t,v,s,b))], 
                                      parallel)[["hessian"]]

  gs <- list()
  gs[["optim"]] <- optim(rep(0,length(c(t,v,s,b))), obj, gr = gra, method="L-BFGS-B")
  gs[["optimp"]] <- optim(rep(0,length(c(t,v,s,b))), objp, gr = grap, method="L-BFGS-B")

  cs <- list()
  cs[["optim"]] <- inlassle_test_Problem_optimize (N, Y, X, D, parallel, tol)
  cs[["optimp"]] <- inlassle_test_Problem_penalize (N, Y, X, D, parallel, tol)
  cs[["optimi"]] <- inlassle_test_Problem_priors (N, Y, X, D, parallel, tol)

  list(gs, cs)
}

.test_ImproperWishart <- function(seed = 1, dim = 4)
{
  set.seed(seed)
  S <- rWishart(1, dim, diag(dim))[,,1]
  eta <- rexp(1,1); gamma <- rexp(dim,1)
  L <- KFAS::ldl(S)
  diag(L) <- log(sqrt(diag(L)))
  l <- L[lower.tri(L,diag=T)]

  make_S <- function(l)
  {
    L <- matrix(0, nrow(S), ncol(S))
    L[lower.tri(L,diag=T)] <- l
    D <- diag(exp(diag(L)))
    diag(L) <- 1
    L %*% D %*% D %*% t(L)
  }

  prior <- function(l)
  {
    S <- make_S(l)
    lp <- 0.5*determinant(S)$modulus
    -lp
  }

  gs <- list()
  gs[["logprob"]] <- prior(l)
  gs[["gradient"]] <- numDeriv::grad(function(l) inlassle_test_ImproperWishart(l)[["logprob"]], l)
  gs[["hessian"]] <- numDeriv::hessian(function(l) inlassle_test_ImproperWishart(l)[["logprob"]], l)

  cs <- inlassle_test_ImproperWishart(l)

  tst <- sapply(1:length(gs), function(i) .near(cs[[i]], gs[[i]]))
  tst
}

.test_Gaussian <- function(seed = 1, dim = 4)
{
  set.seed(seed)
  x <- rnorm(dim)
  mu <- rnorm(dim); Lambda <- rWishart(1, dim, diag(dim))[,,1]

  prior <- function(x, mu, Lambda)
  {
    lp <- -0.5 * t(x-mu) %*% Lambda %*% (x-mu)
    -lp
  }

  gs <- list()
  gs[["logprob"]] <- prior(x, mu, Lambda)
  gs[["gradient"]] <- numDeriv::grad(function(x) inlassle_test_Gaussian(x, mu, Lambda)[["logprob"]], x)
  gs[["hessian"]] <- numDeriv::hessian(function(x) inlassle_test_Gaussian(x, mu, Lambda)[["logprob"]], x)

  cs <- inlassle_test_Gaussian(x, mu, Lambda)

  tst <- sapply(1:length(gs), function(i) .near(cs[[i]], gs[[i]]))
  tst
  # TODO test for normalizing constant
}

.test_LKJ <- function(seed = 1, dim = 4)
{
  set.seed(seed)
  S <- rWishart(1, dim, diag(dim))[,,1]
  eta <- rexp(1,1); gamma <- rexp(dim,1)
  L <- KFAS::ldl(S)
  diag(L) <- log(sqrt(diag(L)))
  l <- L[lower.tri(L,diag=T)]

  make_S <- function(l, ret.mat=FALSE)
  {
    L <- matrix(0, nrow(S), ncol(S))
    L[lower.tri(L,diag=T)] <- l
    D <- diag(exp(diag(L)))
    diag(L) <- 1
    out <- L %*% D %*% D %*% t(L)
    if (ret.mat) # just to recover original matrix, to check if works
      return (out)
    else # actual transform
    {
      s <- sqrt(diag(out))
      out <- cov2cor(out)
      diag(out) <- s
      return(out[lower.tri(out,diag=T)])
    }
  }

  prior <- function(l, eta, gamma)
  {
    S <- make_S(l, TRUE)
    jac <- numDeriv::jacobian (make_S, l)
    lp <- determinant(cov2cor(S))$modulus*(eta-1) - sum(log(diag(S) + gamma^2)) + determinant(jac)$modulus
    -lp
  }

  gs <- list()
  gs[["logprob"]] <- prior(l, eta, gamma)
  gs[["gradient"]] <- numDeriv::grad(function(l) inlassle_test_LKJ(l, eta, gamma)[["logprob"]], l)
  gs[["hessian"]] <- numDeriv::hessian(function(l) inlassle_test_LKJ(l, eta, gamma)[["logprob"]], l)

  cs <- inlassle_test_LKJ(l, eta, gamma)

  tst <- sapply(1:length(gs), function(i) .near(cs[[i]], gs[[i]]))
  tst
  # TODO test for normalizing constant
}

.test_LogBetaPrime <- function(seed = 1, dim = 4)
{
  set.seed(seed)
  s <- rnorm(dim)
  mu <- rbeta(dim,1,1); delta <- rexp(dim,1)

  s_to_f <- function(s)
    1/(exp(s) + 1)

  prior <- function(s, mu, delta)
  {
    f <- s_to_f(s)
    jac <- numDeriv::jacobian (s_to_f, s)
    lp <- sum(log(f^(mu*delta - 1) * (1-f)^((1-mu)*delta - 1))) + log(abs(det(jac)))
    -lp
  }

  gs <- list()
  gs[["logprob"]] <- prior(s, mu, delta)
  gs[["gradient"]] <- numDeriv::grad(function(s) inlassle_test_LogBetaPrime(s, mu, delta)[["logprob"]], s)
  gs[["hessian"]] <- numDeriv::hessian(function(s) inlassle_test_LogBetaPrime(s, mu, delta)[["logprob"]], s)

  cs <- inlassle_test_LogBetaPrime(s, mu, delta)

  tst <- sapply(1:length(gs), function(i) .near(cs[[i]], gs[[i]]))
  tst
  # TODO test for normalizing constant
}

.manual_hess <- function(state, tol = sqrt(.Machine$double.eps))
{
  tind <- 1:(dim(state$D)[3]+1)
  vind <- (tind[length(tind)]+1):(tind[length(tind)]+ncol(state$X)*(ncol(state$X) + 1)/2)
  sind <- (vind[length(vind)]+1):(vind[length(vind)]+nrow(state$Y))
  bind <- (sind[length(sind)]+1):(sind[length(sind)]+ncol(state$X))

  state$He <- inlassle_test_Problem_likelihood (state$N, state$Y, state$X, state$D, 
                                        state$V[tind], state$V[vind],
                                        state$V[sind], state$V[bind], FALSE)[["hessian"]]
  state$Gr <- inlassle_test_Problem_likelihood (state$N, state$Y, state$X, state$D, 
                                                state$V[tind], state$V[vind],
                                                state$V[sind], state$V[bind], FALSE)[["gradient"]]
  state$V <- state$V - MASS::ginv(state$He, tol=tol) %*% state$Gr
  return (state)
}

.manual_phess <- function(state, tol = sqrt(.Machine$double.eps))
{
  tind <- 1:(dim(state$D)[3]+1)
  vind <- (tind[length(tind)]+1):(tind[length(tind)]+ncol(state$X)*(ncol(state$X) + 1)/2)
  sind <- (vind[length(vind)]+1):(vind[length(vind)]+nrow(state$Y))
  bind <- (sind[length(sind)]+1):(sind[length(sind)]+ncol(state$X))

  state$He <- inlassle_test_Problem_plikelihood (state$N, state$Y, state$X, state$D, 
                                        state$V[tind], state$V[vind],
                                        state$V[sind], state$V[bind], FALSE)[["hessian"]]
  state$Gr <- inlassle_test_Problem_plikelihood (state$N, state$Y, state$X, state$D, 
                                                state$V[tind], state$V[vind],
                                                state$V[sind], state$V[bind], FALSE)[["gradient"]]
  state$V <- state$V - MASS::ginv(state$He, tol = tol) %*% state$Gr
  return (state)
}

#---------------------------------- ResistanceSolver

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
  Rd2 <- hi$resistance_distances(c(1,1))
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
  Rd2 <- hi$resistance_distances(pars)
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
  Rd2 <- hi$resistance_distances(pars)
  dRd2 <- hi$rd_resistance_distances(dif)
  emp <- as.vector(dif) %*% numDeriv::jacobian(hi$resistance_distances, pars)
  list(Rd2, dRd2, emp)
}

#---------------------------------- ResistanceOptim

.test_ResistanceOptim <- function(seed=1, dim=20, p=10, l=100, n=10, miss=0, parallel=FALSE)
{
  # TODO: check dimensions of target
  set.seed(seed)

  # spatial data
  rr <- raster::raster(nrows=dim, ncols=dim, 
                       xmn=0, xmx=1, ymn=0, ymx=1)
  targ <- sample(0:(dim*dim-2), p, replace=FALSE)
  adj <- raster::adjacent(rr, 1:(dim*dim))-1
  spd <- cbind(rnorm(dim*dim), rnorm(dim*dim))

  # resistance distance
  rs <- ResistanceSolver$new(spd, targ, adj, parallel)
  pars <- c(0.5,0.5)
  D <- array(rs$resistance_distances_logit(pars),c(p,p,1))
  C <- cov2cor(inlassle_test_Matern_C (D, 2, 0.1, c(0.5, 1.5)))

  # genetic data
  x <- plogis(c(t(MASS::mvrnorm(l, rep(0,p), C))))
  Y <- matrix(rbinom(p*l, prob=rbeta(p*l,2*x,2*(1-x)), size=n), p, l)
  N <- matrix(n, p, l)
  X <- cbind(rep(1,p),rnorm(p))
  D <- array(as.matrix(dist(cbind(runif(p),runif(p)))), c(p,p,1))
  miss <- sample(c(TRUE,FALSE), prob = c(miss, 1-miss), p * l, replace=TRUE)
  Y[miss] <- N[miss] <- 0

  #
  obj <- function(par) inlassle_test_ResistanceOptim_likelihood (par, N,Y,X,D, spd,targ,adj)$loglik * l
  mod <- function(par) inlassle_test_ResistanceOptim_likelihood (par, N,Y,X,D, spd,targ,adj)$start
  gra <- function(par) inlassle_test_ResistanceOptim_likelihood (par, N,Y,X,D, spd,targ,adj)$gradient * l

  # takeaway: if the spatial variance is very small, the gradient for resistance parameters will vanish.
  # so maybe a good idea to use penalization that forces non-zero parameters

  list(numDeriv::grad(obj, c(0,0)), numDeriv::jacobian(mod, c(0,0)), gra(c(0,0)))
}
