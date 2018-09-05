#sloppy sketch of getting predictive measure from LIPO fits
#Burnham & Anderson 2004 pg 61 KEY RESULT
get_predictive <- function(testset, datapack, outer_par, inner_par, parallel=FALSE)
{
  hi <- ResistanceSolver$new(datapack$spd, datapack$targ, datapack$adj, parallel=parallel)
  D <- array(hi$resistance_distances(outer_par), c(nrow(datapack$Y),nrow(datapack$Y),1))
  field <- run_partial(datapack$Y, datapack$N, datapack$X, datapack$sp, datapack$sm, datapack$spd, datapack$targ, datapack$adj, D, inner_par, parallel)
  field2 <- run_partial(testset$Y, testset$N, datapack$X, datapack$sp, datapack$sm, datapack$spd, datapack$targ, datapack$adj, D, inner_par, parallel)
  out <- bedassle(testset$Y, testset$N, datapack$X, D, field$fixef, matrix(inner_par[3],1,1), field$dispersion, c(1,inner_par[1:2]), parallel)

  list(test_fit=out, training_fit=field, test_fit_alt=field2)
}

get_conductance <- function(testset, datapack, outer_par, inner_par, parallel=FALSE)
{
  hi <- ResistanceSolver$new(datapack$spd, datapack$targ, datapack$adj, parallel=parallel)
  D <- array(hi$resistance_distances(outer_par), c(nrow(datapack$Y),nrow(datapack$Y),1))
  cond <- hi$getConductance()
  field <- run_partial(datapack$Y, datapack$N, datapack$X, datapack$sp, datapack$sm, datapack$spd, datapack$targ, datapack$adj, D, inner_par, parallel)
  blah <- hi$rd_resistance_distances (field$grad_distances[,,1])
  gr_cond <- hi$getGradConductance()

  list(training_fit=field, conductances=cond, grad_conductances=gr_cond)
}

predictive_model <- function(VAR, SVD=TRUE, parallel=FALSE, conductance=FALSE){

  if(SVD)
    runName <- paste0("_SVDmodel_", paste(VAR, collapse="_"), "_")
  else
    runName <- paste0("_SVDDWTmodel_", paste(VAR, collapse="_"), "_")

logname <- paste0(runName, ".log")
outname <- paste0(runName, ".RData")

if(SVD)
  colClass <- rep("NULL", 8)
else
  colClass <- rep("NULL", 24)
colClass[VAR] <- NA

load("r60_datapack.RData")
load("template_small.RData")
load("test_set.RData")

if(SVD)
  spd = scale(as.matrix(read.table("Svd_spatial_data.txt", header=TRUE, colClasses=colClass))[!ocean_small,,drop=FALSE])
else
  spd = scale(as.matrix(read.table("Svd_spatial_data_DWT.txt", header=TRUE, colClasses=colClass))[!ocean_small,,drop=FALSE])

library(inlassle)
datapack <- list(Y = Y, N = N, X = matrix(1,16,1), sp = rep(1,16), sm = rep(1,16),
                                  spd = spd,
                                                   adj = adj_remapped-1, targ = targ_remapped-1)
testset <- list(Y = Ytest, N = Ntest)

load(outname)

if (!conductance)
  out <- get_predictive(testset, datapack, fit$resist[-c(1)], fit$spatial[-c(1)], parallel)
else
  out <- get_conductance(testset, datapack, fit$resist[-c(1)], fit$spatial[-c(1)], parallel)

list(fits=out, optimized=fit)
}

predictive_plot_data <- function(parallel=TRUE)
{
  out <- data.frame()
  for (i in 1:8)
  {
    svd_f <- predictive_model(i, SVD=TRUE, parallel=parallel)
    varname <- paste0("SVD", i)
    out <- rbind(out, data.frame(SVD=varname,scale="Full",test=2*sum(svd_f$fits$test_fit$loglik),training=2*svd_f$fits$training_fit$loglik,test_alt=2*svd_f$fits$test_fit_alt$loglik))

    for (j in 1:7)
    {
      if (j==1) { k=1 } else if (j==2) { k=2} else if (j==3) {k=3} else if (j==4) {k=c(1,2)} else if (j==5) {k=c(1,3)} else if(j==6) {k=c(2,3)} else if(j==7) {k=c(1,2,3)}
      sc <- c("C", "M", "F")
      scalename <- paste(sc[k],collapse="")
      vind <- (i-1)*3 + k
      svd_f <- predictive_model(vind, SVD=FALSE, parallel=parallel)
      out <- rbind(out, data.frame(SVD=varname,scale=scalename,test=2*sum(svd_f$fits$test_fit$loglik),training=2*svd_f$fits$training_fit$loglik,test_alt=2*svd_f$fits$test_fit_alt$loglik))
    }
  }
  out
}

plot_conductance <- function(VAR, SVD=TRUE, parallel=TRUE)
{
  svd_f <- predictive_model(VAR, SVD, conductance=TRUE, parallel=parallel)
  load("template_small.RData")
  cond <- tmplt_small
  grcond <- tmplt_small
  grcond2 <- tmplt_small
  cond[!ocean_small] <- svd_f$fits$conductances;
  grcond[!ocean_small] <- svd_f$fits$grad_conductances;
  grcond2[!ocean_small] <- svd_f$fits$conductances * (1-svd_f$fits$conductances) * svd_f$fits$grad_conductances;
  return(list(cond=cond, grcond=grcond, grcond2=grcond2))
}
