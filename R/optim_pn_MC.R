#' @title Maximizing excursion probabilities by sampling uniformly points in Dopt
#' @description Computes or maximizes excursion probability of points in Dinv based on either a fixed 
#' matrix of simulation points in D or on a Monte-Carlo optimization procedure.
#' @param inv.integration.points Matrix of dimension n x dinv where n is the number of integration
#' points and dinv the number of controlled parameters.
#' @param model The current kriging model. km object.
#' @param T Target threshold.
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters. 
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain. 
#' @param n.optpoints The number of simulation points in Dopt. Multivariate integrals in dimension 
#' \code{n.optpoints} are performed, so no value larger than 10 should be used.
#' @param pop.size The number of Monte-Carlo samples used to maximize the excursion probability. 
#' Does not apply if \code{MC = FALSE}.
#' @param MC Boolean. Is Monte-Carlo optimization required ?
#' @param deterministicMat Applies only if \code{MC = FALSE}. Matrix of dimension (n*p) x d 
#' containing the simulation points used to compute (only once) pn through a multivariate integral.
#' @return A list with the following fields. 
#' (i) integration.points: (n*p) x d matrix containing the simulation points optimizing the excursion probability. 
#' First p rows are associated to integration point 1, and so on.
#' (ii) cov: (n*p) x p matrix containing the n different non-conditional covariance matrix in each 
#' batch of p points selected by the optimization procedure.
#' (iii) mean: Array of size n*p with the kriging means of the simulation points.
#' (iv) pn: Array of size n with the (underestimated) excursion probabilities for each integration point.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' library(KrigInv)
#' library(mnormt)
#' myfun <- branin_robinv
#' d <- 3
#' 
#' set.seed(8)
#' 
#' n0 <- 30
#' T <- 10
#' opt.index <- c(3)
#' inv.index <- c(1,2)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' d.inv <- length(inv.index);d.opt <- length(opt.index)
#' lower.inv <- lower[inv.index];upper.inv <- upper[inv.index]
#' lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
#' 
#' n <- n.points <- 50 # number of integration points
#' inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
#' 
#' p <- n.optpoints <- 2
#' pop.size <- 20
#' \dontrun{
#' result <- optim_pn_MC(inv.integration.points=inv.integration.points,model=model,T=T,
#'                       opt.index=opt.index,inv.index=inv.index,lower=lower,upper=upper,
#'                       n.optpoints=n.optpoints,pop.size=pop.size,MC=TRUE)
#' }
#' # an example with no MC
#' opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints,dim=d.opt))*(upper.opt-lower.opt))
#' fullpoints <- matrix(c(0),nrow=n*p , ncol=d)
#' indices <- expand.grid(c(1:p),c(1:n)) 
#' index1 <- indices[,1]
#' index2 <- indices[,2]
#' fullpoints[,opt.index] <- opt.simulation.points[index1,]
#' fullpoints[,inv.index] <- inv.integration.points[index2,]
#' \dontrun{
#' result2 <- optim_pn_MC(inv.integration.points=inv.integration.points,model=model,T=T,
#'                        opt.index=opt.index,inv.index=inv.index,lower=lower,upper=upper,
#'                        n.optpoints=n.optpoints,pop.size=1,MC=FALSE,
#'                        deterministicMat = fullpoints)
#' }                        
optim_pn_MC <- function(inv.integration.points,model,T,
                        opt.index,inv.index,lower,upper,
                        n.optpoints,pop.size,MC=TRUE,
                        deterministicMat=NULL){
  
  #un parametre de vitesse important
  #valeur optimale: 100
  
  #the.tests <- c(10,20,30,40,50,80,100,150,200,250)
  #res <- rep(0,times=length(the.tests))
  #for(hh in 1:length(the.tests)){
  #covmatrix.size <- the.tests[hh]
  #t1 <- Sys.time()
  covmatrix.size <- 100
  
  #pour chaque point d'integration on veux le batch de n.optpoints optimal en terme de pn
  N <- nrow(inv.integration.points)
  d.opt <- length(opt.index)
  d.inv <- length(inv.index)
  lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
  d <- model@d
  
  simult.batch <- floor(covmatrix.size/n.optpoints)
  batch.size <- simult.batch*n.optpoints
  total.npoints <- N*pop.size*n.optpoints
  n.batch <- ceiling(N*pop.size / simult.batch)
  
  allmeans <- rep(0,times=N*pop.size*n.optpoints)
  allcov <- matrix(c(0),nrow=N*pop.size*n.optpoints,ncol=n.optpoints)
  
  if(MC){
    rand.mat <- matrix(runif(n=total.npoints*d.opt),nrow=d.opt)
    opt.mat <- t(lower.opt + rand.mat * (upper.opt-lower.opt))
    
    predict.points <- matrix(c(0),nrow=total.npoints,ncol=d)
    
    indices <- expand.grid(c(1:(pop.size*n.optpoints)),c(1:N))
    indices[,1] <- c(1:total.npoints)
    index1 <- indices[,1]
    index2 <- indices[,2]
    
    for(alpha in 1:length(opt.index)) predict.points[,opt.index[alpha]] <- opt.mat[index1,alpha]
    for(alpha in 1:length(inv.index)) predict.points[,inv.index[alpha]] <- inv.integration.points[index2,alpha]
    
  }else{
    predict.points <- deterministicMat
  }
  
  for(i in 1:n.batch){
    index.begin <- 1+(i-1)*batch.size
    index.end <- min(i*batch.size,total.npoints)
    thisbatchsize <- 1 + index.end - index.begin
    
    obj <- predict_nobias_km(object=model,newdata=predict.points[c(index.begin:index.end),],type="UK",cov.compute=TRUE)
    allmeans[index.begin:index.end] <- obj$mean
    tmpcov <- obj$cov
    
    for(j in 1:(thisbatchsize/n.optpoints)){
      ind.start <- 1+(j-1)*n.optpoints
      ind.end <- j*n.optpoints
      cov.start <- index.begin + ind.start - 1
      cov.end <- cov.start + n.optpoints - 1
      allcov[cov.start:cov.end,] <- tmpcov[ind.start:ind.end,ind.start:ind.end]
      
    }
  }
  #les predict sont termines
  
  #il y a maintenant N*pop.size appels a faire a pmnorm
  pn.tab <- rep(0,times=N*pop.size)
  tabT <- rep(T,times=n.optpoints)
  for(i in 1:length(pn.tab)){
    ind.start <- 1 + (i-1)*n.optpoints
    ind.end <- i*n.optpoints
    pn.tab[i] <- 1 - pmnorm(x=tabT,mean=allmeans[ind.start:ind.end],varcov=allcov[ind.start:ind.end,],maxpts=100*n.optpoints)[1]
  }
  pn.mat <- matrix(pn.tab,nrow=pop.size)
  
  best.pn <- rep(0,times=N)
  integration.points <- matrix(c(0),nrow=N*n.optpoints,ncol=d)
  covresult <- matrix(c(0),nrow=N*n.optpoints,ncol=n.optpoints)
  meanresult <- rep(0,times=N*n.optpoints)
  
  for(i in 1:N){
    best.pn[i] <- max(pn.mat[,i])
    argbest <- which.max(pn.mat[,i])
    
    ind.first <- 1 + (i-1)*pop.size*n.optpoints + (argbest - 1)*n.optpoints
    ind.last <- ind.first + n.optpoints - 1
    argbatch <- predict.points[ind.first:ind.last,]
    
    integration.points[(1+(i-1)*n.optpoints):(i*n.optpoints),] <- argbatch
    covresult[(1+(i-1)*n.optpoints):(i*n.optpoints),] <- allcov[ind.first:ind.last,]
    meanresult[(1+(i-1)*n.optpoints):(i*n.optpoints)] <- allmeans[ind.first:ind.last]
    
    #verif
    #obj.verif <- predict_nobias_km(object=model,newdata=argbatch,type="UK",cov.compute=TRUE)
    #1 - pmnorm(x=tabT,mean=obj.verif$mean,varcov=obj.verif$cov,maxpts=100*n.optpoints)[1]
  }
  
  return(list(integration.points=integration.points,cov=covresult,mean=meanresult,pn=best.pn))
  
}