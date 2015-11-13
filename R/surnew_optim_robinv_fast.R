#' @title The surnew criterion in robust inversion
#' @description Compute the surnew criterion at a given point x
#' @param x Array of size d: the point where the criterion is evaluated.
#' @param inv.integration.points n x d.inv matrix containing the integration points
#' @param allsimu p x (nsimu*n) matrix where p is the number of simulation points. 
#' The first nsimu columns are the nsimu simulations in p points for the first integration point, 
#' the next nsimu columns are linked to the next integration point and so on.
#' @param allsimucentered p x (nsimu*n) matrix containing the later simulations centered by 
#' substracting the kriging mean of the points where the simulations are performed.
#' @param allsimupoints (n*p) x d matrix containing all the points where simulations are performed. 
#' The first p rows correspond to the simulation points linked to the first integration points, and so on. 
#' @param allprecomp List with 3 fields obtained from a call to the function 
#' \code{precomputeUpdateData} on the points \code{allsimupoints}. Help about \code{precomputeUpdateData} 
#' is given in the KrigInv package.
#' @param allKn.inv (n*p) x p matrix containing n different p x p matrices. Matrix number i is the inverse 
#' of the p x p non-conditional covariance matrix of the p simulation points associated to the 
#' integration point i.
#' @param integration.weights Array of size n containing the weights given to each integration points.
#' @param model The current kriging model. km object.
#' @param T Target threshold.
#' @param new.noise.var Noise variance of the new observations. 
#' Leave to NULL for noiseless functions. For noisy functions, any non zero value 
#' is valid and will give the same result.
#' @param current.sur The current integral of pn (1-pn) where pn is the excursion probability. 
#' This argument is used to floor the value of the criterion since, if a new point is evaluated, the expected 
#' future uncertainty is supposed to be lower than current.sur.
#' @param randmatrix n x nsimu matrix containing independent realizations of a standard gaussian random variable.
#' @param penalty_visited For points which are too close to the already visited points, the 
#' function outputs \code{current.sur + penalty_visited} (i.e. a 'bad' value for the criterion) in order to force exploration. 
#' 
#' @return A scalar: the value of the surnew criterion
#' @details The arguments inv.integration.points, allsimu, allsimucentered, allsimupoints, 
#' allprecomp, allKn.inv, integration.weights can be generated with a single call to the 
#' \code{integration_design_robinv} function.
#' 
#' The code of this function is meant to be as fast as possible. A large part of 
#' the computation of the criterion is performed through a call to an optimized 
#' C++ function. Suggestions to further improve the speed of these computations 
#' are of course welcome.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' 
#' library(KrigInv)
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
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
#' 
#' n <- 20 ;p <- 50; nsimu <- 1000
#' 
#' integcontrol <- list(distrib = "surnew",n.points = n,finaldistrib="surnew",
#'                      n.candidates=50,nsimu=nsimu,n.optpoints = p,
#'                      choose_optpoints=TRUE,n.optpoints.candidates=500)
#' \dontrun{
#' obj <- integration_design_robinv(integcontrol = integcontrol,d=d,lower=lower,upper=upper,
#'                                  opt.index=opt.index,inv.index=inv.index,model=model,T=T)
#' 
#' randmatrix <- matrix(rnorm(nsimu*n),nrow=n)
#' current.sur <- sum(obj$integration.weights*obj$pn*(1-obj$pn))
#' 
#' x <- c(0.841295 , 0.0757517 , 0.7507468)
#' 
#' result <- surnew_optim_robinv_fast(x = x,allsimu = obj$allsimu,
#'                                    inv.integration.points = obj$inv.integration.points,
#'                                    allsimucentered = obj$allsimucentered,
#'                                    allsimupoints = obj$allsimupoints,
#'                                    allprecomp = obj$allprecomp,
#'                                    allKn.inv = obj$allKn.inv,
#'                                    integration.weights = obj$integration.weights,
#'                                    model = model,T = T,new.noise.var = NULL,
#'                                    current.sur = current.sur,randmatrix = randmatrix)
#' 
#' result
#' current.sur
#' }
surnew_optim_robinv_fast <- function(x,
                                     inv.integration.points,
                                     allsimu,allsimucentered,allsimupoints,
                                     allprecomp,allKn.inv,
                                     integration.weights=NULL,
                                     model, T, new.noise.var=NULL,
                                     current.sur,
                                     randmatrix,penalty_visited=0.01){
  
  
  if(!is.null(new.noise.var)){
    if(new.noise.var == 0) new.noise.var <- NULL
  }
  #x is a vector of size d
  d <- model@d
  n <- model@n
  N <- nrow(inv.integration.points)
  
  X.new <- matrix(x,nrow=d)
  tp1 <- c(as.numeric(t(model@X)),x)
  
  #distance between the point and all other points in the DOE
  xx <- X.new[,1]
  tp2<-matrix(tp1-as.numeric(xx),ncol=d,byrow=TRUE)^2
  mysums <- sqrt(rowSums(tp2))
  mysums[n+1] <- Inf #because this one is always equal to zero...
  mindist <- min(mysums)
  
  if ((mindist > 1e-5) || (!is.null(new.noise.var))){
    
    xnew <- matrix(x,nrow=1)
    pred <- predict_nobias_km(object=model,newdata=xnew,type="UK",checkNames=FALSE)
    
    mn <- pred$mean
    sn <- pred$sd
    F.newdata <- pred$F.newdata
    c.newdata <- pred$c
    
    colnames(xnew) <- colnames(model@X)
    nsimu <- ncol(allsimu)/N
    p <-     ncol(allKn.inv)
    
    kn_xnplus1_simupoints <- computeQuickKrigcov(model=model,integration.points=allsimupoints,
                                                 X.new=xnew,
                                                 precalc.data=allprecomp ,
                                                 F.newdata=F.newdata ,c.newdata=c.newdata )
    
    
    kn_xnplus1_simupoints <- matrix(kn_xnplus1_simupoints,nrow=p,ncol=N)
    lambda <- kn_xnplus1_simupoints/sn^2
    
    big.prob  <- fast_suroptinv(N=N,p=p,nsimu=nsimu,
                                allsimu=allsimu,
                                allsimu_centered=allsimucentered,
                                kn=kn_xnplus1_simupoints,
                                allKninv=allKn.inv,
                                lambda=lambda,
                                mn=mn,sn=sn,T=T,
                                randmatrix=randmatrix)
    
    if (is.null(integration.weights)) {crit <- min(current.sur, mean(big.prob))
    }else crit <- min(current.sur, sum(big.prob*integration.weights))
  }else crit <- current.sur + penalty_visited
  
  return(crit)
}