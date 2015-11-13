#' @title Independent simulations in several batches of points 
#' @description Performs conditional or non-conditional simulations in n different batches of p points. 
#' @param object The current kriging model. km object.
#' @param nsim Integer. Number of Gaussian process simulation to perform for each batch of p points.
#' @param nbatch Integer. The number n of batches.
#' @param seed The random seed. For repeatability. 
#' @param newdata Matrix of dimension (n*p) x d containing the n different batches of p points where 
#' the simulations are performed. The number of rows must be a multiple of n.
#' @param cond Boolean. When set to TRUE, the simulations are conditional simulations.
#' @param nugget.sim An optimal additional nugget to the simulations. 
#' @param checkNames Boolean. See the documentation of the \code{predict} function of the DiceKriging package.
#' @param type See the documentation of the \code{predict} function of the DiceKriging package.
#' @param pn.only Boolean. When set to TRUE, only the conditional simulations are returned, and other auxiliary results are not returned.
#' @return A list with the following fields. (i) allsimu: (n*p) x nsim matrix. The first p rows correspond 
#' to the nsim simulations (each simulation being a column) associated to batch 1. Next p rows are 
#' for batch 2 and so on. (ii) allKn.inv: (n*p) x p matrix containing n different p x p matrices. Matrix number i is the inverse 
#' of the p x p non-conditional covariance matrix of the p simulation points associated to the 
#' batch i. (iii) allmn: n*p array with the kriging means of all n*p points stored in newdata.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
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
#' d.inv <- length(inv.index);d.opt <- length(opt.index)
#' lower.inv <- lower[inv.index];upper.inv <- upper[inv.index]
#' lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
#' 
#' p <- n.optpoints <- 40
#' nsimu <- 1000
#' 
#' n <- n.points <- 50 # number of integration points
#' inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
#' opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints,dim=d.opt))*(upper.opt-lower.opt))
#' 
#' allsimupoints <- matrix(c(0) , nrow=(n*p) , ncol = d)
#' 
#' for(i in 1:n){
#'   # deal with integration point i
#'   index.first <- 1+(i-1)*p
#'   index.last <- i*p
#'   
#'   inv.point <- inv.integration.points[i,]
#'   simu_points <- matrix(c(0),nrow=p,ncol=d)
#'   
#'   for(alpha in 1:length(opt.index)) simu_points[,opt.index[alpha]] <- opt.simulation.points[,alpha]
#'   for(alpha in 1:length(inv.index)) simu_points[,inv.index[alpha]] <- as.numeric(inv.point[alpha])
#'   
#'   allsimupoints[index.first:index.last ,] <- simu_points
#' }
#' \dontrun{
#' result <- simulate_multi_km(object=model,nsim=nsimu,nbatch = n, newdata=allsimupoints,
#'                             cond=TRUE,checkNames=FALSE,type="UK",seed=NULL,
#'                             pn.only=FALSE)
#' }
simulate_multi_km <- function(object, nsim = 1, nbatch = 1, seed = NULL, newdata = NULL,
                              cond = FALSE, nugget.sim = 0, checkNames = TRUE, type="UK",
                              pn.only=FALSE) {
  
  
  if (is.null(newdata)) {
    return(0)
  } else {
    
    #nbatch <- length(newdata)
    p <- nrow(newdata)/nbatch
    d <- object@d
    
    if (!identical(ncol(newdata), d)) stop("newdata must have the same numbers of columns than the experimental design")
  }
  
  allsimu <- matrix( c(0), nrow = nbatch * p, ncol = nsim )
  allKn.inv <- matrix( c(0) , nrow = nbatch * p, ncol = p )
  allmn <- rep(0,times=nbatch * p)
  
  if(!is.null(seed)) set.seed(seed)
  
  for(i in 1:nbatch){
    
    index.first <- 1+(i-1)*p
    index.last  <- i*p
    
    mypred <- predict(object=object,newdata=newdata[c(index.first:index.last),],type=type,cov.compute=TRUE,se.compute=FALSE,checkNames=FALSE)
    T.cond <- chol(mypred$cov + diag(nugget.sim, p, p))
    
    myrandoms <- rnorm(p*nsim)
    
    
    white.noise <- matrix(myrandoms, p, nsim)
    y.rand.cond <- crossprod(T.cond, white.noise)
    
    y <- matrix(mypred$mean, p, nsim) + y.rand.cond
    
    allsimu[c(index.first:index.last),] <- y
    if(!pn.only){
      allKn.inv[c(index.first:index.last),] <- chol2inv(T.cond)
      allmn[c(index.first:index.last)] <- mypred$mean
    }
  }
  
  return(list(allsimu=allsimu,allKn.inv=allKn.inv,allmn=allmn))
  
}