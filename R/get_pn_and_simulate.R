#' @title Excursion probability through conditional simulations
#' @description Compute excursion probabilities for integration points using 
#' conditional simulations in a set of (possibly: well-chosen) simulation points.
#' @param inv.integration.points Matrix of dimension n x dinv where n is the number of integration
#' points and dinv the number of controlled parameters.
#' @param opt.simulation.points Matrix with dopt columns where 
#' dopt the number of nuisance parameters. Each row is a of simulation points or 
#' candidate simulation point. 
#' @param nsimu Integer. Number of conditional Gaussian process simulation to perform per integration points.
#' @param model The current kriging model. km object.
#' @param T Target threshold.
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters. 
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain. 
#' @param seed The random seed. For repeatability. 
#' @param lowmemory Boolean. When set to TRUE some precomputations on kriging update are not performed.
#' @param subsample_simupoints Boolean. When set to TRUE, p simulation points will be selected among 
#' opt.simulation.points.
#' @param p Integer. Applies only if \code{subsample_simupoints = TRUE}. In that case, this is the number
#' of simulation points which is subsampled (without replacement) among the candidates opt.simulation.points.
#' @param unscale Boolean. When set to TRUE, scaling is applied to the simulation points.
#' @param pn.only Boolean. When set to TRUE, only the excursion probability is returned.
#' @return A list with the following fields. (i) pn: array of size n containing the excursion probabilities 
#' of each integration points. For each point, these excursion probabilities are based on nsimu 
#' conditional simulations. (ii) allsimu: p x (nsimu*n) matrix containing the conditional simulations.
#' The first nsimu columns are the nsimu simulations in p points for the first integration point, 
#' the next nsimu columns are linked to the next integration point and so on. 
#' (iii) allsimucentered: p x (nsimu*n) matrix containing the later simulations centered by 
#' substracting the kriging mean of the points where the simulations are performed.
#' (iv) allsimupoints: (n*p) x d matrix containing all the points where simulations are performed. 
#' The first p rows correspond to the simulation points linked to the first integration points, and so on.
#' In case \code{subsample_simupoints = TRUE}, these points are 'well-chosen' among the candidates.
#' (v) inv.integration.points: n x d.inv matrix containing the integration points passed with the 
#' argument inv.integration.points.
#' (vi) allprecomp: when the argument \code{lowmemory} is set to \code{FALSE}, the function 
#' \code{precomputeUpdateData} from the KrigInv package is called on the (possibly: 'well chosen') 
#' simulation points. The result is a list with several fields which help to compute quickly 
#' updated kriging means and variances. See, the documentation in the KrigInv 
#' package.
#' (vii) allKn.inv: (n*p) x p matrix containing n different p x p matrices. Matrix number i is the inverse 
#' of the p x p non-conditional covariance matrix of the p simulation points associated to the 
#' integration point i. 
#' (viii) allmn: n*p array with the kriging means of all n*p simulation points stored in field (iv).
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
#' n.optpoints <- 40
#' n.optpoints.candidates <- 200
#' nsimu <- 1000
#' 
#' n.points <- 50 # number of integration points
#' inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
#' opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints.candidates,dim=d.opt))*(upper.opt-lower.opt))
#' 
#' \dontrun{
#' result <- get_pn_and_simulate(inv.integration.points=inv.integration.points,
#'    opt.simulation.points=opt.simulation.points,nsimu=nsimu,model=model,T=T,
#'    opt.index=opt.index,inv.index=inv.index,lower=lower,upper=upper,
#'    lowmemory = TRUE,subsample_simupoints = TRUE,unscale=FALSE,
#'    p = n.optpoints)
#' }
get_pn_and_simulate <- function(inv.integration.points,
                                opt.simulation.points,
                                nsimu,
                                model,T,
                                opt.index,inv.index,
                                lower,upper,seed=NULL,lowmemory=FALSE,
                                subsample_simupoints = FALSE,p=NULL,
                                unscale=FALSE,
                                pn.only=FALSE){
  
  n <- nrow(inv.integration.points)
  d.opt <- length(opt.index)
  d.inv <- length(inv.index)
  lower.opt <- lower[opt.index];
  upper.opt <- upper[opt.index]
  d <- model@d
  n.optpoints <- nrow(opt.simulation.points)
  pn <- rep(0,times=n)
  
  if(!subsample_simupoints) p <- n.optpoints
  if(subsample_simupoints){
    if(is.null(p)){
      print("error in get_pn_and_simulate. Specify how many simulation points you want")
      return(0)
    }
  }
  
  allsimupoints <- matrix(c(0) , nrow=(n*p) , ncol = d)
  topt.simulation.points <- t(opt.simulation.points)
  
  # each integration point is dealt with separately
  for(i in 1:n){
    # deal with integration point i
    index.first <- 1+(i-1)*p
    index.last <- i*p
    
    inv.point <- inv.integration.points[i,]
    
    # need to build simu_points (in the space of dimension d)
    simu_points <- matrix(c(0),nrow=n.optpoints,ncol=d)
    
    for(alpha in 1:length(opt.index)) simu_points[,opt.index[alpha]] <- opt.simulation.points[,alpha]
    for(alpha in 1:length(inv.index)) simu_points[,inv.index[alpha]] <- as.numeric(inv.point[alpha])
    
    # choosing simulation points among candidates (when the option is activated)
    if(subsample_simupoints){
      
      index <- picksimulationpoints_weightedPI(model = model,p = p,fullpoints = simu_points,
                                               opt.simulation.points = opt.simulation.points,
                                               topt.simulation.points = topt.simulation.points,
                                               opt.index=opt.index,
                                               T = T,unscale=unscale)
      

      simu_points <- simu_points[index,]
      
    }
    
    allsimupoints[index.first:index.last ,] <- simu_points
  }
  
  # We are out of the loop now. 
  # The allsimupoints matrix (p*n x d) is built
  
  allprecomp <- NULL
  if(!lowmemory) allprecomp <- precomputeUpdateData(model=model,integration.points=allsimupoints)
  
  big.object <- simulate_multi_km(object=model,nsim=nsimu,nbatch = n, newdata=allsimupoints,
                                  cond=TRUE,checkNames=FALSE,type="UK",seed=seed,
                                  pn.only=pn.only)
  
  allsimu <- big.object$allsimu
  
  if(!pn.only){
    allKn.inv <- big.object$allKn.inv
    allmn <- big.object$allmn
    allsimucentered <- allsimu - allmn
  }else{
    allKn.inv <- allmn <- allsimucentered <- NULL
  }
  
  # re-arrange data differently, for faster computations later on
  allsimucentered_test <- allsimu_test <- matrix(c(0),nrow=p,ncol=(nsimu*n))
  if(!pn.only){
    for(i in 1:n){
      index.first <- 1+(i-1)*nsimu
      index.last <- i*nsimu
      index.first2 <- 1+(i-1)*p
      index.last2 <- i*p
      allsimucentered_test[,index.first:index.last] <- allsimucentered[index.first2:index.last2,]
      allsimu_test[,index.first:index.last] <- allsimu[index.first2:index.last2,]
    }
  }
  
  for(i in 1:n){
    index.first <- 1+(i-1)*p
    index.last <- i*p
    
    pn[i] <- sum(colSums(allsimu[index.first:index.last,]>T)>0)/nsimu
  }
  
  return(list(pn=pn,allsimu=allsimu_test,allsimucentered=allsimucentered_test,
              allsimupoints=allsimupoints,
              inv.integration.points=inv.integration.points,
              allprecomp=allprecomp,
              allKn.inv=allKn.inv,
              allmn=allmn)
  )
  
}