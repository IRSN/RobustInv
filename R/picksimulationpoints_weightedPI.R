#' @title Choosing simulation points to maximize an exceedance probability
#' @description Selects (without replacement) p points among some candidates in order to maximize
#' the exceedance probability over a threshold T of a Gaussian process simulated in these p points. 
#' The algorithm consist in choosing sequentially the points with highest individual exceedance 
#' probability and to use a penalty which depends on the distance to the points already selected.
#' @param model The current kriging model. km object.
#' @param p Number of simulation points to select, without replacement.
#' @param fullpoints A matrix with d columns containing all the candidate points.
#' @param opt.simulation.points A matrix with d.opt columns, where d.opt is the number of nuisance 
#' parameters, containing the projection of fullpoints on the space Dopt of nuisance parameters.
#' @param topt.simulation.points The transpose of opt.simulation.points. Used to improve speed.
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param T Target threshold.
#' @param type String. The type of kriging used here. Recommended value is \code{"UK"}.
#' @param unscale Boolean. Used to apply scaling on the simulation points. This does have an influence 
#' on the computed distances between the points.
#' @return An array of size p containing the indices of the selected points. If we call this output 
#' \code{result}, then the final simulation points are \code{fullpoints[result,]}.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' library(KrigInv)
#' library(randtoolbox)
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
#' 
#' n <- n.points <- 50 # number of integration points
#' n.optpoints.candidates <- 500
#' inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
#' opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints.candidates,dim=d.opt))*(upper.opt-lower.opt))
#' 
#' # build the simulation points associated to the first integration point only
#' fullpoints <- matrix(c(0), nrow= n.optpoints.candidates,ncol=d)
#' fullpoints[,inv.index[1]] <-  inv.integration.points[1,1]
#' fullpoints[,inv.index[2]] <-  inv.integration.points[1,2]
#' fullpoints[,opt.index] <- opt.simulation.points
#' 
#' result <- picksimulationpoints_weightedPI(model = model,p = p,fullpoints = fullpoints,
#'                                          opt.simulation.points = opt.simulation.points,
#'                                          opt.index=opt.index,T = T,unscale=FALSE)
#' 
#' 
#' ############################################################
#' # an example with scaling
#' 
#' library(KrigInv)
#' myfun <- function(x){ return(-1*branin_robinv(x) - 50*sin(min(100,1/x[3]))  ) }
#' d <- 3
#' 
#' set.seed(8)
#' 
#' n0 <- 60
#' T <- -60
#' opt.index <- c(3)
#' inv.index <- c(1,2)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' d.inv <- length(inv.index);d.opt <- length(opt.index)
#' lower.inv <- lower[inv.index];upper.inv <- upper[inv.index]
#' lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- apply(X = design,FUN = myfun,MARGIN = 1) 
#' knots.number <- c(0,0,3)
#' knots <- generate_knots(knots.number = knots.number , d = d)
#' 
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",
#' scaling = TRUE,knots=knots)
#' model@@covariance@@eta
#' 
#' p <- n.optpoints <- 40
#' n <- n.points <- 50 # number of integration points
#' n.optpoints.candidates <- 500
#' inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
#' opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints.candidates,dim=d.opt))*(upper.opt-lower.opt))
#' 
#' # more candidates in regions where the gp moves
#' opt.simulation.points <- unscalingFun(mat = opt.simulation.points,model = model,indices = opt.index,standardize = TRUE,lower=lower,upper=upper)
#' 
#' fullpoints <- matrix(c(0), nrow= n.optpoints.candidates,ncol=d)
#' fullpoints[,inv.index[1]] <-  inv.integration.points[1,1]
#' fullpoints[,inv.index[2]] <-  inv.integration.points[1,2]
#' fullpoints[,opt.index] <- opt.simulation.points
#' 
#' result <- picksimulationpoints_weightedPI(model = model,p = p,fullpoints = fullpoints,opt.index = opt.index,
#'                                           opt.simulation.points = opt.simulation.points,
#'                                           T = T,unscale=TRUE)
#' 
#' plot(opt.simulation.points[result,])
picksimulationpoints_weightedPI <- function(model,p,fullpoints,opt.simulation.points,
                                            topt.simulation.points=NULL,opt.index,T,type="UK",
                                            unscale=FALSE){
  
  pred1 <- predict(object = model,newdata = fullpoints,type = type,se.compute = TRUE,cov.compute = FALSE,checkNames = FALSE,light.return = TRUE)
  if(is.null(topt.simulation.points)) topt.simulation.points <- t(opt.simulation.points)
  
  # sequential sampling with a probability of exceedance (PI) penalized by the distance to the points 
  # already chosen
  
  mn <- pred1$mean;sn <- pred1$sd;tmp <- (mn-T)/sn
  PI_all <- pnorm(tmp) # this contains all the probabilities of exceedance
  index <- which.max(PI_all) 
  
  d.opt <- ncol(opt.simulation.points)
  tp1 <- matrix(as.numeric(topt.simulation.points) , nrow=d.opt)
  knots_exist <- (class(model@covariance)=="covAffineScaling" | class(model@covariance)=="covScaling")
  if(knots_exist & unscale){
    # computes the coordinates of tp1 in the rescaled space
    jj=0
    for (j in 1:model@d){
      if (any(opt.index==j)){
        jj = jj+1
        d.opt.j = colnames(model@X)[j]
        tp1[jj,]=scalingFun1d(tp1[jj,], knots=model@covariance@knots[[d.opt.j]],eta=model@covariance@eta[[d.opt.j]])
      }
    }
  }
  
  xmat <- as.numeric(tp1[,index]) # we start with the point with largest PI
  allindex <- c(index)
  
  tp2 <- matrix(tp1-xmat,ncol=d.opt,byrow=TRUE)^2
  mindist <- sqrt(rowSums(tp2))
  
  for(i in 2:p){
    # update distance between potential points and already chosen points
    weighted.PI <- PI_all*mindist
    index <- which.max(weighted.PI)
    xmat.i <- as.numeric(tp1[,index])
    xmat <- c(xmat,xmat.i)
    allindex <- c(allindex,index)
    
    # distance between the j^th point and all other points
    tp2 <- matrix(tp1-xmat.i,ncol=d.opt,byrow=TRUE)^2
    mindist <- pmin(mindist, sqrt(rowSums(tp2)) )
  }
  
  unique_allindex <- unique(allindex)
  if(length(unique_allindex)!=p){
    tmpp <- c(allindex,1:length(mn))
    allindex <- unique(tmpp)[1:p]
  }
  
  return(allindex)
}
