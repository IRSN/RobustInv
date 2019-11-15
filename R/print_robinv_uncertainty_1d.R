#' @title Excursion probability and maximum kriging mean plot in 1d
#' @description Computes and can show two plots for functions with only one controlled parameters 
#' and one or several nuisance parameters. One plot is the maximum of the kriging mean 
#' (taken w.r.t. the nuisance parameters) computed on a grid of controlled parameter values. 
#' The second plot is the excursion probability, pn, which is 
#' computed using conditional simulations. 
#' @param model The current kriging model. km object.
#' @param T Target threshold.
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain. 
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters.
#' @param control A list with fields that will control how the different quantities involved are computed. 
#' \code{control$resolution} is the number of points (in the space of the controlled parameters) where the computations 
#' are performed. A grid in dimension 1 of size resolution is used. The field n.optpoints controls the number of points 
#' (in the space of the nuisance parameters) taken to compute the maximum of the kriging mean, or to do conditional 
#' simulations. For the computation of the excursion probability, the fields n.optpoints.candidates, 
#' choose_optpoints, nsimu, unscale.opt.simulation.points are used in a similar way than in the \code{integration_design_robinv} function. See 
#' help on the integcontrol argument in that function. 
#' @param onlycompute Boolean. When FALSE, no plot is performed, but the maximum of the kriging mean (resp. the excursion probability pn) 
#' is still computed if \code{maxmkplot=TRUE} (resp. pnplot = TRUE).
#' @param lwd Line width for the different plots involved in this function.
#' @param cex.main Title size for the different plots involved in this function.
#' @param cex.axis Axis label size for the different plots involved in this function.
#' @param cex.lab Label size for the different plots involved in this function.
#' @param cex.points Point size for the maximum kriging mean plot. Useless if \code{maxmkplot=FALSE}.
#' @param pch.points Point pch for the maximum kriging mean plot. Useless if \code{maxmkplot=FALSE}.
#' @param color.up Color of the points where there is threshold exceedance. Useless if \code{maxmkplot=FALSE}.
#' @param color.down Color of the points where there is no threshold exceedance. Useless if \code{maxmkplot=FALSE}.
#' @param maxmkplot Boolean. When TRUE, the maximum of the kriging mean (taken w.r.t. the nuisance parameters) is 
#' computed. It is also ploted if \code{onlycompute=FALSE}.
#' @param xlab1 x axis label for the maximum of the kriging mean plot.
#' @param ylab1 y axis label for the maximum of the kriging mean plot.
#' @param main1 Title of the maximum of the kriging mean plot.
#' @param pnplot Boolean. When TRUE, the excursion probability function, pn, is computed. 
#' It is also ploted if \code{onlycompute=FALSE}.
#' @param xlab2 x axis label for the excursion probability plot.
#' @param ylab2 y axis label for the excursion probability plot.
#' @param main2 Title of the excursion probability plot.
#' @param newpar Boolean. When TRUE, the par() function is called. Usefull only if the two plots available in this function are both performed.
#' @return A list containing the important computed quantities: 
#' (i) all.points: the value of the controlled parameters, 
#' (ii) maxmk: array of the same size than all.points containing the maximum of the kriging mean taken w.r.t. the nuisance 
#' parameters, 
#' (iii) pn: array of the same size than all.points containing the excursion probability of the considered points (all.points), 
#' (iv) uncertainty: scalar equal to \code{mean(pn*(1-pn))} giving a measure of the current global uncertainty on the excursion set, 
#' (v) colors.transluded: the colors used to plot the points on the maximum of kriging mean plot. 
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' 
#' library(KrigInv)
#' library(randtoolbox)
#' myfun <- function(x) return(-1 * branin_robinv(x))
#' d <- 3
#' 
#' set.seed(8)
#' 
#' n0 <- 30
#' T <- -10
#' opt.index <- c(2,3)
#' inv.index <- c(1)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
#' control <- list(resolution = 100, n.optpoints = 300)
#' \dontrun{
#' print_robinv_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
#'                             opt.index = opt.index,inv.index = inv.index,
#'                             control = control)
#' }
#' ########################################
#' # A more complicated example with scaling
#' library(KrigInv)
#' myfun <- function(x){ return(-1*branin_robinv(x) - 50*sin(min(100,1/x[3]))  ) }
#' d <- 3
#' 
#' set.seed(8)
#' 
#' n0 <- 60
#' T <- 40
#' opt.index <- c(2,3)
#' inv.index <- c(1)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- apply(X = design,FUN = myfun,MARGIN = 1) 
#' knots.number <- c(0,3,3)
#' knots <- generate_knots(knots.number = knots.number , d = d)
#' 
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",scaling = TRUE,knots=knots)
#' # have a look at model@@covariance@@eta
#' 
#' control <- list(resolution = 100, n.optpoints = 300, unscale.opt.simulation.points=TRUE)
#' 
#' \dontrun{
#' print_robinv_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
#'                             opt.index = opt.index,inv.index = inv.index,
#'                             control = control)
#' }
print_robinv_uncertainty_1d <- function(model,T,lower,upper,opt.index,inv.index,control=NULL,onlycompute=FALSE,
                                        lwd=1,cex.main=1,cex.axis=1,cex.lab=1,cex.points=1,pch.points = 17,
                                        color.up = "red",color.down="blue",
                                        maxmkplot=TRUE,xlab1=NULL,ylab1=NULL,main1=NULL,
                                        pnplot=TRUE,xlab2=NULL,ylab2=NULL,main2=NULL,
                                        newpar=TRUE){
  
  if(length(inv.index)!=1) print("Error in print_robinv_uncertainty_1d: the dimension of the inversion domain is not 1.")
  
  d.inv <- 1;d.opt <- length(opt.index)
  lower.xinv <- lower[inv.index]; upper.xinv <- upper[inv.index]
  lower.opt  <- lower[opt.index]; upper.opt  <- upper[opt.index]
  
  the.points <- model@X
  mynames <- colnames(the.points)
  the.resp <- model@y
  d <- model@d;n <- model@n
  
  absissa <- the.points[,inv.index] # cordinates of the evaluated points on Dinv
  resolution <- control$resolution
  if(is.null(resolution)) resolution <- 100
  n.inv.points <- resolution
  all.points <-  matrix(seq(from=lower.xinv,to=upper.xinv,length=resolution),ncol=1)
  inv.pred.points <- all.points
  n.optpoints <- control$n.optpoints;if(is.null(n.optpoints)) n.optpoints <- 250*d.opt
  
  if(pnplot & maxmkplot & newpar) par(mfrow=c(1,2))
  maxmk <- pn <- colors.transluded <- NULL
  
  if(maxmkplot){
    # Plot 1: maximum of kriging mean, with the points.
    
    # builds the color of of each point:
    # Red = exceedance, Blue = no-exceedance
    # Transparency: opaque = close to T, transparent = far from T
    color.tab <- rep(color.down,times=n)
    color.tab[the.resp > T] <- color.up

    # builds the "alpha" array with is used later on to translude the colors.
    if((T > 0) & (min(the.resp)>=0)){
      percentage.to.T <- abs(the.resp - T)/T
      percentage.to.T[percentage.to.T>1] <- 1
      alpha.tab  <- 1-percentage.to.T
    }else{
      therange <- max(the.resp) - min(the.resp)
      if(therange==0) therange <- 1
    
      dist.to.T <- abs(the.resp - T)/therange
      dist.to.T[dist.to.T>1] <- 1
      alpha.tab  <- 1-dist.to.T
    }
    colors.transluded <- translude(color=color.tab,alpha=alpha.tab)
  
    opt.pred.points <- t(lower.opt + t(sobol(n=n.optpoints,dim=d.opt))*(upper.opt-lower.opt))
    if(d.opt==1) opt.pred.points <- matrix(opt.pred.points,ncol=1)
    
    # unscale opt.pred.points:
    unscale.opt.simulation.points <- control$unscale.opt.simulation.points
    if(is.null(unscale.opt.simulation.points)) unscale.opt.simulation.points <- TRUE
    if ((class(model@covariance)=="covAffineScaling" | class(model@covariance)=="covScaling") & unscale.opt.simulation.points) {
      opt.pred.points <- unscalingFun(mat = opt.pred.points,model = model,indices = opt.index,standardize = TRUE,lower = lower,upper=upper)
    }
    
    indices <- expand.grid(c(1:n.optpoints),c(1:resolution)) 
    index1 <- indices[,1] ; index2 <- indices[,2]
    prediction.points <- matrix(c(0),nrow= n.optpoints*resolution, ncol = d)
    prediction.points[,opt.index] <- as.matrix(opt.pred.points[index1,])
    prediction.points[,inv.index] <- as.matrix(inv.pred.points[index2,])
  
    pred <- predict(object=model,newdata=prediction.points,type="UK",se.compute=TRUE,cov.compute=FALSE,light.return=TRUE, checkNames = FALSE)
    mk <- matrix(pred$mean,nrow = n.optpoints)
    u95<- matrix(pred$upper95,nrow = n.optpoints)
    l95<- matrix(pred$lower95,nrow = n.optpoints)
    
    maxmk <- apply(X = mk,FUN = max,MARGIN = 2)
    maxu95 <- apply(X = u95,FUN = max,MARGIN = 2)
    maxl95 <- apply(X = l95,FUN = max,MARGIN = 2)
    ylim <- c(min(maxl95),max(maxu95))
    
      
    if(is.null(xlab1)) xlab1 <- mynames[inv.index]
    if(is.null(ylab1)) ylab1 <- "kriging mean"
    if(is.null(main1)) main1 <- ""
    if(!onlycompute){
      plot( x =   as.numeric(inv.pred.points) , y = maxmk , type="l" , lwd=lwd,xlab=xlab1,ylab=ylab1 , main=main1 , cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis , ylim=ylim)
      lines(x =   as.numeric(inv.pred.points) , y = maxu95 , lwd=lwd , lty = 2)
      lines(x =   as.numeric(inv.pred.points) , y = maxl95 , lwd=lwd , lty = 2)
      abline(h=T,lwd=lwd, lty=2 )
      points(x=the.points[,inv.index] , y = the.resp , col = colors.transluded,pch = pch.points, cex = cex.points)
    }
  }
  
  if(pnplot){
    # a plot of approximated pn
    
    n.optpoints.candidates <- control$n.optpoints.candidates;if(is.null(n.optpoints.candidates)) n.optpoints.candidates <- n.optpoints
    opt.simulation.points <- control$opt.simulation.points
    if(is.null(opt.simulation.points)) opt.simulation.points <- t(lower.opt + t(sobol(n=n.optpoints.candidates,dim=d.opt))*(upper.opt-lower.opt))
    
    if(d.opt==1) opt.simulation.points <- matrix(opt.simulation.points,ncol=1)
    
    # unscale opt.simulation.points:
    unscale.opt.simulation.points <- control$unscale.opt.simulation.points
    if(is.null(unscale.opt.simulation.points)) unscale.opt.simulation.points <- TRUE
    if ((class(model@covariance)=="covAffineScaling" | class(model@covariance)=="covScaling") & unscale.opt.simulation.points) {
      opt.simulation.points <- unscalingFun(mat = opt.simulation.points,model = model,indices = opt.index,standardize = TRUE,lower = lower,upper=upper)
    }
    
    subsample_simupoints <- control$choose_optpoints;if(is.null(subsample_simupoints)) subsample_simupoints <- FALSE
    
    nsimu <- control$nsimu
    if(is.null(nsimu)) nsimu <- 500
      
    max.points <- 500
    n.batch.invpoints <- ceiling(n.inv.points/max.points)
    pn.result <- rep(0,times=n.inv.points)
      
    for(i in 1:n.batch.invpoints){
      index.first <- 1+(i-1)*max.points
      index.last <- min(i*max.points,n.inv.points)
      points.i <- all.points[index.first:index.last,]
      points.i <- matrix(points.i,ncol=1)
        
      res.i <- get_pn_and_simulate(inv.integration.points=points.i , opt.simulation.points=opt.simulation.points,pn.only = TRUE,
                                   nsimu=nsimu,model=model,T=T,opt.index=opt.index,inv.index=inv.index,lower=lower,
                                   upper=upper,lowmemory=TRUE,subsample_simupoints = subsample_simupoints,p= n.optpoints)
        
      pn.result[index.first:index.last] <- res.i$pn
    }
      
    pn <- pn.result
    if(is.null(xlab2)) xlab2 <- mynames[inv.index]
    if(is.null(ylab2)) ylab2 <- "pn - lower bound"
    if(is.null(main1)) main2 <- ""
    
    if(!onlycompute) plot( x =   as.numeric(inv.pred.points) , y = pn , type="l" , lwd=lwd,xlab=xlab2,ylab=ylab2 , main=main2 , cex.main = cex.main, cex.lab = cex.lab, cex.axis = cex.axis)
      
  }
  
  return( list( all.points=all.points, maxmk= maxmk , pn=pn, uncertainty=mean(pn*(1-pn)) , colors.transluded = colors.transluded  ) )
}
