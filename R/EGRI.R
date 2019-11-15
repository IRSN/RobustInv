#' @title Efficient Global Robust Inversion
#' 
#' @description Sequential design of experiments for robust inversion. 
#' 
#' @param T Target threshold.
#' @param model The current kriging model. km object.
#' @param method Criterion used for choosing observations. 
#' Currently, only the \code{"surnew"} criterion is available. 
#' @param fun Objective function.
#' @param iter Number of iterations. At each iteration, a batch of batchsize point is evaluated in parallel.
#' @param batchsize The size of the batch of points evaluated at each iteration.  
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters. 
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain. 
#' @param new.noise.var Optional scalar value of the noise variance of the new observations. 
#' @param integcontrol The integcontrol argument used in the \code{integration_design_robinv} function. 
#' This applies for the following criterion: \code{"surnew"}. See the help of the 
#' \code{integration_design_robinv} function. 
#' @param optimcontrol A list used to set up the optimization procedure of the chosen sampling criterion. 
#' For the \code{"surnew"} criterion, see the help of the \code{max_surnew_robinv} function. 
#' @param kmcontrol Optional list representing the control variables for the re-estimation of the kriging model once new points are sampled. 
#' The items are the same as in the \code{km} function of the DiceKriging package.
#' @param ... Other arguments of the objective function fun.
#' @return A list with the following fields. 
#' (i) par: The added observations (iter*batchsize) x d matrix,
#' (ii) value: The value of the function fun at the added observations, 
#' (iii) lastmodel: The current (last) kriging model of km class, 
#' (iv) lastvalue: The value of the criterion at the last added point,
#' (v) allvalues: If an optimization on a discrete set of points is chosen, the value of the criterion at all these points, for the last iteration,
#' (vi) lastintegration.param: For debug. The last value returned by the integration_design_robinv function, if applicable.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' 
#' library(KrigInv)
#' myfun <- function(x) return(-1 * branin_robinv(x))
#' d <- 4
#' 
#' set.seed(8)
#' 
#' # an example with scaling
#' n0 <- 40
#' T <- -10
#' opt.index <- c(3,4)
#' inv.index <- c(1,2)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' 
#' knots.number <- c(2,2,2,2)
#' knots <- generate_knots(knots.number = knots.number , d = d)
#' 
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",scaling = TRUE,knots = knots)
#' 
#' integcontrol <- list(distrib = "surnew",n.points = 100,finaldistrib="surnew",n.candidates=300,
#'                      nsimu=100,n.optpoints = 50,choose_optpoints=TRUE,n.optpoints.candidates=200)
#' 
#' optimcontrol <- list(method = "genoud", pop.size = 400, max.generations = 4, wait.generation = 1)
#' \dontrun{
#' obj <- EGRI(T = T,method="surnew",model = model,fun = myfun,iter = 2,
#'             batchsize = 2,opt.index = opt.index,inv.index = inv.index,
#'             integcontrol=integcontrol,optimcontrol=optimcontrol)
#' }
EGRI <- function(T, model, method=NULL, fun, iter=1, batchsize=1,
                 opt.index,inv.index,lower=NULL, upper=NULL, 
                 new.noise.var=NULL,integcontrol=NULL,
                 optimcontrol=NULL, kmcontrol=NULL,...) {
    
  n <- nrow(model@X)
  d <- model@d
    
  if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
  if (length(model@penalty==0)) kmcontrol$penalty <- NULL 
  if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
  if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
  if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
  if (is.null(kmcontrol$cov.reestim)) kmcontrol$cov.reestim <- model@param.estim
  if (is.null(kmcontrol$trend.reestim)) kmcontrol$trend.reestim <- model@param.estim
  if (is.null(kmcontrol$nugget.reestim)) kmcontrol$nugget.reestim <- model@param.estim
  if (is.null(optimcontrol$optim.option)) optimcontrol$optim.option <- 2
  if (is.null(method)) method <- "surnew"
  if (is.null(lower)) lower <- rep(0,times=d)
  if (is.null(upper)) upper <- rep(1,times=d)
  integration.param <- NULL
  
  # so far, we force the criterion to be surnew :
  method <- "surnew"
    
  for (i in 1:iter) {
    
    if(method== "surnew"){
       
      if(i==1){
        if(is.null(integcontrol)) integcontrol$distrib <- "surnew"
        if(!is.null(integcontrol)){
          integcontrol$finaldistrib <- integcontrol$distrib
          integcontrol$distrib <- "surnew"
        }
      }
    
      tmp_model <- model
      result <- matrix(c(0),nrow=batchsize,ncol=d)
      
      for(the_iter in 1:batchsize){    
        
        integration.param <- integration_design_robinv(integcontrol=integcontrol,d=d,lower=lower,upper=upper,
                                                         opt.index=opt.index,inv.index=inv.index,
                                                         model=tmp_model,T=T,seed=8)
          
        oEGOI <- max_surnew_robinv(lower=lower, upper=upper, optimcontrol=optimcontrol, opt.index=opt.index,inv.index=inv.index,
                                  integration.param=integration.param,T=T, model=tmp_model, new.noise.var=new.noise.var)
      
        Xnext <- oEGOI$par
        result[the_iter,] <- as.numeric(Xnext)
        # need to set a lie
        if(the_iter < batchsize){
          apred <- predict(object = model,newdata = Xnext,type = "UK",checkNames = FALSE,light.return = TRUE,se.compute = TRUE)
          the_lie <- apred$mean+3*apred$sd
          tmp_model <- update(object = tmp_model,newX = Xnext,newy = the_lie,newX.alreadyExist = FALSE,cov.reestim = FALSE,trend.reestim = FALSE,nugget.reestim = FALSE,newnoise.var = new.noise.var,kmcontrol = kmcontrol)
        }
      }
    }
                
    print("New points"); print(result)
    X.new <- result; y.new <- rep(0,times=nrow(X.new))
    for (i in 1:nrow(X.new)) y.new[i] <- fun(X.new[i,],...)
        
    model <- update(object=model,newX=X.new,newy=y.new,cov.reestim=kmcontrol$cov.reestim,trend.reestim=kmcontrol$trend.reestim,nugget.reestim=kmcontrol$nugget.reestim,  newnoise.var=rep(new.noise.var,times=batchsize),kmcontrol=kmcontrol)
  }
  
  return(list(
        par=model@X[(n+1):(n+iter*batchsize),, drop=FALSE], 
        value=model@y[(n+1):(n+iter*batchsize),, drop=FALSE], 
        lastmodel=model,
        lastvalue=oEGOI$value,
        allvalues=oEGOI$allvalues,
        lastintegration.param=integration.param
        )
      )
}
