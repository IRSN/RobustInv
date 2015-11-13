#' @title Constructing a set of integration and simulation points in robust inversion
#'
#' @description Constructing a set of integration and simulation points in robust inversion
#'   
#' @param integcontrol An important list with the following possible fields. 
#' \code{distrib}: This field needs to be set to either \code{"sur"} or \code{"surnew"} 
#' depending on the sampling criterion that is used. 
#' 
#' (I) If \code{integcontrol$distrib = "surnew"} (default) the 
#' user can set the field \code{n.points} (default: 20) which is the number of integration points. The points 
#' live in the space of controlled parameters, called Dinv here. 
#' The distribution of these n.points points 
#' is controlled through the \code{finaldistrib} field. Possible values are 
#' \code{"spec"} to set up manually the integration points through the \code{inv.integration.points} 
#' field, \code{"MC"} for 
#' uniformly distributed points over Dinv, \code{"sobol"} (default) to use the sobol sequence and \code{"surnew"} 
#' to use importance sampling based on a density proportional to pn (1-pn) with pn being the 
#' excursion probability of the integration points. In that last case, extra fields are required: 
#' \code{n.candidates} is the number of points which are candidate to become integration points (default: n.points*3). 
#' The important sampling scheme will choose (with replacement) n.points integration points 
#' among the candidates. The distribution of the candidates is set with the field 
#' \code{init.distrib} which can be either \code{"MC"} for uniform points, 
#' \code{"sobol"} (default) to use the sobol 
#' sequence or \code{"spec"} to set up manually the matrix of candidates 
#' with the \code{inv.integration.points} field. 
#' When scaling is used in the kriging model, the candidate points can be unscaled by setting the 
#' \code{unscale.inv.integration.points} 
#' to \code{TRUE} (default: \code{FALSE}).
#' Regarding the nuisance parameters, the number of simulation points used, e.g., to compute 
#' pn is set with the field \code{n.optpoints} (default: 100). There is the option to choose these points 
#' among candidates by setting the field \code{choose_optpoints} to \code{TRUE} 
#' (default: \code{FALSE}).
#' In that case, the number of candidates is set with the field \code{n.optpoints.candidates} 
#' (default: 10*n.optpoints). 
#' The field \code{opt.simulation.points.distrib} refers to the distribution of the 
#' simulation points if \code{integcontrol$choose_optpoints=FALSE} and to the 
#' candidates otherwise. Possible values are \code{"spec"} for points set up manually with the
#' \code{opt.simulation.points} field and \code{"sobol"} (default) for the sobol sequence. 
#' When the sobol sequence is used, there is the option to add the vertices of the input 
#' domain (of the nuisance parameters) with the field \code{addvertices} (default: \code{TRUE}). 
#' When scaling is used in the kriging model, the simulation points (or the candidates) 
#' can be unscaled by setting the \code{unscale.opt.simulation.points} field 
#' to \code{TRUE} (default: \code{TRUE}). 
#' The number of simulation of Gaussian process sample paths per integration points 
#' is controlled with the field 
#' \code{nsimu} (default: 1000).
#' Finally, when multiple calls are performed, the \code{seed} argument can be overwritten 
#' by setting the \code{seed} field.
#' 
#' (II) The option \code{integcontrol$distrib = "sur"} needs to be used only when the 
#' sur criterion (and not surnew) is used. In that case, the fields \code{n.points} and 
#' \code{finaldistrib} work similarly except that the importance sampling option is set 
#' with \code{integcontrol$finaldistrib="sur"} and not \code{"surnew"}. In that case, 
#' the field \code{init.distrib} works as in (I). Since there is no Gaussian process simulation here, 
#' the excursion probability computation is set with the field \code{pnmethod}. 
#' Possible values are \code{"fixed"} and \code{"optimMC"} (recommended). The number 
#' of MC trials to get pn (through a multivariate normal cdf call) is set with the 
#' field \code{pop.size}. Finally, the discretization in the space of nuisance parameters 
#' is also obtained with the field \code{n.optpoints} (maximum value: 10).
#' @param d Dimension (i.e. number of scalar input variables) of the objective function.
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain. 
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters. 
#' We must have \code{d = length(opt.index) + length(inv.index)}.
#' @param model The current kriging model. km object.
#' @param T  Target threshold.
#' @param min.prob Minimum probability for a candidate integration point to be picked.
#' @param seed The random seed. For repeatability. 
#' @return Depeding on the fields in the \code{integcontrol} list, the return is different. 
#' See the description on the integcontrol argument and the details.
#' @details (I) If the \code{"surnew"} criterion is used, then the output is a list with 
#' the following fields:
#' \cr
#' (i) inv.integration.points: n x dinv matrix containing the integration points. 
#' n is the number of integration points and dinv the dimension of the space of controlled parameters.
#' \cr
#' (ii) integration.weights: array of size n containing the weights given to each integration points. 
#' \cr
#' (iii) pn: array of size n containing the current excursion probabilities for each integration points. These
#' excursion probabilities are based on nsimu conditional Gaussian process simulations. 
#' \cr
#' (iv) allsimu: p x (nsimu*n) matrix where p is the number of simulation points. 
#' The first nsimu columns are the nsimu simulations in p points for the first integration point, 
#' the next nsimu columns are linked to the next integration point and so on.
#' \cr
#' (v) allsimucentered: p x (nsimu*n) matrix containing the later simulations centered by 
#' substracting the kriging mean of the points where the simulations are performed.
#' \cr
#' (vi) allsimupoints: (n*p) x d matrix containing all the points where simulations are performed. 
#' The first p rows correspond to the simulation points linked to the first integration points, and so on.
#' \cr
#' (vii) allprecomp: a list with 3 fields obtained from a call to the function 
#' \code{precomputeUpdateData} on the points \code{allsimupoints}. Help about \code{precomputeUpdateData} 
#' is given in the KrigInv package.
#' \cr
#' (viii) allKn.inv: (n*p) x p matrix containing n different p x p matrices. Matrix number i is the inverse 
#' of the p x p non-conditional covariance matrix of the p simulation points associated to the 
#' integration point i. 
#' \cr
#' (ix) allmn: n*p array with the kriging means of all n*p simulation points stored in \code{allsimupoints}.
#' \cr
#' (II) If the \code{"sur"} criterion is used, then the output is still a list, which however has different fields 
#' than in (I):
#' \cr
#' (i) integration.points : (p*n) x d matrix containing all the simulation points. The first p points are 
#' associated to the first integration point and so on.
#' \cr
#' (ii) cov: see description of the field allKn.inv (viii) in (I).
#' \cr
#' (iii) mean: see description of the field allmn (ix) in (I).
#' \cr
#' (iv) pn: array of size n containing the current excursion probabilities for each integration points.
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
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- myfun(design)
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
#' 
#' integcontrol <- list(distrib = "surnew",n.points = 20,finaldistrib="surnew",
#'                      n.candidates=50,nsimu=1000,n.optpoints = 50,
#'                      choose_optpoints=TRUE,n.optpoints.candidates=500)
#' \dontrun{
#' obj <- integration_design_robinv(integcontrol = integcontrol,d=d,lower=lower,upper=upper,
#'                                  opt.index=opt.index,inv.index=inv.index,model=model,T=T)                                  
#' }
#' #####################################
#' # An example with scaling
#' library(KrigInv)
#' myfun <- function(x){ return(-1*branin_robinv(x) - 50*sin(min(100,1/x[3]))  ) }
#' d <- 3
#'
#' set.seed(8)
#' 
#' n0 <- 60
#' T <- 30
#' opt.index <- c(3)
#' inv.index <- c(1,2)
#' lower <- rep(0,times=d)
#' upper <- rep(1,times=d)
#' 
#' design <- matrix(runif(d*n0),nrow=n0)
#' response <- apply(X = design,FUN = myfun,MARGIN = 1) 
#' knots.number <- c(0,0,3)
#' knots <- generate_knots(knots.number = knots.number , d = d)
#' 
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",scaling = TRUE,knots=knots)
#' 
#' integcontrol <- list(distrib = "surnew",n.points = 20,finaldistrib="surnew",
#'                      n.candidates=200,nsimu=100,n.optpoints = 50,
#'                      choose_optpoints=TRUE,n.optpoints.candidates=500,
#'                      unscale.opt.simulation.points=TRUE)
#' \dontrun{
#' obj <- integration_design_robinv(integcontrol = integcontrol,d=d,lower=lower,upper=upper,
#'                                  opt.index=opt.index,inv.index=inv.index,model=model,T=T)
#' }                                
integration_design_robinv <- function(integcontrol=NULL,d=NULL,lower,upper,
                                      opt.index,inv.index,
                                      model,T,
                                      min.prob=0.001,seed=NULL){
  
  result<-NULL
  if(!is.null(seed)) set.seed(seed)
  if(is.null(d)) d <- length(lower)
  if (length(lower) != length(upper) ) print("integration_design_robinv: 'lower' and 'upper' must have the same length")
  if(is.null(integcontrol$distrib)) integcontrol$distrib <- "surnew"
  
  if(integcontrol$distrib=="sur"){
    # only for the "sur" criterion.
    n.optpoints <- integcontrol$n.optpoints
    if(is.null(n.optpoints)) n.optpoints <- 10
    if(n.optpoints > 10) n.optpoints <- 10 # maximum value because we cannot go in dimension higher than 20
    
    finaldistrib <- integcontrol$finaldistrib
    if(is.null(finaldistrib)) finaldistrib <- "sobol" # other possible value is "sur"
    
    lower.inv <- lower[inv.index];upper.inv <- upper[inv.index]
    lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
    d.inv <- length(inv.index);d.opt <- length(opt.index)
    n.points <- integcontrol$n.points
    if(is.null(n.points)) n.points <- 20
    
    if(finaldistrib=="MC") inv.integration.points    <- t(lower.inv + t(matrix(runif(d.inv*n.points),ncol=d.inv))*(upper.inv-lower.inv))
    if(finaldistrib=="sobol") inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
    if(finaldistrib=="spec") inv.integration.points <- integcontrol$inv.integration.points
    
    # We now deal with the importance sampling case
    if(finaldistrib=="sur"){
      # construction of the candidates
      n.candidates <- integcontrol$n.candidates
      if(is.null(n.candidates)) n.candidates <- n.points*5
      init.distrib <- integcontrol$init.distrib
      if(is.null(init.distrib)) init.distrib <- "sobol"
      
      if(init.distrib=="sobol") initial.inv.integration.points <- t(lower.inv+t(sobol(n=n.candidates,dim=d.inv))*(upper.inv-lower.inv))
      if(init.distrib=="MC") initial.inv.integration.points <-    t(lower.inv+t(matrix(runif(d.inv*n.candidates),ncol=d.inv))*(upper.inv-lower.inv))
      if(init.distrib=="spec") initial.inv.integration.points <- integcontrol$inv.integration.points
      
      inv.integration.points <- initial.inv.integration.points
    }
    if(d.inv==1) inv.integration.points <- matrix(inv.integration.points,ncol=1)
    
    pnmethod <- integcontrol$pnmethod
    if(is.null(pnmethod)) pnmethod <- "fixed" # other possible value is "optimMC"
    
    if(pnmethod=="fixed"){
      # in Dopt, some points are generated (and fixed) using the sobol sequence (option 1)
      # or uniformly (option 2)
      # This option is computationally cheap; but is expected to give excursion probabilities
      # that might be widely underestimated.
      N <- nrow(inv.integration.points)
      integration.points <- matrix(c(0),nrow=N*n.optpoints,ncol=d)
      if(is.null(integcontrol$pnmethod.fixed.option)) integcontrol$pnmethod.fixed.option <- 1
      
      if(integcontrol$pnmethod.fixed.option==1){
        opt.integration.points <- t(lower.opt + t(sobol(n = n.optpoints, dim = d.opt))*(upper.opt - lower.opt))
        if(d.opt==1) opt.integration.points <- matrix(opt.integration.points,ncol=1)
        
        indices <- expand.grid(c(1:n.optpoints),c(1:N)) 
        index1 <- indices[,1]
        index2 <- indices[,2]
        integration.points[,opt.index] <- as.matrix(opt.integration.points[index1,])
        integration.points[,inv.index] <- as.matrix(inv.integration.points[index2,])
      }else{
        opt.integration.points <- t(lower.opt + t(matrix(runif(n.optpoints*N*d.opt),ncol=d.opt))*(upper.opt - lower.opt))
        if(d.opt==1) opt.integration.points <- matrix(opt.integration.points,ncol=1)
        indices <- expand.grid(c(1:n.optpoints),c(1:N))
        indices[,1] <- c(1:(n.optpoints*N))
        index1 <- indices[,1]
        index2 <- indices[,2]
        integration.points[,opt.index] <- as.matrix(opt.integration.points[index1,])
        integration.points[,inv.index] <- as.matrix(inv.integration.points[index2,])
      }
      
      #calculation of all cov matrices, means and pn
      res <- optim_pn_MC(inv.integration.points=inv.integration.points,model=model,T=T,
                         opt.index=opt.index,inv.index=inv.index,lower=lower,upper=upper,
                         n.optpoints=n.optpoints,pop.size=1,
                         MC=FALSE,deterministicMat=integration.points)
      
    }
    if(pnmethod=="optimMC"){
      pop.size <- integcontrol$pop.size
      if(is.null(pop.size)) pop.size <- 50
      
      res <- optim_pn_MC(inv.integration.points=inv.integration.points,
                         model=model,T=T,
                         opt.index=opt.index,inv.index=inv.index,
                         lower=lower,upper=upper,
                         n.optpoints=n.optpoints,pop.size=pop.size,
                         MC=TRUE)
      
      
    }
    
    if(finaldistrib!="sur") return(res)
    if(finaldistrib=="sur"){
      # some importance sampling
      pn <- res$pn;Tau.n <- pn*(1-pn);Tau.n.sum <- sum(Tau.n)
      if(Tau.n.sum==0) Tau.n.sum <- 1
      prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/n.candidates);prob.n <- prob.n/sum(prob.n)
      weight.n <- 1/(prob.n*n.candidates*n.points)
      prob.n.copy <-c(0,prob.n);prob.n.cum <-cumsum(prob.n.copy)
      
      my.indices <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
      
      pn.my.indices <- pn[my.indices]
      integration.weights <- weight.n[my.indices]
      
      new.mean <- rep(0,times=n.optpoints*n.points)
      old.mean <- res$mean
      new.cov <- matrix(c(0),nrow=n.optpoints*n.points,ncol=n.optpoints)
      old.cov <- res$cov
      new.intpoints <- matrix(c(0),nrow=n.optpoints*n.points,ncol=d)
      old.intpoints <- res$integration.points
      for(i in 1:n.points){
        index.i <- my.indices[i]
        new.mean[(1+(i-1)*n.optpoints):(i*n.optpoints)] <- old.mean[(1+(index.i-1)*n.optpoints):(index.i*n.optpoints)]
        new.cov[(1+(i-1)*n.optpoints):(i*n.optpoints),] <- old.cov[(1+(index.i-1)*n.optpoints):(index.i*n.optpoints),]
        new.intpoints[(1+(i-1)*n.optpoints):(i*n.optpoints),] <- old.intpoints[(1+(index.i-1)*n.optpoints):(index.i*n.optpoints),]
      }
      result <- list(integration.points=new.intpoints,
                     integration.weights=integration.weights,
                     mean=new.mean,
                     cov=new.cov,
                     pn=pn.my.indices)
      return(result)
      
    }
  }
  
  if(integcontrol$distrib=="surnew"){
    ## only for the "surnew" criterion.
    
    n.optpoints <- integcontrol$n.optpoints
    choose_optpoints <- integcontrol$choose_optpoints
    n.optpoints.candidates <- integcontrol$n.optpoints.candidates
    nsimu <- integcontrol$nsimu
    unscale.inv.integration.points <- integcontrol$unscale.inv.integration.points
    
    if(is.null(nsimu)) nsimu <- 1000
    if(is.null(n.optpoints)) n.optpoints <- 100
    if(n.optpoints > 1000) n.optpoints <- 1000
    if(is.null(choose_optpoints)) choose_optpoints <- FALSE
    if(choose_optpoints){if(is.null(n.optpoints.candidates)) n.optpoints.candidates <- min(1000,10*n.optpoints)}
    if(is.null(n.optpoints.candidates)) n.optpoints.candidates <- n.optpoints
    if(is.null(unscale.inv.integration.points)) unscale.inv.integration.points <- FALSE
    
    finaldistrib <- integcontrol$finaldistrib
    if(is.null(finaldistrib)) finaldistrib <- "sobol" #other possible values are "MC" and "surnew"
    
    lower.inv <- lower[inv.index];upper.inv <- upper[inv.index]
    lower.opt <- lower[opt.index];upper.opt <- upper[opt.index]
    d.inv <- length(inv.index);d.opt <- length(opt.index)
    n.points <- integcontrol$n.points
    if(is.null(n.points)) n.points <- 20
    
    if(finaldistrib=="MC") inv.integration.points    <- t(lower.inv + t(matrix(runif(d.inv*n.points),ncol=d.inv))*(upper.inv-lower.inv))
    if(finaldistrib=="sobol") inv.integration.points <- t(lower.inv + t(sobol(n=n.points,dim=d.inv))*(upper.inv-lower.inv))
    if(finaldistrib=="spec") inv.integration.points <- integcontrol$inv.integration.points
    # l'option points d'integration choisis a l'avance est traitee
    # on traite donc maintenant le cas importance sampling
    
    if(finaldistrib=="surnew"){
      # importance sampling will be perfomed later in the function
      # Construction of the points which are 'candidates' to become integration points.
      
      n.candidates <- integcontrol$n.candidates
      if(is.null(n.candidates)) n.candidates <- n.points*3
      init.distrib <- integcontrol$init.distrib
      if(is.null(init.distrib)) init.distrib <- "sobol"
      
      if(init.distrib=="sobol") initial.inv.integration.points <- t(lower.inv+t(sobol(n=n.candidates,dim=d.inv))*(upper.inv-lower.inv))
      if(init.distrib=="MC") initial.inv.integration.points <-    t(lower.inv+t(matrix(runif(d.inv*n.candidates),ncol=d.inv))*(upper.inv-lower.inv))
      if(init.distrib=="spec") initial.inv.integration.points <- integcontrol$inv.integration.points
      inv.integration.points <- initial.inv.integration.points
    }  
    
    if ((class(model@covariance)=="covAffineScaling" | class(model@covariance)=="covScaling") & unscale.inv.integration.points) {
      # unscaling of integration points (or candidates)
      inv.integration.points <- unscalingFun(mat = inv.integration.points,model = model,indices = inv.index,standardize = TRUE,lower = lower,upper = upper)
    }
    
    if(d.inv==1) inv.integration.points <- matrix(inv.integration.points,ncol=1)
    
    # opt.simulation.points now
    
    opt.simulation.points.distrib <- integcontrol$opt.simulation.points.distrib
    if(is.null(opt.simulation.points.distrib)) opt.simulation.points.distrib <- "sobol"
    if(opt.simulation.points.distrib=="spec") opt.simulation.points <- integcontrol$opt.simulation.points
    
    if(opt.simulation.points.distrib=="sobol"){
      addvertices <- integcontrol$addvertices
      if(is.null(addvertices)) addvertices <- TRUE
      if(addvertices)  opt.simulation.points <- rbind(allvertices(d = d.opt,lower = lower.opt,upper = upper.opt) ,  t(lower.opt + t(sobol(n= n.optpoints.candidates-2^d.opt,dim=d.opt))*(upper.opt-lower.opt)) )
      if(!addvertices) opt.simulation.points <-t(lower.opt + t(sobol(n= n.optpoints.candidates,dim=d.opt))*(upper.opt-lower.opt))
    } 
    
    unscale.opt.simulation.points <- integcontrol$unscale.opt.simulation.points
    if(is.null(unscale.opt.simulation.points)) unscale.opt.simulation.points <- TRUE
    
    if ((class(model@covariance)=="covAffineScaling" | class(model@covariance)=="covScaling") & unscale.opt.simulation.points) {
      # unscaling of candidate simulation points
      opt.simulation.points <- unscalingFun(mat = opt.simulation.points,model = model,indices = opt.index,standardize = TRUE,lower = lower,upper = upper)
    }
    if(d.opt==1) opt.simulation.points <- matrix(opt.simulation.points,ncol=1)
    
    
    # Getting pn
    
    if(!is.null(integcontrol$seed)) seed <- integcontrol$seed
    
    res <- get_pn_and_simulate(inv.integration.points=inv.integration.points,
                               opt.simulation.points=opt.simulation.points,
                               nsimu=nsimu,
                               model=model,T=T,
                               opt.index=opt.index,inv.index=inv.index,
                               lower=lower,upper=upper,seed=seed,
                               lowmemory = (finaldistrib=="surnew"),
                               subsample_simupoints = choose_optpoints,
                               unscale=unscale.opt.simulation.points,
                               p = n.optpoints)
    
    if(finaldistrib!="surnew") return(res)
    if(finaldistrib=="surnew"){
      #some importance sampling
      pn <- res$pn;Tau.n <- pn*(1-pn);Tau.n.sum <- sum(Tau.n)
      if(Tau.n.sum==0) Tau.n.sum <- 1
      prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/n.candidates);
      prob.n <- prob.n/sum(prob.n)
      weight.n <- 1/(prob.n*n.candidates*n.points)
      prob.n.copy <-c(0,prob.n);
      prob.n.cum <-cumsum(prob.n.copy)
      
      myrandoms <- runif(integcontrol$n.points)
        
      my.indices <- findInterval(myrandoms,prob.n.cum,all.inside=TRUE)
      
      pn.my.indices <- pn[my.indices]
      integration.weights <- weight.n[my.indices]
      inv.integration.points.my.indices <- inv.integration.points[my.indices,]
      if(d.inv==1) inv.integration.points.my.indices <- matrix(inv.integration.points.my.indices,ncol=1)
      
      
      testobj <- get_pn_and_simulate(inv.integration.points=inv.integration.points.my.indices,
                                     opt.simulation.points=opt.simulation.points,
                                     nsimu=2,
                                     model=model,T=T,
                                     opt.index=opt.index,inv.index=inv.index,
                                     lower=lower,upper=upper,
                                     subsample_simupoints = choose_optpoints,
                                     lowmemory = FALSE,
                                     p = n.optpoints)
      
      allprecomp <- testobj$allprecomp
      allKn.inv <- testobj$allKn.inv
      allsimupoints <- testobj$allsimupoints
      allmn <- testobj$allmn
      
      allsimu <- allsimucentered <- matrix( c(0),  nrow=n.optpoints , ncol=n.points*nsimu)
      
      for(i in 1:n.points){
        index.i <- my.indices[i]
        
        indice.first <- 1+(index.i-1)*nsimu
        indice.last <- index.i*nsimu
        indice.first2 <- 1+(i-1)*nsimu
        indice.last2 <- i*nsimu
        
        allsimu[,indice.first2:indice.last2] <- res$allsimu[,indice.first:indice.last]
        allsimucentered[,indice.first2:indice.last2] <- res$allsimucentered[,indice.first:indice.last]
      }
      result <- list(inv.integration.points=inv.integration.points.my.indices,
                     integration.weights=integration.weights,
                     pn=pn.my.indices,
                     allsimu=allsimu,
                     allsimucentered=allsimucentered,
                     allsimupoints=allsimupoints,
                     allprecomp=allprecomp,
                     allKn.inv=allKn.inv,
                     allmn=allmn
      )
      return(result)
      
    }
  }
  
}