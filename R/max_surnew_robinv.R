#' @title Optimization of the surnew criterion
#' @description Minimizes the surnew criterion using either discrete or genetic optimization
#' @param lower Array of size d. Lower bound of the input domain.
#' @param upper Array of size d. Upper bound of the input domain.
#' @param optimcontrol A list with the following fields. method : can be either
#' equal to \code{"discrete"} for discrete optimization or \code{"genoud"} (default)
#' for optimization using the \code{genoud} package.
#'
#' (I) If \code{optimcontrol$method="discrete"} the user can set the field
#' optim.points which is a matrix with d columns containing all the points
#' tested to find the optimum. If not set we use 100*d random points generated
#' uniformly.
#'
#' (II) If \code{optimcontrol$method="genoud"} the user can set the arguments of the
#' genoud algorithm. pop.size: population size of each generation (default: 50*d),
#' max.generation: maximum number of generations used to find the optimum (default: 2*d),
#' wait.generation: the algorithm stops if no improvement is done for wait.generation
#' generations (default: 1). Other parameters of the genoud algorithm can be set, namely
#' BFGSburnin, parinit, unif.seed, int.seed, P1, P2, ..., until P9. See the documentation
#' of the genoud function.
#' @param opt.index Array with integers corresponding to the indices of the nuisance parameters.
#' @param inv.index Array with integers corresponding to the indices of the controlled parameters.
#' @param integration.param An object obtained from a call to the
#' \code{integration_design_robinv}.
#' @param T Target threshold.
#' @param model The current kriging model. km object.
#' @param new.noise.var Noise variance of the new observations.
#' Leave to NULL for noiseless functions. For noisy functions, any non zero value
#' is valid and will give the same result.
#' @return A list with the following fields.
#' (i) par: optimizer of the sampling criterion.
#' (ii) value: minimum of the sampling criterion.
#' (iii) allvalues: when a discrete optimization is performed, this contains all the values of the criterion for the tested points.
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
#'
#' # one try with discrete optimization:
#' optimcontrol <- list(method="discrete")
#'
#' result <- max_surnew_robinv(lower = lower,upper = upper,optimcontrol = optimcontrol,
#'                             opt.index = opt.index,inv.index = inv.index,
#'                             integration.param = obj,T = T,model = model)
#'
#' result$par
#' result$value
#'
#' # one try with genoud optimization:
#' optimcontrol <- list(method="pso",pop.size = 200,max.generations=3)
#'
#' result2 <- max_surnew_robinv(lower = lower,upper = upper,optimcontrol = optimcontrol,
#'                              opt.index = opt.index,inv.index = inv.index,
#'                              integration.param = obj,T = T,model = model)
#'
#' result2$par
#' result2$value
#' }
max_surnew_robinv <- function(lower, upper, optimcontrol=NULL,
                              opt.index,inv.index,
                              integration.param,
                              T, model, new.noise.var=NULL){

    if(is.null(optimcontrol$method)) optimcontrol$method <- "pso"

    d <- model@d
    inv.integration.points <- as.matrix(integration.param$inv.integration.points)
    pn <- integration.param$pn
    integration.weights <- integration.param$integration.weights

    if(is.null(integration.weights)) current.sur <- mean(pn*(1-pn))
    if(!is.null(integration.weights)) current.sur <- sum(integration.weights*pn*(1-pn))

    allsimu <- integration.param$allsimu
    allsimucentered <- integration.param$allsimucentered
    allsimupoints <- integration.param$allsimupoints
    allprecomp <- integration.param$allprecomp
    allKn.inv <- integration.param$allKn.inv

    fun.optim <- surnew_optim_robinv_fast

    p <- ncol(allKn.inv)
    N <- nrow(allKn.inv)/p
    nsimu <- ncol(allsimu)/N
    if(is.null(optimcontrol$randmatrix)){
        randmatrix <- matrix(rnorm(nsimu*N),nrow=N)
    }else{
        randmatrix <- optimcontrol$randmatrix
    }

    #discrete Optimisation
    if(optimcontrol$method=="discrete_nopar"){

        if (is.null(optimcontrol$optim.points)){
            n.discrete.points <- d*100
            optimcontrol$optim.points <- t(lower + t(matrix(runif(d*n.discrete.points),ncol=d)) * (upper - lower))
        }
        optim.points <- optimcontrol$optim.points
        optim.points <- data.frame(optim.points)
        colnames(optim.points) <- colnames(model@X)
        all.crit <- seq(1,nrow(optim.points))

        for (i in 1:nrow(optim.points)){
            all.crit[i] <- fun.optim(x=t(optim.points[i,]),
                                     inv.integration.points=inv.integration.points,
                                     integration.weights=integration.weights,
                                     allsimu=allsimu,allsimucentered=allsimucentered,allsimupoints=allsimupoints,
                                     allprecomp=allprecomp,allKn.inv=allKn.inv,
                                     T=T, model=model, new.noise.var=new.noise.var,
                                     current.sur=current.sur,randmatrix=randmatrix,penalty_visited=0)
        }

        ibest <- which.min(all.crit)
        other.points <- as.numeric(optim.points[ibest,])

        o <- list(3)
        o$par <- other.points;o$par <- t(matrix(o$par,nrow=d)); colnames(o$par) <- colnames(model@X)
        o$value <- min(all.crit); o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
        o$allvalues <- all.crit
        return(list(par=o$par, value=o$value,allvalues=o$allvalues))
    }

    #discrete Optimisation
    if(optimcontrol$method=="discrete"){

        if (!getDoParRegistered()) {
            warning("No %dopar% backend registered. You should provide a %dopar% backend (using registerDo*(...) )")
            if (isTRUE(any(rownames(installed.packages())=="doParallel"))) {
                warning("... registering doParallel.")
                library(doParallel)
                registerDoParallel(makePSOCKcluster(rep("localhost",4)),cores=4)
            } else {
                warning("... no doParallel backend avaialble, using sequential instead.")
                registerDoSEQ()
            }
        }
        print(paste0("Using parallel eval of 'surnew_optim_robinv_fast' (",getDoParWorkers()," ",getDoParName()," workers)"))

        if (is.null(optimcontrol$optim.points)){
            n.discrete.points <- d*500
            optimcontrol$optim.points <- t(lower + t(matrix(runif(d*n.discrete.points),ncol=d)) * (upper - lower))
        }
        optim.points <- optimcontrol$optim.points
        optim.points <- data.frame(optim.points)
        colnames(optim.points) <- colnames(model@X)

        all.crit <- unlist(foreach(i=1:nrow(optim.points),.packages = "KrigInv") %dopar% {
            fun.optim(x=t(optim.points[i,]),
                      inv.integration.points=inv.integration.points,
                      integration.weights=integration.weights,
                      allsimu=allsimu,allsimucentered=allsimucentered,allsimupoints=allsimupoints,
                      allprecomp=allprecomp,allKn.inv=allKn.inv,
                      T=T, model=model, new.noise.var=new.noise.var,
                      current.sur=current.sur,randmatrix=randmatrix,penalty_visited=0)
        })

        ibest <- which.min(all.crit)
        other.points <- as.numeric(optim.points[ibest,])

        o <- list(3)
        o$par <- other.points;o$par <- t(matrix(o$par,nrow=d)); colnames(o$par) <- colnames(model@X)
        o$value <- min(all.crit); o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)
        o$allvalues <- all.crit
        return(list(par=o$par, value=o$value,allvalues=o$allvalues))
    }

    #Optimization with Genoud
    if(optimcontrol$method=="genoud"){

        if (is.null(optimcontrol$pop.size))  optimcontrol$pop.size <- 50*d
        if (is.null(optimcontrol$max.generations))  optimcontrol$max.generations <- 2*d
        if (is.null(optimcontrol$wait.generations))  optimcontrol$wait.generations <- 1
        if (is.null(optimcontrol$BFGSburnin)) optimcontrol$BFGSburnin <- 2
        if (is.null(optimcontrol$parinit))  optimcontrol$parinit <- NULL
        if (is.null(optimcontrol$unif.seed))  optimcontrol$unif.seed <- 1
        if (is.null(optimcontrol$int.seed))  optimcontrol$int.seed <- 1

        #mutations
        if (is.null(optimcontrol$P1)) optimcontrol$P1<-0#50
        if (is.null(optimcontrol$P2)) optimcontrol$P2<-0#50
        if (is.null(optimcontrol$P3)) optimcontrol$P3<-0#50
        if (is.null(optimcontrol$P4)) optimcontrol$P4<-0#50
        if (is.null(optimcontrol$P5)) optimcontrol$P5<-50
        if (is.null(optimcontrol$P6)) optimcontrol$P6<-50#50
        if (is.null(optimcontrol$P7)) optimcontrol$P7<-50
        if (is.null(optimcontrol$P8)) optimcontrol$P8<-50
        if (is.null(optimcontrol$P9)) optimcontrol$P9<-0


        #one unique optimisation in dimension d
        domaine <- cbind(lower, upper)

        o <- genoud(fn=fun.optim, nvars=d, max=FALSE, pop.size=optimcontrol$pop.size,
                    max.generations=optimcontrol$max.generations,wait.generations=optimcontrol$wait.generations,
                    hard.generation.limit=TRUE, starting.values=optimcontrol$parinit, MemoryMatrix=TRUE,
                    Domains=domaine, default.domains=10, solution.tolerance=0.000000001,
                    boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
                    data.type.int=FALSE, hessian=FALSE, unif.seed=optimcontrol$unif.seed,
                    int.seed=optimcontrol$int.seed,print.level=0, share.type=0, instance.number=0,
                    output.path="stdout", output.append=FALSE, project.path=NULL,
                    P1=optimcontrol$P1, P2=optimcontrol$P2, P3=optimcontrol$P3,
                    P4=optimcontrol$P4, P5=optimcontrol$P5, P6=optimcontrol$P6,
                    P7=optimcontrol$P7, P8=optimcontrol$P8, P9=optimcontrol$P9,
                    P9mix=NULL, BFGSburnin=optimcontrol$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
                    cluster=FALSE, balance=FALSE, debug=FALSE,
                    model=model, T=T,
                    allsimu=allsimu,allsimucentered=allsimucentered,allsimupoints=allsimupoints,
                    allprecomp=allprecomp,allKn.inv=allKn.inv,
                    inv.integration.points=inv.integration.points,
                    integration.weights=integration.weights,
                    new.noise.var=new.noise.var,current.sur=current.sur,randmatrix=randmatrix)

        o$par <- t(matrix(o$par,nrow=d)); colnames(o$par) <- colnames(model@X)
        o$value <- as.matrix(o$value); colnames(o$value) <- colnames(model@y)

        return(list(par=o$par, value=o$value))
    }

    #Optimization with PSO
    if(optimcontrol$method=="pso"){

        if (is.null(optimcontrol$pop.size))  optimcontrol$pop.size <- 10*d
        if (is.null(optimcontrol$max.generations))  optimcontrol$max.generations <- 10*d
        if (is.null(optimcontrol$wait.generations))  optimcontrol$wait.generations <- 1
        if (is.null(optimcontrol$parinit))  optimcontrol$parinit <- NULL

        if (!getDoParRegistered()) {
            warning("No %dopar% backend registered. You should provide a %dopar% backend (using registerDo*(...) )")
            if (isTRUE(any(rownames(installed.packages())=="doParallel"))) {
                warning("... registering doParallel.")
                library(doParallel)
                registerDoParallel(makePSOCKcluster(rep("localhost",4)),cores=4)
            } else {
                warning("... no doParallel backend avaialble, using sequential instead.")
                registerDoSEQ()
            }
        }
        print(paste0("Using parallel eval of 'surnew_optim_robinv_fast' (",getDoParWorkers()," ",getDoParName()," workers)"))

        vecfun.optim = function(x,
                                inv.integration.points,integration.weights,
                                allsimu,allsimucentered, allsimupoints,
                                allprecomp, allKn.inv,
                                T, model, new.noise.var,
                                current.sur,randmatrix) {

            if (!is.matrix(x)) x = matrix(x,ncol=model@d)
            unlist(foreach(i=1:nrow(x),.packages = "KrigInv") %dopar% {
                fun.optim(x=t(x[i,]),
                          inv.integration.points=inv.integration.points,
                          integration.weights=integration.weights,
                          allsimu=allsimu,allsimucentered=allsimucentered,allsimupoints=allsimupoints,
                          allprecomp=allprecomp,allKn.inv=allKn.inv,
                          T=T, model=model, new.noise.var=new.noise.var,
                          current.sur=current.sur,randmatrix=randmatrix,penalty_visited=0)
            })
        }

        if (is.null(optimcontrol$parinit)) optimcontrol$parinit=matrix((upper+lower)/2,ncol=d)

        o <- psoptim(par = optimcontrol$parinit,
                          fn = vecfun.optim,
                          lower = lower,upper=upper,
                          control = list(maxit=optimcontrol$max.generations,s=optimcontrol$pop.size,vectorize=TRUE),
                          model=model, T=T,
                          allsimu=allsimu,allsimucentered=allsimucentered,allsimupoints=allsimupoints,
                          allprecomp=allprecomp,allKn.inv=allKn.inv,
                          inv.integration.points=inv.integration.points,
                          integration.weights=integration.weights,
                          new.noise.var=new.noise.var,current.sur=current.sur,randmatrix=randmatrix)

        return(list(par=o$par, value=o$value))
    }

}