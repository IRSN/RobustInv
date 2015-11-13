#' @title Unscaling function applied to all the dimensions
#'
#' @description Unscaling the coordinates of points in dimension d
#'   
#' @param mat Matrix with values to be unscaled
#' @param model km object containing knots and eta's. If provided, the arguments allknots and alleta are ignored.
#' @param allknots List of arrays obtained from the field \code{model@@covariance@@knots} of a \code{km} object
#' @param alleta List of arrays obtained from the field \code{model@@covariance@@eta} of a km object
#' @param indices Array containing the indices of the dimensions to rescale.
#' @param standardize If the initial values in the columns of mat are not in [0,S], 
#' with S the integral of the piecewise linear function equal to eta at points knots, 
#' then there is the possibility to rescale the values y by indicating in which interval they are.
#' @param lower If \code{standardize=TRUE}, this is an array of lower bounds.
#' @param upper If \code{standardize=TRUE}, this is an array of upper bounds.
#' @return A rescaled matrix with the same size.
#' @details There are two possible ways to use this function. If the number of columns of mat 
#' is not equal to \code{d = model@@d} (i.e. the number of dimensions), we assume that mat 
#' is allready a submatrix of some bigger matrix with d columns. In that case, all the columns 
#' of mat are rescaled an so \code{length(indices)} Must be equal to the number of columns of mat. 
#' The second way to use the function is to have a matrix mat with the same number of columns as d. 
#' In that case only some (or all) columns of mat will be rescaled, depending on the values of indices. 
#' If indices is missing, all the columns are rescaled.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' library(DiceKriging)
#' myfun <- function(x) return(-1 * branin_robinv(x))
#' d <- 4
#' 
#' set.seed(8)
#' 
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
#' knots.number <- c(3,0,2,2)
#' knots <- generate_knots(knots.number = knots.number , d = d)
#' 
#' model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",scaling = TRUE,knots = knots)
#' allknots <- model@@covariance@@knots
#' alleta <- model@@covariance@@eta
#' 
#' # generates a sequence of points in dimension 4:
#' myrands <- matrix( runif(2000),ncol=4  )
#' 
#' result <- unscalingFun(mat = myrands[,1:2], allknots =  allknots,alleta = alleta,indices = c(1,2),standardize = TRUE , lower = lower,upper = upper)
#' 
#' result2 <- unscalingFun(mat = myrands,      allknots =  allknots,alleta = alleta,indices = c(1,2),standardize = TRUE , lower = lower,upper = upper)
#' 
#' # result and result2 are the same
unscalingFun <- function(mat, model = NULL, allknots = NULL, alleta = NULL, indices=NULL , 
                         standardize = FALSE , lower = NULL, upper = NULL){
  
  if(!is.null(model)){
    allknots <- model@covariance@knots
    mynames <- names(allknots)
    alleta <- model@covariance@eta[mynames]
  }
  if(is.null(allknots)) return(mat) # do nothing
  
  d <- length(allknots)
  if(is.null(indices)) indices <- 1:d
  if(is.null(lower)) lower <- rep(0,times=d)
  if(is.null(upper)) upper <- rep(1,times=d)
  
  d2 <- ncol(mat)
  result <- mat
  
  if(d2 < d){
    # only a few columns are passed, in that case 
    # length(indices) must be equal to d2
    if(length(indices) != d2) return(mat) # do nothing
    for(i in 1:d2) result[,i] <- unscalingFun1d(y = mat[,i] , knots = allknots[[indices[i]]] , 
                                                eta =  alleta[[indices[i]]] , standardize = standardize, 
                                                lower=lower[indices[i]],upper=upper[indices[i]])
  }else{
    # case where all the colums are passed
    # so only the ones present in indices will be modified
    j <- 1
    for(i in indices){
      result[,i] <- unscalingFun1d(y = mat[,i] , knots = allknots[[indices[j]]] , 
                                   eta =  alleta[[indices[j]]] , standardize = standardize, 
                                   lower=lower[indices[j]],upper=upper[indices[j]])
      j <- j+1
    }
  }
  
  return(result)
}