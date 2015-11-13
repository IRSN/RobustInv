#' @title Branin function with additional nuisance parameters
#'
#' @description Similar to the branin function of the DiceKriging package, but with additional
#' mute variables, serving as nuisance parameters.
#' 
#' @param x vector of size 3 or more, or matrix with 3 columns or more.
#' @return a scalar
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' library(DiceKriging)
#' 
#' branin_robinv(c(0.5,0.4,0.6))
#' 
#' xx <- matrix(runif(50),nrow=10)
#' branin_robinv(xx)
branin_robinv <- function(x) {
  if(is.null(dim(x))){
    return(branin(x[1:2]))
  }else{
    return(  apply( X = x[,1:2] , FUN = branin , MARGIN = 1  ) )
  }
}
