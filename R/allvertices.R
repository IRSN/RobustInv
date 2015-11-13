#' @title Vertices of an hypercube
#'
#' @description Builds the 2^d vertices of an hypercube
#' 
#' @param d dimension of the hypercube
#' @param lower Array of size d containing the lower bound of the domain in each dimension.
#' @param upper Array of size d containing the upper bound of the domain in each dimension.
#' @return A matrix with 2^d rows and d columns containing all the vertices.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' allvertices(d=3,lower=c(0,-10,10),upper=c(100,0,1000))
allvertices <- function(d,lower=NULL,upper=NULL){
  if(is.null(lower)) lower = rep(0,times=d)
  if(is.null(upper)) upper = rep(1,times=d)
  
  e <- c(lower[1],upper[1])
  if(d==1) return(matrix(e,ncol=1))
  if(d >1){
    for(i in 2:d) e <- rbind(cbind(e,lower[i]),cbind(e,upper[i]))
  }
  return(e)
}
