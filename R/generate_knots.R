#' @title Generating the knots used in a km call
#' @description Generates automatically the object used as an argument of the km function of the DiceKriging package
#' @param knots.number Array of size d containing integers. 
#' \code{knots.number[i]} contains the number of knots in dimension i. Values of 1 or 0 are considered as 
#' no knots. 
#' @param d Integer. Number of dimensions. This is also equal to \code{length(knots.number)}
#' @param lower Array of size d containing the lower bounds of the input domain.
#' @param upper Array of size d containing the upper bounds of the input domain.
#' @return A list with d fields. Field i is an array of size \code{knots.number[i]} when there is 2 knots 
#' or more, and is 0 otherwise. 
#' @details The knots are placed at the boundary of the domain, so, if we work in the unit hypercube then 
#' the array \code{result[[i]]} will be \code{c(0,1)} if  \code{knots.number[i] = 2}, 
#' \code{c(0,0.5,1)} if  \code{knots.number[i] = 3}, \code{c(0,1/3,2/3,1)} if  \code{knots.number[i] = 4} 
#' and so on. 
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' 
#' knots <- c(5,3,1,0,2,4)
#' 
#' result <- generate_knots(knots.number=knots,d=6)
#' 
generate_knots <- function(knots.number=NULL,d,lower=NULL,upper=NULL){
  
  if(is.null(lower)) lower <- rep(0,times=d)
  if(is.null(upper)) upper <- rep(1,times=d)
  
  if(is.null(knots.number)) return(NULL)
  
  if(length(knots.number)==1){
    if(knots.number>1){
      knots.number <- rep(knots.number,times=d)
    }else{
      return(NULL)
    }
  }
  
  if(length(knots.number)!=d){
    print("Error in function generate_knots. The size of the vector knots.number needs to be equal to d")
    return(NULL)
  }
  
  knots.number <- pmax(1,knots.number) # 2 knots at least per dimension
  
  thelist <- list()
  for(i in 1:d){
    thelist[[i]] <- seq(from=lower[i],to=upper[i],length=knots.number[i])
  }
  return(thelist)
}