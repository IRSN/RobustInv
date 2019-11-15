#' @title Unscaling in one dimension
#'
#' @description Unscaling in one dimension using the knots of a km object
#'   
#' @param y Array with values to be unscaled
#' @param knots Array obtained from the field \code{model@@covariance@@knots} of a \code{km} object
#' @param eta Array obtained from the field \code{model@@covariance@@eta} of a km object
#' @param standardize If the initial values y are not in [0,S], 
#' with S the integral of the piecewise linear function equal to eta at points knots, 
#' then there is the possibility to rescale the values y by indicating in which interval they are.
#' @param lower If \code{standardize=TRUE}, this is the lower bound of the interval where the yi's are.
#' @param upper If \code{standardize=TRUE}, this is the upper bound of the interval where the yi's are.
#' @return Array with size \code{length(y)} containing the scaled uncoordinates.
#' @export
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}
#' @examples 
#' library(DiceKriging)
#' knots <- c(1,2,3)
#' eta <- c(2,1,4)
#' t <- seq(from = 1, to = 3, length = 101)
#' scaled_t <- scalingFun1d(x = t,knots = knots,eta = eta)
#' 
#' result <- unscalingFun1d(scaled_t,knots=knots,eta=eta)
#' # now result is equal to t !
#' 
#' # an example of unscaling of uniformly distributed points to have more points in regions where eta is large.
#' myrands <- matrix( runif(2000),ncol=2  )
#' 
#' knots <- c(0,0.5,1)
#' eta <- c(5,1,5) # large on the bounds, low in the middle
#' res1 <- unscalingFun1d(y = myrands[,1] , knots=knots , eta = eta , standardize = TRUE, lower=0, upper = 1)
#' res2 <- unscalingFun1d(y = myrands[,2] , knots=knots , eta = eta , standardize = TRUE, lower=0, upper = 1)
#' 
#' plot(x=res1 , y=res2, type="p") # more points in the corners
unscalingFun1d <- function(y,knots,eta,standardize=FALSE,lower=NULL,upper=NULL){
  
  nknots <- length(knots)
  if(nknots < 2) return(y) # do nothing
  
  # first step: build the vector S
  dknots <- diff(knots)
  avg_height <- (eta[-1] + eta[-nknots])/2
  
  Si <- dknots*avg_height
  S <- c(0,cumsum(Si))
  # S contains the value of scalingFun1d at the different knots.
  
  if(standardize){
    # the values y are in the interval [lower,upper] and not in [ 0, S[nknots] ]
    # so we want to standardize them
    y <- (y - lower) * S[nknots]/(upper - lower)
  }
  
  interval <- findInterval(y,vec = S,all.inside = TRUE)
  
  slope <- (eta[interval+1] - eta[interval])/(knots[interval+1] - knots[interval])
  etai <- eta[interval]
  currentS <- S[interval]
  
  delta <- etai*etai - 2*slope*(currentS-y)
  solution <- (-etai + sqrt(delta) )/slope + knots[interval]
  
  slope0 <- (slope==0)
  solution[slope0] <- ((y - currentS)/etai + knots[interval])[slope0]
  
  return(solution)
}