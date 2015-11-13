#' @title 3d test function computing a k-effective in nuclear safety
#'
#' @description 3d test function computing a k-effective in nuclear safety
#' 
#' @param x Array of size 3 or matrix with 3 columns. All input variables are 
#' in the interval [0,1].
#' @return If x is a matrix with n rows and 3 columns, the 
#' function is evaluated on each row and returns an array of size n. If x is an array, 
#' the function returns a scalar (the k-effective).
#' @details This function is important in robust inversion since the first parameter 
#' of the function is considered as a nuisance parameter while the second and third 
#' parameters are controlled.
#' @export
#' @author Gregory Caplin \email{gregory.caplin@@irsn.fr}
#' @examples 
#' keff(runif(3))
keff <- function(x){
  data(b2_l_kinf_cx_hx)
  
  #cx=seq(f=0.035,t=18.9197,l=700)
  kinf <- function(cx) approx(x=b2_l_kinf_cx_hx$cx, y=b2_l_kinf_cx_hx$kinf,xout=cx)$y
  b2 <-   function(cx) approx(x=b2_l_kinf_cx_hx$cx, y=b2_l_kinf_cx_hx$b2  ,xout=cx)$y
  l <-    function(cx) approx(x=b2_l_kinf_cx_hx$cx, y=b2_l_kinf_cx_hx$l   ,xout=cx)$y

  tmp <- function(m,d,cx){
    return(kinf(cx) / (1 + (kinf(cx)-1) / b2(cx) * 
                        (2.405^2/(d/2 + l(cx))^2 + 
                          pi^2/(4*m/(pi*d^2*cx) + 2*l(cx))^2) 
                       )
           )
  }
  
  tmp2 <- function(y) return(tmp(y[1],y[2],y[3]))
  
  if(is.null(dim(x))){
    arg1 <- 500+x[3]*(45000-500)
    arg2 <- 7.8+x[2]*(25-7.8)
    arg3 <- 0.035+x[1]*(18.9197-0.035)
  
    return(tmp(arg1,arg2,arg3))
  }else{
    # x is a matrix, we apply the function to each row
    mat <- cbind( 500+x[,3]*(45000-500) , 7.8+x[,2]*(25-7.8) , 0.035+x[,1]*(18.9197-0.035) )
    return(apply(X = mat,FUN = tmp2,MARGIN = 1))
  }  
}