#' @title Using transparence to modify colors
#' @description Modifies a color by playing on its transparence. Low alpha make the color more transparent 
#' while larger alpha make it opaque.
#' @param colors A color or an array where each element is a color.
#' @param alpha Scalar taking values between 0 (fully transparent) and 1 (opaque).
#' @return An array with the colors transluded with opacity alpha
#' @export
#' @author Yann Richet \email{yann.richet@@irsn.fr}
#' @examples 
#' library("grid")
#' onecolor <- "blue"
#'
#' shades.count <- 10
#'
#' alpha_tab <- seq(from = 0, to = 1, length=shades.count)
#' result <- rep("",times = shades.count)
#' 
#' for(i in 1:shades.count) result[i] <- translude(colors=onecolor,alpha_tab[i])
#' plotCol(result,txt.col=result)
translude <- function(colors, alpha = 0.6) {
  
  alpha <- rep(alpha, length.out = length(colors))
  rgb <- as.matrix(col2rgb(colors)/255)
  colors2 <- rgb(red = rgb["red", ],
                 green = rgb["green", ],
                 blue = rgb["blue", ],
                 alpha = alpha)
  
  return(colors2)
}
