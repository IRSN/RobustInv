% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_robinv_uncertainty_1d.R
\name{print_robinv_uncertainty_1d}
\alias{print_robinv_uncertainty_1d}
\title{Excursion probability and maximum kriging mean plot in 1d}
\usage{
print_robinv_uncertainty_1d(model, T, lower, upper, opt.index, inv.index,
  control = NULL, onlycompute = FALSE, lwd = 1, cex.main = 1,
  cex.axis = 1, cex.lab = 1, cex.points = 1, pch.points = 17,
  color.up = "red", color.down = "blue", maxmkplot = TRUE,
  xlab1 = NULL, ylab1 = NULL, main1 = NULL, pnplot = TRUE,
  xlab2 = NULL, ylab2 = NULL, main2 = NULL, newpar = TRUE)
}
\arguments{
\item{model}{The current kriging model. km object.}

\item{T}{Target threshold.}

\item{lower}{Array of size d. Lower bound of the input domain.}

\item{upper}{Array of size d. Upper bound of the input domain.}

\item{opt.index}{Array with integers corresponding to the indices of the nuisance parameters.}

\item{inv.index}{Array with integers corresponding to the indices of the controlled parameters.}

\item{control}{A list with fields that will control how the different quantities involved are computed. 
\code{control$resolution} is the number of points (in the space of the controlled parameters) where the computations 
are performed. A grid in dimension 1 of size resolution is used. The field n.optpoints controls the number of points 
(in the space of the nuisance parameters) taken to compute the maximum of the kriging mean, or to do conditional 
simulations. For the computation of the excursion probability, the fields n.optpoints.candidates, 
choose_optpoints, nsimu, unscale.opt.simulation.points are used in a similar way than in the \code{integration_design_robinv} function. See 
help on the integcontrol argument in that function.}

\item{onlycompute}{Boolean. When FALSE, no plot is performed, but the maximum of the kriging mean (resp. the excursion probability pn) 
is still computed if \code{maxmkplot=TRUE} (resp. pnplot = TRUE).}

\item{lwd}{Line width for the different plots involved in this function.}

\item{cex.main}{Title size for the different plots involved in this function.}

\item{cex.axis}{Axis label size for the different plots involved in this function.}

\item{cex.lab}{Label size for the different plots involved in this function.}

\item{cex.points}{Point size for the maximum kriging mean plot. Useless if \code{maxmkplot=FALSE}.}

\item{pch.points}{Point pch for the maximum kriging mean plot. Useless if \code{maxmkplot=FALSE}.}

\item{color.up}{Color of the points where there is threshold exceedance. Useless if \code{maxmkplot=FALSE}.}

\item{color.down}{Color of the points where there is no threshold exceedance. Useless if \code{maxmkplot=FALSE}.}

\item{maxmkplot}{Boolean. When TRUE, the maximum of the kriging mean (taken w.r.t. the nuisance parameters) is 
computed. It is also ploted if \code{onlycompute=FALSE}.}

\item{xlab1}{x axis label for the maximum of the kriging mean plot.}

\item{ylab1}{y axis label for the maximum of the kriging mean plot.}

\item{main1}{Title of the maximum of the kriging mean plot.}

\item{pnplot}{Boolean. When TRUE, the excursion probability function, pn, is computed. 
It is also ploted if \code{onlycompute=FALSE}.}

\item{xlab2}{x axis label for the excursion probability plot.}

\item{ylab2}{y axis label for the excursion probability plot.}

\item{main2}{Title of the excursion probability plot.}

\item{newpar}{Boolean. When TRUE, the par() function is called. Usefull only if the two plots available in this function are both performed.}
}
\value{
A list containing the important computed quantities: 
(i) all.points: the value of the controlled parameters, 
(ii) maxmk: array of the same size than all.points containing the maximum of the kriging mean taken w.r.t. the nuisance 
parameters, 
(iii) pn: array of the same size than all.points containing the excursion probability of the considered points (all.points), 
(iv) uncertainty: scalar equal to \code{mean(pn*(1-pn))} giving a measure of the current global uncertainty on the excursion set, 
(v) colors.transluded: the colors used to plot the points on the maximum of kriging mean plot.
}
\description{
Computes and can show two plots for functions with only one controlled parameters 
and one or several nuisance parameters. One plot is the maximum of the kriging mean 
(taken w.r.t. the nuisance parameters) computed on a grid of controlled parameter values. 
The second plot is the excursion probability, pn, which is 
computed using conditional simulations.
}
\examples{

library(KrigInv)
library(randtoolbox)
myfun <- function(x) return(-1 * branin_robinv(x))
d <- 3

set.seed(8)

n0 <- 30
T <- -10
opt.index <- c(2,3)
inv.index <- c(1)
lower <- rep(0,times=d)
upper <- rep(1,times=d)

design <- matrix(runif(d*n0),nrow=n0)
response <- myfun(design)
model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2")
control <- list(resolution = 100, n.optpoints = 300)
\dontrun{
print_robinv_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
                            opt.index = opt.index,inv.index = inv.index,
                            control = control)
}
########################################
# A more complicated example with scaling
library(KrigInv)
myfun <- function(x){ return(-1*branin_robinv(x) - 50*sin(min(100,1/x[3]))  ) }
d <- 3

set.seed(8)

n0 <- 60
T <- 40
opt.index <- c(2,3)
inv.index <- c(1)
lower <- rep(0,times=d)
upper <- rep(1,times=d)

design <- matrix(runif(d*n0),nrow=n0)
response <- apply(X = design,FUN = myfun,MARGIN = 1) 
knots.number <- c(0,3,3)
knots <- generate_knots(knots.number = knots.number , d = d)

model <- km(formula = ~1,design = design,response = response,covtype = "matern3_2",scaling = TRUE,knots=knots)
# have a look at model@covariance@eta

control <- list(resolution = 100, n.optpoints = 300, unscale.opt.simulation.points=TRUE)

\dontrun{
print_robinv_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
                            opt.index = opt.index,inv.index = inv.index,
                            control = control)
}
}
\author{
Clement Chevalier \email{clement.chevalier@unine.ch}
}
