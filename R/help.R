#'RobustInv - Robust Inversion of Expensive Black-box functions
#'
#'
#' @description The robust inversion problem in computer experiments can be seen as a generalization 
#' of the inversion problem tackled in the KrigInv package. It applies to functions with inputs that can 
#' be classified in two categories: controlled parameters and nuisance parameters. 
#' The input domain D, of dimension d can thus be written as a tensor product D = Dinv x Dopt, where 
#' Dinv is the space of the d.inv controlled parameters and Dopt is the space of the d.opt nuisance 
#' parameters. 
#' 
#' The goal is thus to build a sequential design of experiment which aims at identifying the set
#' Gamma* = \{ xinv in Dinv : f(xinv,xopt) < T for all xopt in Dopt \}. 
#' In short, we are interested in the configuration of controlled parameters such that the system at 
#' hand f remains 'safe' (i.e. below a target threshold T) for all possible values of the nuisance 
#' parameters. 
#' 
#' The suffix 'inv' is often used in this package when we refer to the controlled parameters. 
#' This is due to the fact that the inversion is actually performed in the set of controlled parameters. 
#' The suffix 'opt' is often used when we refer to the nuisance parameters. This is due to the fact that 
#' the set Gamma* that we aim at identifying can be rewritten:
#' Gamma* = \{ xinv in Dinv : max_xopt f(xinv,xopt) < T \}; meaning that some kind of optimization 
#' is performed with respect to the nuisance parameters. When xinv is fixed the optimizer 
#' corresponds to the most penalysing (or most dangerous) value of xopt, 
#' leading to the highest response f(xinv,xopt).
#' 
#' This package shares many similarities with the KrigInv package and a good understanding of KrigInv 
#' is necessary to use it properly. The sampling criteria used to build the sequential design 
#' of experiments are detailed in Clement Chevalier's PhD manuscript.
#' @author Clement Chevalier \email{clement.chevalier@@unine.ch}, 
#' Yann Richet \email{yann.richet@@irsn.fr}, 
#' Gregory Caplin \email{gregory.caplin@@irsn.fr}
#' @docType package
#' @name RobustInv
#' @importFrom Rcpp evalCpp
#' @useDynLib RobustInv
#' @aliases RobustInv RobustInv-package
NULL