#' @title Sure Independence Screening
#' @name fscreen
#' @description Performs the sure independence screening using correlation learning
#' between each predictor and the response.
#' @param X A data matrix or object of potential predictors. Can be of any dimension.
#' @param Y A vector of response. This implementation only supports a univariate
#' response.
#' @param method The screening procedure based on a particular correlation
#' measure. Four methods are available: "SIS", "DC-SIS", "ECD-SIS", and "ECCFIC-SIS".
#' The method names are case insensitive meaning one can either type DC-SIS or dc-sis
#' for the same method. See the 'Details' section for more information.
#' @param est The estimation procedure for ECD or ECCFIC. Use 'slice' for the
#' slicing approach and 'kernel' for the kernel estimation procedure.
#' @param ns The number of slices to use when the slicing estimation procedure
#' is specified.
#'
#' @return \code{fscreen} returns a list object containing the following:
#' \describe{
#' \item{utility}{ A vector of the marginal correlations between each predictor/feature
#' and the response Y. Each value is a measure of the estimated strength of
#' association with Y.}
#' \item{rank}{A vector of the position for each predictor/feature in terms of their
#' predictive power. Ranks are assigned from 1 to the number of columns in the feature
#' matrix X, with 1 indicating the best predictive performance.
#' The lower the rank the more active the predictor at that index is in explaining/predicting Y.}
#' }
#' @section Screening methods:
#' \itemize{
#' \item{SIS: Sure Independence Screening based on the Pearson correlation coefficient}
#' \item{DC-SIS: Distance Correlation Sure Independence Screening. Uses dcor() from
#' the energy package.}
#' \item{ECD-SIS: Expectation of the Conditional Difference Sure Independence Screening.
#'  When 'est' is set default the slicing method is used, otherwise, the kernel
#'  estimation method can be specified via the 'est' parameter.}
#'  \item{ECCFIC-SIS: The default method uses the kernel regression estimation.
#'  Otherwise, the slicing procedure can be specified through the 'est' parameter.}
#' }
# Feature screening function
#'@export
fscreen <- function(X, Y, method="ECD-SIS", est="default", ns=NULL) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  # X <- data1$X
  # Y <- data1$Y
  # method <- "ecdsis"
  #  method <- "eccficsis"
  # method = "SIS"

  method <- toupper(method)
  if(method=="ECDSIS") method <- "ECD-SIS" # defaults to slicing
  if(method=="ECCFICSIS") method <- "ECCFIC-SIS"
  if(method=="DCSIS") method <- "DC-SIS"


  # Choose the correlation type function
  if(method=="ECD-SIS") cor_fun <- eccficm::ecd
  if(method=="ECCFIC-SIS") cor_fun <- eccficm::eccfic
  if(method=="DC-SIS") cor_fun <- energy::dcor
  if(method=="SIS") cor_fun <- stats::cor

  #computing distance correlation between response Y and each predictor

  if((method != "SIS") & (method != "DC-SIS")) {
    if(est=="default") {

      omega_k <- function(Xk) cor_fun(Xk, Y)[[2]] # correlation type statistic is the 2nd returned value for ECCFIC measures

    } else {
      # ns <- ifelse((est=="slice")&!is.null(ns), ns, 5) # use 5 slices

      omega_k <- function(Xk) cor_fun(Xk, Y, est, ns)[[2]] # correlation type statistic is the 2nd returned value for ECCFIC measures
    }

  } else {
    omega_k <- function(Xk) cor_fun(Xk, Y)
  }

  # Compute the marginal strength of association
  cor <- abs(apply(X,2,omega_k))

  # Get ranks
  ordered_inds <- sort.list(importance, decreasing = T, method = "shell")
  rank <- match(importance, importance[ordered_inds]) # rank in order of estimated association strengths;
  # these ranks are the positions for each predictor
  # in order of predictive power /significance
  return(list(importance=importance, rank=rank))

}
