#' eccficm: Expected Conditional Characteristic Function-based
#' Independence Criterion measures
#'
#' A package for computing new correlation/covariance type measures
#' for testing independence between two random vectors using the expected
#' conditional characteristic function-based independence criterion (ECCFIC)
#' methods developed by \insertCite{yin2019new}{eccficm} and
#' \insertCite{ke2019expected}{eccficm}.
#'
#' The eccficm package provides four main functions (
#' ecd, eccfic, ecd.test and eccfic.test) related to the measuring and
#' testing of independence, and two additional functions - one for simulating data
#' and the other for screening important predictors.
#'
#' @section Available Functions:
#' \itemize{
#' \item{\code{\link{ecd}} and \code{\link{eccfic}}: compute the
#' covariance and correlation type statistics.}
#'
#' \item{\code{\link{ecd.test}} and  \code{\link{eccfic.test}}: perform a
#' permutation test of independence based on the ECD and ECCFIC correlation
#' and covariance type statistics, respectively.}
#'
#' \item{\code{\link{generateData}}: generates optional data for feature screening based on
#' some specific models.}
#'
#' \item{\code{\link{fscreen}}: Performs the sure independence feature screening using
#' correlation learning between each predictor and the response.}
#'}
#'
#' @docType package
#' @name eccficm
NULL
