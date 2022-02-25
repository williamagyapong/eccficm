#' @title Data Generation Function
#' @name generateData
#' @description Generates data for demonstrating the performance of Sure Independence Screening (SIS)
#' procedures.
#' @param n The number of observations or sample size.
#' @param p Dimension of the data matrix of predictors to be generated.
#' @param rho Correlation coefficient between predictors.
#' @param model A numeric value representing the model for the response.
#' See the 'Details' section for available options.
#' @param sigma The variance of the error term. Defaults to 1.
#' @param penalize A boolean indicating whether the penalizing constants should
#' be used or not.
#' @param beta A vector of model beta coefficients. Length must correspond to the
#' number of beta coefficients in the model. It is only used when the fixed coefficient
#' models numbered 5 - 8 are desired.
#' be used or not. These constants are meant to challenge the screening procedures.
#'@return \code{generateData()} returns a list object of the following:
#'\describe{
#'\item{X}{A data matrix or object of predictors of dimension p.}
#'\item{Y}{A vector of response. This implementation only generates a univariate
#' response.}
#'}
#'@details
#' The response is simulated based on models provided in example 3.1 of
#' \insertCite{li2012feature}{eccficm}. The available models are
#' numbered 1 - 4, from first to last in the order as they appear in the paper.
#' @author William Agyapong and Qingcong Yuan (\email{yuanq3@@miamioh.edu})
#'
#' @references \insertRef{li2012feature}{eccficm}


# Adapted from VariableScreening package
#'@export
generateData <- function(n=200,
                          p=2000,
                          rho=0,
                          model= 1,
                          sigma = 1,
                          penalize=T,
                          beta = c(2.01,  3.98, -2.51, -1.79)) {
  # if ((p<30)|(p>100000)) {stop("Please select a number p of predictors between 30 and 100000.")}
  if ((rho<0)|(rho>=1)) {stop("Please select a rho parameter in the interval [0,1).")}

  # Simulate the predictors
  if (rho != 0) {
    X <- matrix(NA,n,p)
    X[,1] <- rnorm(n)
    for (i in 2:p)  {
      X[,i]=rho*X[,i-1]+sqrt(1-rho^2)*rnorm(n)
    }
  } else {
    # We can save computational time by just treating all the X as independent.
    X <- matrix(rnorm(n*p),n,p)
  }

  # Simulate Y
  # Models based on Example 1 in section 3 of Li, Zhong & Zhu (2012)
  consts <- c(2,0.5,3,2) # predefined model constants given in the paper to
  #challenge the feature screening procedures

  a <- 4*log(n)/sqrt(n)
  e <- rnorm(n,0,sigma)  # Residuals for Y
  if (model == 1) {
    Z <- rnorm(4)
    U <- rbinom(4, size = 1, prob = 0.4)
    if(penalize) {coeff <- c(2,0.5,3,2) * ((-1)^U)*(a + abs(Z))} else {coeff<-((-1)^U)*(a + abs(Z))}
    Y = coeff[1]*X[,1] + coeff[2]*X[,2] + coeff[3]*(X[,12]<0) + coeff[4]*X[,22] + e

  } else if (model==2){
    # set.seed(22)
    Z <- rnorm(3)
    U <- rbinom(3, size = 1, prob = 0.4)
    if(penalize) {coeff <- c(2,3,2) * ((-1)^U)*(a + abs(Z))} else {coeff<-((-1)^U)*(a + abs(Z))}
    Y = coeff[1]*X[,1]*X[,2] + coeff[2]*(X[,12]<0) + coeff[3]*X[,22] + e

  } else if (model==3){
    Z <- rnorm(2)
    U <- rbinom(2, size = 1, prob = 0.4)

    if(penalize) { coeff <- c(2,3) * ((-1)^U)*(a + abs(Z))} else {coeff<-((-1)^U)*(a + abs(Z))}
    Y = coeff[1]*X[,1]*X[,2] + coeff[2]*(X[,12]<0)*X[,22] + e

  } else if (model==4){
    Z <- rnorm(3)
    U <- rbinom(3, size = 1, prob = 0.4)
    if(penalize) {coeff <- c(2,0.5,3) * ((-1)^U)*(a + abs(Z))}else{coeff<-((-1)^U)*(a + abs(Z))}
    Y = coeff[1]*X[,1] + coeff[2]*X[,2] + coeff[3]*(X[,12]<0) + exp(2*abs(X[,22]))* e
    # message(round(coeff, 2))
  }
  else if (model==5){
    # same as model 1 except that constant betas generated
    # from the underlying distribution with seed set to 22 are used.
    betas <- c(2.01,  3.98, -2.51, -1.79)# c(2.010730,  3.983774, -2.506417, -1.791405)
    if(penalize) {coeff <- c(2,0.5,3,2) * betas}else{coeff<-betas}
    Y = coeff[1]*X[,1] + coeff[2]*X[,2] + coeff[3]*(X[,12]<0) + coeff[4]*X[,22] + e
    # message(round(coeff, 2))
  }
  else if (model==6){
    # same as model 4 except that constant betas generated
    # from the underlying distribution with seed set to 22 are used.
    betas <- c(2.01,  3.98, -2.51) # c(2.010730,  3.983774, -2.506417, -1.791405)
    if(penalize) {coeff <- c(2,0.5,3) * betas}else{coeff<-betas}
    Y = coeff[1]*X[,1] + coeff[2]*X[,2] + coeff[3]*(X[,12]<0) + exp(2*abs(X[,22]))* e

  }else if (model==8){
    Z <- rnorm(3)
    U <- rbinom(3, size = 1, prob = 0.4)
    if(penalize) {coeff <- c(2,0.5,3) * ((-1)^U)*(a + abs(Z))}else{coeff<-((-1)^U)*(a + abs(Z))}
    Y = coeff[1]*X[,1] + coeff[2]*X[,2] + coeff[3]*(X[,12]<0) + exp(2*abs(X[,22]))* e
    # message(round(coeff, 2))
  }
  else {
    stop("Please choose a model number in (1,2,3,4)")
  }

  return(list(X=X, Y=Y))
}
