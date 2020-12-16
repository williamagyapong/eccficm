#------------------------------------------------------------------------------
# This function depends on the sliceY() function available in the sliceY.R file
#______________________________________________________________________________

# @title ECD slicing method
# @description Implements ECD slicing method
# @param x data of first sample serving as the response
# @param y data of second sample conditioned on
# @param ns number of slices
# @param index exponent on Euclidean distance, in (0,2]

ecd.slice <-
function(x, y, ns=NULL,bw=NULL, index=1.0) {

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  q <- ncol(y)

  # validate inputs
  if(q>1)  stop("Sorry, this current implementation does not support a multivariate Y. Check out for updates.")

  if (n != m) stop("Sample sizes for x and y must agree")

  if (index <= 0 || index > 2) {
    warning("Index must be in (0,2], using default index=1")
    index=1.0 # reset index to default
  }

  x0=x[order(y),] # sort X according to the order of Y
  xdist=(as.matrix(dist(x0)))^index
  term1=sum(xdist)/(n^2)  # first term in eqn (9)/the same as the term in eqn (10)

  # computing Cn(X|X)
  CnXX <- CnXY <- sqrt(term1) # a measure analogous to variance of X (not very sure)

  if(!identical(x,y)) {
    # computing Cn(X|Y)
    if(!(is.factor(y) | is.character(y)) & !is.null(ns)) {
      # slice the continuous Y variable
      ysliced <- sliceY(y, ns)
    } else {
      ysliced <- y # use the originally categorical Y
    }
    sizes <- tabulate(factor(as.character(ysliced)))
    ns <- length(sizes)
    diagind=matrix(0,nrow=m,ncol=m)
    for(j in 1:ns){
      if(j==1) {
        beg <- j
      }
      else {
        beg <- (sum(sizes[1:(j-1)])+1)
      }
      end <- sum(sizes[1:j])
      diagind[beg:end,beg:end] <- matrix(rep(1/(sizes[j]),sizes[j]*sizes[j]),sizes[j])
    }

    term2 <- sum(xdist*diagind)/n # second term in eqn (9)
    CnXY <- sqrt(term1 - term2) # Cn(X|Y), a measure analogous to covariance between X and Y (not very sure)
  }

  rhoc <- CnXY/CnXX
  return(list(ecdCov=CnXY, ecdCor=rhoc, ecdVarX=CnXX))

}
