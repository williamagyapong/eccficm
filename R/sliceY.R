#-------------------------------------------------------------
# A function for slicing a continuous Y into finite categories
#_____________________________________________________________
sliceY <-
function(y, ns) {
  # y is the random vector we condition on and must be in a matrix form
  # ns is the number of slices/levels
  if(!('matrix'%in%class(y))) y <- as.matrix(y)
  n <- nrow(y) # number of observations in y

  if(ns > n) stop("Number of slices(ns) cannot be more than the number of observations in Y")

  ysliced <- as.matrix(ceiling(rank(y,ties.method = "min")*ns/n))

  # check whether number of obs per slice are enough/too many for good result
  max.obs <- max(as.data.frame(table(ysliced))$Freq) # maximum number of observations per slice
  upper.bound <- n/2
  if(max.obs < 5 | max.obs > upper.bound) {
    warning(paste0("Number of slices is either too small or too large for good result.
    Choose the number of slices such that the number of data points in each slice is greater
    than or equal to 5 but not more than ", upper.bound))
  }

  return(ysliced)
}
