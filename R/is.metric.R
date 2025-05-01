#' Metricity
#' 
#' Checks whether the input dissimilarity matrix is metric (i.e. all triplets fulfill the triangle inequality).
#' 
#' @param x Either an object of class \code{trajectories},  a symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecological states.
#' @param tol Tolerance value for metricity
#' 
#' @return A boolean indicating metric property
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @encoding UTF-8
#' @export
is.metric<-function(x, tol=0.0001) {
  if(inherits(x, "trajectories")) {
    return(.ismetricC(as.matrix(x$d), tol))
  } else if(inherits(x, "dist")) {
    return(.ismetricC(as.matrix(x), tol))
  } else if(inherits(x, "matrix")) {
    return(.ismetricC(x, tol))
  }
  stop("'x' should be of class `trajectories`, `dist` or `matrix`")
}
