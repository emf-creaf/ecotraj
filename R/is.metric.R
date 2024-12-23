#' Metricity
#' 
#' Checks whether the input dissimilarity matrix is metric (i.e. all triplets fulfill the triangle inequality).
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states (see details).
#' @param tol Tolerance value for metricity
#' 
#' @return a boolean indicating metric property
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @encoding UTF-8
#' @export
#' @keywords internal
is.metric<-function(d, tol=0.0001) {
  return(.ismetricC(as.matrix(d), tol))
}
