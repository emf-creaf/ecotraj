#' Trajectory definition
#' 
#' Defines data structures for trajectory analysis
#'
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states..
#' @param sites A character vector indicating the site corresponding to each ecosystem state (other types are converted to character).
#' @param surveys An integer vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @param times A numeric vector indicating survey times (if missing, survey times are made equal to surveys).
#'
#' @returns An object (list) of class \code{trajectories} with the following elements:
#' \itemize{
#'    \item{\code{d}: An object of class \code{\link{dist}} containing relationships between ecosystem states}
#'    \item{\code{metadata}: A dataframe describing trajectory states, with the following columns:
#'      \itemize{
#'        \item{\code{sites}: A character vector indicating the site corresponding to each ecosystem state.}
#'        \item{\code{surveys}: An integer vector indicating the survey corresponding to each ecosystem state.}
#'        \item{\code{times}: A numeric vector indicating survey times.}
#'      } 
#'    }
#' }
#' 
#' @aliases trajectories
#' @seealso \code{\link{subsetTrajectories}}
#' @export
#' @examples 
#' #Description of sites and surveys
#' sites <- c("1","1","1","2","2","2")
#' surveys <- c(1,2,3,1,2,3)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=6, ncol=2)
#' xy[2,2]<-1
#' xy[3,2]<-2
#' xy[4:6,1] <- 0.5
#' xy[4:6,2] <- xy[1:3,2]
#' xy[6,1]<-1
#' 
#' d <- dist(xy)
#' 
#' # Defines trajectories
#' x <- defineTrajectories(d, sites, surveys)
#' x
defineTrajectories<-function(d, sites, surveys = NULL, times = NULL) {
  if(!inherits(d,"dist") && !inherits(d, "matrix")) {
    stop("'d' should be of class `dist` or `matrix`")
  } else {
    if(inherits(d, "matrix")){
      if(nrow(d)!=ncol(d)) stop("Number of rows should be equal to number of columns in 'x'")
      d <- as.dist(d)
    } 
  }
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) {
    if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
    if(!is.numeric(surveys)) stop("'surveys' should be a numeric vector of integers")
    if(!all(as.integer(surveys)==surveys)) stop("'surveys' should be a numeric vector of integers")
  } else {
    surveys <- rep(NA, length(sites))
    for(s in unique(sites)) {
      surveys[sites==s] <- 1:sum(sites==s)
    }
  }
  if(!is.null(times)) {
    if(length(sites)!=length(times)) stop("'sites' and 'times' need to be of the same length")
    if(!is.numeric(times)) stop("'times' should be a numeric vector")
  } else {
    times <- surveys
  }
  
  df <- data.frame(sites = as.character(sites),
                   surveys = as.integer(surveys),
                   times = as.numeric(times))
  
  l <- list(d = d, metadata = df)
  class(l) <- c("trajectories", "list")
  return(l)
}
