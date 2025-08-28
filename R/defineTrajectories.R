#' Trajectory definition
#' 
#' Defines data structures for trajectory analysis
#'
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecological states..
#' @param sites A character vector indicating the ecological entity (site, individual, community) corresponding to each ecological state (other types are converted to character).
#' @param surveys An integer vector indicating the survey corresponding to each ecological state (only necessary when surveys are not in order and \code{times} is not provided).
#' @param times A numeric vector indicating survey times.
#'
#' @returns An object (list) of class \code{trajectories} with the following elements:
#' \itemize{
#'    \item{\code{d}: An object of class \code{\link{dist}} containing relationships between ecological states}
#'    \item{\code{metadata}: A data frame describing trajectory states, with the following columns:
#'      \itemize{
#'        \item{\code{sites}: A character vector indicating the ecological entity corresponding to each ecological state.}
#'        \item{\code{surveys}: An integer vector indicating the survey corresponding to each ecological state.}
#'        \item{\code{times}: A numeric vector indicating survey times.}
#'      } 
#'    }
#' }
#' 
#' @details
#' If \code{surveys} is not provided, but \code{times} is available, surveys will be taken as the order of times. Otherwise, \code{surveys} will be assumed to be in order for all the occurrences of the same value of \code{sites}.
#' If \code{times} is not provided, then it is made equal to \code{surveys}.
#'
#' @aliases trajectories
#' @seealso \code{\link{subsetTrajectories}}
#' @export
#' @examples 
#' #Description of entities (sites) and surveys
#' entities <- c("1","1","1","2","2","2")
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
#' x <- defineTrajectories(d, entities, surveys)
#' x
defineTrajectories<-function(d, sites, surveys = NULL, times = NULL) {
  # Check inputs
  if(!inherits(d,"dist") && !inherits(d, "matrix")) {
    stop("'d' should be of class `dist` or `matrix`")
  } else {
    if(inherits(d, "matrix")){
      if(nrow(d)!=ncol(d)) stop("Number of rows should be equal to number of columns in 'x'")
      d <- as.dist(d)
    } 
  }
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(times)) {
    if(length(sites)!=length(times)) stop("'sites' and 'times' need to be of the same length")
    if(!is.numeric(times)) stop("'times' should be a numeric vector")
  }
  if(!is.null(surveys)) {
    if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
    if(!is.numeric(surveys)) stop("'surveys' should be a numeric vector of integers")
    if(!all(as.integer(surveys)==surveys)) stop("'surveys' should be a numeric vector of integers")
  }
  
  # Imputation of missing surveys and times
  if(is.null(surveys)) {
    if(!is.null(times)) { # If times is provided make surveys be the order of times
      surveys <- rep(NA, length(sites))
      for(s in unique(sites)) {
        surveys[sites==s] <- order(times[sites==s])
      }
    } else { # Otherwise take surveys as the order of sites provided
      surveys <- rep(NA, length(sites))
      for(s in unique(sites)) {
        surveys[sites==s] <- 1:sum(sites==s)
      }
    }
  }
  if(is.null(times)) { # If times is not provided, make times equal to surveys
    times <- surveys
  }
  
  df <- data.frame(sites = as.character(sites),
                   surveys = as.integer(surveys),
                   times = as.numeric(times))
  
  l <- list(d = d, metadata = df)
  class(l) <- c("trajectories", "list")
  return(l)
}
