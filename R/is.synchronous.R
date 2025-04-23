#' Synchronicity in trajectory observations
#' 
#' Checks whether trajectories are synchronous, meaning that observation times are equal
#'
#' @param x An object of class \code{trajectories} (or its children subclasses \code{fd.trajectories} or \code{cycles})
#'
#' @returns A boolean indicating whether trajectories are synchronous
#' 
#' @seealso \code{\link{defineTrajectories}}
#' @export
#'
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
#' #Synchronous trajectories
#' x1 <- defineTrajectories(dist(xy), sites, surveys)
#' is.synchronous(x1)
#' 
#' # Non synchronous trajectories
#' x2 <- defineTrajectories(dist(xy[1:5,]), sites[1:5], surveys[1:5])
#' is.synchronous(x2)
is.synchronous<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  times <- x$metadata$times
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  unique_sites <- unique(sites)
  unique_times <- sort(unique(times))
  for(i in 1:length(unique_sites)) {
    unique_times_i <- sort(unique(times[sites==unique_sites[i]]))
    if(length(unique_times_i)< length(unique_times)) return(FALSE)
    if(!all(unique_times==unique_times_i)) return(FALSE)
  }
  return(TRUE)
}