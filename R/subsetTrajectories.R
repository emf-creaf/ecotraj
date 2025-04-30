#' Trajectory subsetting
#' 
#' Subsets data structures for trajectory analysis
#'
#' @param x An object of class \code{trajectories} (or its children subclasses \code{fd.trajectories} or \code{cycles})
#' @param site_selection A character vector indicating the subset of entity (site) trajectories to be selected (if NULL, all sites are included).  
#' @param subtrajectory_selection A character vector indicating the subset of cycles or fixed date trajectories to be selected (only used when \code{x} is of class \code{fd.trajectories} or \code{cycles}).
#' @param survey_selection An integer vector indicating the subset of surveys to be included (if NULL, all surveys are included).
#'
#' @returns An object (list) of class \code{trajectories} (or its children subclasses \code{fd.trajectories} or \code{cycles}), depending on the input.
#' @details
#' When using function \code{subsetTrajectories} on cycles or fixed-date trajectories then the parameter \code{site_selection} applies to sites 
#' (hence allows selecting multiple cycles or fixed-date trajectories). Specific cycles or fixed-date trajectories can be selected using \code{trajectory_selection}.
#' 
#' @seealso \code{\link{defineTrajectories}}, \code{\link{trajectoryCyclical}}
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
#' 
#' # Extracts (subset) second trajectory
#' x_2 <- subsetTrajectories(x, "2")
#' x_2
subsetTrajectories<-function(x, site_selection = NULL, subtrajectory_selection = NULL, survey_selection = NULL) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  sites <- x$metadata$sites
  selection <- rep(TRUE, length(sites))
  if(!is.null(site_selection)) {
    if(!inherits(site_selection, "character")) stop("`site_selection` must be a character vector")
    if(!all(site_selection %in% sites)) stop("At least one element in `site_selection` is not a trajectory site.")
    selection <- selection & (sites %in% site_selection)
  }
  if(!is.null(subtrajectory_selection)) {
    if(!inherits(subtrajectory_selection, "character")) stop("`subtrajectory_selection` must be a character vector")
    if(inherits(x, "fd.trajectories")) {
      subtrajectories <- x$metadata$fdT
    } else if(inherits(x, "cycles")) {
      subtrajectories <- x$metadata$cycles
    } else if(inherits(x, "sections")) {
      subtrajectories <- x$metadata$sections
    } else {
      subtrajectories <- NULL
      warning("Parameter `subtrajectory_selection` was ignored.")
    }
    if(!is.null(subtrajectories)) {
      if(!all(subtrajectory_selection %in% subtrajectories)) stop("At least one element in `subtrajectory_selection` is not a cycle/fixed-date trajectory.")
      selection <- selection & (subtrajectories %in% subtrajectory_selection)
    }
  }
  if(!is.null(survey_selection)) {
    if(!is.numeric(survey_selection)) stop("`survey_selection` must be a numeric (integer) vector")
    selection <- selection & (x$metadata$surveys %in% as.integer(survey_selection))
  }
  msel <- as.matrix(x$d)[selection, selection]
  rownames(msel) <- 1:nrow(msel)
  colnames(msel) <- 1:ncol(msel)
  dsel <- as.dist(msel)
  dfsel <- x$metadata[selection,]
  row.names(dfsel) <- 1:nrow(dfsel)
  #Reorder surveys
  for(si in site_selection) {
    if(inherits(x, "fd.trajectories")) {
      sel <- which(dfsel$fdT == si) 
    } else if(inherits(x, "cycles")) {
      sel <- which(dfsel$cycles == si) 
    } else if(inherits(x, "sections")) {
      sel <- which(dfsel$sections == si) 
    } else {
      sel <- which(dfsel$sites == si) 
    }
    dfsel$surveys[sel] <- order(dfsel$surveys[sel])
  }
  l <- list(d = dsel, metadata = dfsel)
  class(l) <- class(x)
  return(l)
}
