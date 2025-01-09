#' Trajectory projection
#' 
#' Performs an projection of a set of target points onto a specified trajectory and returns the distance to the trajectory (i.e. rejection) and the relative position of the projection point within the trajectory.
#'
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states (see details).
#' @param target An integer vector of the ecosystem states to be projected.
#' @param trajectory An integer vector of the ecosystem states conforming the trajectory onto which target states are to be projected.
#' @param tol Numerical tolerance value to determine that projection of a point lies within the trajectory.
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#'
#' @returns 
#' A data frame with the following columns:
#' \itemize{
#'   \item{\code{distanceToTrajectory}: Distances to the trajectory, i.e. rejection. If there is no orthogonal projection the distance corresponds to the minimum distance to the trajectory.}
#'   \item{\code{segment}: Segment that includes the projected point or the closest state.}
#'   \item{\code{relativeSegmentPosition}: Relative position of the projected point within the segment, i.e. values from 0 to 1 with 0 representing the start of the segment and 1 representing its end.}
#'   \item{\code{relativeTrajectoryPosition}: Relative position of the projected point within the trajectory, i.e. values from 0 to 1 with 0 representing the start of the trajectory and 1 representing its end.}
#' }
#' @export
trajectoryProjection<-function(d, target, trajectory, tol = 0.000001, add=TRUE) {
  if(length(trajectory)<2) stop("Trajectory needs to include at least two states")
  if(length(trajectory)!=length(unique(trajectory))) stop("Trajectory states must be different")
  dmat <- as.matrix(d)
  npoints <- length(target)
  nsteps <- length(trajectory) -1
  #Distance between target points and trajectory points
  d2ref <- dmat[target, trajectory, drop=FALSE]
  #Distance between trajectory steps
  dsteps <- diag(dmat[trajectory[1:(length(trajectory)-1)], trajectory[2:length(trajectory)]])
  #Cumulative distance between steps
  dstepcum <- rep(0,nsteps+1)
  if(nsteps>1) {
    for(i in 2:nsteps) {
      dstepcum[i] <- dstepcum[i-1]+dsteps[i-1]
    }
  }
  dstepcum[nsteps+1] <- sum(dsteps)
  
  projA1 <- matrix(NA, nrow=npoints, ncol = nsteps)
  whichstep <- rep(NA, npoints)
  dgrad <- rep(NA, npoints)
  posgradseg <- rep(NA, npoints)
  posgradtraj <- rep(NA, npoints)
  
  for(i in 1:npoints) {
    for(j in 1:nsteps) {
      p <- .distanceToSegmentC(dsteps[j], d2ref[i, j], d2ref[i, j+1], add)
      if((p[1]>-tol) & (p[2]>-tol)) {
        projA1[i,j] <- p[1]
        if(is.na(dgrad[i])) {
          dgrad[i] <- p[3]
          whichstep[i] <- j
          posgradseg[i] <- p[1]/dsteps[j]
        } else {
          if(p[3]<dgrad[i]) { # Replace if closer
            dgrad[i] <- p[3]
            whichstep[i] <- j
            posgradseg[i] <- p[1]/dsteps[j]
          }
        }
      }
    }
    if(!is.na(whichstep[i])) { # If we determined a closest segment determine trajectory relative postion 
      dg <- dstepcum[whichstep[i]]+projA1[i,whichstep[i]]
      posgradtraj[i] <- dg/sum(dsteps)
    }
  }
  # Assemble output
  res <- data.frame(distanceToTrajectory=dgrad, segment = whichstep, relativeSegmentPosition = posgradseg, relativeTrajectoryPosition = posgradtraj)
  row.names(res)<-row.names(d2ref)
  return(res)
}



