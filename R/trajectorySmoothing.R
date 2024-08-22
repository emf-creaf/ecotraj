#' Trajectory smoothing
#'
#' Trajectory smoothing using a Gaussian kernel
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites  A vector indicating the site corresponding to each ecosystem state.
#' @param surveys  A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @param survey_times A vector indicating the survey time for all surveys (if \code{NULL}, time between consecutive surveys is considered to be one)
#' @param kernel_scale Scale of the Gaussian kernel, related to survey times
#' @param fixed_endpoints A logical flag to force keeping the location of trajectory endpoints unmodified
#'
#' @return An object of class \code{\link{dist}} 
#' @export
#'
trajectorySmoothing<-function(d, sites, surveys = NULL, survey_times = NULL, kernel_scale = 1, fixed_endpoints = TRUE) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  
  dmat = as.matrix(d)
  
  umat_smooth = matrix(0, length(sites), length(sites))
  for(i1 in 1:length(sites)) {
    if(!is.null(surveys)) {
      surveys_i1 <- surveys
      surveys_i1[sites!=sites[i1]] <- 0
    } else {
      surveys_i1 <- cumsum(sites==sites[i1])
    }
    x1 <- surveys_i1[i1]
    is_endpoint <- (x1 == 1) || (x1 == sum(sites==sites[i1]))
    if(fixed_endpoints && is_endpoint) { # Keep point
      umat_smooth[i1, i1] <- 1.0
    } else { # Apply kernel
      # Translate to time axis if provided
      if(!is.null(survey_times)) {
        x1 <- survey_times[x1]
      }
      for(i2 in 1:length(sites)) {
        if(sites[i2]!=sites[i1]) { # Keep different trajectories separated
          umat_smooth[i1, i2] <- 0
        } else {
          x2 <- surveys_i1[i2]
          # Translate to time axis if provided
          if(!is.null(survey_times)) {
            x2 <- survey_times[x2]
          }
          umat_smooth[i1, i2] <- exp(-1*((x1-x2)^2)/kernel_scale)
        }
      }
      # Normalize to one
      umat_smooth[i1, ] <- umat_smooth[i1, ]/sum(umat_smooth[i1, ])
    }
  }
  dmat_smooth <- .distanceBetweenClusters(dmat, umat_smooth)
  return(as.dist(dmat_smooth))
}