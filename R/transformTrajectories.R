#' Transform trajectories
#' 
#' The following functions are provided to transform trajectories:
#' \itemize{
#' \item{Function \code{smoothTrajectories} performs multivariate smoothing on trajectory data using a Gaussian kernel.}
#' \item{Function \code{centerTrajectories} shifts all trajectories to the center of the multivariate space and returns a modified distance matrix.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name transformTrajectories
#' @aliases centerTrajectories smoothTrajectories
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' 
#' @details 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (2019). 
#' Function \code{centerTrajectories} performs centering of trajectories using matrix algebra as explained in Anderson (2017).
#'
#' @return 
#' A modified object of class \code{\link{trajectories}}, where distance matrix has been transformed.
#' 
#' @author 
#' Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' Nicolas Djeghri, UBO
#' 
#' @references
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
#' Anderson (2017). Permutational Multivariate Analysis of Variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article ID: stat07841.
#' 
#' @seealso \code{\link{trajectoryPlot}} \code{\link{trajectoryMetrics}}
#' 


#' @rdname transformTrajectories
#' 
#' @param survey_times A vector indicating the survey time for all surveys (if \code{NULL}, time between consecutive surveys is considered to be one)
#' @param kernel_scale Scale of the Gaussian kernel, related to survey times
#' @param fixed_endpoints A logical flag to force keeping the location of trajectory endpoints unmodified
#' @export
smoothTrajectories<-function(x, survey_times = NULL, kernel_scale = 1, fixed_endpoints = TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")

  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] <- sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  
  dmat <- as.matrix(d)
  
  umat_smooth <- matrix(0, length(sites), length(sites))
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
          # Gaussian kernel
          umat_smooth[i1, i2] <- exp(-1*((x1-x2)^2)/(2*(kernel_scale^2)))
        }
      }
      # Normalize to one
      umat_smooth[i1, ] <- umat_smooth[i1, ]/sum(umat_smooth[i1, ])
    }
  }
  # Replace distance matrix with distances between clusters
  x$d <- .distanceBetweenClusters(dmat, umat_smooth)
  return(x)
}

#' @rdname transformTrajectories
#' @param exclude An integer vector indicating sites that are excluded from trajectory centroid computation
#' @export
centerTrajectories<-function(x, exclude = integer(0)) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
  if(length(exclude)>0) {
    if(!is.numeric(exclude)) stop("`exclude` needs to be an integer vector")
    exclude <- as.integer(exclude)
    if(max(exclude)>length(sites)) stop("`exclude` contains values outside range")
    if(min(exclude)<1) stop("`exclude` contains values outside range")
    s_non_excluded <- sites[-exclude] 
    for(s in unique(sites)) {
      if(sum(s_non_excluded==s)==0) stop("`exclude` cannot include all sites of a trajectory")
    }
  }
  
  Dmat <-as.matrix(d)
  
  # Anderson (2017). Permutational Multivariate Analysis of Variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article ID: stat07841.
  Amat <- (-0.5)*(Dmat^2)
  n <- nrow(Dmat)
  #Identity matrix  
  I <- diag(n)
  #Centering matrix
  One <- matrix(1, n, n)
  Cmat <- I - (One/n)
  #Gower matrix
  G = Cmat %*% Amat %*% Cmat
  #model matrix
  df <- data.frame(a = factor(sites))
  M <- model.matrix(~a,df, contrasts = list(a = "contr.helmert"))
  if(length(exclude)>0) M[exclude,] <- 0
  #Projection matrix
  H <- M%*%MASS::ginv(t(M)%*%M)%*%t(M)
  if(length(exclude)>0) {
    non_exclude <- (1:length(sites))[-exclude]
    #Copy projection values from non-excluded site of the trajectory that the external site belongs to
    for(i in 1:length(exclude)) {
      s <- sites[exclude[i]]
      copy_from <- which((sites==s)&((1:n)%in%non_exclude))[1]
      H[,exclude[i]] <- H[,copy_from] # Copies values for centroid removal
    }
  }
  #Residual G matrix (when there are no excluded sites, the H matrix is symmetrical)
  R <- (I-t(H))%*%G%*%(I-H)
  #Backtransform to distances
  dcent<-matrix(0,n,n)
  for(i in 1:n) {
    for(j in i:n) {
      dsq <- (R[i,i]-2*R[i,j]+R[j,j])
      if(dsq > 0) {
        dcent[i,j] = sqrt(dsq) #truncate negative squared distances
        dcent[j,i] = dcent[i,j]
      }
    }
  }
  # Replace distance matrix by dcent
  x$d <- as.dist(dcent)
  return(x)
}
