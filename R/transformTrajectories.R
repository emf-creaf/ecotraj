.gowerCentered<-function(D) {
  Dmat <- as.matrix(D)
  Amat <- (-0.5)*(Dmat^2)
  n <- nrow(Dmat)
  #Identity matrix  
  I <- diag(n)
  #Centering matrix
  One <- matrix(1, n, n)
  Cmat <- I - (One/n)
  #Gower matrix
  G = Cmat %*% Amat %*% Cmat
  return(G)
}
#' Transform trajectories
#' 
#' The following functions are provided to transform trajectories:
#' \itemize{
#' \item{Function \code{smoothTrajectories} performs multivariate smoothing on trajectory data using a Gaussian kernel.}
#' \item{Function \code{centerTrajectories} shifts all trajectories to the center of the multivariate space and returns a modified distance matrix.}
#' \item{Function \code{interpolateTrajectories} relocates trajectory ecological states to those corresponding to input times, via interpolation.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name transformTrajectories
#' @aliases centerTrajectories smoothTrajectories interpolateTrajectories
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' 
#' @details 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (2019). 
#' Function \code{centerTrajectories} performs centering of trajectories using matrix algebra as explained in Anderson (2017).
#'
#' @return 
#' A modified object of class \code{\link{trajectories}}, where distance matrix has been transformed. When calling \code{interpolateTrajectories}, also the
#' number of observations and metadata is likely to be affected.
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
  surveys <- x$metadata$surveys
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] <- sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  
  n <- length(sites)
  umat_smooth <- matrix(0, n, n)
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
  x$d <- .distanceBetweenClusters(as.matrix(d), umat_smooth)
  return(x)
}

#' @rdname transformTrajectories
#' @param exclude An integer vector indicating sites that are excluded from trajectory centroid computation. Note: for objects of class \code{\link{cycles}}, \code{external} are excluded by default.
#' @export
centerTrajectories<-function(x, exclude = integer(0)) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  surveys <- x$metadata$surveys
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
    exclude <- unique(c(exclude,which(x$metadata$internal==FALSE)))
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
    exclude <- unique(c(exclude,which(x$metadata$internal==FALSE)))
  } else {
    sites <- x$metadata$sites
  }
  
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
  n <- nrow(Dmat)
  I <- diag(n)
  
  # Anderson (2017). Permutational Multivariate Analysis of Variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article ID: stat07841.
  G <- .gowerCentered(d)
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


#' @rdname transformTrajectories
#' @param times A numeric vector indicating new observation times for trajectories. Values should be comprised between time limits of the original trajectories.
#' @export
interpolateTrajectories<-function(x, times) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  if(!inherits(times, "numeric") && !inherits(times, "integer")) stop("'times' should be a numeric vector")
  if(length(unique(times)) < length(times)) stop("Time values must be unique")
  if(length(times)<0) stop("'times' should be of length > 1")
  times <- sort(as.numeric(times))
    
  dmat <- as.matrix(x$d)
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }  
  times_old <- x$metadata$times
  
  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] <- sum(sites==siteIDs[i])
  if(sum(nsurveysite<2)>0) stop("All sites need to be surveyed at least twice")
  t_start <- tapply(x$metadata$times,sites,min)[siteIDs]
  t_end <- tapply(x$metadata$times,sites,max)[siteIDs]
  if(any(t_start > times[1])) stop("Initial time beyond trajectory common temporal range")
  if(any(t_end < times[length(times)])) stop("Final time beyond trajectory common temporal range")
  
  times_final <- times
  times <- sort(unique(c(times, times_old)))
  survey_selection <- which(times %in% times_final)
  
  n_new <- nsite*length(times)
  sites_new <- as.character(gl(nsite, length(times), labels = siteIDs))
  times_new <- rep(times, nsite)
  dmat_new <- matrix(NA, nrow = n_new, ncol = n_new)
  rownames(dmat_new) <- 1:n_new
  colnames(dmat_new) <- 1:n_new
  is_old <- rep(FALSE, n_new)
  for(i in 1:n_new) {
    site_i <- sites_new[i]
    time_i <- times_new[i]
    if(time_i %in% times_old[sites==site_i]) is_old[i] <- TRUE
    dmat_new[i,i] <- 0
  }
  # Distances between existing points
  for(i in which(is_old)) {
    site_i <- sites_new[i]
    time_i <- times_new[i]
    r_i <- which((times_old ==time_i) & (sites == site_i))
    for(j in which(is_old)) {
      site_j <- sites_new[j]
      time_j <- times_new[j]
      c_j <- which((times_old ==time_j) & (sites == site_j))
      dmat_new[i,j] <- dmat[r_i, c_j]
    }
  }
  # Distances between new points and existing ones
  for(i in which(!is_old)) {
    site_i <- sites_new[i]
    time_i <- times_new[i]
    bef_i <- max(which(times_old[sites == site_i] < time_i))
    aft_i <- min(which(times_old[sites == site_i] > time_i))
    time_bef_i <- times_old[sites == site_i][bef_i]
    time_aft_i <- times_old[sites == site_i][aft_i]
    p <- (time_aft_i - time_i)/(time_aft_i - time_bef_i)
    r_bef_i <- which((times_old ==time_bef_i) & (sites == site_i))
    r_aft_i <- which((times_old ==time_aft_i) & (sites == site_i))
    dref <- dmat[r_bef_i, r_aft_i]
    for(j in which(is_old)) {
      site_j <- sites_new[j]
      time_j <- times_new[j]
      c_j <- which((times_old ==time_j) & (sites == site_j))
      d1 <- dmat[r_bef_i, c_j]
      d2 <- dmat[r_aft_i, c_j]
      # print(c(time_i, bef_i, aft_i, r_bef_i, r_aft_i, i, j, dref, d1, d2, p))
      d_int <- .distanceToInterpolatedC(dref, d1, d2, p)
      dmat_new[i,j] <- d_int
      dmat_new[j,i] <- d_int
    }
  }
  
  # Distances between new points
  for(i in which(!is_old)) {
    site_i <- sites_new[i]
    time_i <- times_new[i]
    bef_i <- max(which(times_old[sites == site_i] < time_i))
    aft_i <- min(which(times_old[sites == site_i] > time_i))
    time_bef_i <- times_old[sites == site_i][bef_i]
    time_aft_i <- times_old[sites == site_i][aft_i]
    p <- (time_aft_i - time_i)/(time_aft_i - time_bef_i)
    r_bef_i <- which((times_old ==time_bef_i) & (sites == site_i))
    r_aft_i <- which((times_old ==time_aft_i) & (sites == site_i))
    dref <- dmat[r_bef_i, r_aft_i]
    r_bef_new_i <- which((times_new ==time_bef_i) & (sites_new == site_i))
    r_aft_new_i <- which((times_new ==time_aft_i) & (sites_new == site_i))
    other_new <- which(!is_old)
    other_new <- other_new[other_new>i]
    for(j in other_new) {
      site_j <- sites_new[j]
      time_j <- times_new[j]
      c_j_new <- which((times_new ==time_j) & (sites_new == site_j))
      d1 <- dmat_new[r_bef_new_i, c_j_new]
      d2 <- dmat_new[r_aft_new_i, c_j_new]
      # print(c(time_i, bef_i, aft_i, r_bef_i, r_aft_i, i, j, dref, d1, d2, p))
      d_int <- .distanceToInterpolatedC(dref, d1, d2, p)
      dmat_new[i,j] <- d_int
      dmat_new[j,i] <- d_int
    }
  }
  sel <- times_new %in% times_final
  return(defineTrajectories(as.dist(dmat_new[sel,sel]), sites = sites_new[sel], times = times_new[sel]))
}