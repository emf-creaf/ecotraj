#' Trajectory comparison
#'
#' Functions to compare pairs of trajectories or trajectory segments.
#' \itemize{
#' \item{Function \code{segmentDistances} calculates the distance between pairs of trajectory segments.}
#' \item{Function \code{trajectoryDistances} calculates the distance between pairs of trajectories.}
#' \item{Function \code{trajectoryConvergence} performs the Mann-Kendall trend test on the distances between trajectories (symmetric test) or the distance between points of one trajectory to the other.}
#' }
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param distance.type The type of distance index to be calculated (see section Details).
#' @param symmetrization Function used to obtain a symmetric distance, so that DSPD(T1,T2) = DSPD(T2,T1) (e.g., \code{mean} or \code{min}). If \code{symmetrization = NULL} then the symmetrization is not conducted and the output dissimilarity matrix is not symmetric. 
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#'
#' @details 
#' Ecological Trajectory Analysis (ETA) is a framework to analyze dynamics of ecosystems described as trajectories in a chosen space of multivariate resemblance (De \enc{Cáceres}{Caceres} et al. 2019).
#' ETA takes trajectories as objects to be analyzed and compared geometrically. 
#' 
#' The input distance matrix \code{d} should ideally be metric. That is, all subsets of distance triplets should fulfill the triangle inequality (see utility function \code{\link{is.metric}}). 
#' All ETA functions that require metricity include a parameter '\code{add}', which by default is TRUE, meaning that whenever the triangle inequality is broken the minimum constant required to fulfill it is added to the three distances.
#' If such local (an hence, inconsistent across triplets) corrections are not desired, users should find another way modify \code{d} to achieve metricity, such as PCoA, metric MDS or non-metric MDS (see vignette 'Introduction to Ecological Trajectory Analysis'). 
#' If parameter '\code{add}' is set to FALSE and problems of triangle inequality exist, ETA functions may provide missing values in some cases where they should not.
#' 
#' The resemblance between trajectories is done by adapting concepts and procedures used for the analysis of trajectories in space (i.e. movement data) (Besse et al. 2016).   
#' 
#' Parameter \code{distance.type} is the type of distance index to be calculated which for function \code{segmentDistances} has the following options (Besse et al. 2016, De \enc{Cáceres}{Caceres} et al. 2019:
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two segments.}
#'     \item{\code{directed-segment}: Directed segment distance (default).}
#'     \item{\code{PPA}: Perpendicular-parallel-angle distance.}
#'   }
#'  In the case of function \code{trajectoryDistances} the following values are possible (De \enc{Cáceres}{Caceres} et al. 2019):
#'  \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two trajectories.}
#'     \item{\code{SPD}: Segment path distance.}
#'     \item{\code{DSPD}: Directed segment path distance (default).}
#'  }
#'  
#' @returns 
#' Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories (if \code{symmetrization = NULL} then the object returned is of class \code{matrix}). 
#' 
#' Function \code{segmentDistances} list with the following elements:
#' \itemize{
#'   \item{\code{Dseg}: Distance matrix between segments.}
#'   \item{\code{Dini}: Distance matrix between initial points of segments.}
#'   \item{\code{Dfin}: Distance matrix between final points of segments.}
#'   \item{\code{Dinifin}: Distance matrix between initial points of one segment and the final point of the other.}
#'   \item{\code{Dfinini}: Distance matrix between final points of one segment and the initial point of the other.}
#' }
#' Function \code{trajectoryConvergence} returns a list with two elements:
#' \itemize{
#'   \item{\code{tau}: A matrix with the statistic (Mann-Kendall's tau) of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the statistic of the test of the row trajectory approaching the column trajectory.}
#'   \item{\code{p.value}: A matrix with the p-value of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the p-value indicates the test of the row trajectory approaching the column trajectory.}
#' }
#' 
#' @references
#' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
#' 
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). 
#' Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
#' @seealso \code{\link{trajectoryMetrics}}, \code{\link{trajectoryPlot}}, \code{\link{transformTrajectories}}, \code{\link{trajectoryProjection}}
#' 
#' @examples 
#' #Description of sites and surveys
#' sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
#' surveys <- c(1,2,3,4,1,2,3,4,1,2,3,4)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=12, ncol=2)
#' xy[2,2]<-1
#' xy[3,2]<-2
#' xy[4,2]<-3
#' xy[5:6,2] <- xy[1:2,2]
#' xy[7,2]<-1.5
#' xy[8,2]<-2.0
#' xy[5:6,1] <- 0.25
#' xy[7,1]<-0.5
#' xy[8,1]<-1.0
#' xy[9:10,1] <- xy[5:6,1]+0.25
#' xy[11,1] <- 1.0
#' xy[12,1] <-1.5
#' xy[9:10,2] <- xy[5:6,2]
#' xy[11:12,2]<-c(1.25,1.0)
#'   
#' #Draw trajectories
#' trajectoryPlot(xy, sites, surveys,  
#'                traj.colors = c("black","red", "blue"), lwd = 2)
#' 
#' #Distance matrix
#' d <- dist(xy)
#' d
#'   
#' #Trajectory data
#' x <- defineTrajectories(d, sites, surveys)
#' 
#' #Distances between trajectory segments
#' segmentDistances(x, distance.type = "Hausdorff")
#' segmentDistances(x, distance.type = "directed-segment")
#' 
#' #Distances between trajectories
#' trajectoryDistances(x, distance.type = "Hausdorff")
#' trajectoryDistances(x, distance.type = "DSPD")
#'   
#' #Trajectory convergence/divergence
#' trajectoryConvergence(x)
#' 
#' @name trajectoryComparison
#' @export
segmentDistances<-function(x, distance.type ="directed-segment", add = TRUE) {
  distance.type <- match.arg(distance.type, c("directed-segment", "Hausdorff", "PPA"))
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
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat <- as.matrix(d)
  n <- nrow(dmat)
  nseg <- sum(nsurveysite)-nsite
  segnames <- character(nseg)
  cnt <- 1
  for(i in 1:nsite) {
    if(!is.null(surveys)) {
      surv <- surveys[sites==siteIDs[i]]
      surv <- sort(surv) #Surveys may not be in order
    }
    else surv <- 1:nsurveysite[i]
    for(j in 1:(nsurveysite[i]-1)) {
      segnames[cnt] <- paste0(siteIDs[i],"[",surv[j],"-",surv[j+1],"]")
      cnt <- cnt+1
    }
  }
  dsegmat <- matrix(0, nseg, nseg)
  rownames(dsegmat) <-segnames
  colnames(dsegmat) <-segnames
  dinisegmat <- dsegmat
  dfinsegmat <- dsegmat
  dinifinsegmat <- dsegmat
  
  os1 <- 1
  for(i1 in 1:nsite) {
    ind_surv1 <- which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 <- ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-1)) {
      os2 <- 1
      for(i2 in 1:nsite) {
        ind_surv2 <- which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 <- ind_surv2[order(surveys[sites==siteIDs[i2]])]
        for(s2 in 1:(nsurveysite[i2]-1)) {
          #Select submatrix from dmat
          ind12 <- c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])
          # print(ind12)
          # print(c(os1, os2))
          dmat12 <- dmat[c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1]),
                         c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])]
          dsegmat[os1,os2] <- .twoSegmentDistanceC(dmat12, type=distance.type, add)
          dsegmat[os2,os1] <- dsegmat[os1,os2]
          dinisegmat[os2,os1] <- dinisegmat[os1,os2]<-dmat[ind_surv1[s1],ind_surv2[s2]]
          dfinsegmat[os2,os1] <- dfinsegmat[os1,os2]<-dmat[ind_surv1[s1+1],ind_surv2[s2+1]]
          dinifinsegmat[os1,os2] <- dmat[ind_surv1[s1],ind_surv2[s2+1]]
          dinifinsegmat[os2,os1] <- dmat[ind_surv1[s1+1],ind_surv2[s2]]
          os2 <- os2+1
        }
      }
      os1 <- os1+1
    }
  }
  
  return(list(Dseg = as.dist(dsegmat), Dini=as.dist(dinisegmat), Dfin = as.dist(dfinsegmat),
              Dinifin=dinifinsegmat))
}

#' @rdname trajectoryComparison
#' @export
trajectoryDistances<-function(x, distance.type="DSPD", symmetrization = "mean" , add=TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  distance.type <- match.arg(distance.type, c("DSPD", "SPD", "Hausdorff"))
  
  d <- x$d
  surveys <- x$metadata$surveys
  # This allows treating fixed date trajectories as sites for plotting purposes
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
  
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  n <- nrow(as.matrix(d))
  nseg <- sum(nsurveysite)-nsite
  
  #Init output
  dtraj <- matrix(0, nrow=nsite, ncol = nsite)
  rownames(dtraj) <- siteIDs
  colnames(dtraj) <- siteIDs
  if(distance.type=="DSPD"){
    lsd <- segmentDistances(x, distance.type="directed-segment", add)
    dsegmat <- as.matrix(lsd$Dseg)
    for(i1 in 1:nsite) {
      for(i2 in 1:nsite) {
        dt12 = 0
        for(s1 in 1:(nsurveysite[i1]-1)) {
          dt12ivec = numeric(0)
          iseg1 = sum(nsurveysite[1:i1]-1)-(nsurveysite[i1]-1)+s1
          for(s2 in 1:(nsurveysite[i2]-1)) {
            iseg2 = sum(nsurveysite[1:i2]-1)-(nsurveysite[i2]-1)+s2
            dt12ivec = c(dt12ivec, dsegmat[iseg1, iseg2])
          }
          dt12 = dt12 + min(dt12ivec)
        }
        dt12 = dt12/(nsurveysite[i1]-1) #Average of distances between segments of T1 and trajectory T2
        dt21 = 0 
        for(s2 in 1:(nsurveysite[i2]-1)) {
          dt21ivec = numeric(0)
          iseg2 = sum(nsurveysite[1:i2]-1)-(nsurveysite[i2]-1)+s2
          for(s1 in 1:(nsurveysite[i1]-1)) {
            iseg1 = sum(nsurveysite[1:i1]-1)-(nsurveysite[i1]-1)+s1
            dt21ivec = c(dt21ivec, dsegmat[iseg1, iseg2])
          }
          dt21 = dt21 + min(dt21ivec)
        }
        dt21 = dt21/(nsurveysite[i2]-1) #Average of distances between segments of T2 and trajectory T1
        
        if(!is.null(symmetrization)) {
          dtraj[i1,i2] = do.call(symmetrization, list(c(dt12,dt21))) #Symmetrization
          dtraj[i2,i1] = dtraj[i1,i2]
        } else {
          dtraj[i1,i2] = dt12
          dtraj[i2,i1] = dt21
        }
      }
    }
    
  } 
  else if(distance.type=="SPD") {
    dmat = as.matrix(d)
    for(i1 in 1:nsite) {
      ind_surv1 = which(sites==siteIDs[i1])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        dt12 = 0
        for(p1 in 1:nsurveysite[i1]) {
          dt12ivec = numeric(0)
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12ivec = c(dt12ivec, .distanceToSegmentC(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2], add)[3])
          }
          dt12 = dt12 + min(dt12ivec)
        }
        dt12 = dt12/nsurveysite[i1] #Average of distances between points of T1 and trajectory T2
        dt21 = 0 
        for(p2 in 1:nsurveysite[i2]) {
          dt21ivec = numeric(0)
          ip2 = ind_surv2[p2]
          for(s1 in 1:(nsurveysite[i1]-1)) {
            ipi1 = ind_surv1[s1] #initial point
            ipe1 = ind_surv1[s1+1] #end point
            dt21ivec = c(dt21ivec, .distanceToSegmentC(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1], add)[3])
          }
          dt21 = dt21 + min(dt21ivec)
        }
        dt21 = dt21/nsurveysite[i2] #Average of distances between points of T2 and trajectory T1
        
        if(!is.null(symmetrization)) {
          dtraj[i1,i2] = (dt12+dt21)/2 #Symmetrization
          dtraj[i2,i1] = dtraj[i1,i2]
        } else {
          dtraj[i1,i2] = dt12
          dtraj[i2,i1] = dt21
        }
      }
    }
  }
  else if(distance.type=="Hausdorff") {
    dmat = as.matrix(d)
    for(i1 in 1:nsite) {
      ind_surv1 = which(sites==siteIDs[i1])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        dt12 = 0
        dt12vec = numeric(0)
        for(p1 in 1:nsurveysite[i1]) {
          ip1 = ind_surv1[p1]
          for(s2 in 1:(nsurveysite[i2]-1)) {
            ipi2 = ind_surv2[s2] #initial point
            ipe2 = ind_surv2[s2+1] #end point
            dt12vec = c(dt12vec, .distanceToSegmentC(dmat[ipi2,ipe2], dmat[ip1, ipi2], dmat[ip1,ipe2], add)[3])
          }
        }
        dt12 = max(dt12vec) #Maximum of distances between points of T1 and segments of T2
        dt21 = 0 
        dt21vec = numeric(0)
        for(p2 in 1:nsurveysite[i2]) {
          ip2 = ind_surv2[p2]
          for(s1 in 1:(nsurveysite[i1]-1)) {
            ipi1 = ind_surv1[s1] #initial point
            ipe1 = ind_surv1[s1+1] #end point
            dt21vec = c(dt21vec, .distanceToSegmentC(dmat[ipi1,ipe1], dmat[ip2, ipi1], dmat[ip2,ipe1], add)[3])
          }
        }
        dt21 = max(dt21vec) #Maximum of distances between points of T2 and segments of T1
        
        dtraj[i1,i2] = max(dt12, dt21) #maximum of maximums
        dtraj[i2,i1] = dtraj[i1,i2]
      }
    }
  } 
  else stop("Wrong distance type")
  if(!is.null(symmetrization)) return(as.dist(dtraj))
  return(dtraj)
}

#' @rdname trajectoryComparison
#' @param symmetric A logical flag to indicate a symmetric convergence comparison of trajectories.
#' @export
trajectoryConvergence<-function(x, symmetric = FALSE, add=TRUE){
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  surveys <- x$metadata$surveys
  # This allows treating fixed date trajectories as sites for plotting purposes
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
  n <- nrow(as.matrix(d))
  
  #Init output
  tau <- matrix(NA, nrow=nsite, ncol = nsite)
  rownames(tau) <- siteIDs
  colnames(tau) <- siteIDs
  p.value <- tau
  dmat <- as.matrix(d)
  for(i1 in 1:(nsite-1)) {
    ind_surv1 <- which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 <- ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(i2 in (i1+1):nsite) {
      ind_surv2 <- which(sites==siteIDs[i2])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv2 <- ind_surv2[order(surveys[sites==siteIDs[i2]])]
      if(!symmetric) {
        trajectory <- ind_surv2
        target <- ind_surv1
        trajProj <- trajectoryProjection(d,target, trajectory, add=add)
        dT <- trajProj$distanceToTrajectory
        mk.test <- MannKendall(dT)
        tau[i1,i2] <- mk.test$tau
        p.value[i1,i2] <- mk.test$sl
        trajectory <- ind_surv1
        target <- ind_surv2
        trajProj <- trajectoryProjection(d,target, trajectory, add=add)
        dT <- trajProj$distanceToTrajectory
        mk.test <- MannKendall(dT)
        tau[i2,i1] <- mk.test$tau
        p.value[i2,i1] <- mk.test$sl
      } 
      else {
        if(length(ind_surv1)==length(ind_surv2)) {
          dT <- numeric(length(ind_surv1))
          for(j in 1:length(ind_surv1)) dT[j] = dmat[ind_surv1[j], ind_surv2[j]]
          mk.test <- MannKendall(dT)
          tau[i1,i2] <- mk.test$tau
          p.value[i1,i2] <- mk.test$sl
          tau[i2,i1] <- mk.test$tau
          p.value[i2,i1] <- mk.test$sl
        } else {
          warning(paste0("sites ",i1, " and ",i2," do not have the same number of surveys."))
        }
      }
    }
  }
  return(list(tau = tau, p.value = p.value))
}
