#' Metrics for Ecological Trajectory Analysis
#' 
#' Ecological Trajectory Analysis (ETA) is a framework to analyze dynamics of ecosystems described as trajectories in a chosen space of multivariate resemblance (De \enc{Cáceres}{Caceres} et al. 2019).
#' ETA takes trajectories as objects to be analyzed and compared geometrically. 
#' 
#' Given a distance matrix between ecosystem states, the set of functions that provide ETA metrics are:
#' \itemize{
#' \item{Functions \code{segmentDistances} and \code{trajectoryDistances} calculate the distance between pairs of directed segments and ecosystem trajectories, respectively.}
#' \item{Function \code{trajectoryLengths} calculates lengths of directed segments and total path lengths of trajectories.}
#' \item{Function \code{trajectoryLengths2D} calculates lengths of directed segments and total path lengths of trajectories from 2D coordinates given as input.} 
#' \item{Function \code{trajectoryAngles} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryAngles2D} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryProjection} performs an orthogonal projection of a set of target points onto a specified trajectory and returns the distance to the trajectory (i.e. rejection) and the relative position of the projection point within the trajectory.}
#' \item{Function \code{trajectoryConvergence} performs the Mann-Kendall trend test on the distances between trajectories (symmetric test) or the distance between points of one trajectory to the other.}
#' \item{Function \code{trajectoryDirectionality} returns (for each trajectory) a statistic that measures directionality of the whole trajectory.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name trajectorymetrics
#' @aliases segmentDistances trajectoryDistances trajectoryLengths trajectoryLengths2D trajectoryAngles trajectoryAngles2D
#'          trajectoryProjection trajectoryConvergence trajectoryDirectionality 
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states (see details).
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param surveys A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @param distance.type 
#' The type of distance index to be calculated (Besse et al. 2016; De Cáceres et al. submitted). For \code{segmentDistances} the available indices are:
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two segments.}
#'     \item{\code{directed-segment}: Directed segment distance (default).}
#'     \item{\code{PPA}: Perpendicular-parallel-angle distance.}
#'   }
#' whereas for \code{trajectoryDistances} the available indices are:
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two trajectories.}
#'     \item{\code{SPD}: Segment path distance.}
#'     \item{\code{DSPD}: Directed segment path distance (default).}
#'   }
#' @param symmetrization Function used to obtain a symmetric distance, so that DSPD(T1,T2) = DSPD(T2,T1) (e.g., \code{mean} or \code{min}). If \code{symmetrization = NULL} then the symmetrization is not conducted and the output dissimilarity matrix is not symmetric. 
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#' @param verbose Provides console output informing about process (useful for large dataset).
#' 
#' @details 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (2019). 
#' The input distance matrix \code{d} should ideally be metric. That is, all subsets of distance triplets should fulfill the triangle inequality (see utility function \code{\link{is.metric}}). 
#' All ETA functions that require metricity include a parameter '\code{add}', which by default is TRUE, meaning that whenever the triangle inequality is broken the minimum constant required to fulfill it is added to the three distances.
#' If such local (an hence, inconsistent across triplets) corrections are not desired, users should find another way modify \code{d} to achieve metricity, such as PCoA, metric MDS or non-metric MDS (see vignette 'Introduction to Ecological Trajectory Analysis'). 
#' If parameter '\code{add}' is set to FALSE and problems of triangle inequality exist, ETA functions may provide missing values in some cases where they should not.
#' 
#' The resemblance between trajectories is done by adapting concepts and procedures used for the analysis of trajectories in space (i.e. movement data) (Besse et al. 2016).   
#' 
#' Function \code{trajectoryAngles} calculates angles between consecutive segments in degrees. For each pair of segments, the angle between the two is defined on the plane that contains the two segments, and measures the change in direction (in degrees) from one segment to the other. 
#' Angles are always positive, with zero values indicating segments that are in a straight line, and values equal to 180 degrees for segments that are in opposite directions. If \code{all = TRUE}
#' angles are calculated between the segments corresponding to all ordered triplets. Alternatively, if \code{relativeToInitial = TRUE} angles are calculated for each segment with respect to the initial survey.
#' 
#' Function \code{trajectoryAngles2D} calculates angles between consecutive segments in degrees from 2D coordinates given as input. For each pair of segments, the angle between the two is defined on the plane that contains the two segments, and measures the change in direction (in degrees) from one segment to the other. 
#' Angles are always positive (O to 360), with zero values indicating segments that are in a straight line, and values equal to 180 degrees for segments that are in opposite directions. 
#' If \code{all = TRUE} angles are calculated between the segments corresponding to all ordered triplets. Alternatively, if \code{relativeToInitial = TRUE} angles are calculated for each segment with respect to the initial survey.
#' If \code{betweenSegments = TRUE} angles are calculated between segments of trajectory, otherwise, If \code{betweenSegments = FALSE}, angles are calculated considering Y axis as the North (0°).
#' 
#' @return Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories (if \code{symmetrization = NULL} then the object returned is of class \code{matrix}). 
#' 
#' Function \code{trajectorySegments} returns a list with the following elements:
#' \itemize{
#'   \item{\code{Dseg}: Distance matrix between segments.}
#'   \item{\code{Dini}: Distance matrix between initial points of segments.}
#'   \item{\code{Dfin}: Distance matrix between final points of segments.}
#'   \item{\code{Dinifin}: Distance matrix between initial points of one segment and the final point of the other.}
#'   \item{\code{Dfinini}: Distance matrix between final points of one segment and the initial point of the other.}
#' }
#' 
#' Function \code{trajectoryLengths} returns a data frame with the length of each segment on each trajectory and the total length of all trajectories. 
#' If \code{relativeToInitial = TRUE} lengths are calculated between the initial survey and all the other surveys.
#' If \code{all = TRUE} lengths are calculated for all segments.
#' 
#' Function \code{trajectoryLengths2D} returns a data frame with the length of each segment on each trajectory and the total length of all trajectories. 
#' If \code{relativeToInitial = TRUE} lengths are calculated between the initial survey and all the other surveys.
#' If \code{all = TRUE} lengths are calculated for all segments.
#' 
#' Function \code{trajectoryAngles} returns a data frame with angle values on each trajectory. If \code{stats=TRUE}, then the mean, standard deviation and mean resultant length of those angles are also returned. 
#' 
#' Function \code{trajectoryAngles2D} returns a data frame with angle values on each trajectory. If \code{betweenSegments=TRUE}, then angles are calculated between trajectory segments, alternatively, If \code{betweenSegments=FALSE}, angles are calculated considering Y axis as the North (0°).
#' 
#' Function \code{trajectoryProjection} returns a data frame with the following columns:
#' \itemize{
#'   \item{\code{distanceToTrajectory}: Distances to the trajectory, i.e. rejection (\code{NA} for target points whose projection is outside the trajectory).}
#'   \item{\code{segment}: Segment that includes the projected point (\code{NA} for target points whose projection is outside the trajectory).}
#'   \item{\code{relativeSegmentPosition}: Relative position of the projected point within the segment, i.e. values from 0 to 1 with 0 representing the start of the segment and 1 representing its end (\code{NA} for target points whose projection is outside the trajectory).}
#'   \item{\code{relativeTrajectoryPosition}: Relative position of the projected point within the trajectory, i.e. values from 0 to 1 with 0 representing the start of the trajectory and 1 representing its end (\code{NA} for target points whose projection is outside the trajectory).}
#' }
#' 
#' Function \code{trajectoryConvergence} returns a list with two elements:
#' \itemize{
#'   \item{\code{tau}: A matrix with the statistic (Mann-Kendall's tau) of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the statistic of the test of the row trajectory approaching the column trajectory.}
#'   \item{\code{p.value}: A matrix with the p-value of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the p-value indicates the test of the row trajectory approaching the column trajectory.}
#' }
#' 
#' Function \code{trajectoryDirectionality} returns a vector with directionality values (one per trajectory).
#' 
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @author Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la Baie de Saint-Brieuc
#' 
#' @references
#' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
#' 
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). 
#' Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
#' @seealso \code{\link{trajectoryplots}}, \code{\link{trajectoryutils}} 
#' 
#' @examples 
#' #Description of sites and surveys
#' sites = c(1,1,1,2,2,2)
#' surveys=c(1,2,3,1,2,3)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=6, ncol=2)
#' xy[2,2]<-1
#' xy[3,2]<-2
#' xy[4:6,1] <- 0.5
#' xy[4:6,2] <- xy[1:3,2]
#' xy[6,1]<-1
#'   
#' #Draw trajectories
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)
#' 
#' #Distance matrix
#' d = dist(xy)
#' d
#'   
#' trajectoryLengths(d, sites, surveys)
#' trajectoryLengths2D(xy, sites, surveys)
#' trajectoryAngles(d, sites, surveys)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = TRUE)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = FALSE)
#' segmentDistances(d, sites, surveys)$Dseg
#' trajectoryDistances(d, sites, surveys, distance.type = "Hausdorff")
#' trajectoryDistances(d, sites, surveys, distance.type = "DSPD")
#'   
#'   
#' #Should give the same results if surveys are not in order 
#' #(here we switch surveys for site 2)
#' temp = xy[5,]
#' xy[5,] = xy[6,]
#' xy[6,] = temp
#' surveys[5] = 3
#' surveys[6] = 2
#'   
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)   
#' trajectoryLengths(dist(xy), sites, surveys)
#' trajectoryLengths2D(xy, sites, surveys)
#' segmentDistances(dist(xy), sites, surveys)$Dseg
#' trajectoryAngles(dist(xy), sites, surveys)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = TRUE)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = FALSE)
#' trajectoryDistances(dist(xy), sites, surveys, distance.type = "Hausdorff")
#' trajectoryDistances(dist(xy), sites, surveys, distance.type = "DSPD")
#'  
#' @export
segmentDistances<-function(d, sites, surveys=NULL, distance.type ="directed-segment", add = TRUE, verbose=FALSE) {
  distance.type <- match.arg(distance.type, c("directed-segment", "Hausdorff", "PPA"))
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat = as.matrix(d)
  n = nrow(dmat)
  nseg = sum(nsurveysite)-nsite
  segnames = character(nseg)
  cnt=1
  for(i in 1:nsite) {
    if(!is.null(surveys)) {
      surv = surveys[sites==siteIDs[i]]
      surv = sort(surv) #Surveys may not be in order
    }
    else surv = 1:nsurveysite[i]
    for(j in 1:(nsurveysite[i]-1)) {
      segnames[cnt] = paste0(siteIDs[i],"[",surv[j],"-",surv[j+1],"]")
      cnt = cnt+1
    }
  }
  dsegmat = matrix(0, nseg, nseg)
  rownames(dsegmat) =segnames
  colnames(dsegmat) =segnames
  dinisegmat = dsegmat
  dfinsegmat = dsegmat
  dinifinsegmat = dsegmat

  os1 = 1
  if(verbose) {
    cat("\nCalculating segment distances...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-1)) {
      os2 = 1
      for(i2 in 1:nsite) {
        ind_surv2 = which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
        for(s2 in 1:(nsurveysite[i2]-1)) {
          # os2 = sum(nsurveysite[1:(i2-1)]-1)+s2 #Output index of segment 2
          #Select submatrix from dmat
          ind12 = c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])
          # print(ind12)
          # print(c(os1, os2))
          dmat12 = dmat[c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1]),
                        c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])]
          dsegmat[os1,os2] <- .twoSegmentDistanceC(dmat12, type=distance.type, add)
          dsegmat[os2,os1] <- dsegmat[os1,os2]
          dinisegmat[os2,os1] <- dinisegmat[os1,os2]<-dmat[ind_surv1[s1],ind_surv2[s2]]
          dfinsegmat[os2,os1] <- dfinsegmat[os1,os2]<-dmat[ind_surv1[s1+1],ind_surv2[s2+1]]
          dinifinsegmat[os1,os2] <- dmat[ind_surv1[s1],ind_surv2[s2+1]]
          dinifinsegmat[os2,os1] <- dmat[ind_surv1[s1+1],ind_surv2[s2]]
          os2 = os2+1
        }
      }
      os1 = os1+1
    }
  }

  return(list(Dseg = as.dist(dsegmat), Dini=as.dist(dinisegmat), Dfin = as.dist(dfinsegmat),
              Dinifin=dinifinsegmat))
}

#' @rdname trajectorymetrics
#' @export
trajectoryDistances<-function(d, sites, surveys=NULL, distance.type="DSPD", symmetrization = "mean" , add=TRUE, verbose=FALSE) {
  distance.type <- match.arg(distance.type, c("DSPD", "SPD", "Hausdorff"))
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  n = nrow(as.matrix(d))
  nseg = sum(nsurveysite)-nsite
  
  #Init output
  dtraj = matrix(0, nrow=nsite, ncol = nsite)
  rownames(dtraj) = siteIDs
  colnames(dtraj) = siteIDs
  if(distance.type=="DSPD"){
    lsd = segmentDistances(d,sites, surveys,distance.type="directed-segment", add, verbose)
    dsegmat = as.matrix(lsd$Dseg)
    if(verbose) {
      cat("\nCalculating trajectory distances...\n")
      tb = txtProgressBar(1, nsite, style=3)
    }
    for(i1 in 1:nsite) {
      if(verbose) setTxtProgressBar(tb, i1)
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

#' @rdname trajectorymetrics
#' @param relativeToInitial Flag to indicate that lengths or angles should be calculated with respect to initial survey.
#' @param all Flag to indicate that lengths or angles are desired for all segments or for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, length or angles are calculated according to relativeToInitial flag.
#' @export
trajectoryLengths<-function(d, sites, surveys=NULL, relativeToInitial = FALSE, all=FALSE, verbose= FALSE) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  surveyIDs<-unique(surveys)
  nsurvey<-length(surveyIDs)
  
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat = as.matrix(d)
  n = nrow(dmat)

  maxnsurveys = max(nsurveysite)

  if(!all) {
    lengths = as.data.frame(matrix(NA, nrow=nsite, ncol=maxnsurveys))
    row.names(lengths)<-siteIDs
    if(relativeToInitial) names(lengths)<-c(paste0("Lt1_t",as.character(2:(maxnsurveys))),"Trajectory")
    else names(lengths)<-c(paste0("S",as.character(1:(maxnsurveys-1))),"Trajectory")
    if(verbose) {
      cat("\nCalculating trajectory lengths...\n")
      tb = txtProgressBar(1, nsite, style=3)
    }
    for(i1 in 1:nsite) {
      if(verbose) setTxtProgressBar(tb, i1)
      ind_surv1 = which(sites==siteIDs[i1])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
      for(s1 in 1:(nsurveysite[i1]-1)) {
        if(relativeToInitial) lengths[i1,s1] = dmat[ind_surv1[1], ind_surv1[s1+1]]
        else  lengths[i1,s1] = dmat[ind_surv1[s1], ind_surv1[s1+1]]
      }
      lengths[i1, maxnsurveys] = sum(lengths[i1,], na.rm=T)
    }
  } else{ 
    
    if(nsite==1) {
    vectord<-as.vector(d)
    lengths<-as.data.frame(matrix(NA, nrow=nsite, ncol=((nsurvey*(nsurvey-1))/2)))
    lengths[1,]<-vectord
    listsurvey<-c(1:nsurvey)#creating column names
    tsurvey<-c()
    for (i in 1: nsurvey){
      tsurvey<-c(tsurvey,paste0("Lt",listsurvey[i]))
    }
    comb<-combn(tsurvey, 2)
    tsurvey<-c(paste0(comb[1,], comb[2,]))
    colnames(lengths)<-c(tsurvey)
    rownames(lengths)<-c(siteIDs)
        
  }else{
 #vector to indicate line for selection in the matrix
    seqline<-c(nsurvey-1, nsurvey, nsurvey)
    for(i in 1:(nsite-1)){
      seqline<-c(seqline, seqline[length(seqline)-2]+nsurvey,seqline[length(seqline)-1]+nsurvey,seqline[length(seqline)-1]+nsurvey)
    }
 #vector to indicate column for selection in the matrix
    seqcolumn<-c(1,1,2)
    for(i in 1:(nsite-1)){
      seqcolumn<-c(seqcolumn, seqcolumn[length(seqcolumn)-2]+nsurvey,seqcolumn[length(seqcolumn)-2]+nsurvey,seqcolumn[length(seqcolumn)]+nsurvey)
    }
#create a vector off all lengths using seqline and seqcolumn to oper selection in the matrix
alllengths<-c()
for(i in 1:(length(seqline))){
  alllengths<-c(alllengths, dmat[seqline[i], seqcolumn[i]])
}

lengths<-as.data.frame(matrix(NA, nrow=nsite, ncol=((nsurvey*(nsurvey-1))/2)))
seqlength<-c(seq(0,length(alllengths),length(alllengths)/nsite))
for(i in 1:nsite){
  lengths[i,]<-alllengths[(seqlength[i]+1):(seqlength[i+1])]
}
listsurvey<-c(1:nsurvey)#creating column names
tsurvey<-c()
for (i in 1: nsurvey){
  tsurvey<-c(tsurvey,paste0("Lt",listsurvey[i]))
}
comb<-combn(tsurvey, 2)
tsurvey<-c(paste0(comb[1,], comb[2,]))
colnames(lengths)<-c(tsurvey)
rownames(lengths)<-c(siteIDs)
  }}

  return(lengths)
}

#' @rdname trajectorymetrics
#' @param xy Matrix with 2D coordinates in a Cartesian space (typically an ordination of ecosystem states).
#' @export
trajectoryLengths2D<-function(xy,sites,surveys, relativeToInitial=FALSE, all=FALSE, verbose = FALSE) {
  
  #order inputs by sites and surveys
xy_temp<-as.data.frame(xy)
xy_temp$sites<-sites
xy_temp$surveys<-surveys
xy_temp<-xy_temp[order(xy_temp$sites,xy_temp$surveys),]
xy<-xy_temp[,1:2]
sites<-c(xy_temp$sites)
surveys<-c(xy_temp$surveys)
  
  siteIDs = unique(sites)
  surveyIDs<-unique(surveys)
  nsite<-length(siteIDs)
  nsurvey<-length(surveyIDs)
  if(nsite!=nrow(xy)/nsurvey) stop("'sites' needs to be of length equal in xy")
  if(nrow(xy)!=nsurvey*nsite) stop("All sites need to be surveyed at least twice")
  #prep all
  D<-dist(xy)
  return(trajectoryLengths(D,sites,surveys,relativeToInitial = relativeToInitial, all = all, verbose = verbose))
}

#' @rdname trajectorymetrics
#' @param all A flag to indicate that angles are desired for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, angles are calculated for consecutive segments only.
#' @param stats A flag to indicate that circular statistics are desired (mean, standard deviation and mean resultant length, i.e. rho)
#' @export
trajectoryAngles<-function(d, sites, surveys=NULL, all = FALSE, relativeToInitial = FALSE, stats = TRUE, add=TRUE, verbose= FALSE) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat = as.matrix(d)
  n = nrow(dmat)
  
  maxnsurveys = max(nsurveysite)
  if(!all) angles = matrix(NA, nrow=nsite, ncol=maxnsurveys+1)
  else {
    angles = matrix(NA, nrow=nsite, ncol=choose(maxnsurveys,3)+3)
  }
  if(verbose) {
    cat("\nCalculating trajectory angles...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    if(!all) {
      for(s1 in 1:(nsurveysite[i1]-2)) {
        if(relativeToInitial) {
          d12 = dmat[ind_surv1[1], ind_surv1[s1 + 1]]
          d23 = dmat[ind_surv1[s1 + 1], ind_surv1[s1 +2]]
          d13 = dmat[ind_surv1[1], ind_surv1[s1 + 2]]
        } else {
          d12 = dmat[ind_surv1[s1], ind_surv1[s1+1]]
          d23 = dmat[ind_surv1[s1+1], ind_surv1[s1+2]]
          d13 = dmat[ind_surv1[s1], ind_surv1[s1+2]]
        }
        angles[i1, s1] = .angleConsecutiveC(d12,d23,d13, add)
        # cat(paste(i1,s1,":", d12,d23,d13,.angleConsecutiveC(d12,d23,d13, TRUE),"\n"))
      }
      x <- angles[i1,1:(nsurveysite[i1]-2)] 
      angles[i1, ncol(angles)-2] = .meancircular(x, na.rm=TRUE)
      angles[i1, ncol(angles)-1] = .sdcircular(x, na.rm=TRUE)
      angles[i1, ncol(angles)] = .rhocircular(x, na.rm=TRUE)
    } else {
      cs = combn(length(ind_surv1),3)
      dsub = dmat[ind_surv1, ind_surv1]
      for(s in 1:ncol(cs)) {
        d12 = dsub[cs[1,s],cs[2,s]]
        d23 = dsub[cs[2,s],cs[3,s]]
        d13 = dsub[cs[1,s],cs[3,s]]
        angles[i1, s] = .angleConsecutiveC(d12,d23,d13, add)
      }
      x <- angles[i1,]
      angles[i1, ncol(angles)-2] = .meancircular(x, na.rm=TRUE)
      angles[i1, ncol(angles)-1] = .sdcircular(x, na.rm=TRUE)
      angles[i1, ncol(angles)] = .rhocircular(x, na.rm=TRUE)
    }
  }
  angles = as.data.frame(angles)
  row.names(angles)<-siteIDs
  if(!all) {
    if(relativeToInitial) {
      names(angles) <- c(paste0("t", rep("1",maxnsurveys -2), "-S", 
                                as.character(2:(maxnsurveys - 1))), 
                         "mean", "sd", "rho")    
    } else {
      names(angles)<-c(paste0("S",as.character(1:(maxnsurveys-2)),"-S",
                              as.character(2:(maxnsurveys-1))),
                       "mean", "sd", "rho")
    }
  } else {
    names(angles)<-c(paste0("A",as.character(1:(ncol(angles)-3))),"mean", "sd", "rho")
  }
  if(!stats) angles = angles[,1:(ncol(angles)-3), drop=FALSE]
  return(angles)
}

#' @rdname trajectorymetrics
#' @param betweenSegments Flag to indicate that angles should be calculated between trajectory segments or with respect to X axis.
#' @export
trajectoryAngles2D<-function(xy,sites,surveys,relativeToInitial=FALSE, betweenSegments=TRUE) {
  
  xy_temp<-as.data.frame(xy)
  xy_temp$sites<-sites
  xy_temp$surveys<-surveys
  xy_temp<-xy_temp[order(xy_temp$sites,xy_temp$surveys),]
  xy<-xy_temp[,1:2]
  sites<-c(xy_temp$sites)
  surveys<-c(xy_temp$surveys)
  
  siteIDs<-unique(sites)
  surveyIDs<-unique(surveys)
  nsite<-length(siteIDs)
  nsurvey<-length(surveyIDs)
  nsurveysite<-numeric(nsite)
  maxnsurveys<-length(surveys)
  
  dx<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey))
  x<-xy[,1]
  seq1<-c(nsurvey)
  for (i in 1:(nsite-1)){
    seq1<-c(seq1,seq1[i]+nsurvey)
  }
  seq2<-c(1)
  for (i in 1:(nsite-1)){
    seq2<-c(seq2,seq2[i]+nsurvey)
  }
  
  for (i in 1:nsite){
    dx[i,]<-x[seq2[i]:seq1[i]]
  }
  
  dy<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey))
  y<-xy[,2]
  for (i in 1:nsite){
    dy[i,]<-y[seq2[i]:seq1[i]]
  }
  
  if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  if(nrow(xy)!=nsite*nsurvey) stop("nrow(xy) need to be equal to 'number of sites * number of surveys'")
  
  
  
  #x modification
  dxmod<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
  for(i in 1:ncol(dxmod)){
    dxmod[,i]<-dx[,i]-dx[,i+1]
  }
  #y modification
  dymod<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
  for(i in 1:ncol(dymod)){
    dymod[,i]<-dy[,i]-dy[,i+1]
  }
  #dfxy
  dmod<-as.data.frame(cbind(dxmod,dymod))
  
  if(!betweenSegments){
    #angle alpha measurement
    Angle_alpha_temp<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    for(i in 1:ncol(Angle_alpha_temp)){
      Angle_alpha_temp[,i]<-apply(dmod[c(i,i+nsurvey-1)], 1, function(irow) {
        atan(irow[2]/irow[1])
      })
    }
    Angle_alpha_temp<-Angle_alpha_temp*(180/pi)
    
    
    dxy<-as.data.frame (cbind(dx,dy))
    Angle_alpha<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    for (i in 1:ncol(Angle_alpha_temp)) {
      Angle_alpha[,i]<-as.numeric(c(ifelse(dxy[,i]==dxy[,i+1] & dxy[,i+(nsurvey)]<dxy[,i+(nsurvey+1)],0,
                                           ifelse(dxy[,i]==dxy[,i+1] & dxy[,i+(nsurvey)]>dxy[,i+(nsurvey+1)],180,
                                                  ifelse(dxy[,i]<dxy[,i+1] & dxy[,i+(nsurvey)]==dxy[,i+(nsurvey+1)],90,
                                                         ifelse(dxy[,i]>dxy[,i+1] & dxy[,i+(nsurvey)]==dxy[,i+(nsurvey+1)],270,
                                                                ifelse(dxy[,i] < dxy[,i+1] & dxy[,i+(nsurvey)]< dxy[,i+(nsurvey+1)],90-Angle_alpha_temp[,i],
                                                                       ifelse(dxy[,i]< dxy[,i+1] & dxy[,i+(nsurvey)] > dxy[,i+(nsurvey+1)],90-Angle_alpha_temp[,i],
                                                                              ifelse(dxy[,i]>dxy[,i+1] & dxy[,i+(nsurvey)] > dxy[,i+(nsurvey+1)],270-Angle_alpha_temp[,i],
                                                                                     ifelse(dxy[,i]>dxy[,i+1] & dxy[,i+(nsurvey)] < dxy[,i+(nsurvey+1)],270-Angle_alpha_temp[,i],"ERROR"))))))))))
    }
    colnames(Angle_alpha) <- c(paste0("Axis2", "-t", 
                                      as.character(1:(nsurvey-1))))  
    angles_out = Angle_alpha
  }
  else if(!relativeToInitial){
    #angle Theta measurement
    #position according to first trajectory of each triplet
    dxy<-as.data.frame (cbind(dx,dy))
    pos<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:(nsurvey-2)){
      pos[,i]<-c((dxy[,i+1]-dxy[,i])*(dxy[,nsurvey+i+2]-dxy[,nsurvey+i])-(dxy[,nsurvey+i+1]-dxy[,nsurvey+i])*(dxy[,i+2]-dxy[,i]))
    }
    
    #consecutive segement lenght
    Scons<-as.data.frame(trajectoryLengths2D(xy,sites,surveys, relativeToInitial=FALSE))
    Scons<-Scons[,-c(ncol(Scons)-1,ncol(Scons))]
    
    #calculation of length t to t+2: 1to3,2to4....
    distn2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:(nsurvey-2)){
      distn2[,i]<-sqrt(((dxy[,i+2]-dxy[,i])^2)+((dxy[,nsurvey+i+2]-dxy[,nsurvey+i])^2))
    }
    
    #recovering or distancing pattern for each consecutive triplet
    RDTcons<-distn2-Scons
    
    # Angle theta temp (0-180°)
    Angle_theta_temp<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    xvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    xvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    for(i in 1:(nsurvey-2)){
      xvector1[,i] = dxy[,i] - dxy[,i+1]
      yvector1[,i] = dxy[,i+nsurvey] - dxy[,i+nsurvey+1]
      xvector2 [,i] = dxy[,i+2] - dxy[,i+1]
      yvector2 [,i] = dxy[,i+nsurvey+2] - dxy[,i+nsurvey+1]
      
      num = (xvector1[,i] * xvector2[,i] + yvector1[,i] * yvector2[,i])
      den = sqrt(xvector1[,i]^2 + yvector1[,i]^2) * sqrt(xvector2[,i]^2 +yvector2[,i]^2)
      Angle_theta_temp[,i] = (360 * acos(num/den))/(2 * pi)
    }
    
    # Angle theta (0-360°)
    Angle_theta<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:ncol(Angle_theta_temp)){
      Angle_theta[,i]<-c(ifelse(Angle_theta_temp[,i]==180,0,
                                ifelse(Angle_theta_temp[,i]==0,180,
                                       ifelse(pos[,i]<0 & RDTcons[,i]<0,180-Angle_theta_temp[,i],
                                              ifelse(pos[,i]>0 & RDTcons[,i]<0 ,360-(180-Angle_theta_temp[,i]),
                                                     ifelse(pos[,i]<0 & RDTcons[,i]>0 & Angle_theta_temp[,i]<90,180-Angle_theta_temp[,i],
                                                            ifelse(pos[,i]<0 & RDTcons[,i]>0 & Angle_theta_temp[,i]>90,180-Angle_theta_temp[,i], 
                                                                   ifelse(pos[,i]>0 & RDTcons[,i]>0 & Angle_theta_temp[,i]<90 ,270-(90-Angle_theta_temp[,i]),
                                                                          ifelse(pos[,i]>0 & RDTcons[,i]>0 & Angle_theta_temp[,i]>90 ,270+(Angle_theta_temp[,i]-90),"ERROR")))))))))
      
      
    }
    
    colnames(Angle_theta) <- c(paste0("t", as.character(1:(nsurvey-2)), "-t", 
                                      as.character(2:(nsurvey-1))))    
    angles_out = Angle_theta
  } else {
    #angle Theta temp measurement
    #position according to first trajectory of each triplet
    
    dxy<-as.data.frame (cbind(dx,dy))
    pos<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:(nsurvey-2)){
      pos[,i]<-c((dxy[,i+1]-dxy[,i])*(dxy[,nsurvey+i+2]-dxy[,nsurvey+i])-(dxy[,nsurvey+i+1]-dxy[,nsurvey+i])*(dxy[,i+2]-dxy[,i]))
    }
    
    #consecutive segment length
    Scons<-as.data.frame(trajectoryLengths2D(xy,sites,surveys, relativeToInitial=FALSE))
    Scons<-Scons[,-c(ncol(Scons)-1,ncol(Scons))]
    
    #calculation of length t to t+2: 1to3,2to4....
    distn2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:(nsurvey-2)){
      distn2[,i]<-sqrt(((dxy[,i+2]-dxy[,i])^2)+((dxy[,nsurvey+i+2]-dxy[,nsurvey+i])^2))
    }
    
    #recovering or distancing pattern for each consecutive triplet
    RDTcons<-distn2-Scons
    
    # Angle theta temp (0-180°)
    Angle_theta_temp<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    xvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    xvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    for(i in 1:(nsurvey-2)){
      xvector1[,i] = dxy[,i] - dxy[,i+1]
      yvector1[,i] = dxy[,i+nsurvey] - dxy[,i+nsurvey+1]
      xvector2 [,i] = dxy[,i+2] - dxy[,i+1]
      yvector2 [,i] = dxy[,i+nsurvey+2] - dxy[,i+nsurvey+1]
      
      num = (xvector1[,i] * xvector2[,i] + yvector1[,i] * yvector2[,i])
      den = sqrt(xvector1[,i]^2 + yvector1[,i]^2) * sqrt(xvector2[,i]^2 +yvector2[,i]^2)
      Angle_theta_temp[,i] = (360 * acos(num/den))/(2 * pi)
    }
    
    
    #Angle_omega measurement
    #position according to first trajectory of each triplet
    pos<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    for (i in 1:(nsurvey-2)){
      pos[,i]<-c((dxy[,2]-dxy[,1])*(dxy[,nsurvey+i+2]-dxy[,nsurvey+1])-(dxy[,nsurvey+2]-dxy[,nsurvey+1])*(dxy[,i+2]-dxy[,1]))
    }
    
    #Angle_omega_temp(0-180°)
    Angle_omega_temp<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-2))
    xvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector1<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    xvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    yvector2<-as.data.frame(matrix(NA, nrow=nsite, ncol=nsurvey-1))
    for(i in 1:(nsurvey-2)){
      xvector1[,i] = dxy[,1] - dxy[,2]
      yvector1[,i] = dxy[,1+nsurvey] - dxy[,2+nsurvey]
      xvector2 [,i] = dxy[,i+2] - dxy[,2]
      yvector2 [,i] = dxy[,i+nsurvey+2] - dxy[,2+nsurvey]
      
      num = (xvector1[,i] * xvector2[,i] + yvector1[,i] * yvector2[,i])
      den = sqrt(xvector1[,i]^2 + yvector1[,i]^2) * sqrt(xvector2[,i]^2 +yvector2[,i]^2)
      Angle_omega_temp[,i] = (360 * acos(num/den))/(2 * pi)
    }
    
    #Angle_omega(0-360°)
    Angle_omega<-as.data.frame(matrix(NA, nrow=nrow(Angle_omega_temp), ncol=ncol(Angle_omega_temp)-1))
    for (i in 1:ncol(Angle_omega_temp)){
      Angle_omega[,i]<-c(ifelse(Angle_theta_temp[,i]==180,0,
                                ifelse(Angle_theta_temp[,i]==0,180,
                                       ifelse(pos[,i]<0,180-Angle_omega_temp[,i],
                                              ifelse(pos[,i]>0,360-(180-Angle_omega_temp[,i]),"ERROR")))))
    }
    
    colnames(Angle_omega) <- c(paste0("S1", "-t", 
                                      as.character(2:(nsurvey-1))))
    angles_out = Angle_omega
  }
  return(angles_out)
}

#' @rdname trajectorymetrics
#' @param target An integer vector of the ecosystem states to be projected.
#' @param trajectory An integer vector of the trajectory onto which target states are to be projected.
#' @param tol Numerical tolerance value to determine that projection of a point lies within the trajectory.
#' @export
trajectoryProjection<-function(d, target, trajectory, tol = 0.000001, add=TRUE) {
  if(length(trajectory)<2) stop("Trajectory needs to include at least two states")
  if(length(trajectory)!=length(unique(trajectory))) stop("Trajectory states must be different")
  dmat = as.matrix(d)
  npoints = length(target)
  nsteps = length(trajectory) -1
  #Distance betwen target points and trajectory points
  d2ref = dmat[target, trajectory, drop=FALSE]
  #Distance between trajectory steps
  dsteps = diag(dmat[trajectory[1:(length(trajectory)-1)], trajectory[2:length(trajectory)]])
  #Cumulative distance between steps
  dstepcum = rep(0,nsteps+1)
  if(nsteps>1) {
    for(i in 2:nsteps) {
      dstepcum[i] = dstepcum[i-1]+dsteps[i-1]
    }
  }
  dstepcum[nsteps+1] = sum(dsteps)
  
  projH = matrix(NA, nrow=npoints, ncol = nsteps)
  projA1 = matrix(NA, nrow=npoints, ncol = nsteps)
  projA2 = matrix(NA, nrow=npoints, ncol = nsteps)
  whichstep = rep(NA, npoints)
  dgrad = rep(NA, npoints)
  posgradseg = rep(NA, npoints)
  posgradtraj = rep(NA, npoints)
  
  for(i in 1:npoints) {
    for(j in 1:nsteps) {
      p <-.projectionC(dsteps[j], d2ref[i, j], d2ref[i, j+1], add)
      if((!is.na(p[3])) & (p[1]>-tol) & (p[2]>-tol)) {
        projA1[i,j] = p[1]
        projA2[i,j] = p[2]
        projH[i,j] = p[3]
        if(is.na(dgrad[i])) {
          dgrad[i] = p[3]
          whichstep[i] = j
          posgradseg[i] = p[1]/dsteps[j]
        } else {
          if(p[3]<dgrad[i]) {
            dgrad[i] = p[3]
            whichstep[i] = j
            posgradseg[i] = p[1]/dsteps[j]
          }
        }
      }
    }
    if(!is.na(whichstep[i])) {
      dg = dstepcum[whichstep[i]]+projA1[i,whichstep[i]]
      posgradtraj[i] = dg/sum(dsteps)
    }
  }
  res = data.frame(distanceToTrajectory=dgrad, segment = whichstep, relativeSegmentPosition = posgradseg, relativeTrajectoryPosition = posgradtraj)
  row.names(res)<-row.names(d2ref)
  return(res)
}


#' @rdname trajectorymetrics
#' @param symmetric A logical flag to indicate a symmetric convergence comparison of trajectories.
#' @export
trajectoryConvergence<-function(d, sites, surveys = NULL, symmetric = FALSE, add=TRUE, verbose = FALSE){
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  n = nrow(as.matrix(d))

  #Init output
  tau = matrix(NA, nrow=nsite, ncol = nsite)
  rownames(tau) = siteIDs
  colnames(tau) = siteIDs
  p.value = tau
  dmat = as.matrix(d)
  if(verbose) {
    cat("\nCalculating trajectory convergence...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }
  for(i1 in 1:(nsite-1)) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(i2 in (i1+1):nsite) {
      ind_surv2 = which(sites==siteIDs[i2])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv2 = ind_surv2[order(surveys[sites==siteIDs[i2]])]
      if(!symmetric) {
        trajectory = ind_surv2
        target = ind_surv1
        trajProj = trajectoryProjection(d,target, trajectory, add=add)
        dT = trajProj$distanceToTrajectory
        mk.test = MannKendall(dT)
        tau[i1,i2] = mk.test$tau
        p.value[i1,i2] = mk.test$sl
        trajectory = ind_surv1
        target = ind_surv2
        trajProj = trajectoryProjection(d,target, trajectory, add=add)
        dT = trajProj$distanceToTrajectory
        mk.test = MannKendall(dT)
        tau[i2,i1] = mk.test$tau
        p.value[i2,i1] = mk.test$sl
      } 
      else {
        if(length(ind_surv1)==length(ind_surv2)) {
          dT = numeric(length(ind_surv1))
          for(j in 1:length(ind_surv1)) dT[j] = dmat[ind_surv1[j], ind_surv2[j]]
          mk.test = MannKendall(dT)
          tau[i1,i2] = mk.test$tau
          p.value[i1,i2] = mk.test$sl
          tau[i2,i1] = mk.test$tau
          p.value[i2,i1] = mk.test$sl
        } else {
          warning(paste0("sites ",i1, " and ",i2," do not have the same number of surveys."))
        }
      }
    }
  }
  return(list(tau = tau, p.value = p.value))
}


#' @rdname trajectorymetrics
#' @export
trajectoryDirectionality<-function(d, sites, surveys = NULL, add=TRUE, verbose = FALSE) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<3)>0) stop("All sites need to be surveyed at least three times")
  
  dmat = as.matrix(d)
  #Init output
  dir = rep(NA, nsite)
  names(dir) = siteIDs
  if(verbose) {
    cat("\nAssessing trajectory directionality...\n")
    tb = txtProgressBar(1, nsite, style=3)
  }

  for(i1 in 1:nsite) {
    if(verbose) setTxtProgressBar(tb, i1)
    ind_surv1 = which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 = ind_surv1[order(surveys[sites==siteIDs[i1]])]
    dsub = dmat[ind_surv1, ind_surv1]
    n = length(ind_surv1)
    den = 0
    num = 0
    if(n>2) {
      for(i in 1:(n-2)) {
        for(j in (i+1):(n-1)) {
          for(k in (j+1):n) {
            da = dsub[i,j]
            db = dsub[j,k]
            dab = dsub[i,k]
            theta = .angleConsecutiveC(da,db,dab, add)
            if(!is.na(theta)) {
              den = den + (da + db)
              num = num + (da + db)*((180-theta)/180)
            }
          }
        }
      }
      dir[i1] = num/den
    }
  }
  return(dir)
}
