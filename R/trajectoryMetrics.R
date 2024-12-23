#' Metrics for Ecological Trajectory Analysis
#' 
#' Ecological Trajectory Analysis (ETA) is a framework to analyze dynamics of ecosystems described as trajectories in a chosen space of multivariate resemblance (De \enc{Cáceres}{Caceres} et al. 2019).
#' ETA takes trajectories as objects to be analyzed and compared geometrically. 
#' 
#' Given a distance matrix between ecosystem states, the set of functions that provide ETA metrics are:
#' \itemize{
#' \item{Function \code{trajectoryDistances} calculates the distance between pairs of trajectories.}
#' \item{Function \code{trajectoryLengths} calculates lengths of directed segments and total path lengths of trajectories.}
#' \item{Function \code{trajectoryLengths2D} calculates lengths of directed segments and total path lengths of trajectories from 2D coordinates given as input.} 
#' \item{Function \code{trajectoryAngles} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryAngles2D} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryConvergence} performs the Mann-Kendall trend test on the distances between trajectories (symmetric test) or the distance between points of one trajectory to the other.}
#' \item{Function \code{trajectoryDirectionality} returns (for each trajectory) a statistic that measures directionality of the whole trajectory.}
#' \item{Function \code{trajectoryVariability} returns (for each trajectory) a statistic that measures the variability between the states included in the trajectory.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name trajectoryMetrics
#' @aliases trajectoryDistances trajectoryLengths trajectoryLengths2D trajectoryAngles trajectoryAngles2D
#'          trajectoryConvergence trajectoryDirectionality trajectoryVariability
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param distance.type The type of distance index to be calculated (Besse et al. 2016; De Cáceres et al. 2019):
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two trajectories.}
#'     \item{\code{SPD}: Segment path distance.}
#'     \item{\code{DSPD}: Directed segment path distance (default).}
#'   }
#' @param symmetrization Function used to obtain a symmetric distance, so that DSPD(T1,T2) = DSPD(T2,T1) (e.g., \code{mean} or \code{min}). If \code{symmetrization = NULL} then the symmetrization is not conducted and the output dissimilarity matrix is not symmetric. 
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
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
#' @return
#' Function \code{trajectoryDistances} returns an object of class \code{\link{dist}} containing the distances between trajectories (if \code{symmetrization = NULL} then the object returned is of class \code{matrix}). 
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
#' Function \code{trajectoryConvergence} returns a list with two elements:
#' \itemize{
#'   \item{\code{tau}: A matrix with the statistic (Mann-Kendall's tau) of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the statistic of the test of the row trajectory approaching the column trajectory.}
#'   \item{\code{p.value}: A matrix with the p-value of the convergence/divergence test between trajectories. If \code{symmetric=TRUE} then the matrix is square. Otherwise the p-value indicates the test of the row trajectory approaching the column trajectory.}
#' }
#' 
#' Function \code{trajectoryDirectionality} returns a vector with directionality values (one per trajectory).
#' 
#' Function \code{trajectoryVariability} returns a vector with total variability values (one per trajectory).
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
#' @seealso \code{\link{segmentDistances}}, \code{\link{trajectoryPlot}}, \code{\link{transformTrajectories}}, \code{\link{trajectoryProjection}}
#' 
#' @examples 
#' #Description of sites and surveys
#' sites <- c("1","1","1","2","2","2")
#' surveys <- c(1, 2, 3, 1, 2, 3)
#'   
#' #Raw data table
#' xy <- matrix(0, nrow=6, ncol=2)
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
#' d <- dist(xy)
#' d
#'   
#' #Trajectory data
#' x <- defineTrajectories(d, sites, surveys)
#' 
#' #Trajectory lengths
#' trajectoryLengths(x)
#' trajectoryLengths2D(xy, sites, surveys)
#' 
#' #Trajectory angles
#' trajectoryAngles(x)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = TRUE)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = FALSE)
#' 
#' #Distances between trajectories
#' trajectoryDistances(x, distance.type = "Hausdorff")
#' trajectoryDistances(x, distance.type = "DSPD")
#'   
#' #Should give the same results if surveys are not in order 
#' #(here we switch surveys for site 2)
#' temp = xy[5,]
#' xy[5,] = xy[6,]
#' xy[6,] = temp
#' surveys[5] = 3
#' surveys[6] = 2
#' d <- dist(xy)
#' x <- defineTrajectories(d, sites, surveys)
#'   
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)   
#' trajectoryLengths(x)
#' trajectoryLengths2D(xy, sites, surveys)
#' trajectoryAngles(x)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = TRUE)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = FALSE)
#' trajectoryDistances(x, distance.type = "Hausdorff")
#' trajectoryDistances(x, distance.type = "DSPD")
#'  
#' @export
trajectoryDistances<-function(x, distance.type="DSPD", symmetrization = "mean" , add=TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  distance.type <- match.arg(distance.type, c("DSPD", "SPD", "Hausdorff"))
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
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

#' @rdname trajectoryMetrics
#' @param relativeToInitial Flag to indicate that lengths or angles should be calculated with respect to initial survey.
#' @param all Flag to indicate that lengths or angles are desired for all segments or for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, length or angles are calculated according to relativeToInitial flag.
#' @export
trajectoryLengths<-function(x, relativeToInitial = FALSE, all=FALSE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
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
    for(i1 in 1:nsite) {
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

#' @rdname trajectoryMetrics
#' @param xy Matrix with 2D coordinates in a Cartesian space (typically an ordination of ecosystem states).
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param surveys A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @export
trajectoryLengths2D<-function(xy, sites, surveys = NULL, relativeToInitial=FALSE, all=FALSE) {
  if(length(sites)!=nrow(xy)) stop("'sites' needs to be of length equal to the number of rows in xy")
  if(!is.null(surveys)) {
    if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  } else {
    surveys <- rep(NA, length(sites))
    for(s in unique(sites)) {
      surveys[sites==s] <- 1:sum(sites==s)
    }
  }
  
  #order inputs by sites and surveys
  xy_temp<-as.data.frame(xy)
  xy_temp$sites<-sites
  xy_temp$surveys<-surveys
  xy_temp<-xy_temp[order(xy_temp$sites,xy_temp$surveys),]
  xy<-xy_temp[,1:2]
  sites<-c(xy_temp$sites)
  surveys<-c(xy_temp$surveys)
  
  siteIDs <- unique(sites)
  surveyIDs<-unique(surveys)
  nsite<-length(siteIDs)
  nsurvey<-length(surveyIDs)
  if(nsite!=nrow(xy)/nsurvey) stop("'sites' needs to be of length equal in xy")
  if(nrow(xy)!=nsurvey*nsite) stop("All sites need to be surveyed at least twice")
  
  #prep all
  x <- defineTrajectories(dist(xy),sites,surveys)
  return(trajectoryLengths(x,relativeToInitial = relativeToInitial, all = all))
}

#' @rdname trajectoryMetrics
#' @param all A flag to indicate that angles are desired for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, angles are calculated for consecutive segments only.
#' @param stats A flag to indicate that circular statistics are desired (mean, standard deviation and mean resultant length, i.e. rho)
#' @export
trajectoryAngles<-function(x, all = FALSE, relativeToInitial = FALSE, stats = TRUE, add=TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat <- as.matrix(d)
  n <- nrow(dmat)
  
  maxnsurveys <- max(nsurveysite)
  if(!all) angles <- matrix(NA, nrow=nsite, ncol=maxnsurveys+1)
  else {
    angles <- matrix(NA, nrow=nsite, ncol=choose(maxnsurveys,3)+3)
  }
  for(i1 in 1:nsite) {
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

#' @rdname trajectoryMetrics
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



#' @rdname trajectoryMetrics
#' @param symmetric A logical flag to indicate a symmetric convergence comparison of trajectories.
#' @export
trajectoryConvergence<-function(x, symmetric = FALSE, add=TRUE){
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
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


#' @rdname trajectoryMetrics
#' @export
trajectoryDirectionality <- function(x, add=TRUE) {
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
  #Init output
  dir <- rep(NA, nsite)
  names(dir) <- siteIDs

  for(i1 in 1:nsite) {
    ind_surv1 <- which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 <- ind_surv1[order(surveys[sites==siteIDs[i1]])]
    dsub <- dmat[ind_surv1, ind_surv1]
    n <- length(ind_surv1)
    den <- 0
    num <- 0
    if(n>2) {
      for(i in 1:(n-2)) {
        for(j in (i+1):(n-1)) {
          for(k in (j+1):n) {
            da <- dsub[i,j]
            db <- dsub[j,k]
            dab <- dsub[i,k]
            theta <- .angleConsecutiveC(da,db,dab, add)
            if(!is.na(theta)) {
              den <- den + (da + db)
              num <- num + (da + db)*((180-theta)/180)
            }
          }
        }
      }
      dir[i1] <- num/den
    }
  }
  return(dir)
}


#' @rdname trajectoryMetrics
#' @export
trajectoryVariability<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<2)>0) stop("All sites need to be surveyed at least twice")
  
  dmat <- as.matrix(d)
  var <- rep(NA, nsite)
  names(var) <- siteIDs
  for(i in 1:nsite) {
    ind_surv = which(sites==siteIDs[i])
    dsub <- dmat[ind_surv, ind_surv]
    r <- ncol(dsub)
    var[i] <- sum(as.vector(as.dist(dsub))^2)/(r^2)
  }
  return(var)
}
