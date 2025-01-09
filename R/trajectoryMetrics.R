#' Trajectory metrics
#' 
#' Set of functions to estimate metrics describing individual trajectories. Given input trajectory data, the set of functions that provide ETA metrics are:
#' \itemize{
#' \item{Function \code{trajectoryLengths} calculates lengths of directed segments and total path lengths of trajectories.}
#' \item{Function \code{trajectoryLengths2D} calculates lengths of directed segments and total path lengths of trajectories from 2D coordinates given as input.} 
#' \item{Function \code{trajectorySpeeds} calculates speeds of directed segments and total path speed of trajectories.}
#' \item{Function \code{trajectorySpeeds2D} calculates speeds of directed segments and total path speed of trajectories from 2D coordinates given as input.} 
#' \item{Function \code{trajectoryAngles} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryAngles2D} calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points.}
#' \item{Function \code{trajectoryDirectionality} calculates (for each trajectory) a statistic that measures directionality of the whole trajectory.}
#' \item{Function \code{trajectoryVariability} calculates (for each trajectory) a statistic that measures the variability between the states included in the trajectory.}
#' \item{Function \code{trajectoryMetrics} evaluates several trajectory metrics at once.}
#' \item{Function \code{trajectoryWindowMetrics} evaluates several trajectory metrics on subtrajectories defined using moving windows.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name trajectoryMetrics
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param relativeToInitial Flag to indicate that lengths or angles should be calculated with respect to initial survey.
#' @param all Flag to indicate that lengths or angles are desired for all segments or for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, length or angles are calculated according to relativeToInitial flag.
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
#' Function \code{trajectoryAngles} calculates angles between consecutive segments in degrees. For each pair of segments, the angle between the two is defined on the plane that contains the two segments, and measures the change in direction (in degrees) from one segment to the other. 
#' Angles are always positive, with zero values indicating segments that are in a straight line, and values equal to 180 degrees for segments that are in opposite directions. If \code{all = TRUE}
#' angles are calculated between the segments corresponding to all ordered triplets. Alternatively, if \code{relativeToInitial = TRUE} angles are calculated for each segment with respect to the initial survey.
#' 
#' Function \code{trajectoryAngles2D} calculates angles between consecutive segments in degrees from 2D coordinates given as input. For each pair of segments, the angle between the two is defined on the plane that contains the two segments, and measures the change in direction (in degrees) from one segment to the other. 
#' Angles are always positive (O to 360), with zero values indicating segments that are in a straight line, and values equal to 180 degrees for segments that are in opposite directions. 
#' If \code{all = TRUE} angles are calculated between the segments corresponding to all ordered triplets. Alternatively, if \code{relativeToInitial = TRUE} angles are calculated for each segment with respect to the initial survey.
#' If \code{betweenSegments = TRUE} angles are calculated between segments of trajectory, otherwise, If \code{betweenSegments = FALSE}, angles are calculated considering Y axis as the North (0°).
#' 
#' Function \code{trajectoryDirectionality} evaluates the directionality metric proposed in De \enc{Cáceres}{Caceres} et al (2019). If \code{nperm} is supplied, then the function
#' performs a permutational test to evaluate the significance of directionality, where the null hypothesis entails a random order of surveys within each trajectory. The p-value corresponds to the proportion of
#' permutations with a directional value equal or larger than the observed.
#'
#' @return
#' 
#' Functions \code{trajectoryLengths} and  \code{trajectoryLengths2D} return a data frame with the length of each segment on each trajectory and the total length of all trajectories. 
#' If \code{relativeToInitial = TRUE} lengths are calculated between the initial survey and all the other surveys.
#' If \code{all = TRUE} lengths are calculated for all segments. 
#' 
#' Functions \code{trajectorySpeeds} and \code{trajectorySpeeds2D} return a data frame with the speed of each segment on each trajectory and the total speeds of all trajectories. Units depend on the units of distance matrix and the units of \code{times} of the input trajectory data.
#'
#' Function \code{trajectoryAngles} returns a data frame with angle values on each trajectory. If \code{stats=TRUE}, then the mean, standard deviation and mean resultant length of those angles are also returned. 
#' 
#' Function \code{trajectoryAngles2D} returns a data frame with angle values on each trajectory. If \code{betweenSegments=TRUE}, then angles are calculated between trajectory segments, alternatively, If \code{betweenSegments=FALSE}, angles are calculated considering Y axis as the North (0°).
#' 
#' Function \code{trajectoryDirectionality} returns a vector with directionality values (one per trajectory). If \code{nperm} is not missing, the function returns a data frame
#' with a column of directional values and a column of p-values corresponding to the result of the permutational test.
#' 
#' Function \code{trajectoryVariability} returns a vector with total variability values (one per trajectory).
#' 
#' Function \code{trajectoryMetrics} returns a data frame where rows are trajectories and columns are different trajectory metrics.
#' 
#' Function \code{trajectoryWindowMetrics} returns a data frame where rows are midpoints over trajectories and columns correspond to different trajectory metrics.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @author Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la Baie de Saint-Brieuc
#' @author Nicolas Djeghri, UBO
#' 
#' @references
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). 
#' Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
#' @seealso \code{\link{trajectoryComparison}}, \code{\link{trajectoryPlot}}, \code{\link{transformTrajectories}}, \code{\link{cycleMetrics}}
#' 
#' @examples 
#' #Description of sites and surveys
#' sites <- c("1","1","1","2","2","2")
#' surveys <- c(1, 2, 3, 1, 2, 3)
#' times <- c(0, 1.5, 3, 0, 1.5, 3)
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
#' x <- defineTrajectories(d, sites, surveys, times)
#' 
#' #Trajectory lengths
#' trajectoryLengths(x)
#' trajectoryLengths2D(xy, sites, surveys)
#' 
#' #Trajectory speeds
#' trajectorySpeeds(x)
#' trajectorySpeeds2D(xy, sites, surveys, times)
#'
#' #Trajectory angles
#' trajectoryAngles(x)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = TRUE)
#' trajectoryAngles2D(xy, sites, surveys, betweenSegments = FALSE)
#' 
#' #Several metrics at once
#' trajectoryMetrics(x)  
#'  
#' @export
trajectoryLengths<-function(x, relativeToInitial = FALSE, all=FALSE) {
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
  surveyIDs<-unique(surveys)
  nsurvey<-length(surveyIDs)
  
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  
  dmat <- as.matrix(d)
  n <- nrow(dmat)

  maxnsurveys <- max(nsurveysite)

  if(!all) {
    lengths = as.data.frame(matrix(NA, nrow=nsite, ncol=maxnsurveys))
    row.names(lengths)<-siteIDs
    if(relativeToInitial) names(lengths)<-c(paste0("Lt1_t",as.character(2:(maxnsurveys))),"Path")
    else names(lengths)<-c(paste0("S",as.character(1:(maxnsurveys-1))),"Path")
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
  return(trajectoryLengths(defineTrajectories(dist(xy),sites,surveys),
                           relativeToInitial = relativeToInitial, 
                           all = all))
}

#' @rdname trajectoryMetrics
#' @export
trajectorySpeeds<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  d <- x$d
  surveys <- x$metadata$surveys
  times <- x$metadata$times
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
  surveyIDs<-unique(surveys)
  nsurvey<-length(surveyIDs)
  
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  maxnsurveys <- max(nsurveysite)
  
  tl <- trajectoryLengths(x, relativeToInitial = FALSE, all = FALSE)
  speeds <- tl
  for(i in 1:nsite) {
    ind_surv <- which(sites==siteIDs[i])
    #Surveys may not be in order
    ind_surv <- ind_surv[order(surveys[sites==siteIDs[i]])]
    for(s in 1:(nsurveysite[i]-1)) {
      speeds[i,s] <- speeds[i,s]/(times[ind_surv[s+1]] - times[ind_surv[s]])
    }
    speeds[i, maxnsurveys] <- speeds[i, maxnsurveys]/(times[ind_surv[length(ind_surv)]] - times[ind_surv[1]])
  }
  return(speeds)
}

#' @rdname trajectoryMetrics
#' @param times A numeric vector indicating the time corresponding to each ecosystem state.
#' @export
trajectorySpeeds2D<-function(xy, sites, surveys = NULL, times = NULL) {
  if(length(sites)!=nrow(xy)) stop("'sites' needs to be of length equal to the number of rows in xy")
  return(trajectorySpeeds(defineTrajectories(dist(xy),sites,surveys,times)))
}

#' @rdname trajectoryMetrics
#' @param all A flag to indicate that angles are desired for all triangles (i.e. all pairs of segments) in the trajectory. If FALSE, angles are calculated for consecutive segments only.
#' @param stats A flag to indicate that circular statistics are desired (mean, standard deviation and mean resultant length, i.e. rho)
#' @export
trajectoryAngles<-function(x, all = FALSE, relativeToInitial = FALSE, stats = TRUE, add=TRUE) {
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
#' @param nperm The number of permutations to be used in the directionality test.
#' @export
trajectoryDirectionality <- function(x, add=TRUE, nperm = NA) {
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
  
  dir_traj <- function(d, sites, surveys, add = TRUE) {
    siteIDs <- unique(sites)
    nsite <- length(siteIDs)
    dmat <- as.matrix(d)
    #Init output
    dir <- rep(NA, nsite)
    for(s in 1:nsite) {
      ind_surv <- which(sites==siteIDs[s])
      #Surveys may not be in order
      if(!is.null(surveys)) ind_surv <- ind_surv[order(surveys[sites==siteIDs[s]])]
      dsub <- dmat[ind_surv, ind_surv]
      n <- length(ind_surv)
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
        dir[s] <- num/den
      }
    }
    return(dir)
  }
  
  if(is.na(nperm)) {
    dir <- dir_traj(d, sites, surveys, add)
    names(dir) <- siteIDs
    return(dir)
  } else {
    if(!is.numeric(nperm)) stop("'nperm' should be an integer")
    if(length(nperm)>1) stop("'nperm' should be a single number")
    if(nperm<1) stop("'nperm' should be larger than 0")
    nperm <- as.integer(nperm)
    dir <- dir_traj(d, sites, surveys, add)
    dirp <- matrix(NA, nperm, nsite)
    for(p in 1:nperm) {
      surveysp <- surveys
      for(i in 1:nsite) surveysp[sites==siteIDs[i]] <- sample(surveysp[sites==siteIDs[i]])
      dirp[p,] <- dir_traj(d,sites,surveysp, add)
    }
    p.value <- (colSums(sweep(dirp,2,dir, "-")>=0)+1)/(nperm+1)
    return(data.frame(directionality = dir, p.value = p.value, row.names = siteIDs))
  }
}

#' @rdname trajectoryMetrics
#' @export
trajectoryVariability<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    selec <- x$metadata$internal==TRUE
    x$metadata <- x$metadata[selec,]
    x$d <- as.dist(as.matrix(x$d)[selec,selec])
    
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    selec <- x$metadata$internal==TRUE
    x$metadata <- x$metadata[selec,]
    x$d <- as.dist(as.matrix(x$d)[selec,selec])
    
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  d <- x$d
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

#' @rdname trajectoryMetrics
#' @export
trajectoryMetrics <- function(x, add = TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
    warning("Function cycleMetrics() may be more appropriate for cycles")
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  siteIDs <- unique(sites)
  df <-  data.frame(trajectory = siteIDs, site = NA,
                    n = NA, t_start = NA, t_end = NA, 
                    length = NA, mean_speed = NA, mean_angle = NA,
                    directionality = NA, variability = NA)
  if (inherits(x, "cycles") || inherits(x, "fd.trajectories") || inherits(x, "sections")){
    for (i in 1:length(siteIDs)){
      df$site[i] <-unique(x$metadata$sites[sites==siteIDs[i]])
    }
  } else {
    df$site <- NULL
  }
  for(i in 1:length(siteIDs)) {
    df$n[i] <- sum(sites==siteIDs[i])
  }
  df$t_start <- tapply(x$metadata$times,sites,min)
  df$t_end <- tapply(x$metadata$times,sites,max)
  df$length <- trajectoryLengths(x)$Path
  df$mean_speed <- trajectorySpeeds(x)$Path
  df$mean_angle <- trajectoryAngles(x, add = add)$mean
  df$directionality <- trajectoryDirectionality(x, add = add)
  df$variability <- trajectoryVariability(x)
  return(df)
}
#' @rdname trajectoryMetrics
#' @param bandwidth Bandwidth of the moving windows (in units of surveys or times, depending on \code{type})
#' @param type A string, either "surveys" or "times", indicating how windows are defined.
#' @export
trajectoryWindowMetrics <- function(x, bandwidth, type = "surveys", add = TRUE) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  match.arg(type, c("surveys", "times"))
  if(!is.numeric(bandwidth)) stop("'bandwidth' should be an numeric value")
  if(length(bandwidth)>1) stop("'bandwidth' should be a single number")
  if(type=="surveys") {
    if(bandwidth<1) stop("'bandwidth' should be at least 1")
    bandwidth <- as.integer(bandwidth)
  } else {
    if(bandwidth<=0) stop("'bandwidth' should be larger than 0")
  }
  
  surveys <- x$metadata$surveys
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
  
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(sum(nsurveysite<2)>0) stop("All sites need to be surveyed at least twice")
  
  df <-  data.frame(trajectory = sites, site = NA, midpoint = surveys,
                    t_start = NA,t_end = NA, n = NA, 
                    length = NA, mean_speed = NA, mean_angle = NA,
                    directionality = NA, variability = NA)
  if (inherits(x, "cycles") || inherits(x, "fd.trajectories") || inherits(x, "sections")){
    for (i in 1:length(sites)){
      df$site[i] <-unique(x$metadata$sites[sites==sites[i]])
    }
  } else {
    df$site <- NULL
  }
  for(i in 1:length(sites)) {
    surveys_i <- surveys[sites==sites[i]]
    times_i  <- times[sites==sites[i]]
    if(type=="surveys") {
      surveys_window <- surveys_i[abs(surveys_i - surveys[i])<=bandwidth]
    } else {
      surveys_window <- surveys_i[abs(times_i - times[i])<=bandwidth]
    }
    if (inherits(x, "cycles")|inherits(x, "fd.trajectories")|inherits(x, "sections")){
      x_i <- subsetTrajectories(x = x, subtrajectory_selection = sites[i], survey_selection = surveys_window)
    }else{
      x_i <- subsetTrajectories(x = x, site_selection = sites[i], survey_selection = surveys_window)
    }
    df$t_start[i] <- min(x_i$metadata$times)
    df$t_end[i] <- max(x_i$metadata$times)
    df$n[i] <- length(surveys_window)
    df$length[i] <- trajectoryLengths(x_i)$Path
    df$mean_speed[i] <- trajectorySpeeds(x_i)$Path
    if(length(surveys_window)>2) df$mean_angle[i] <- trajectoryAngles(x_i, add = add)$mean
    if(length(surveys_window)>2) df$directionality[i] <- trajectoryDirectionality(x_i, add = add)
    df$variability[i] <- trajectoryVariability(x_i)
  }
  return(df)
}
