#' Functions for building Trajectory Sections
#'
#' Trajectory sections are flexible way to cut longer trajectories. They are presently used chiefly in building cycles for cyclical ecological trajectory analysis (CETA) but might have other applications.
#' 
#' Trajectory sections functions:
#' \itemize{
#' \item{Function \code{extractTrajectorySections} reformats an object of class \code{\link{trajectories}} describing one or more trajectories into another object of class \code{\link{trajectories}} describing specified trajectory sections. Trajectory sections represent a way to subset trajectories flexibly. Cycles (see \code{\link{extractCycles}}) are a particular case of trajectory sections.}
#' \item{Function \code{interpolateEcolStates} compute interpolated ecological states and the new distance matrix associated (used in extractTrajectorySections).}
#' }
#' 
#' 
#' @encoding UTF-8
#' @name trajectorySections
#' @aliases extractTrajectorySections interpolateEcolStates sections
#' 
#' @details
#' Trajectory sections can be obtained using \code{extractTrajectorySections}. Trajectory sections allow to cut a longer trajectory into parts for further analyses. Cycles are specical case of trajectory sections.
#' A trajectory section TS(Traj,(tstart, BCstart),(tend, BCend)) is defined by the trajectory (Traj) it is obtained from, by an start and end times (tstart and tend) and start and end boundary conditions (BCstart, BCend).
#' The function extractTrajectorySections builds trajectory sections as a function of its arguments \code{Traj}, \code{tstart}, \code{tend}, \code{BCstart}, \code{BCend}.
#' 
#' Function \code{interpolateEcolStates} is called within \code{extractTrajectorySections} to interpolate ecological states when \code{tstart} and or \code{tend} do not have an associated measured ecological state within matrix \code{d}.
#' 
#' IMPORTANT: Trajectory sections comprises both "internal" and "external" ecological states (indicated in vector \code{internal}, see the output of function \code{extractTrajectorySections}).
#' "external" ecological states need a specific treatment in some calculations and for some operations within ETA, namely:
#' \itemize{
#'  \item{Centering, where external ecological states must be excluded from computation but included nonetheless in the procedure. This is automatically handled by function \code{\link{centerTrajectories}}.}
#'  \item{Trajectory variability, where only internal ecological states must be taken in account. This is handled automatically by function \code{\link{trajectoryVariability}}.}
#' }
#' Special care must also be taken when processing the data through principal coordinate analysis as external ecological states are effectively duplicated or interpolated in the output of \code{extractTrajectorySections}.
#' 
#' 
#' @return 
#' Function \code{extractTrajectorySections} returns the base information needed to describe trajectory sections. Its outputs are meant to be used as inputs for other ETA functions in order to obtain desired metrics. Importantly, within trajectory sections, ecological states can be considered "internal" or "external" and may necessitate special treatment (see details). Function \code{extractTrajectorySections} returns an object of class \code{\link{sections}} containing:
#' \itemize{
#'  \item{\code{d}: an object of class \code{\link{dist}}, the new distance matrix describing the trajectory sections. Ecological states may be duplicated in this matrix if trajectory sections overlap. As compared to the input matrix, \code{d} may also present deletions of ecological states that do not belong to any trajectory sections.}
#'  \item{\code{metadata}: an object of class \code{\link{data.frame}} describing the ecological states in \code{d} with columns:
#'    \itemize{
#'      \item{\code{sites}: the sites associated to each ecological states.}
#'      \item{\code{sections}: the names of the trajectory sections each ecological states belongs to.}
#'      \item{\code{surveys}: renumbering of the surveys to describe individual trajectory sections.}
#'      \item{\code{times}: the times associated to each ecological states.}
#'      \item{\code{internal}: a boolean vector with \code{TRUE} indicating "internal" ecological states whereas \code{FALSE} indicates "external" ecological states. This has important implications for the use of \code{extractTrajectorySections} outputs (see details).}
#'      }
#'    }
#'  \item{\code{interpolationInfo}: an output that only appear if ecological states have been interpolated. It is used internally by plotting functions (see \code{cyclePCoA}) but is not intended to be of interest to the end user.}
#' }
#' 
#' Function \code{interpolateEcolStates} returns an object of class \code{\link{dist}} including the desired interpolated ecological states.
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @rdname trajectorySections
#' @param x An object of class \code{\link{trajectories}} describing a cyclical trajectory.
#' @param Traj A vector of length equal to the number of desired trajectory sections indicating the trajectories from which trajectory sections must be build (see details).
#' @param tstart A vector of start times for each of the desired trajectory sections (see details).
#' @param tend A vector of end times for each of the desired trajectory sections (see details).
#' @param BCstart A vector of start boundary conditions (either \code{"internal"} or \code{"external"}) for each of the desired trajectory sections (see details).
#' @param BCend A vector of end boundary conditions (either \code{"internal"} or \code{"external"}) for each of the desired trajectory sections (see details). 
#' @param namesTS An optional vector giving a name for each of the desired trajectory sections (by default trajectory sections are simply numbered). 
#' 
#' @export
#' @keywords internal
extractTrajectorySections <- function(x,
                                      Traj,
                                      tstart,tend,BCstart,BCend,
                                      namesTS=1:length(Traj))
{
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  sites <- x$metadata$sites
  times <- x$metadata$times
  
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times))
    stop("The lengths of sites and times must corespond to the dimension of d")
  if (any(BCstart%in%c("internal","external")==FALSE))
    stop("BCstart and BCend can only take values 'internal' and 'external'")
  if (any(BCend%in%c("internal","external")==FALSE))
    stop("BCstart and BCend can only take values 'internal' and 'external'")
  
  #those will contain the sites, and times of all ecological states including the interpolated ones
  sitesTS <- sites
  timesTS <- times

  #this one will contain the interpolation coefficients
  intCoef <- rep(NA,length(sites))
  
  #first, check if there is ecological states to interpolate
  ToInterpolate <- integer(0)
  interpolated <- rep(FALSE,length(sites))
  for (i in unique(Traj)){
    timesTraj <- times[which(sites==i)]
    timesInt <- c(tstart[which(Traj==i)],tend[which(Traj==i)])
    
    if (any(is.na(cut(timesInt,range(timesTraj),include.lowest=TRUE)))) 
      stop(paste("tstart and/or tend out of bounds for trajectory",i))
    
    timesInt <- timesInt[timesInt%in%timesTraj==FALSE]
    timesInt <- unique(timesInt)#this prevents computing several times the same interpolated ecological state
    if (length(timesInt)>0){
      for (j in 1:length(timesInt)){
        truc <- timesTraj-timesInt[j]
        truc[truc>0] <- NA
        A <- timesTraj[which(truc==max(truc,na.rm=TRUE))]
        
        truc <- timesTraj-timesInt[j]
        truc[truc<0] <- NA
        B <- timesTraj[which(truc==min(truc,na.rm=TRUE))]
        
        newline <- c(intersect(which(times==A),which(sites==i)),
                     intersect(which(times==B),which(sites==i)),
                     (timesInt[j]-A)/(B-A))
        ToInterpolate <- rbind(ToInterpolate,newline)
      }
      sitesTS <- c(sitesTS,rep(i,length(timesInt)))
      timesTS <- c(timesTS,timesInt)
      interpolated <- c(interpolated,rep(TRUE,length(timesInt)))
      intCoef <- c(intCoef,ToInterpolate[,3])
    }
  }
  dTS <- d
  if(length(ToInterpolate)>0){
    dTS <- interpolateEcolStates(dTS,ToInterpolate)
  }
  dTS <- as.matrix(dTS)
  
  #Now build distances matrices and associated tags describing the trajectory sections
  #prepare the list that will store everything
  TS <- list()
  
  #and the element that will go into TS$metadata
  sections <- integer(0)
  selec <- integer(0)
  
  for (i in 1:length(Traj)){
    selection <- intersect(
      intersect(which(timesTS>=tstart[i]),
                which(timesTS<=tend[i])),
      which(sitesTS==Traj[i]))
    
    selec <- c(selec,selection)
    
    #Find out which are the external ecological states
    internali <- rep(TRUE,length(selection))
    #boundary conditions
    if (BCstart[i]=="external"){
      internali[timesTS[selection]==tstart[i]] <- FALSE
    }
    if (BCend[i]=="external"){
      internali[timesTS[selection]==tend[i]] <- FALSE
    }
    #and consider anything interpolated as external
    internali[interpolated[selection]] <- FALSE
    
    #starting filling TS
    if (i==1){
      internal <- internali
    }else{
      internal <- c(internal,internali)
    }
    sections <- c(sections,rep(namesTS[i],length(selection)))
  }
  #get the corresponding sites, and times
  sites <- sitesTS[selec]
  times <- timesTS[selec]
  
  #renumber the surveys
  surveys <- integer(0)
  for (i in unique(sections)){
    surveys <- c(surveys,order(times[sections==i]))
  }
  
  #Filling TS
  TS$d <- as.dist(dTS[selec,selec])
  TS$metadata <- data.frame(sites,sections,surveys,times,internal)
  
  #adding interpolation information if appropriate
  if (sum(interpolated)>0){
    #interpolated <- interpolated[selec]
    intCoef <- intCoef[selec]
    TS$interpolationInfo <- intCoef
  }
  class(TS) <- c("sections","trajectories", "list")
  return(TS)
}


#' @rdname trajectorySections
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param ToInterpolate a matrix with three columns: 1) the positions of ecosystem states A in d; 2) the positions of ecosystem states B in d; 3) an interpolation coefficient (i.e. at what proportion of directed segment AB the interpolate ecosystem state X needs to be).
#' @noRd
interpolateEcolStates <- function(d,ToInterpolate)
{
  if (any(c(ToInterpolate[,3]>1),ToInterpolate[,3]<0)) 
    stop("Interpolation coefficients in ToInterpolate need to be between 0 and 1")
  
  if (any(c(ToInterpolate[,1]>dim(as.matrix(d))[1]),ToInterpolate[,2]>dim(as.matrix(d))[1])) 
    stop("Wrong indexing, ToInterpolate has values superior to the size of the distance matrix")
  
  dInt <- as.matrix(d)
  
  for (i in 1:dim(ToInterpolate)[1]){
    A <- ToInterpolate[i,1]
    B <- ToInterpolate[i,2]
    int <- ToInterpolate[i,3]
    
    AB <- dInt[A,B]
    if (AB == 0){#this is to cover the case where A and B are confounded
      CX <- dInt[A,]
    }else{#otherwise, use trigonometry
      AC <- dInt[A,]
      BC <- dInt[B,]
      AX <- int*AB
      
      #add a small value to 0 AC distances so that trigonometry can still be computed
      AC[which(AC==0)] <- 10^-6
      
      alpha <- (AB^2+AC^2-BC^2)/(2*AC*AB)
      alpha[alpha>=1] <- 1#If I'm not mistaken, this forces the places where triangle inequality is not respected to respected it (MIQUEL I'M NOT SUR THIS IS THE BEST WAY TO DO IT!!!)
      alpha <- acos(alpha)
      CX <- sqrt(AX^2+AC^2-2*AX*AC*cos(alpha))
      CX[A] <- AX
      CX[B] <- (1-int)*AB
    }
    
    dInt <- cbind(dInt,CX)
    dInt <- rbind(dInt,c(CX,0))
    
    rownames(dInt)[dim(dInt)[1]] <- paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
    colnames(dInt)[dim(dInt)[2]] <- paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
  }
  dInt <- as.dist(dInt)
  return(dInt)
}
