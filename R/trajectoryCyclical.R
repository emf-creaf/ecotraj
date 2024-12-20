#' Functions for Cyclical Ecological Trajectory Analysis
#'
#' The Cyclical extension of Ecological Trajectory Analysis (CETA) aims at allowing ETA to describe ecological trajectories presenting cyclical dynamics such as seasonal or day/night cycles. We call such trajectories "cyclical".
#' CETA operates by subdividing cyclical trajectories into two types of sub-trajectories of interest: cycles and fixed-date trajectories.
#' \itemize{
#' \item{Cycles are sub-trajectories joining the ecological states belonging to the same cycle.}
#' \item{Fixed-date trajectories are sub-trajectories joining the ecological states of the same date in different cycles (e.g. in a multi-annual cyclical trajectory with seasonality, a fixed-date trajectory might join all the ecological states associated with the January months of the different years).}
#' }
#' We recommend reading the vignette on CETA prior to use it.The CETA functions provided here achieve one of two goals:
#' 1) Reformatting data to analyze either cycles or fixed-date trajectories. The reformatted data can then be fed into existing ETA functions to obtain desired metrics (although special care need to be taken with cycles, see details).
#' 2) Providing new metrics relevant to cycles complementing other ETA functions.
#' 
#' CETA functions:
#' \itemize{
#' \item{Function \code{extractCycles} reformats a dataset describing one or more cyclical trajectories for the analysis of their cycles.}
#' \item{Function \code{extractFixedDateTrajectories} reformats a dataset describing one or more cyclical trajectories for the analysis of their fixed-date trajectories.}
#' \item{Function \code{cycleConvexity} computes the "convexity" of the cycles embedded in one or more cyclical trajectories.}
#' \item{Function \code{cycleShift} computes the cyclical shifts (i.e. advances and delays) that can be obtain from one or more cyclical trajectories.}
#' }
#' 
#' @encoding UTF-8
#' @name trajectoryCyclical
#' @aliases extractCycles extractFixedDateTrajectories cycleConvexity cycleShift
#' 
#' @details
#' CETA is a little more time-explicit than the rest of ETA. Hence the parameter \code{surveys} having only an ordinal meaning in other ETA functions is replaced by \code{times}.
#' CETA also distinguishes between times and dates. Times represent linear time whereas dates represent circular time (e.g. the month of year). Dates are circular variables, coming back to zero when reaching their maximum value \code{cycleDuration} corresponding to the duration of a cycle.
#' In CETA, dates are by default assumed to be \code{times} modulo \code{cycleDuration}. This should fit many applications but if this is not the case (i.e. if there is an offset between times and dates), dates can be specified. \code{dates} however need to remain compatible with \code{times} and \code{cycleDuration} (i.e. (times modulo cycleDuration) - (dates modulo cycleDuration) needs to be a constant).
#' 
#' IMPORTANT: Cycles within CETA comprises both \code{"internal"} and \code{"external"} ecological states (see the output of function \code{extractCycles}). This distinction is a solution to what we call the "December-to-January segment problem". Taking the example of a monthly resolved multi-annual time series, a way to make cycles would be to take the set of ecological states representing months from January to December of each year. However, this omits the segment linking December of year Y to January of year Y+1. However, including this segments means having two Janury months in the same cycle.
#' The proposed solution in CETA (in the case of this specific example) is to set the January month of year Y+1 as \code{"external"}. \code{"external"} ecological states MUST BE EXCLUDED from the calculation of some metrics and for some operations within ETA namely:
#' \itemize{
#'  \item{centering where external ecological states must be excluded from computation but included nonetheless in the procedure (see function \code{centerTrajectories})}
#'  \item{cycle convexity (readily handled within the function \code{cycleConvexity}, see below)}
#'  \item{trajectory variability (see the CETA vignette for the use of the function \code{\link{trajectoryVariability}} in the CETA context)}
#'  \item{Visualization through principal coordinate analysis of the cycles. The dedicated function \code{\link{cyclePCoA}} must be preferred over \code{\link{trajectoryPCoA}} (see details below).}
#' }
#' 
#' 
#' As a general rule the outputs of \code{extractCycles} should be used as inputs in other, non-CETA function (e.g. \code{trajectoryDistances}), taking care of removing external ecological states when appropriate. There is two important exceptions to that rule: the functions \code{cycleConvexity} and \code{cycleShift}. Instead, the inputs of these two functions should parallel the inputs of \code{extractCycles} in a given analysis. For \code{cycleConvexity}, this is because smoothness uses angles obtained from the whole cyclical trajectory, and not only the cycles. For \code{cycleShift}, this is because cyclical shifts are not obtained with respect to a particular set of cycles. The function instead compute the most adapted set of cycles to obtain the metric.
#' 
#' Visualization of cycles through principal coordinate analysis is handled by the function \code{cyclePCoA} ADD DETAILS HERE!!!
#' 
#' Note: Function \code{cycleShift} is computation intensive for large datasets, it may not execute immediately.
#' 
#' Further information and detailed examples of the use of CETA functions can be found in the associated vignette.
#' 
#' @return 
#' Function \code{extractCycles} returns the base information needed to describe cycles. Its outputs are meant to be used as inputs for other ETA functions in order to obtain desired metrics. Importantly, within cycles, ecological states can be considered \code{"internal"} or \code{"external"}. Some operations and metrics within ETA use all ecological states whereas others use only \code{"internal"} ones (see details). Function \code{extractCycles} returns a list containing:
#' \itemize{
#'  \item{\code{d}: an object of class \code{\link{dist}}, the new distance matrix describing the cycles. To take in account ecological states that are both the end of a cycle and the start of another,\code{d} contains duplications. As compared to the input matrix, \code{d} may present deletions of ecological states that do not belong to any cycles (e.g. due to \code{minEcolStates}))}
#'  \item{\code{metadata}: an object of class \code{\link{data.frame}} describing the ecological states in \code{d} with columns:
#'    \itemize{
#'      \item{\code{sites}: the sites associated to  each ecological states.}
#'      \item{\code{Cycles}: the names of the cycle each ecological states belongs to. The cycle name is built by combining the site name with C1, C2, C3... in chronological order.}
#'      \item{\code{times}: the times associated to each ecological states.}
#'      \item{\code{internal}: a boolean vector with \code{TRUE} indicating "internal" ecological states whereas \code{FALSE} indicates "external" ecological states. This has important implications for the use of \code{extractCycles} outputs (see details).}
#'      \item{\code{dates}: the dates associated to each ecological states.}
#'      }
#'    }
#'  \item{\code{interpolationInfo}: an output that only appear if ecological states have been interpolated. It is used by plotting functions (see \code{cyclePCoA}) but is not intended to be of interest to the end user.}
#' }
#' 
#' Function \code{extractFixedDateTrajectories} returns the base information needed to describe fixed-date trajectories. Its outputs are meant to be used as inputs for other ETA functions in order to obtain desired metrics. Unlike cycles,fixed-date trajectories do not include \code{"internal"} or \code{"external"} ecological states and can be treated as any trajectories. Function \code{extractFixedDateTrajectories} returns a list containing:
#' \itemize{
#'  \item{\code{d}: an object of class \code{\link{dist}}, the new distance matrix describing the fixed-date trajectories. As compared to the input matrix, \code{d} may present deletions of ecological states that do not belong to any fixed-date trajectories (e.g. due to \code{minEcolStates}))}
#'  \item{\code{metadata}: an object of class \code{\link{data.frame}} describing the ecological states in \code{d} with columns:
#'    \itemize{
#'      \item{\code{sites}: the sites to  each ecological states.}
#'      \item{\code{fdT}: the names of the fixed-date trajectory each ecological states belongs to. The fixed-date trajectory name is built by combining the site name with "fdT" and the name of the fixed date (from \code{namesFixedDate}).}
#'      \item{\code{times}: the times associated to each ecological states.}
#'      \item{\code{dates}: the dates associated to each ecological states.}
#'      }
#'    }
#' }
#' 
#' Function \code{cycleConvexity} returns the a vector containing values between 0 and 1 describing the convexity of cycles. Importantly, outputs of \code{extractCycles} should not be used as inputs for \code{cycleConvexity} (see details).
#'
#' Function \code{cycleShift} returns an object of class \code{\link{data.frame}} describing cyclical shifts (i.e. advances and delays). Importantly, outputs of \code{extractCycles} should not be used as inputs for \code{cycleShift} (see details). The columns of the \code{\link{data.frame}} are:
#' \itemize{
#'      \item{\code{site}: the site for which each cycle shift has been computed.}
#'      \item{\code{dateCS}: the date for which a cycle shift has been computed.}
#'      \item{\code{timeCS}: the time of the ecological state for which a cycle shift has been computed (i.e. the time associated to the projected ecological state).}
#'      \item{\code{timeRef}: the time associated to the reference ecological state.}
#'      \item{\code{timeScale}: the time difference between the reference and the projected ecological state.}
#'      \item{\code{cyclicalShift}: the cyclical shift computed (an advance if positive, a delay if negative) in the same units as the times input.}
#' }
#' 
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#'
#' @rdname trajectoryCyclical
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param cycleDuration A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{cycleDuration}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo cycleDuration (see details).
#' @param startdate An optional value indicating at which date the cycles must begin. Must be in the same units as times. Defaults to \code{min(dates)}.
#' @param externalBoundary An optional string, either \code{"end"} or \code{"start"}, indicating whether the start or end of the cycles must be considered "external". Defaults to \code{"end"}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
extractCycles <- function(d,sites,times,cycleDuration,dates=times%%cycleDuration,startdate=min(dates),externalBoundary="end",minEcolStates=3)
{
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times)|length(sites)!=length(dates))
    stop("The lengths of sites, times, and dates must corespond to the dimension of d")
  
  check <- (times%%cycleDuration-(dates%%cycleDuration))%%cycleDuration
  if (any(round(check-check[1],12)!=0))
    stop("provided times and dates are not compatible given cycleDuration")
  
  #Goal of the function: reshape its inputs to give them to extractTrajectorySections
  Traj <- integer(0)
  tstart <- integer(0)
  tend <- integer(0)
  namesCycles <- integer(0)
  
  #This loops build Traj, tstart and tend that will be fed into extractTrajectorySections
  for (i in unique(sites)){
    truc <- data.frame(times[sites==i],dates[sites==i])
    colnames(truc) <- c("times","dates")
    truc <- truc[order(truc$times),]
    
    tstarti <- truc$times[1]+((startdate+cycleDuration)-truc$dates[1])%%cycleDuration
    tstarti <- seq(from=tstarti,to=max(truc$times),by=cycleDuration)
    
    tendi <- tstarti[2:length(tstarti)]
    tstarti <- tstarti[1:(length(tstarti)-1)]
    
    #take in account the minimum number of ecological states required
    if (externalBoundary=="end"){
      selec <- table(cut(truc$times,c(tstarti,tendi[length(tendi)]),right=F))>=minEcolStates
      
    }else if (externalBoundary=="start"){
      selec <- table(cut(truc$times,c(tstarti,tendi[length(tendi)])))>=minEcolStates
      
    }else{
      stop("externalBoundary string invalid, it can only be 'end' or 'start'")
    }
    
    tstarti <- tstarti[selec]
    tendi <- tendi[selec]
    
    tstart <- c(tstart,tstarti)
    tend <- c(tend,tendi)
    Traj <- c(Traj,rep(i,length(tstarti)))
    
    namesCycles <- c(namesCycles,paste(i,paste("C",1:length(tstarti),sep="")))
  }
  
  #End of the loop, we then add the vectors for the boundary conditions
  if (externalBoundary=="end"){
    BCstart <- rep("internal",length(tstart))
    BCend <- rep("external",length(tend))
  }else if (externalBoundary=="start"){
    BCstart <- rep("external",length(tstart))
    BCend <- rep("internal",length(tend))
  }
  
  #Now feed this into extractTrajectorySections
  output <- extractTrajectorySections(d,sites,times,Traj,tstart,tend,BCstart,BCend,namesTS=namesCycles)
  
  #Add dates to the output
  offset <- tapply((times%%cycleDuration-dates%%cycleDuration),sites,min)#computing a potential offset between times and dates (the min function here is just to give a singular output, the values should be the same for a given site)
  
  dates <- rep(NA,nrow(output$metadata))
  
  for (i in unique(sites)){
    dates[output$metadata$sites==i] <- ((output$metadata$times[output$metadata$sites==i]%%cycleDuration)-offset[i])%%cycleDuration
  }
  output$metadata <- data.frame(output$metadata,dates)
  
  #change the name of the "TrajSec" column in output$metadata to "Cycles"
  colnames(output$metadata)[2] <- "Cycles"
  
  return(output)
}

#' @rdname trajectoryCyclical
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param cycleDuration A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{cycleDuration}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo cycleDuration (see details).
#' @param fixedDate An optional vector of dates for which fixed-date trajectories must be computed. Defaults to \code{unique(dates)}, resulting in returning all possible fixed-date trajectories.
#' @param namesFixedDate An optional vector of names associated to each \code{fixedDate}. Defaults to \code{round(fixedDate,2)}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a fixed-date trajectory. Fixed-date trajectories comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 2.
#' @export 
extractFixedDateTrajectories <- function (d,sites,times,cycleDuration,dates=times%%cycleDuration,fixedDate=sort(unique(dates%%cycleDuration)),namesFixedDate=as.character(round(fixedDate,2)),minEcolStates=2)
{
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times)|length(sites)!=length(dates))
    stop("The lengths of sites, times, and dates must corespond to the dimension of d")
  
  check <- (times%%cycleDuration-(dates%%cycleDuration))%%cycleDuration
  if (any(round(check-check[1],12)!=0))
    stop("provided times and dates are not compatible given cycleDuration")
  
  if (length(fixedDate)!=length(namesFixedDate))
    stop("fixedDate and namesFixedDate must have the same length")
  
  #build a convenient "metadata" to find what's needed easily to compute stuff on fixed dates trajectories
  fdT <- rep(NA,length(sites))
  metadata <- data.frame(sites,fdT,times,dates)
  for (i in unique(sites)){
    for (j in 1:length(fixedDate)){
      selec <- (sites==i)&(dates==fixedDate[j])
      if (sum(selec)>=minEcolStates){
        metadata$fdT[selec] <- paste(i,"fdT",namesFixedDate[j])
      }
    }
  }
  #remove the lines not belonging to any fixedDate_trajectories
  selec <- is.na(metadata$fdT)==F
  
  d <- as.matrix(d)
  d <- d[selec,selec]
  d <- as.dist(d)
  
  metadata <- metadata[selec,]
  
  #build and return the output
  output <- list()
  output$d <- d
  output$metadata <- metadata
  
  return(output)
}

#' @rdname trajectoryCyclical
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param cycleDuration A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{cycleDuration}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo cycleDuration (see details).
#' @param startdate An optional value indicating at which date the cycles must begin. Must be in the same units as times. Defaults to \code{min(dates)}.
#' @param externalBoundary An optional string, either \code{"end"} or \code{"start"}, indicating whether the start or end of the cycles must be considered "external". Defaults to \code{"end"}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
cycleConvexity <- function (d,sites,times,cycleDuration,dates=times%%cycleDuration,startdate=min(dates),externalBoundary="end",minEcolStates=3)
{
  Cycles <- extractCycles(d=d,sites=sites,times=times,dates=dates,cycleDuration=cycleDuration,startdate=startdate,externalBoundary=externalBoundary,minEcolStates=minEcolStates)
  Angles <- trajectoryAngles(d,sites,surveys=times)
  
  Convexity <- integer(0)
  
  for (i in unique(sites)){
    Anglesi <- Angles[i,]
    Anglesi <- Anglesi[1:(length(Anglesi)-3)]
    Anglesi <- Anglesi[is.na(Anglesi)==F]
    
    timesi <- times[sites==i]
    timesi <- timesi[2:(length(timesi)-1)]
    
    bidule <- data.frame(timesi,Anglesi)
    bidule <- bidule[order(bidule$timesi),]
    
    timesCyclesi <- Cycles$metadata$times[Cycles$metadata$sites==i&Cycles$metadata$internal]
    Cyclesi <- Cycles$metadata$Cycles[Cycles$metadata$sites==i&Cycles$metadata$internal]
    anglesCyclesi <- rep(NA,length(Cyclesi))
    truc <- data.frame(Cyclesi,timesCyclesi,anglesCyclesi)
    truc <- truc[order(truc$timesCyclesi),]
    
    toreplace <- truc$timesCyclesi%in%bidule$timesi
    replacewith <- bidule$timesi%in%truc$timesCyclesi[toreplace]
    truc$anglesCyclesi[toreplace] <- bidule$Anglesi[replacewith]
    
    
    Convexity <- c(Convexity,360/tapply(truc$anglesCyclesi,truc$Cyclesi,sum))
  }
  Convexity <- Convexity[unique(Cycles$metadata$Cycles)]
  return(Convexity)
}

#' @rdname trajectoryCyclical
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param cycleDuration A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{cycleDuration}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo cycleDuration (see details).
#' @param datesCS An optional vector indicating the dates for which a cyclical shift must be computed. Default to \code{unique(dates)} resulting in the computation of all possible cyclical shifts.
#' @param centering An optional boolean. Should the cycles be centered before computing cyclical shifts? Defaults to \code{TRUE}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
cycleShift <- function (d,sites,times,cycleDuration,dates=times%%cycleDuration,datesCS=sort(unique(dates%%cycleDuration)),centering=TRUE,minEcolStates=3)#add xCS and DeltaCS at some point to allow more targeted computations!
{
  Output <- integer(0)#this will contain the final output
  
  #for each date in datesCS, we compute all possible non-overlapping Cycles starting and ending cycleDuration/2 before and after the date.
  #those cycles will become the cycles for which a cyclical shift will be computed and the reference cycles
  
  for (i in datesCS){
    Cycles <- extractCycles(d=d,sites=sites,times=times,dates=dates,cycleDuration=cycleDuration,startdate=(i-cycleDuration/2)%%cycleDuration,externalBoundary="end",minEcolStates)
    
    #optional (but advised) centering (only done on complete cycles, the internal only cycles are not needed for this application)
    if (centering==T){
      Cycles$d <- centerTrajectories(d=Cycles$d,sites=Cycles$metadata$Cycles,exclude=which(Cycles$metadata$internal==FALSE))
    }
    
    
    for (j in unique(sites)){#This loop goes through the sites (we only compare cyclical shifts from the same sites)
      siteTag <- which(Cycles$metadata$sites==j)
      dateTag <- which(Cycles$metadata$dates==i)
      dateTag <- intersect(siteTag,dateTag)
      
      AllRefCycles <- unique(Cycles$metadata$Cycles[siteTag])
      AllRefCycles <- AllRefCycles[1:(length(AllRefCycles)-1)]
      for (k in AllRefCycles){#This loop goes through all the Cycles that will be used as reference.
        CrefTag <- which(Cycles$metadata$Cycles==k)
        CrefTag <- CrefTag[order(Cycles$metadata$times[CrefTag])]#this line re-orders the cycle according to its times to be sure its properly represented in trajectoryProjection below
        drefTag <- dateTag[dateTag%in%CrefTag]#a tag for the reference date
        dateTag <- dateTag[dateTag>max(CrefTag)]#This line ensures that the date (dateTag) which will be compared to the reference are posterior in time
        
        if (length(drefTag)==1&length(dateTag)>0){#a condition to test if the dates we want to compare exist (i.e. do they have associated ecological states?)
          #Note: there is no need to reorder the cycles as extractCycles always output them in chronological order for a given site
          
          #Find out onto which segment of Cref the ecological states of interest are projected
          projCS <- trajectoryProjection(d=Cycles$d,
                                      target=dateTag,
                                      trajectory=CrefTag)
          
          segmentsTag <- cbind(CrefTag[projCS$segment],CrefTag[projCS$segment+1],projCS$relativeSegmentPosition)
          
          #get the "times" of the projection of the ecological states
          timesProj <- Cycles$metadata$times[segmentsTag[,1]]+
            ((Cycles$metadata$times[segmentsTag[,2]]-Cycles$metadata$times[segmentsTag[,1]])*segmentsTag[,3])
          
          #get the Cyclical shifts:
          cyclicalShift <- timesProj-Cycles$metadata$times[drefTag]
          
          #Add some "metadata" to accompany Dt
          timeRef <- rep(Cycles$metadata$times[drefTag],length(dateTag))
          timeCS <- Cycles$metadata$times[dateTag]
          timeScale <- timeCS-timeRef
          dateCS <- rep(i,length(timeRef))
          site <- rep(j,length(timeRef))
          
          PartialOutput <- data.frame(site,dateCS,timeCS,timeRef,timeScale,cyclicalShift)
          Output <- rbind(Output,PartialOutput)
        }
      }
    }
  }
  return(Output)
}

