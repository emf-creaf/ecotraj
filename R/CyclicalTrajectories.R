#' Functions for Cyclical Ecological Trajectory Analysis
#'
#' The Cyclical extension of Ecological Trajectory Analysis (CETA) aims at allowing ETA to describe ecological trajectories presenting cyclical dynamics such as seasonal or day/night cycles.
#' CETA operates by subdividing a trajectory presenting cyclical dynamics into two types of sub-trajectories of interest: cycles and fixed-dates trajectories
#' 
#' CETA functions achieve one of two goals:
#' 1) Reformating data to analyse either cycles or fixed-dates trajectories. The reformated data can then be fed into existing ETA functions to obtain desired metrics.
#' 2) Providing new metrics relevant to cycles complementing other ETA functions
#' 
#' CETA functions:
#' \itemize{
#' \item{Function \code{cycleBuild} reformates a dataset describing one or more trajectories with cyclical dynamics for the analysis of their cycles.}
#' \item{Function \code{fdtrajBuild} reformates a dataset describing one or more trajectories with cyclical dynamics for the analysis of their fixed-dates trajectories.}
#' \item{Function \code{cycleSmoothness} computes the "smoothness" of the cycles embedded in one or more cyclical trajectories.}
#' \item{Function \code{cycleShift} computes the cyclical shifts (i.e. advances and delays) that can be obtain from one or more cyclical trajectories.}
#' \item{Function \code{trajectorysectionBuild} reformates a dataset describing one or more trajectories into specified trajectory sections. Trajectory sections represent a way to subset trajectories flexibly. Cycles are a particular case of trajectory sections.}
#' \item{Function \code{interpolateEcolStates} compute interpolated ecological states and the new distance matrix associated (used in trajectorysectionBuild).}
#' }
#' 
#' 
#' @encoding UTF-8
#' @name cyclicaltrajectories
#' @aliases cycleBuild fdtrajBuild cycleSmoothness cycleShift trajectorysectionBuild interpolateEcolStates
#' 
#' @details
#' CETA is a little more time-expicit than the rest of ETA. Hence the parameters "surveys" having only an ordinal meaning in other ETA functions is replaced by "times".
#' CETA also distinguishes between times and dates. Times represent linear time whereas dates represent circular time (e.g. the month of year). Dates are circular variables, coming back to zero when reaching their maximum value DurC corresponding to the duration of a cycle.
#' In CETA, dates are by default assumed to be times modulo DurC. This should fit many applications but if this is not the case (i.e. if there is an offset between times and dates), dates can be specified. Dates however need to remain compatible with times and DurC (i.e. (times modulo DurC) - (dates modulo DurC) needs to be a constant).
#' 
#' Trajectory sections can be obtained using trajectorysectionBuild. Trajectory sections allow to cut a longer trajectory into parts for further analyses. Cycles are specical case of trajectory sections.
#' A trajectory section TS(Traj,(tstart, BCstart),(tend, BCend)) is defined by the trajectory (Traj) it is obtained from, by an start and end times (tstart and tend) and start and end boundary conditions (BCstart, BCend).
#' The function trajectorysectionBuild builds trajectory sections as a function of its arguments \code{Traj}, \code{tstart}, \code{tend}, \code{BCstart}, \code{BCend}.
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#'
#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param DurC A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{DurC}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo DurC (see details).
#' @param startdate An optional value indicating at which date the cycles must begin. Must be in the same units as times. Defaults to \code{min(dates)}.
#' @param extBound An optional value, either \code{"end"} or \code{"start"}, indicating whether the start or end of the cycles must be considered "external". Defaults to \code{"end"}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
cycleBuild <- function(d,sites,times,DurC,dates=times%%DurC,startdate=min(dates),extBound="end",minEcolStates=3)
{
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times)|length(sites)!=length(dates))
    stop("The lengths of sites, times, and dates must corespond to the dimension of d")
  
  check=(times%%DurC-(dates%%DurC))%%DurC
  if (any(round(check-check[1],12)!=0))
    stop("provided times and dates are not compatible given cycle duration DurC")
  
  #Goal of the function: reshape its inputs to give them to trajectorysectionBuild
  Traj=integer(0)
  tstart=integer(0)
  tend=integer(0)
  namesCycles=integer(0)
  
  #This loops build Traj, tstart and tend that will be fed into trajectorysectionBuild
  for (i in unique(sites)){
    truc=data.frame(times[sites==i],dates[sites==i])
    colnames(truc)=c("times","dates")
    truc=truc[order(truc$times),]
    
    tstarti=truc$times[1]+((startdate+DurC)-truc$dates[1])%%DurC
    tstarti=seq(from=tstarti,to=max(truc$times),by=DurC)
    
    tendi=tstarti[2:length(tstarti)]
    tstarti=tstarti[1:(length(tstarti)-1)]
    
    #take in account the minimum number of ecological states required
    if (extBound=="end"){
      selec=table(cut(truc$times,c(tstarti,tendi[length(tendi)]),right=F))>=minEcolStates
      
    }else if (extBound=="start"){
      selec=table(cut(truc$times,c(tstarti,tendi[length(tendi)])))>=minEcolStates
      
    }else{
      stop("extBound value invalid, it can only be 'end' or 'start'")
    }
    
    tstarti=tstarti[selec]
    tendi=tendi[selec]
    
    tstart=c(tstart,tstarti)
    tend=c(tend,tendi)
    Traj=c(Traj,rep(i,length(tstarti)))
    
    namesCycles=c(namesCycles,paste(i,paste("C",1:length(tstarti),sep="")))
  }
  
  #End of the loop, we then add the vectors for the boundary conditions
  if (extBound=="end"){
    BCstart=rep("internal",length(tstart))
    BCend=rep("external",length(tend))
  }else if (extBound=="start"){
    BCstart=rep("external",length(tstart))
    BCend=rep("internal",length(tend))
  }
  
  #Now feed this into trajectorysectionBuild
  output=trajectorysectionBuild(d,sites,times,Traj,tstart,tend,BCstart,BCend,namesTS=namesCycles)
  
  #Add dates to the output
  offset=tapply((times%%DurC-dates%%DurC),sites,min)#computing a potential offset between times and dates (the min function here is just to give a singular output, the values should be the same for a given site)
  
  dates=rep(NA,nrow(output$metadata))
  
  for (i in unique(sites)){
    dates[output$metadata$sites==i]=((output$metadata$times[output$metadata$sites==i]%%DurC)-offset[i])%%DurC
  }
  output$metadata=data.frame(output$metadata,dates)
  
  #change the name of the "TrajSec" column in output$metadata to "Cycles"
  colnames(output$metadata)[2]="Cycles"
  
  return(output)
}

#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param DurC A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{DurC}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo DurC (see details).
#' @param fixed_dates An optional vector of dates for which fixed-dates trajectories must be computed. Defaults to \code{unique(dates)}, resulting in returning all possible fixed-dates trajectories.
#' @param names_fixed_dates An optional vector of names associated to each \code{fixed_dates}. Defaults to round(fixed_dates,2).
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a fixed-date trajectory. Fixed-dates trajectories comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 2.
#' @export 
fdtrajBuild <- function (d,sites,times,DurC,dates=times%%DurC,fixed_dates=unique(dates%%DurC),names_fixed_dates=as.character(round(fixed_dates,2)),minEcolStates=2)
{
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times)|length(sites)!=length(dates))
    stop("The lengths of sites, times, and dates must corespond to the dimension of d")
  
  check=(times%%DurC-(dates%%DurC))%%DurC
  if (any(round(check-check[1],12)!=0))
    stop("provided times and dates are not compatible given cycle duration DurC")
  
  if (length(fixed_dates)!=length(names_fixed_dates))
    stop("fixed_dates and names_fixed_dates must have the same length")
  
  #build a convenient "metadata" to find what's needed easily to compute stuff on fixed dates trajectories
  fdT=rep(NA,length(sites))
  metadata=data.frame(sites,fdT,times,dates)
  for (i in unique(sites)){
    for (j in 1:length(fixed_dates)){
      selec=(sites==i)&(dates==fixed_dates[j])
      if (sum(selec)>=minEcolStates){
        metadata$fdT[selec]=paste(i,"fdT",names_fixed_dates[j])
      }
    }
  }
  #remove the lines not belonging to any fixed_dates_trajectories
  selec=is.na(metadata$fdT)==F
  
  d=as.matrix(d)
  d=d[selec,selec]
  d=dist(d)
  
  metadata=metadata[selec,]
  
  #build and return the output
  output=list()
  output$d<-d
  output$metadata<-metadata
  
  return(output)
}

#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param DurC A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{DurC}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo DurC (see details).
#' @param datesCS An optional vector indicating the dates for which a cyclical shift must be computed. Default to unique(dates) resulting in the computation of all possible cyclical shifts.
#' @param centering An optional boolean. Should the cycles be centered before computing cyclical shifts? Defaults to \code{T}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
cycleShift <- function (d,sites,times,DurC,dates=times%%DurC,datesCS=unique(dates),centering=T,minEcolStates=3)#add xCS and DeltaCS at some point to allow more targeted computations!
{
  Output=integer(0)#this will contain the final output
  
  #for each date in datesCS, we compute all possible non-overlapping Cycles starting and ending DurC/2 before and after the date.
  #those cycles will become the cycles for which a cyclical shift will be computed and the reference cycles
  
  for (i in datesCS){
    Cycles=cycleBuild(d,sites,times,dates,DurC,startdate=(i-DurC/2)%%DurC,extBound="end",minEcolStates)
    
    #optional (but advised) centering (only done on complete cycles, the internal only cycles are not needed for this application)
    if (centering==T){
      Cycles$d=centerTrajectories(d=Cycles$d,sites=Cycles$metadata$Cycles,exclude=which(Cycles$metadata$IntExt=="external"))
    }
    
    
    for (j in unique(sites)){#This loop goes through the sites (we only compare cyclical shifts from the same sites)
      siteTag <- which(Cycles$metadata$sites==j)
      dateTag <- which(Cycles$metadata$dates==i)
      dateTag <- intersect(siteTag,dateTag)
      
      AllRefCycles=unique(Cycles$metadata$Cycles[siteTag])
      AllRefCycles=AllRefCycles[1:(length(AllRefCycles)-1)]
      for (k in AllRefCycles){#This loop goes through all the Cycles that will be used as reference.
        CrefTag <- which(Cycles$metadata$Cycles==k)
        CrefTag <- CrefTag[order(Cycles$metadata$times[CrefTag])]#this line re-orders the cycle according to its times to be sure its properly represented in trajectoryProjection below
        drefTag <- dateTag[dateTag%in%CrefTag]#a tag for the reference date
        dateTag <- dateTag[dateTag>max(CrefTag)]#This line ensures that the date (dateTag) which will be compared to the reference are posterior in time
        
        if (length(drefTag)==1&length(dateTag)>0){#a condition to test if the dates we want to compare exist (i.e. do they have associated ecological states?)
          #Note: there is no need to reorder the cycles as cycleBuild always output them in chronological order for a given site
          
          #Find out onto which segment of Cref the ecological states of interest are projected
          projCS=trajectoryProjection(d=Cycles$d,
                                      target=dateTag,
                                      trajectory=CrefTag)
          
          segmentsTag <- cbind(CrefTag[projCS$segment],CrefTag[projCS$segment+1],projCS$relativeSegmentPosition)
          
          #get the "times" of the projection of the ecological states
          timesProj=Cycles$metadata$times[segmentsTag[,1]]+
            ((Cycles$metadata$times[segmentsTag[,2]]-Cycles$metadata$times[segmentsTag[,1]])*segmentsTag[,3])
          
          #get the Cyclical shifts:
          Cyclical_Shift=timesProj-Cycles$metadata$times[drefTag]
          
          #Add some "metadata" to accompany Dt
          timeRef=rep(Cycles$metadata$times[drefTag],length(dateTag))
          timeCS=Cycles$metadata$times[dateTag]
          timescale=timeCS-timeRef
          dateCS=rep(i,length(timeRef))
          site=rep(j,length(timeRef))
          
          PartialOutput=data.frame(site,dateCS,timeRef,timeCS,timescale,Cyclical_Shift)
          Output=rbind(Output,PartialOutput)
          
        }
      }
    }
  }
  return(Output)
}

#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param DurC A value indicating the duration of a cycle. Must be in the same units as times.
#' @param dates An optional vector indicating the dates (< \code{DurC}) corresponding to each ecosystem state. Must be in the same units as times. Defaults to times modulo DurC (see details).
#' @param startdate An optional value indicating at which date the cycles must begin. Must be in the same units as times. Defaults to \code{min(dates)}.
#' @param extBound An optional value, either \code{"end"} or \code{"start"}, indicating whether the start or end of the cycles must be considered "external". Defaults to \code{"end"}.
#' @param minEcolStates An optional integer indicating the minimum number of ecological states to return a cycle. Cycle comprising less ecological states than minEcolStates are discarded and do not appear in the output of the function. Defaults to 3.
#' @export
cycleSmoothness <- function (d,sites,times,DurC,dates=times%%DurC,startdate=min(dates),extBound="end",minEcolStates=3)
{
  Cycles=cycleBuild(d,sites,times,dates,DurC,startdate,extBound,minEcolStates)
  Angles=trajectoryAngles(d,sites,surveys=times)
  
  SC=integer(0)
  
  for (i in unique(sites)){
    Anglesi=Angles[i,]
    Anglesi=Anglesi[1:(length(Anglesi)-3)]
    Anglesi=Anglesi[is.na(Anglesi)==F]
    
    timesi=times[sites==i]
    timesi=timesi[2:(length(timesi)-1)]
    
    bidule=data.frame(timesi,Anglesi)
    bidule=bidule[order(bidule$timesi),]
    
    timesCyclesi=Cycles$metadata$times[Cycles$metadata$sites==i&Cycles$metadata$IntExt=="internal"]
    Cyclesi=Cycles$metadata$Cycles[Cycles$metadata$sites==i&Cycles$metadata$IntExt=="internal"]
    anglesCyclesi=rep(NA,length(Cyclesi))
    truc=data.frame(Cyclesi,timesCyclesi,anglesCyclesi)
    truc=truc[order(truc$timesCyclesi),]
    
    toreplace=truc$timesCyclesi%in%bidule$timesi
    replacewith=bidule$timesi%in%truc$timesCyclesi[toreplace]
    truc$anglesCyclesi [toreplace]=bidule$Anglesi[replacewith]
    
    
    SC=c(SC,360/tapply(truc$anglesCyclesi,truc$Cyclesi,sum))
  }
  return(SC)
}


#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param Traj A vector of length equal to the number of desired trajectory sections indicating the trajectories from which trajectory sections must be build (see details).
#' @param tstart A vector of start times for each of the desired trajectory sections (see details).
#' @param tend A vector of end times for each of the desired trajectory sections (see details).
#' @param BCstart A vector of start boundary conditions (either \code{"internal"} or \code{"external"}) for each of the desired trajectory sections (see details).
#' @param BCend A vector of end boundary conditions (either \code{"internal"} or \code{"external"}) for each of the desired trajectory sections (see details). 
#' @param namesTS An optional vector giving a name for each of the desired trajectory sections (by default trajectory sections are simply numbered). 
#' @export
trajectorysectionBuild <- function(d,sites,times,Traj,tstart,tend,BCstart,BCend,namesTS=1:length(Traj))
{
  if (nrow(as.matrix(d))!=length(sites)|length(sites)!=length(times))
    stop("The lengths of sites and times must corespond to the dimension of d")
  if (any(BCstart%in%c("internal","external")==F))
    stop("BCstart and BCend can only take values 'internal' and 'external'")
  if (any(BCend%in%c("internal","external")==F))
    stop("BCstart and BCend can only take values 'internal' and 'external'")
  
  #those two will contain the sites and times of all ecological states including the interpolated ones
  sitesTS=sites
  timesTS=times
  
  #first, check if there is ecological states to interpolate
  ToInterpolate=integer(0)
  interpolated=rep("N",length(sites))
  for (i in unique(Traj)){
    timesTraj=times[which(sites==i)]
    timesInt=c(tstart[which(Traj==i)],tend[which(Traj==i)])
    
    if (any(is.na(cut(timesInt,range(timesTraj),include.lowest=T)))) 
      stop(paste("tstart and/or tend out of bounds for trajectory",i))
    
    timesInt=timesInt[timesInt%in%timesTraj==F]
    if (length(timesInt)>0){
      for (j in 1:length(timesInt)){
        truc=timesTraj-timesInt[j]
        truc[truc>0]=NA
        A=timesTraj[which(truc==max(truc,na.rm=T))]
        
        truc=timesTraj-timesInt[j]
        truc[truc<0]=NA
        B=timesTraj[which(truc==min(truc,na.rm=T))]
        
        ToInterpolate=rbind(ToInterpolate,c(
          intersect(which(times==A),which(sites==i)),
          intersect(which(times==B),which(sites==i)),
          (timesInt[j]-A)/(B-A)))
      }
      sitesTS=c(sitesTS,rep(i,length(timesInt)))
      timesTS=c(timesTS,timesInt)
      interpolated=c(interpolated,rep("Y",length(timesInt)))
    }
  }
  dTS=d
  if(length(ToInterpolate)>0){
    dTS=interpolateEcolStates(dTS,ToInterpolate)
  }
  dTS=as.matrix(dTS)
  
  #Now build distances matrices and associated tags describing the trajectory sections
  
  #prepare the list that will store everything
  TS=list()
  TS$d=integer(0)
  TS$metadata=integer(0)
  
  #and the element that will go into TS$metadata
  IntExt=integer(0)
  TrajSec=integer(0)
  
  selec=integer(0)
  
  for (i in 1:length(Traj)){
    selection=intersect(
      intersect(which(timesTS>=tstart[i]),
                which(timesTS<=tend[i])),
      which(sitesTS==Traj[i]))
    
    selec=c(selec,selection)
    
    #Find out which are the external ecological states
    IntExti=rep("internal",length(selection))
    #boundary conditions
    if (BCstart[i]=="external"){
      IntExti[timesTS[selection]==tstart[i]]="external"
    }
    if (BCend[i]=="external"){
      IntExti[timesTS[selection]==tend[i]]="external"
    }
    #and consider anything interpolated as external
    IntExti[interpolated[selection]=="Y"]="external"
    
    #starting filling TS
    IntExt=c(IntExt,IntExti)
    TrajSec=c(TrajSec,rep(namesTS[i],length(selection)))
  }
  #get the corresponding sites and times
  sites=sitesTS[selec]
  times=timesTS[selec]
  
  #Filling TS
  TS$d=as.dist(dTS[selec,selec])
  TS$metadata=data.frame(sites,TrajSec,times,IntExt)
  
  return(TS)
}


#MIQUEL, THIS ONE WILL NEED AMENDMENTS, FOR NOW, IT WILL ONLY WORK IF THE TRIANGLE INEQUALITY IS VERIFIED (-ISH), OTHERWISE I DON'T KNOW TOO MUCH...
#' @rdname cyclicaltrajectories
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param ToInterpolate a matrix with three columns: 1) the positions of ecosystem states A in d; 2) the positions of ecosystem states B in d; 3) an interpolation coefficient (i.e. at what proportion of directed segment AB the interpolate ecosystem state X needs to be).
#' @export
interpolateEcolStates <- function(d,ToInterpolate)
{
  if (any(c(ToInterpolate[,3]>1),ToInterpolate[,3]<0)) 
    stop("Interpolation coefficients in ToInterpolate need to be between 0 and 1")
  
  if (any(c(ToInterpolate[,1]>dim(as.matrix(d))[1]),ToInterpolate[,2]>dim(as.matrix(d))[1])) 
    stop("Wrong indexing, ToInterpolate has values superior to the size of the distance matrix")
  
  dInt=as.matrix(d)
  
  for (i in 1:dim(ToInterpolate)[1]){
    A=ToInterpolate[i,1]
    B=ToInterpolate[i,2]
    int=ToInterpolate[i,3]
    
    AB=dInt[A,B]
    if (AB == 0){#this is to cover the case where A and B are confounded
      CX=dInt[A,]
    }else{#otherwise, use trigonometry
      AC=dInt[A,]
      BC=dInt[B,]
      AX=int*AB
      
      #add a small value to 0 AC distances so that trigonometry can still be computed
      AC[which(AC==0)]=10^-6
      
      alpha=(AB^2+AC^2-BC^2)/(2*AC*AB)
      alpha[alpha>=1]=1#If I'm not mistaken, this forces the places where triangle inequality is not respected to respected it (MIQUEL I'M NOT SUR THIS IS THE BEST WAY TO DO IT!!!)
      alpha=acos(alpha)
      CX=sqrt(AX^2+AC^2-2*AX*AC*cos(alpha))
      CX[A]=AX
      CX[B]=(1-int)*AB
    }
    
    dInt=cbind(dInt,CX)
    dInt=rbind(dInt,c(CX,0))
    
    rownames(dInt)[dim(dInt)[1]]=paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
    colnames(dInt)[dim(dInt)[2]]=paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
  }
  dInt=as.dist(dInt)
  return(dInt)
}
