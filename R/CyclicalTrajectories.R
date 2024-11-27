#' Utility functions for Cyclical Ecological Trajectory Analysis
#'
#' @author Nicolas Djeghri, UBO
#' 
#' #' The set following set of utility functions are provided:
#' Three functions are available here:
#' InterpolateEcolStates
#' BuildTrajectorySections
#' BuildCycles


# InterpolateEcolStates: Function for interpolating ecological states:-------------------
#I think of this one more as something internal to the package, not normally advertised to the user
#MIQUEL, THIS ONE WILL NEED AMENDMENTS, FOR NOW, IT WILL ONLY WORK IF THE TRIANGLE INEQUALITY IS VERIFIED, OTHERWISE I DON'T KNOW TOO MUCH...
#
#The idea of the function is to interpolate an ecological state X between ecological states A and B and computing the distances to all other ecological states in the distance matrix (noted C). This is done usinng the cosines law in triangles ABC and AXC. The function takes as input the distance matrix d and another matrix "ToInterpolate" with: first column: ecological states A, second column: ecological states B, third column: an "interpolation coefficient" (i.e. at what proportion of directed segment AB X needs to be). The function returns a new distance matrix that includes the interpolated ecological states.

InterpolateEcolStates <- function(d,ToInterpolate)
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
    AC=dInt[A,]
    BC=dInt[B,]
    AX=int*AB
    alpha=acos((AB^2+AC^2-BC^2)/(2*AC*AB))
    CX=sqrt(AX^2+AC^2-2*AX*AC*cos(alpha))
    CX[A]=AX
    CX[B]=(1-int)*AB
    
    dInt=cbind(dInt,CX)
    dInt=rbind(dInt,c(CX,0))
    
    rownames(dInt)[dim(dInt)[1]]=paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
    colnames(dInt)[dim(dInt)[2]]=paste(ToInterpolate[i,1],"-",ToInterpolate[i,2]," interpolated",sep="")
  }
  dInt=as.dist(dInt)
  return(dInt)
}

#Little test:
As=c(1,3,5)
Bs=c(2,4,6)
Ints=c(0.5,0.7,0.2)
MatrixTest=cbind(As,Bs,Ints)

Newdists=InterpolateEcolStates(ToyDist,MatrixTest)

pcoa=pcoa(Newdists)

plot(pcoa$vectors,las=1,col=c("red","red","blue","blue","green3","green3","black","black","black","black","red","blue","green3"),pch=c(16,16,16,16,16,16,1,1,1,1,17,17,17))
text(pcoa$vectors,rownames(pcoa$vectors),pos=1)


# BuildTrajectorySections: Function to generate trajectory sections-----------------------
#
#This function generates, or more accurately formats, the data so that it describes trajectory sections from existing trajectories and is easy to use for existing ETA functions (mostly, centering will need a specific one though...). The idea is that the user describes which trajectory section he/she wants to that function, and then the functions gives the data tags and so on in the shape that is good for the use of other functions.
#
#This one seems particularly cumbersome to me but I don't really see a way around.
#
#Uses: InterpolateEcolStates
#Note: The function allows to specify overlapping trajectory sections. This generates duplication of ecological states. Analyses such as PCoA should not be done on the outputs of this function but on the input distance matrix instead!!!


BuildTrajectorySections <- function(d,sites,times,Traj,tstart,tend,BCstart,BCend,namesTS)
{
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
    dTS=InterpolateEcolStates(dTS,ToInterpolate)
  }
  dTS=as.matrix(dTS)
  
  #Now build distances matrices and associated tags describing the trajectory sections
  
  #prepare the list that will store everything
  TS=list()
  #containing the complete trajectory sections
  TS$CompleteTS$d=integer(0)
  TS$CompleteTS$sites=integer(0)
  TS$CompleteTS$TrajSec=integer(0)
  TS$CompleteTS$times=integer(0)
  TS$CompleteTS$surveys=integer(0)
  TS$CompleteTS$IntExt=integer(0)
  
  #containing only the internal ecological states
  TS$InternalOnly$d=integer(0)
  TS$InternalOnly$sites=integer(0)
  TS$InternalOnly$TrajSec=integer(0)
  TS$InternalOnly$times=integer(0)
  TS$InternalOnly$surveys=integer(0)
  
  Allselec=integer(0)
  InternalOnlyselec=integer(0)
  
  for (i in 1:length(Traj)){
    selection=intersect(
      intersect(which(timesTS>=tstart[i]),
                which(timesTS<=tend[i])),
      which(sitesTS==Traj[i]))
    
    Allselec=c(Allselec,selection)
    
    #Find out which are the external ecological states
    IntExt=rep("internal",length(selection))
    #boundary conditions
    if (BCstart[i]=="external"){
      IntExt[timesTS[selection]==tstart[i]]="external"
    }
    if (BCend[i]=="external"){
      IntExt[timesTS[selection]==tend[i]]="external"
    }
    #and consider anything interpolated as external
    IntExt[interpolated[selection]=="Y"]="external"
    
    #starting filling TS
    TS$CompleteTS$IntExt=c(TS$CompleteTS$IntExt,IntExt)
    TS$CompleteTS$TrajSec=c(TS$CompleteTS$TrajSec,rep(namesTS[i],length(selection)))
    TS$CompleteTS$surveys=c(TS$CompleteTS$surveys,rank(timesTS[selection]))
    
    InternalOnlyselec=c(InternalOnlyselec,selection[IntExt=="internal"])
    
    TS$InternalOnly$TrajSec=c(TS$InternalOnly$TrajSec,rep(namesTS[i],length(selection[IntExt=="internal"])))
    TS$InternalOnly$surveys=c(TS$InternalOnly$surveys,rank(timesTS[selection[IntExt=="internal"]]))
  }
  #Finish filling TS
  #completeTS
  TS$CompleteTS$d=as.dist(dTS[Allselec,Allselec])
  TS$CompleteTS$sites=sitesTS[Allselec]
  TS$CompleteTS$times=timesTS[Allselec]
  
  #only internal TS
  TS$InternalOnly$d=as.dist(dTS[InternalOnlyselec,InternalOnlyselec])
  TS$InternalOnly$sites=sitesTS[InternalOnlyselec]
  TS$InternalOnly$times=timesTS[InternalOnlyselec]
  
  return(TS)
}

# BuildCycles: Function to build all (or most) possible consecutive, non-overlapping cycles from one or more cyclical trajectory----------------
#
#This function is an easier implementation of BuildTrajectorySections specifically designed for cyclical trajectories
#
#Uses: BuildTrajectorySections

BuildCycles <- function(d,sites,times,dates,startdate,DurC,extBound="end",minEcolStates=0)
{
  if (any((times%%DurC)-(dates%%DurC)!=0))
    stop("provided times and dates are not compatible given cycle duration DurC")
  
  #Goal of the function: reshape its inputs to give them to BuildTrajectorySections
  Traj=integer(0)
  tstart=integer(0)
  tend=integer(0)
  namesCycles=integer(0)
  
  #This loops build Traj, tstart and tend that will be fed into BuildTrajectorySections
  for (i in unique(sites)){
    truc=data.frame(times[sites==i],dates[sites==i])
    colnames(truc)=c("times","dates")
    truc=truc[order(truc$times),]
    
    firsttime=truc$times[min(which(truc$dates>=startdate%%DurC))]
    tstarti=firsttime-(truc$dates[truc$times==firsttime]-startdate%%DurC)
    
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
  
  #Now feed this into BuildTrajectorySections
  output=BuildTrajectorySections(d,sites,times,Traj,tstart,tend,BCstart,BCend,namesTS=namesCycles)
  return(output)
}
