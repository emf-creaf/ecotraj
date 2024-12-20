#' Cyclical trajectory plots
#' 
#' Set of plotting functions for Cyclical Ecological Trajectory Analysis:
#' 
#' \itemize{
#'  \item{Function \code{cyclePCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#'  \item{Function \code{fixedDateTrajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#' }
#' 
#' @encoding UTF-8
#' @name trajectoryCyclicalPlots
#' @aliases cyclePCoA fixedDateTrajectoryPCoA
#' 
#' CETA plotting functions:
#' \itemize{
#'  \item{Function \code{cyclePCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws cycles in the ordination scatterplot. Sister function of \code{trajectoryPCoA} adapted to cycles.}
#'  \item{Function \code{fixedDateTrajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws fixed-date trajectories in the ordination scatterplot. Sister function of \code{trajectoryPCoA} adapted to fixed date-trajectories.}
#' }
#' 
#' @details
#' The functions \code{cyclePCoA} and \code{fixedDateTrajectoryPCoA} give adapted graphical representation of cycles and fixed-date trajectories respectively using principal coordinate analysis (PCoA, see \code{\link{cmdscale}}).  
#' Function \code{cyclePCoA} handles external and potential interpolated ecological states so that they are correctly taken in account in PCoA (i.e. avoiding duplication, and reducing the influence of interpolated ecological states as much as possible). In case of centered cycles, the influence of these point will grow as they will not correspond to duplications anymore.
#' In case of centered cycles, the intended use is to set the parameter \code{centered} to \code{TRUE}.  
#' Function \code{fixedDateTrajectoryPCoA} has no flexibility for the colors as it defaults to a circular color palette adapted to representation of cyclical processes.
#' 
#' 
#' @return 
#' Functions \code{cyclePCoA} and \code{fixedDateTrajectoryPCoA} return the results of calling of \code{\link{cmdscale}}.
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' 
#' @seealso \code{\link{trajectoryCyclical}}, \code{\link{cmdscale}}
#' 
#' 
#' @rdname trajectoryCyclicalPlots
#' @param inputCyclePCoA The full output of function \code{\link{extractCycles}}.
#' @param centered Boolean. Have the cycles been centered? Default to FALSE.
#' @param sites.colors The colors applied to the different sites. The cycles will be distinguished (old to recent) by increasingly lighter tones of the provided colors.
#' @param print.names A boolean flag to indicate whether the names of cycles or fixed-date trajectories should be printed.
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
cyclePCoA <- function (inputCyclePCoA,centered=FALSE,sites.colors=NULL,print.names=FALSE,axes=c(1,2),...){
  if (is.list(inputCyclePCoA)==FALSE)
    stop ("cyclePCoA takes as main input the whole output of function extractCycles")
  if (is.null(inputCyclePCoA$metadata$Cycles))
    stop ("cyclePCoA takes as main input the whole output of function extractCycles")
  
  if (is.null(sites.colors)==FALSE){
    if (length(sites.colors)!=length(unique(inputCyclePCoA$metadata$sites)))
      stop ("sites.colors must have the same number of colors as the number of sites")
    names(sites.colors) <- unique(inputCyclePCoA$metadata$sites)
  }
  
  #first isolate the subset of the data on which the PCoA must be performed
  if (centered){
    D <- inputCyclePCoA$d
    PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
    
    metadataD <- inputCyclePCoA$metadata
    
  }else{
    selec <- integer(0)
    #this loop will first isolate the non-duplicated external ecological states
    #then it will add the internal ecological states (and removing possible overlap)
    for (i in unique(inputCyclePCoA$metadata$sites)){
      sitei <- inputCyclePCoA$metadata$sites==i
      timesi <- inputCyclePCoA$metadata$times[sitei]
      inti <- (inputCyclePCoA$metadata$internal)[sitei]
      
      #In there we will have the non-duplicated external ecological states
      nonDuplExt <- table(timesi[inti==FALSE])
      nonDuplExt <-as.numeric((names(nonDuplExt[nonDuplExt==1])))
      
      #Here we ensure the external ecological states do not correspond to already existing internal ecological states
      TimesExtToKeepi <- setdiff(nonDuplExt,timesi[inti])
      selec <- c(selec,which(inputCyclePCoA$metadata$times%in%TimesExtToKeepi
                             &
                             inputCyclePCoA$metadata$sites==i))
    }
    selec <- sort(c(selec,which(inputCyclePCoA$metadata$internal)))#and add all internal ecological states
    
    D <- as.dist(as.matrix(inputCyclePCoA$d)[selec,selec])
    
    PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
    
    metadataD <- inputCyclePCoA$metadata[selec,]
  }
  
  #Plotting
  x <- PCoA$points[,axes[1]]
  y <- PCoA$points[,axes[2]]
  
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*PCoA$eig[axes[1]]/sum(PCoA$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*PCoA$eig[axes[2]]/sum(PCoA$eig)),"%)"))
  
  for (i in unique(metadataD$sites)){
    sitei <- metadataD$sites==i
    cyclesi <- unique(metadataD$Cycles[sitei])
    
    #add a (simple) color ramp for cycle
    if (is.null(sites.colors)==FALSE){
      colorCycles <- rep(sites.colors[i],length(cyclesi))
      colorCycles <- rgb(t(col2rgb(colorCycles)/255)+t(c(1,1,1)-col2rgb(colorCycles)/255)*seq(0,0.5,length.out=length(cyclesi)))
      names(colorCycles) <- cyclesi
    }else{
      colorCycles <- rep("black",length(cyclesi))
      colorCycles <- rgb(t(col2rgb(colorCycles)/255)+t(c(1,1,1)-col2rgb(colorCycles)/255)*seq(0,0.5,length.out=length(cyclesi)))
      names(colorCycles) <- cyclesi
    }
    
    for (j in cyclesi){
      cyclejinput <- inputCyclePCoA$metadata$Cycles==j
      
      timesjinput <- inputCyclePCoA$metadata$times[cyclejinput]
      
      selec <- metadataD$times%in%timesjinput&metadataD$sites==i
      
      timesj <- metadataD$times[selec]
      xarrows <- x[selec][order(timesj)]
      yarrows <- y[selec][order(timesj)]
      
      #main arrows
      arrows(x0=xarrows[1:(length(xarrows)-1)],y0=yarrows[1:(length(yarrows)-1)],
             x1=xarrows[2:length(xarrows)],y1=yarrows[2:length(yarrows)],
             col=colorCycles[j],...)
      
      #potentially print cycle names
      if (print.names==TRUE){
        text(x=xarrows[1],y=yarrows[1],j,col=colorCycles[j])
      }
      #add the interpolated ecological states if any
      if (is.null(inputCyclePCoA$interpolationInfo)==FALSE){
        #addition may be made at the start and/or at the end of the cycle
        #for the start:
        if ((min(timesjinput)%in%timesj)==FALSE){
          selecInt <- inputCyclePCoA$metadata$times==min(timesjinput)&inputCyclePCoA$metadata$Cycles==j
          IntCoef <- inputCyclePCoA$interpolationInfo[selecInt]
          
          timeprevious <- max(metadataD$times[metadataD$times<min(timesjinput)&metadataD$sites==i])
          EcolStatePrevious <- metadataD$times==timeprevious&metadataD$sites==i
          
          xprevious <- x[EcolStatePrevious]
          yprevious <- y[EcolStatePrevious]
          
          x1Int <- xarrows[1]
          y1Int <- yarrows[1]
          
          x0Int <- xprevious+(x1Int-xprevious)*IntCoef
          y0Int <- yprevious+(y1Int-yprevious)*IntCoef
          
          arrows(x0=x0Int,y0=y0Int,x1=x1Int,y1=y1Int,col=colorCycles[j],...)
        }
        #for the end:
        if ((max(timesjinput)%in%timesj)==FALSE){
          selecInt <- inputCyclePCoA$metadata$times==max(timesjinput)&inputCyclePCoA$metadata$Cycles==j
          IntCoef <- inputCyclePCoA$interpolationInfo[selecInt]
          
          timenext <- min(metadataD$times[metadataD$times>max(timesjinput)&metadataD$sites==i])
          EcolStateNext <- metadataD$times==timenext&metadataD$sites==i
          
          xnext <- x[EcolStateNext]
          ynext <- y[EcolStateNext]
          
          x0Int <- xarrows[length(xarrows)]
          y0Int <- yarrows[length(yarrows)]
          
          x1Int <- x0Int+(xnext-x0Int)*IntCoef
          y1Int <- y0Int+(ynext-y0Int)*IntCoef
          
          segments(x0=x0Int,y0=y0Int,x1=x1Int,y1=y1Int,col=colorCycles[j],...)
        }
      }
    }
  }
  invisible(PCoA)
}


#' @rdname trajectoryCyclicalPlots
#' @param inputFixedDateTrajectoryPCoA The full output of function \code{\link{extractFixedDateTrajectories}}.
#' @param sites.lty The line type for the different sites (see \code{\link{par}}, \code{"lty"}). The fixed-date trajectories will be distinguished by a default circular color palette.
#' @param print.names A boolean flag to indicate whether the names of cycles or fixed-date trajectories should be printed.
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
fixedDateTrajectoryPCoA <- function (inputFixedDateTrajectoryPCoA,sites.lty=NULL,print.names=FALSE,axes=c(1,2),...){
  if (is.list(inputFixedDateTrajectoryPCoA)==FALSE)
    stop ("cyclePCoA takes as main input the whole output of function extractFixedDateTrajectories")
  if (is.null(inputFixedDateTrajectoryPCoA$metadata$fdT))
    stop ("cyclePCoA takes as main input the whole output of function extractFixedDateTrajectories")
  
  if (is.null(sites.lty)==FALSE){
    if (length(sites.lty)!=length(unique(inputFixedDateTrajectoryPCoA$metadata$sites)))
      stop ("sites.lty must have the same number of indices as the number of sites")
    names(sites.lty) <- unique(inputFixedDateTrajectoryPCoA$metadata$sites)
  }else{
    sites.lty <- rep(1,length(unique(inputFixedDateTrajectoryPCoA$metadata$sites)))
    names(sites.lty) <- unique(inputFixedDateTrajectoryPCoA$metadata$sites)
  }
  
  D <- inputFixedDateTrajectoryPCoA$d
  PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
  metadataD <- inputFixedDateTrajectoryPCoA$metadata
  
  #Plotting
  x <- PCoA$points[,axes[1]]
  y <- PCoA$points[,axes[2]]
  
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*PCoA$eig[axes[1]]/sum(PCoA$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*PCoA$eig[axes[2]]/sum(PCoA$eig)),"%)"))
  
  for (i in unique(metadataD$sites)){
    sitei <- metadataD$sites==i
    fdTi <- unique(metadataD$fdT[sitei])
    
    colorsfdTi <- cmocean::cmocean("phase")(length(fdTi))
    names(colorsfdTi) <- fdTi
    
    for (j in fdTi){
      selec <- metadataD$fdT==j
      timesj <- metadataD$times[selec]
      
      xarrows <- x[selec][order(timesj)]
      yarrows <- y[selec][order(timesj)]
      
      arrows(x0=xarrows[1:(length(xarrows)-1)],y0=yarrows[1:(length(yarrows)-1)],
             x1=xarrows[2:length(xarrows)],y1=yarrows[2:length(yarrows)],
             col=colorsfdTi[j],lty=sites.lty[i],...)
      
      if (print.names==TRUE){
        text(x=xarrows[1],y=yarrows[1],j,col=colorsfdTi[j])
      }
    }
  }
  invisible(PCoA)
}

