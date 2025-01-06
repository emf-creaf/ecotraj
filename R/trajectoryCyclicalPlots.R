#' Cyclical trajectory plots
#' 
#' Plotting functions for Cyclical Ecological Trajectory Analysis:
#' 
#' \itemize{
#'  \item{Function \code{cyclePCoA} removes unwanted points (see details) and performs principal coordinates analysis (\code{\link{cmdscale}}) and draws cycles in the ordination scatterplot.}
#'  \item{Function \code{fixedDateTrajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws fixed date trajectories in the ordination scatterplot.}
#' }
#' 
#' @encoding UTF-8
#' @name trajectoryCyclicalPlots
#' @aliases cyclePCoA fixedDateTrajectoryPCoA customCircularPalette
#' 
#' @details
#' The functions \code{cyclePCoA} and \code{fixedDateTrajectoryPCoA} give adapted graphical representation of cycles and fixed-date trajectories respectively using principal coordinate analysis (PCoA, see \code{\link{cmdscale}}).  
#' Function \code{cyclePCoA} handles external and potential interpolated ecological states so that they are correctly taken in account in PCoA (i.e. avoiding duplication, and reducing the influence of interpolated ecological states as much as possible). In case of centered cycles, the influence of these ecological states will grow as they will not correspond to duplications anymore.
#' In case of centered cycles, the intended use is to set the parameter \code{centered} to \code{TRUE}.  
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
#' @rdname trajectoryCyclicalPlots
#' @param x The full output of function \code{\link{extractCycles}} or \code{\link{extractFixedDateTrajectories}} as appropriate, an object of class \code{\link{trajectories}}.
#' @param centered Boolean. Have the cycles been centered? Default to FALSE.
#' @param sites.colors The colors applied to the different sites. The cycles will be distinguished (old to recent) by increasingly lighter tones of the provided colors.
#' @param print.names A boolean flag to indicate whether the names of cycles or fixed-date trajectories should be printed.
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
cyclePCoA <- function (x, 
                       centered=FALSE,
                       sites.colors=NULL,
                       print.names=FALSE,
                       axes=c(1,2), ...)
{
  if (!inherits(x, "cycles"))
    stop ("cyclePCoA takes as main input the whole output of function extractCycles")

  if (is.null(sites.colors)==FALSE){
    if (length(sites.colors)!=length(unique(x$metadata$sites)))
      stop ("sites.colors must have the same number of colors as the number of sites")
    names(sites.colors) <- unique(x$metadata$sites)
  }
  
  #first isolate the subset of the data on which the PCoA must be performed
  if (centered){
    D <- x$d
    PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
    
    metadataD <- x$metadata
    
  } else{
    selec <- integer(0)
    #this loop will first isolate the non-duplicated external ecological states
    #then it will add the internal ecological states (and removing possible overlap)
    for (i in unique(x$metadata$sites)){
      sitei <- x$metadata$sites==i
      timesi <- x$metadata$times[sitei]
      inti <- (x$metadata$internal)[sitei]
      
      #In there we will have the non-duplicated external ecological states
      nonDuplExt <- table(timesi[inti==FALSE])
      nonDuplExt <-as.numeric((names(nonDuplExt[nonDuplExt==1])))
      
      #Here we ensure the external ecological states do not correspond to already existing internal ecological states
      TimesExtToKeepi <- setdiff(nonDuplExt,timesi[inti])
      selec <- c(selec,which(x$metadata$times%in%TimesExtToKeepi
                             &
                             x$metadata$sites==i))
    }
    selec <- sort(c(selec,which(x$metadata$internal)))#and add all internal ecological states
    
    D <- as.dist(as.matrix(x$d)[selec,selec])
    
    PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
    
    metadataD <- x$metadata[selec,]
  }
  
  #Plotting
  xp <- PCoA$points[,axes[1]]
  yp <- PCoA$points[,axes[2]]
  
  plot(xp,yp, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*PCoA$eig[axes[1]]/sum(PCoA$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*PCoA$eig[axes[2]]/sum(PCoA$eig)),"%)"))
  
  for (i in 1:length(unique(metadataD$sites))){
    sitei <- metadataD$sites==metadataD$sites[i]
    cyclesi <- unique(metadataD$cycles[sitei])
    
    #add a (simple) color ramp for cycle
    if (!is.null(sites.colors)){
      colorCycles <- rep(sites.colors[i],length(cyclesi))
      colorCycles <- rgb(t(col2rgb(colorCycles)/255)+t(c(1,1,1)-col2rgb(colorCycles)/255)*seq(0,0.5,length.out=length(cyclesi)))
      names(colorCycles) <- cyclesi
    }else{
      colorCycles <- rep("black",length(cyclesi))
      colorCycles <- rgb(t(col2rgb(colorCycles)/255)+t(c(1,1,1)-col2rgb(colorCycles)/255)*seq(0,0.5,length.out=length(cyclesi)))
      names(colorCycles) <- cyclesi
    }
    
    for (j in cyclesi){
      cyclejinput <- x$metadata$cycles==j
      
      timesjinput <- x$metadata$times[cyclejinput]
      
      selec <- (metadataD$times%in%timesjinput) & sitei
      
      timesj <- metadataD$times[selec]
      xarrows <- xp[selec][order(timesj)]
      yarrows <- yp[selec][order(timesj)]
      
      #main arrows
      arrows(x0=xarrows[1:(length(xarrows)-1)],y0=yarrows[1:(length(yarrows)-1)],
             x1=xarrows[2:length(xarrows)],y1=yarrows[2:length(yarrows)],
             col=colorCycles[j],...)
      
      #potentially print cycle names
      if (print.names==TRUE){
        text(x=xarrows[1],y=yarrows[1],j,col=colorCycles[j])
      }
      #add the interpolated ecological states if any
      if (is.null(x$interpolationInfo)==FALSE){
        #addition may be made at the start and/or at the end of the cycle
        #for the start:
        if ((min(timesjinput)%in%timesj)==FALSE){
          selecInt <- x$metadata$times==min(timesjinput)&x$metadata$cycles==j
          IntCoef <- x$interpolationInfo[selecInt]
          
          timeprevious <- max(metadataD$times[(metadataD$times<min(timesjinput)) & sitei])
          EcolStatePrevious <- (metadataD$times==timeprevious) & sitei
          
          xprevious <- xp[EcolStatePrevious]
          yprevious <- yp[EcolStatePrevious]
          
          x1Int <- xarrows[1]
          y1Int <- yarrows[1]
          
          x0Int <- xprevious+(x1Int-xprevious)*IntCoef
          y0Int <- yprevious+(y1Int-yprevious)*IntCoef
          
          arrows(x0=x0Int,y0=y0Int,x1=x1Int,y1=y1Int,col=colorCycles[j],...)
        }
        #for the end:
        if ((max(timesjinput)%in%timesj)==FALSE){
          selecInt <- (x$metadata$times==max(timesjinput)) & (x$metadata$cycles==j)
          IntCoef <- x$interpolationInfo[selecInt]
          
          timenext <- min(metadataD$times[metadataD$times>max(timesjinput) & (metadataD$sites==i)])
          EcolStateNext <- (metadataD$times==timenext) & (metadataD$sites==i)
          
          xnext <- xp[EcolStateNext]
          ynext <- yp[EcolStateNext]
          
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
#' @param x The full output of function \code{\link{extractCycles}} or \code{\link{extractFixedDateTrajectories}} as appropriate, an object of class \code{\link{trajectories}}.
#' @param fixedDates.colors The colors applied to the different fixed dates trajectories. Defaults to a simple RGB circular color palette.
#' @param sites.lty The line type for the different sites (see \code{\link{par}}, \code{"lty"}).
#' @param print.names A boolean flag to indicate whether the names of cycles or fixed-date trajectories should be printed.
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
fixedDateTrajectoryPCoA <- function (x,
                                     fixedDates.colors=NULL,
                                     sites.lty=NULL,
                                     print.names=FALSE,
                                     axes=c(1,2),...)
{
  if (!inherits(x, "fd.trajectories"))
    stop ("fixedDateTrajectoryPCoA takes as main input the whole output of function extractFixedDateTrajectories")
  
  if (is.null(sites.lty)==FALSE){
    if (length(sites.lty)!=length(unique(x$metadata$sites)))
      stop ("sites.lty must have the same number of indices as the number of sites")
    names(sites.lty) <- unique(x$metadata$sites)
  }else{
    sites.lty <- rep(1,length(unique(x$metadata$sites)))
    names(sites.lty) <- unique(x$metadata$sites)
  }
  
  D <- x$d
  PCoA <- cmdscale(D,eig=TRUE, add=TRUE, k=nrow(as.matrix(D))-1)
  metadataD <- x$metadata
  
  #Plotting
  x <- PCoA$points[,axes[1]]
  y <- PCoA$points[,axes[2]]
  
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*PCoA$eig[axes[1]]/sum(PCoA$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*PCoA$eig[axes[2]]/sum(PCoA$eig)),"%)"))
  
  #Build the color palette:
  if (is.null(fixedDates.colors)){
    fixedDates.colors <- customCircularPalette(metadataD)
  }else{
    if (length(fixedDates.colors)!=length(unique(metadataD$fdT)))
      stop ("fixedDates.colors must have the same number of indices as the number of fixed date trajectories")
    names(fixedDates.colors) <- unique(metadataD$fdT)
  }
  
  
  
  for (i in unique(metadataD$sites)){
    sitei <- metadataD$sites==i
    fdTi <- unique(metadataD$fdT[sitei])
    
    colorsfdTi <- fixedDates.colors[fdTi]
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

#' @rdname trajectoryCyclicalPlots
#' @param x the metadata \code{\link{data.frame}} of an object of class \code{\link{trajectories}}.
#' @noRd
customCircularPalette <- function(x)
{
  #this function builds a simple custom circular color palette for fixed date trajectories plotting
  #first retrieve cycleDuration:
  cycleDuration <- sort(x$times[which(x$dates==min(x$dates))])[2]-sort(x$times[which(x$dates==min(x$dates))])[1]
  
  #prepare the colors:
  coefs.col <- tapply(x$dates,x$fdT,min)/cycleDuration
  #reds:
  reds <- rep(0,length(coefs.col))
  reds[coefs.col<1/3] <- 1-coefs.col[coefs.col<1/3]*3
  reds[coefs.col>2/3] <- (coefs.col[coefs.col>2/3]-2/3)*3
  
  #greens:
  greens <- rep(0,length(coefs.col))
  greens[coefs.col<2/3] <- 1-(coefs.col[coefs.col<2/3]-1/3)*3
  greens[coefs.col<1/3] <- coefs.col[coefs.col<1/3]*3
  
  #blues:
  blues <- rep(0,length(coefs.col))
  blues[coefs.col>1/3] <- (coefs.col[coefs.col>1/3]-1/3)*3
  blues[coefs.col>2/3] <- 1-(coefs.col[coefs.col>2/3]-2/3)*3
  
  #build the palette
  palette <- rgb(red = reds, green = greens, blue = blues)
  names(palette) <- names(coefs.col)
  palette <- palette[unique(x$fdT)]#finish with a little re-ordering
  
  return(palette)
} 