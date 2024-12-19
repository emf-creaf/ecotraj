#' Cyclical trajectory plots
#' 
#' Set of plotting functions for Cyclical Ecological Trajectory Analysis:
#' 
#' \itemize{
#' \item{Function \code{cyclePCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#' \item{Function \code{fixedDateTrajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#' }
#' 
#' @encoding UTF-8
#' @name trajectoryplots
#' @aliases cyclePCoA fixedDateTrajectoryPCoA
#' 
#' CETA plotting functions:
#' \itemize{
#'  \item{Function \code{cyclePCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws cycles in the ordination scatterplot. Sister function of \code{trajectoryPCoA} adapted to cycles.}
#' }
#' 
#' @details
#' Additional details...
#' 
#' 
#' @return 
#' Function \code{cyclePCoA} returns a list similar to the calling \code{\link{cmdscale}}.
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' 
#' @seealso \code{\link{trajectoryCyclical}}}, \code{\link{cmdscale}}
#' 
#' 
#' @rdname trajectoryCyclicalPlots
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states.
#' @param cycles A vector indicating the cycles corresponding to each ecosystem state.
#' @param times A vector indicating the times corresponding to each ecosystem state (equivalent to "surveys" in other ETA function but more time-explicit).
#' @param selection A character vector of cycles, a numeric vector of cycle indices or logical vector of the same length as \code{cycles}, indicating a subset of cycle trajectories to be selected.
#' @param IntExt A vector indicating whether each ecosystem state is "internal" or "external".
#' @param cycle.colors A vector of colors (one per cycle). If \code{selection != NULL} the length of the color vector should be equal to the number of cycles selected.
#' @param time.labels A boolean flag to indicate whether times should be plotted as text next to arrow endpoints
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
cyclePCoA <- function (d,cycles,times,IntExt,selection = NULL,cycle.colors = NULL, axes=c(1,2),time.labels = FALSE, ...)
{
  cycleIDs <- unique(cycles)
  ncycle <- length(cycleIDs)
  
  #Apply cycle selection
  
  if(is.null(selection)){
    selection <- 1:ncycle 
  } else if(is.character(selection)){
    selection <- (cycleIDs %in% selection)
  }
  selIDs <- cycleIDs[selection]
  
  D2 <- as.dist(as.matrix(d)[cycles %in% selIDs, cycles %in% selIDs])
  #cmd_D2 <- cmdscale(D2,eig=TRUE, add=TRUE, k=nrow(as.matrix(D2))-1)
  
  #START OF THE BIG WEIRD PCoA MODIFS... THIS IS SUPPOSED TO RETURN SOME SORT OF "WEIGHTED PCoA". BUT I FEAR IT'S ALL CRAPPY!!!
  cmd_D2 <-list()
  
  Corr <- cmdscale(D2,eig=TRUE, add=TRUE, k=nrow(as.matrix(D2))-1)$ac
  
  dCorr <- as.matrix(D2)+Corr
  diag(dCorr) <- 0
  
  n <- ncol(dCorr)
  A <- -0.5*(dCorr)^2
  
  Int <- as.numeric(IntExt[which(cycles %in% selIDs)]=="internal")
  D <- diag(Int)/sum(Int)
  ones <- rep(1,n)
  I <- diag(1,n)
  
  Q <- I - (ones%*%t(ones))%*%D
  
  Delta1 <- Q%*%A%*%t(Q)
  
  U <- eigen(Delta1)$vectors
  eig <- eigen(Delta1)$values
  points <- U*sqrt(matrix(eig,n,n,byrow=T))
  
  cmd_D2$points <- points
  cmd_D2$eig <- eig
  #END OF THE BIG WEIRD PCoA MODIFS...
  
  
  x <- cmd_D2$points[,axes[1]]
  y <- cmd_D2$points[,axes[2]]
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))
  
  cyclesred <- cycles[cycles %in% selIDs]
  if(!is.null(times)){
    timesred <- times[cycles %in% selIDs]
  } else {
    timesred <- NULL
  }
  #Draw arrows
  for(i in 1:length(selIDs)) {
    ind_time <- which(cyclesred==selIDs[i])
    #Times may not be in order
    if(!is.null(timesred)){
      ind_time <- ind_time[order(timesred[cyclesred==selIDs[i]])]
    } 
    for(t in 1:(length(ind_time)-1)) {
      niini <- ind_time[t]
      nifin <- ind_time[t+1]
      if(!is.null(cycle.colors)){
        arrows(x[niini],y[niini],x[nifin],y[nifin], col = cycle.colors[i])#, ...)
      } else {
        arrows(x[niini],y[niini],x[nifin],y[nifin])#, ...)
      }
      if(time.labels) {
        text(x[niini],y[niini], labels = ifelse(!is.null(timesred), timesred[niini],t), pos = 3)
        if(t==(length(ind_time)-1)) {
          text(x[nifin],y[nifin], labels = ifelse(!is.null(timesred), timesred[nifin],t+1), pos = 3)
        }
      }
    }
  }
  #Return cmdscale result
  invisible(cmd_D2)
}


#' @rdname trajectoryCyclicalPlots
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
fixedDateTrajectoryPCoA <- function (){
  
}

