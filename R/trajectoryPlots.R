#' Trajectory plots
#' 
#' Set of plotting functions for Ecological Trajectory Analysis:
#' 
#' \itemize{
#' \item{Function \code{trajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#' \item{Function \code{trajectoryPlot} draws trajectories in a scatter plot corresponding to the input coordinates.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name trajectoryPlot
#' @aliases trajectoryPCoA trajectoryPlot
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' 
#' @return 
#' Function \code{trajectoryPCoA} returns the result of calling \code{\link{cmdscale}}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @author Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la Baie de Saint-Brieuc
#' 
#' @references
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
#' @seealso \code{\link{trajectoryMetrics}}, \code{\link{transformTrajectories}}, \code{\link{cmdscale}}
#' 
#' @examples 
#' #Description of sites and surveys
#' sites <- c("1","1","1","2","2","2")
#' surveys <- c(1,2,3,1,2,3)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=6, ncol=2)
#' xy[2,2]<-1
#' xy[3,2]<-2
#' xy[4:6,1] <- 0.5
#' xy[4:6,2] <- xy[1:3,2]
#' xy[6,1]<-1
#' 
#' #Define trajectory data
#' x <- defineTrajectories(dist(xy), sites, surveys)
#'   
#' #Draw trajectories using original coordinates
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)
#' 
#' #Draw trajectories in a PCoA
#' trajectoryPCoA(x, 
#'                traj.colors = c("black","red"), lwd = 2)   
#'   
#' #Should give the same results if surveys are not in order 
#' #(here we switch surveys for site 2)
#' temp <- xy[5,]
#' xy[5,] <- xy[6,]
#' xy[6,] <- temp
#' surveys[5] <- 3
#' surveys[6] <- 2
#'   
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)   
#'  
#' x <- defineTrajectories(dist(xy), sites, surveys)
#' trajectoryPCoA(x, 
#'                traj.colors = c("black","red"), lwd = 2)   

#' @rdname trajectoryPlot
#' @param traj.colors A vector of colors (one per site). If \code{selection != NULL} the length of the color vector should be equal to the number of sites selected.
#' @param survey.labels A boolean flag to indicate whether surveys should be plotted as text next to arrow endpoints
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
#' @export
trajectoryPCoA<-function(x, traj.colors = NULL, axes=c(1,2), 
                         survey.labels = FALSE, ...) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")

  d <- x$d
  sites <- x$metadata$sites
  surveys <- x$metadata$surveys
  
  siteIDs <- unique(sites)

  cmd_D2 <- cmdscale(d,eig=TRUE, add=TRUE, k=nrow(as.matrix(d))-1)
  
  x<-cmd_D2$points[,axes[1]]
  y<-cmd_D2$points[,axes[2]]
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))
  
  #Draw arrows
  for(i in 1:length(siteIDs)) {
    ind_surv <- which(sites==siteIDs[i])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv = ind_surv[order(surveys[sites==siteIDs[i]])]
    for(t in 1:(length(ind_surv)-1)) {
      niini <-ind_surv[t]
      nifin <-ind_surv[t+1]
      if(!is.null(traj.colors)) arrows(x[niini],y[niini],x[nifin],y[nifin], col = traj.colors[i], ...)
      else arrows(x[niini],y[niini],x[nifin],y[nifin], ...)
      if(survey.labels) {
        text(x[niini],y[niini], labels = ifelse(!is.null(surveys), surveys[niini],t), pos = 3)
        if(t==(length(ind_surv)-1)) {
          text(x[nifin],y[nifin], labels = ifelse(!is.null(surveys), surveys[nifin],t+1), pos = 3)
        }
      }
    }
  }
  #Return cmdscale result
  invisible(cmd_D2)
}

#' @rdname trajectoryPlot
#' @param coords A data.frame or matrix where rows are ecosystem states and columns are coordinates in an arbitrary space
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param surveys A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @export
trajectoryPlot<-function(coords, sites, surveys = NULL, traj.colors = NULL, axes=c(1,2), 
                         survey.labels = FALSE, ...) {
  if(length(sites)!=nrow(coords)) stop("'sites' needs to be of length equal to the number of rows in 'coords'")
  if(!is.null(surveys)) {
    if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  } else {
    surveys <- rep(NA, length(sites))
    for(s in unique(sites)) {
      surveys[sites==s] <- 1:sum(sites==s)
    }
  }
  
  siteIDs <- unique(sites)

  xp <- coords[, axes[1]]
  yp <- coords[,axes[2]]
  plot(xp,yp, type="n", asp=1, xlab=paste0("Axis ",axes[1]), 
       ylab=paste0("Axis ",axes[2]))
  
  #Draw arrows
  for(i in 1:length(siteIDs)) {
    ind_surv <- which(sites==siteIDs[i])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv = ind_surv[order(surveys[sites==siteIDs[i]])]
    for(t in 1:(length(ind_surv)-1)) {
      niini <-ind_surv[t]
      nifin <-ind_surv[t+1]
      if(!is.null(traj.colors)) arrows(xp[niini],yp[niini],xp[nifin],yp[nifin], col = traj.colors[i], ...)
      else arrows(xp[niini],yp[niini],xp[nifin],yp[nifin], ...)
      if(survey.labels) {
        text(xp[niini],yp[niini], labels = ifelse(!is.null(surveys), surveys[niini],t), pos = 3)
        if(t==(length(ind_surv)-1)) {
          text(xp[nifin],yp[nifin], labels = ifelse(!is.null(surveys), surveys[nifin],t+1), pos = 3)
        }
      }
    }
  }
}
