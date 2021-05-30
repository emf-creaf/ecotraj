#' Trajectory plots
#' 
#' Set of plotting functions for Ecosystem Trajectory Analysis:
#' 
#' \itemize{
#' \item{Function \code{trajectoryPCoA} performs principal coordinates analysis (\code{\link{cmdscale}}) and draws trajectories in the ordination scatterplot.}
#' \item{Function \code{trajectoryPlot} Draws trajectories in a scatterplot corresponding to the input coordinates.}
#' }
#'  
#' 
#' @encoding UTF-8
#' @name trajectoryplots
#' @aliases trajectoryPCoA trajectoryPlot
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of community states (see details).
#' @param sites A vector indicating the site corresponding to each community state.
#' @param surveys A vector indicating the survey corresponding to each community state (only necessary when surveys are not in order).
#' 
#' @details 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (2019). 
#' The input distance matrix \code{d} should ideally be metric. That is, all subsets of distance triplets should fulfill the triangle inequality (see function \code{is.metric}). 
#' All CTA functions that require metricity include a parameter '\code{add}', which by default is TRUE, meaning that whenever the triangle inequality is broken the minimum constant required to fulfill it is added to the three distances.
#' If such local (an hence, inconsistent across triplets) corrections are not desired, users should find another way modify \code{d} to achieve metricity, such as PCoA, metric MDS or non-metric MDS (see CTA vignette). 
#' If parameter '\code{add}' is set to FALSE and problems of triangle inequality exist, CTA functions may provide missing values in some cases where they should not.
#' 
#' The resemblance between trajectories is done by adapting concepts and procedures used for the analysis of trajectories in space (i.e. movement data) (Besse et al. 2016).   
#' 
#' @return 
#' Function \code{trajectoryPCoA} returns the result of calling \code{\link{cmdscale}}.
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @author Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la Baie de Saint-Brieuc
#' 
#' @references
#' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
#' 
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). Trajectory analysis in community ecology. Ecological Monographs.
#' 
#' Anderson (2017). Permutational Multivariate Analysis of Variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article ID: stat07841.
#' 
#' @seealso \code{\link{trajectoryDistances}}, \code{\link{cmdscale}}
#' 
#' @examples 
#' #Description of sites and surveys
#' sites = c(1,1,1,2,2,2)
#' surveys=c(1,2,3,1,2,3)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=6, ncol=2)
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
#'   
#' #Should give the same results if surveys are not in order 
#' #(here we switch surveys for site 2)
#' temp = xy[5,]
#' xy[5,] = xy[6,]
#' xy[6,] = temp
#' surveys[5] = 3
#' surveys[6] = 2
#'   
#' trajectoryPlot(xy, sites, surveys, 
#'                traj.colors = c("black","red"), lwd = 2)   
#'  

#' @rdname trajectoryplots
#' @param traj.colors A vector of colors (one per site). If \code{selection != NULL} the length of the color vector should be equal to the number of sites selected.
#' @param survey.labels A boolean flag to indicate whether surveys should be plotted as text next to arrow endpoints
#' @param axes The pair of principal coordinates to be plotted.
#' @param ... Additional parameters for function \code{\link{arrows}}.
trajectoryPCoA<-function(d, sites, surveys = NULL, selection = NULL, traj.colors = NULL, axes=c(1,2), 
                         survey.labels = FALSE, ...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  
  #Apply site selection
  
  if(is.null(selection)) selection = 1:nsite 
  else {
    if(is.character(selection)) selection = (siteIDs %in% selection)
  }
  selIDs = siteIDs[selection]
  
  D2 =as.dist(as.matrix(d)[sites %in% selIDs, sites %in% selIDs])
  cmd_D2<-cmdscale(D2,eig=TRUE, add=TRUE, k=nrow(as.matrix(D2))-1)
  
  x<-cmd_D2$points[,axes[1]]
  y<-cmd_D2$points[,axes[2]]
  plot(x,y, type="n", asp=1, xlab=paste0("PCoA ",axes[1]," (", round(100*cmd_D2$eig[axes[1]]/sum(cmd_D2$eig)),"%)"), 
       ylab=paste0("PCoA ",axes[2]," (", round(100*cmd_D2$eig[axes[2]]/sum(cmd_D2$eig)),"%)"))
  
  sitesred = sites[sites %in% selIDs]
  if(!is.null(surveys)) surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  #Draw arrows
  for(i in 1:length(selIDs)) {
    ind_surv = which(sitesred==selIDs[i])
    #Surveys may not be in order
    if(!is.null(surveysred)) ind_surv = ind_surv[order(surveysred[sitesred==selIDs[i]])]
    for(t in 1:(length(ind_surv)-1)) {
      niini =ind_surv[t]
      nifin =ind_surv[t+1]
      if(!is.null(traj.colors)) arrows(x[niini],y[niini],x[nifin],y[nifin], col = traj.colors[i], ...)
      else arrows(x[niini],y[niini],x[nifin],y[nifin], ...)
      if(survey.labels) {
        text(x[niini],y[niini], labels = ifelse(!is.null(surveysred), surveysred[niini],t), pos = 3)
        if(t==(length(ind_surv)-1)) {
          text(x[nifin],y[nifin], labels = ifelse(!is.null(surveysred), surveysred[nifin],t+1), pos = 3)
        }
      }
    }
  }
  #Return cmdscale result
  invisible(cmd_D2)
}

#' @rdname trajectoryplots
#' @param x A data.frame or matrix where rows are community states and columns are coordinates in an arbitrary space
trajectoryPlot<-function(x, sites, surveys = NULL, selection = NULL, traj.colors = NULL, axes=c(1,2), 
                         survey.labels = FALSE, ...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  
  #Apply site selection
  
  if(is.null(selection)) selection = 1:nsite 
  else {
    if(is.character(selection)) selection = (siteIDs %in% selection)
  }
  selIDs = siteIDs[selection]
  
  xp =x[sites %in% selIDs, axes[1]]
  yp<-x[sites %in% selIDs,axes[2]]
  plot(xp,yp, type="n", asp=1, xlab=paste0("Axis ",axes[1]), 
       ylab=paste0("Axis ",axes[2]))
  
  sitesred = sites[sites %in% selIDs]
  if(!is.null(surveys)) surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  #Draw arrows
  for(i in 1:length(selIDs)) {
    ind_surv = which(sitesred==selIDs[i])
    #Surveys may not be in order
    if(!is.null(surveysred)) ind_surv = ind_surv[order(surveysred[sitesred==selIDs[i]])]
    for(t in 1:(length(ind_surv)-1)) {
      niini =ind_surv[t]
      nifin =ind_surv[t+1]
      if(!is.null(traj.colors)) arrows(xp[niini],yp[niini],xp[nifin],yp[nifin], col = traj.colors[i], ...)
      else arrows(xp[niini],yp[niini],xp[nifin],yp[nifin], ...)
      if(survey.labels) {
        text(xp[niini],yp[niini], labels = ifelse(!is.null(surveysred), surveysred[niini],t), pos = 3)
        if(t==(length(ind_surv)-1)) {
          text(xp[nifin],yp[nifin], labels = ifelse(!is.null(surveysred), surveysred[nifin],t+1), pos = 3)
        }
      }
    }
  }
}
