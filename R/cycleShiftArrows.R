
#' Displaying cycle shifts
#' 
#' Adds arrows representing cyclical shifts (advances/delays) into convergence/divergence plots created using function \code{\link{trajectoryConvergencePlot}}.
#' 
#' @encoding UTF-8
#' @name cycleShiftArrows
#' 
#' @param cycle.shifts Cyclical shifts computed for each fixed date trajectory plotted by \code{trajectoryConvergencePlot}.
#' @param radius The radius of the circles representing trajectories. Defaults to 1.
#' @param top A string indicating if the top of the plotting area should contain a circle representing a trajectory ("circle"), or should be in between two circles ("between"). Defaults to "between". 
#' @param cycle.shifts.inf.conf Lower confidence intervals for cyclical shifts.
#' @param cycle.shifts.sup.conf Upper confidence intervals for cyclical shifts.
#' @param arrows.length.mult A multiplication coefficient for the arrows representing cyclical shifts. Attempts an automatic adjustment by default (dividing by max(cycle.shifts)).
#' @param arrows.lwd Line width of the arrows. Defaults to 2.
#' 
#' @details
#' This function is meant to be used in conjunction with \code{\link{trajectoryConvergencePlot}}, to study fixed-date trajectories convergence/divergence patterns:
#' \itemize{
#' \item{First, setting \code{pointy = TRUE}, in the call to \code{\link{trajectoryConvergencePlot}}, when the studied trajectories are the outputs of \code{\link{extractFixedDateTrajectories}} allows to suggest cyclicity.}
#' \item{Second, function \code{cycleShiftArrows} allows to represent cyclical shifts by adding arrows to the circles representing trajectories. Clockwise arrows will represent advances, anticlockwise arrows will represent delays. Arrows length is proportional to the cyclical shift provided.
#' Relevant measures of cyclical shifts have to be computed by the user. A variety of methods may be employed but the outputs of \code{\link{cycleShifts}} cannot be used immediately as they do not directly correspond to the fixed-date trajectories. An example of how to compute relevant 
#' measures of cyclical shifts is provided in \code{\link{cycleShifts}} and in the CETA vignette.
#' The arguments \code{radius} and \code{top} in \code{cycleShiftArrows} must match those in the corresponding call to \code{\link{trajectoryConvergencePlot}} for proper display.}
#' }
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @seealso \code{\link{trajectoryConvergencePlot}},\code{\link{cycleShifts}}
#' @examples
#' data("avoca")
#' avoca_D_man <- vegclust::vegdiststruct(avoca_strat, 
#'                                        method ="manhattan", 
#'                                        transform = function(x){log(x+1)})
#' years <- c(1971, 1974, 1978, 1983, 1987, 1993, 1999, 2004, 2009)
#' avoca_times <- years[avoca_surveys]
#' avoca_x <- defineTrajectories(d = avoca_D_man,  
#'                               sites = avoca_sites, 
#'                               times = avoca_times)
#' 
#'
#' #CETA-specific example of convergence/divergence plot with made-up cyclical shifts:
#' trajectoryConvergencePlot(avoca_x,
#'                           type = "both",
#'                           alpha.filter = 0.05,
#'                           half.arrows.size = 2,
#'                           pointy = TRUE)
#' 
#' #Made-up cyclical shifts (matching the number of trajectories in avoca_x)
#' CS <- c(2.5, 1.5, 0.5, 0.5, 3, 1, 1, 2)
#' cycleShiftArrows(CS,
#'                   cycle.shifts.inf.conf = CS - 0.2,
#'                   cycle.shifts.sup.conf = CS + 0.2)
#' @export
cycleShiftArrows <- function (cycle.shifts,
                               radius = 1,
                               top = "between",
                               cycle.shifts.inf.conf = NULL,
                               cycle.shifts.sup.conf = NULL,
                               arrows.length.mult = "auto",
                               arrows.lwd = 2){
  if ((!is.null(cycle.shifts.inf.conf))&&(is.null(cycle.shifts.sup.conf))) stop("If confidence intervals are provided, both superiors and inferiors are needed!")
  if ((is.null(cycle.shifts.inf.conf))&&(!is.null(cycle.shifts.sup.conf))) stop("If confidence intervals are provided, both superiors and inferiors are needed!")
  
  #scaling the arrows
  if(arrows.length.mult == "auto"){
    if ((!is.null(cycle.shifts.inf.conf))&&(!is.null(cycle.shifts.sup.conf))){
      cycle.shifts.inf.conf <- cycle.shifts.inf.conf/(max(cycle.shifts))
      cycle.shifts.sup.conf <- cycle.shifts.sup.conf/(max(cycle.shifts))
    }
    cycle.shifts <- cycle.shifts/(max(cycle.shifts))
    arrows.length.mult <- 1
  }
  
  radius <- radius*0.1 #reduce radius (just to have convenient numbers in function calling)
  
  nTraj <- length(cycle.shifts)
  #Find the angles and centers of the shapes representing the trajectories
  quadrant <- (2*pi)/nTraj
  if (top=="circle"){
    angles <- seq(pi/2,length.out=nTraj,by=-quadrant)
    centersX <- cos(angles)
    centersY <- sin(angles)
  }else if (top=="between"){
    angles <- seq(pi/2-quadrant/2,length.out=nTraj,by=-quadrant)
    centersX <- cos(angles)
    centersY <- sin(angles)
  }else{
    stop("top must be either circle or between and should match the same argument in the parent convergence plot")
  }
  
  baseArrowsX <- centersX*(1+radius)
  baseArrowsY <- centersY*(1+radius)
  if ((!is.null(cycle.shifts.inf.conf))&&(!is.null(cycle.shifts.sup.conf))){
    arrows(x0=baseArrowsY*cycle.shifts.inf.conf*arrows.length.mult+baseArrowsX,
           y0=-baseArrowsX*cycle.shifts.inf.conf*arrows.length.mult+baseArrowsY,
           x1=baseArrowsY*cycle.shifts.sup.conf*arrows.length.mult+baseArrowsX,
           y1=-baseArrowsX*cycle.shifts.sup.conf*arrows.length.mult+baseArrowsY,
           angle=90,length=0.05,col="grey30",code=3,lwd=arrows.lwd-1,xpd=NA)
  }
  arrows(x0=baseArrowsX,y0=baseArrowsY,
         x1=baseArrowsY*cycle.shifts*arrows.length.mult+baseArrowsX,
         y1=-baseArrowsX*cycle.shifts*arrows.length.mult+baseArrowsY,
         lwd=arrows.lwd,length=0.1,xpd=NA)
  points(baseArrowsX,baseArrowsY,pch=16)
}



