#' Summary plot for trajectory convergence and divergence
#' 
#' Provides plots to represent trajectory convergence and divergence tests performed by the function \code{\link{trajectoryConvergence}} or to present the results of Relative Trajectory Movement Assessment (\code{\link{trajectoryRMA}}).
#' 
#' @encoding UTF-8
#' @name trajectoryConvergencePlot
#' 
#' @details
#' Function \code{trajectoryConvergencePlot} provides ways to visualize pairwise convergence and divergence between trajectories.
#' It has two modes of functioning:
#' \itemize{
#'    \item{If \code{x} is of class \code{\link{trajectories}}, the function will display the results of convergence/divergence tests by calls to function \code{\link{trajectoryConvergence}}.}
#'    \item{If \code{x} is of class \code{RTMA},the function will display the results of convergence/divergence tests and dynamic correspondence tests stored in the \code{RTMA} object supplied.}
#' }
#' In the plots, trajectories are represented by circles. The convergence or divergence between pairs of trajectories are represented by links. If convergence tests are symmetric, the links are simple. If the convergence tests are asymmetric, the links are displayed as half arrows pointing from the trajectory converging or diverging towards the trajectory being approached or diverged from.
#' The width and color hue of the links are proportional to the tau statistic of the Mann-Kendall test performed by the \code{\link{trajectoryConvergence}} function. 
#' Function \code{trajectoryConvergencePlot} also offers the possibility to plot both tests at the same time.
#' 
#' If \code{x} is of class \code{RTMA}, \code{trajectoryConvergencePlot} will display both convergence tests as explained above, as well as cases of parallelism recognized in \code{\link{trajectoryRMA}}.
#' \code{Parallel} scenarios are indicated by two full parallel black lines linking two trajectories, while in case of \code{Antiparallel} scenarios one of the lines is dotted.
#' 
#' In addition, see function \code{\link{cycleShiftArrows}} for additional graphical elements to be displayed when conducting CETA.
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @references
#' Djeghri et al. (in preparation) Uncovering the relative movements of ecological trajectories.
#' 
#' @seealso \code{\link{trajectoryConvergence}},\code{\link{trajectoryRMA}}, \code{\link{cycleShiftArrows}}
#' 
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
#' #Raw output with asymmetric convergence test (default)
#' trajectoryConvergencePlot(avoca_x)
#' 
#' #More refined output with both type of tests and only plotting significant 
#' #test results (p-value < 0.05)
#' trajectoryConvergencePlot(avoca_x,
#'                           type = "both",
#'                           alpha.filter = 0.05)
#' 
#' #Much more refined output with nicer colors, bigger half arrows, 
#' #personalized trajectory names, controlling the size of circles representing 
#' #trajectories and border customization.
#' trajectoryConvergencePlot(avoca_x,type = "both",alpha.filter = 0.05,
#'                           half.arrows.size = 1.5, 
#'                           conv.color = "orangered",
#'                           div.color = "dodgerblue",
#'                           radius = 1.2, 
#'                           traj.colors = "black",border = "white",lwd = 2,
#'                           traj.names = LETTERS[1:8],traj.names.colors = "white")
#' #RTMA version.
#' avoca_RTMA <- trajectoryRMA(avoca_x)
#' trajectoryConvergencePlot(avoca_RTMA,
#'                           half.arrows.size = 1.5, 
#'                           conv.color = "orangered",
#'                           div.color = "dodgerblue",
#'                           radius = 1.2, 
#'                           traj.colors = "black",border = "white",lwd = 2,
#'                           traj.names = LETTERS[1:8],traj.names.colors = "white")
#'                           
#' @rdname trajectoryConvergencePlot
#' @param x An object of class \code{\link{trajectories}}. Alternatively, an object of class \code{RTMA} (returned by \code{\link{trajectoryRMA}}).
#' @param type A string indicating the convergence test to be displayed, either \code{"pairwise.asymmetric"}, \code{"pairwise.symmetric"} or \code{"both"} (see \code{\link{trajectoryConvergence}}). Disregarded if \code{inherits(x,"RTMA")}.
#' @param alpha.filter The minimum p-value for a link to be drawn (see \code{\link{trajectoryConvergence}}). Defaults to \code{NULL} (all links drawn). Disregarded if \code{inherits(x,"RTMA")} (the RTMA corrected alpha level is used instead).
#' @param traj.colors The colors for the trajectories (circles). Defaults to \code{"grey"}.
#' @param traj.names The names of trajectories. Defaults to the names provided in \code{x}.
#' @param traj.names.colors The color of the names of trajectories on the circles. Defaults to \code{"black"}.
#' @param ... Additional parameters passed to \code{\link{polygon}} to personalize the circles representing trajectories.
#' @param radius The radius of the circles representing trajectories. Defaults to \code{1}.
#' @param conv.color The color used to mark convergent trajectories. Defaults to \code{"red"}.
#' @param div.color The color used to mark divergent trajectories. Defaults to \code{"blue"}.
#' @param half.arrows.size A multiplication coefficient for the size of the arrow heads when representing asymmetric tests results. Defaults to \code{1}.
#' @param tau.links.transp The transparency of the links representing the tau statistic of the Mann-Kendall test (see \code{\link{trajectoryConvergence}}).
#' @param top A string indicating if the top of the plotting area should contain a circle representing a trajectory (\code{"circle"}), or should be in between two circles (\code{"between"}). Defaults to \code{"between"}.
#' @param pointy Boolean. Should the circles representing trajectories be made pointy (i.e. pointing to the next trajectory)? Useful when trajectories have some order, as in the context of CETA to represent fixed date trajectories (see \code{\link{trajectoryCyclical}}).
#' @param add Passed to function \code{\link{trajectoryConvergence}}. Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality. Disregarded if \code{inherits(x,"RTMA")}.
#' @export
trajectoryConvergencePlot <- function (x,
                                       type = "pairwise.asymmetric",
                                       alpha.filter = NULL,
                                       traj.colors = "grey",
                                       traj.names = NULL,
                                       traj.names.colors = "black",
                                       ...,
                                       radius = 1,
                                       conv.color = "red",
                                       div.color = "blue",
                                       half.arrows.size = 1,
                                       tau.links.transp = 0.3,
                                       top = "between",
                                       pointy = FALSE,
                                       add = TRUE){
  
  widthMult <- 1#this is a multiplication coefficient for the width of the "tau links" it will change to make room to display more links in case "both" is selected
  
  if(inherits(x,"trajectories")){
    if (type=="pairwise.asymmetric"){
      ConvTest <- trajectoryConvergence(x,type="pairwise.asymmetric",add=add)
    }else if (type=="pairwise.symmetric"){
      ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric",add=add)
    }else if (type=="both"){
      ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric",add=add)
      ConvTestAsym <- trajectoryConvergence(x,type="pairwise.asymmetric",add=add)
      widthMult <- 1/3
    }else{
      stop("invalid value for 'type'")
    }
  }else if (inherits(x,"RTMA")){
    ConvTest <- x$symmetric_convergence
    ConvTestAsym <- x$asymmetric_convergence
    widthMult <- 1/3
    type <- "both"
    alpha.filter <- x$parameters["corrected_alpha"]
  }else{
    stop ("'x' should be of class 'trajectory' or 'RTMA'")
  }
  
  
  #Check for appropriate length of the parameter traj.colors
  if (length(traj.colors)==1){
    traj.colors <- rep(traj.colors,ncol(ConvTest$tau))
  }else if (length(traj.colors)!=ncol(ConvTest$tau)){
    stop("traj.colors must have the same length as the number of trajectories or length 1")
  }
  
  #Give trajectory names if not provided
  if (is.null(traj.names)){
    traj.names <- colnames(ConvTest$tau)
  }
  
  nTraj <- ncol(ConvTest$tau)
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
    stop("top must be either circle or between")
  }
  names(angles) <- rownames(ConvTest$tau)
  names(centersX) <- rownames(ConvTest$tau)
  names(centersY) <- rownames(ConvTest$tau)
  
  #Call the plot background
  plot(NA,NA,xlim=c(-1,1),ylim=c(-1,1),asp=1,axes=F,xlab="",ylab="")
  par(xpd=NA)
  radius <- radius*0.1 #reduce radius (just to have convenient numbers in function calling)
  
  #Draw the convergence/divergence rectangles/half arrows
  if (type=="both"){
    for (i in rownames(ConvTest$tau)){
      for (j in colnames(ConvTest$tau)){
        if (i!=j){
          if (is.null(alpha.filter)||(ConvTestAsym$p.value[i,j]<alpha.filter)){
            tau <- ConvTestAsym$tau[i,j]
            tauabs <- abs(tau)
            taucol <- ConvDivColor(tau=tau,conv.color=conv.color,div.color=div.color)
            
            #find the appropriate vector (the length of our future half-arrow) and scale it to unit length
            vecX <- centersX[j]-centersX[i]
            vecY <- centersY[j]-centersY[i]
            vecLength <- sqrt(vecX^2+vecY^2)
            vecX <- vecX/vecLength
            vecY <- vecY/vecLength
            
            #use this to find the coordinate of the corners of the half-arrow
            #the width of the half-arrow will be proportional to tau
            
            #prepare things to move out the asymmetric arrows so that they sit on the symmetric rectangle
            if (is.null(alpha.filter)||(ConvTest$p.value[i,j]<alpha.filter)){
              centersXi <- vecY*radius*abs(ConvTest$tau[i,j])*widthMult+centersX[i]
              centersYi <- -vecX*radius*abs(ConvTest$tau[i,j])*widthMult+centersY[i]
              
              centersXj <- vecY*radius*abs(ConvTest$tau[i,j])*widthMult+centersX[j]
              centersYj <- -vecX*radius*abs(ConvTest$tau[i,j])*widthMult+centersY[j]
            }else{
              centersXi <- centersX[i]
              centersYi <- centersY[i]
              
              centersXj <- centersX[j]
              centersYj <- centersY[j]
            }
            
            
            Xcorner1 <- centersXi
            Ycorner1 <- centersYi
            
            Xcorner2 <- vecY*radius*tauabs*widthMult+centersXi
            Ycorner2 <- -vecX*radius*tauabs*widthMult+centersYi
            
            Xcorner3 <- vecY*radius*tauabs*widthMult+centersXj-(vecX*radius+2*vecX*radius*tauabs*widthMult*half.arrows.size)
            Ycorner3 <- -vecX*radius*tauabs*widthMult+centersYj-(vecY*radius+2*vecY*radius*tauabs*widthMult*half.arrows.size)
            
            Xcorner4 <- 1.5*vecY*radius*tauabs*widthMult*half.arrows.size+centersXj-(vecX*radius+2*vecX*radius*tauabs*widthMult*half.arrows.size)
            Ycorner4 <- -1.5*vecX*radius*tauabs*widthMult*half.arrows.size+centersYj-(vecY*radius+2*vecY*radius*tauabs*widthMult*half.arrows.size)
            
            Xcorner5 <- centersXj-vecX*radius
            Ycorner5 <- centersYj-vecY*radius
            
            #finally draw the half-arrow!
            polygon(x=c(Xcorner1,Xcorner2,Xcorner3,Xcorner4,Xcorner5),
                    y=c(Ycorner1,Ycorner2,Ycorner3,Ycorner4,Ycorner5),
                    col=rgb(t(col2rgb(taucol)/255),alpha=tau.links.transp),border=taucol,lty=2)
          }
        }
      }
    }
  }
  if (type!="pairwise.asymmetric"){
    count <- 1
    for (i in rownames(ConvTest$tau)[1:(nrow(ConvTest$tau)-1)]){
      count <- count+1
      for (j in colnames(ConvTest$tau)[count:ncol(ConvTest$tau)]){
        if (is.null(alpha.filter)||(ConvTest$p.value[i,j]<alpha.filter)){
          tau <- ConvTest$tau[i,j]
          tauabs <- abs(tau)
          taucol <- ConvDivColor(tau=tau,conv.color=conv.color,div.color=div.color)
          
          #find the appropriate vector (the length of our future rectangle) and scale it to unit length
          vecX <- centersX[j]-centersX[i]
          vecY <- centersY[j]-centersY[i]
          vecLength <- sqrt(vecX^2+vecY^2)
          vecX <- vecX/vecLength
          vecY <- vecY/vecLength
          
          #use this to find the coordinate of the corners of the rectangle
          #the width of the rectangle will be proportional to tau
          Xcorner1 <- -vecY*radius*tauabs*widthMult+centersX[i]
          Ycorner1 <- vecX*radius*tauabs*widthMult+centersY[i]
          
          Xcorner2 <- vecY*radius*tauabs*widthMult+centersX[i]
          Ycorner2 <- -vecX*radius*tauabs*widthMult+centersY[i]
          
          Xcorner3 <- vecY*radius*tauabs*widthMult+centersX[j]
          Ycorner3 <- -vecX*radius*tauabs*widthMult+centersY[j]
          
          Xcorner4 <- -vecY*radius*tauabs*widthMult+centersX[j]
          Ycorner4 <- vecX*radius*tauabs*widthMult+centersY[j]
          
          #finally draw the rectangle!
          polygon(x=c(Xcorner1,Xcorner2,Xcorner3,Xcorner4),
                  y=c(Ycorner1,Ycorner2,Ycorner3,Ycorner4),
                  col=rgb(t(col2rgb(taucol)/255),alpha=tau.links.transp),border=taucol)
        }
      }
    }
  }
  if (type=="pairwise.asymmetric"){
    for (i in rownames(ConvTest$tau)){
      for (j in colnames(ConvTest$tau)){
        if (i!=j){
          if (is.null(alpha.filter)||(ConvTest$p.value[i,j]<alpha.filter)){
            tau <- ConvTest$tau[i,j]
            tauabs <- abs(tau)
            taucol <- ConvDivColor(tau=tau,conv.color=conv.color,div.color=div.color)
            
            #find the appropriate vector (the length of our future half-arrow) and scale it to unit length
            vecX <- centersX[j]-centersX[i]
            vecY <- centersY[j]-centersY[i]
            vecLength <- sqrt(vecX^2+vecY^2)
            vecX <- vecX/vecLength
            vecY <- vecY/vecLength
            
            #use this to find the coordinate of the corners of the half-arrow
            #the width of the half-arrow will be proportional to tau
            Xcorner1 <- centersX[i]
            Ycorner1 <- centersY[i]
            
            Xcorner2 <- vecY*radius*tauabs*widthMult+centersX[i]
            Ycorner2 <- -vecX*radius*tauabs*widthMult+centersY[i]
            
            Xcorner3 <- vecY*radius*tauabs*widthMult+centersX[j]-(vecX*radius+2*vecX*radius*tauabs*widthMult*half.arrows.size)
            Ycorner3 <- -vecX*radius*tauabs*widthMult+centersY[j]-(vecY*radius+2*vecY*radius*tauabs*widthMult*half.arrows.size)
            
            Xcorner4 <- 1.5*vecY*radius*tauabs*widthMult*half.arrows.size+centersX[j]-(vecX*radius+2*vecX*radius*tauabs*widthMult*half.arrows.size)
            Ycorner4 <- -1.5*vecX*radius*tauabs*widthMult*half.arrows.size+centersY[j]-(vecY*radius+2*vecY*radius*tauabs*widthMult*half.arrows.size)
            
            Xcorner5 <- centersX[j]-vecX*radius
            Ycorner5 <- centersY[j]-vecY*radius
            
            #finally draw the half-arrow!
            polygon(x=c(Xcorner1,Xcorner2,Xcorner3,Xcorner4,Xcorner5),
                    y=c(Ycorner1,Ycorner2,Ycorner3,Ycorner4,Ycorner5),
                    col=rgb(t(col2rgb(taucol)/255),alpha=tau.links.transp),border=taucol)
          }
        }
      }
    }
  }
  #add the parallelism results for RTMA
  if (inherits(x,"RTMA")){
    count <- 1
    for (i in rownames(ConvTest$tau)[1:(nrow(ConvTest$tau)-1)]){
      count <- count+1
      for (j in colnames(ConvTest$tau)[count:ncol(ConvTest$tau)]){
        if ((ConvTest$p.value[i,j]>alpha.filter)&
            (ConvTestAsym$p.value[i,j]>alpha.filter)&
            (ConvTestAsym$p.value[j,i]>alpha.filter)&
            (x$correspondence[j,i]<=alpha.filter)){
          
          vecX <- centersX[j]-centersX[i]
          vecY <- centersY[j]-centersY[i]
          vecLength <- sqrt(vecX^2+vecY^2)
          vecX <- vecX/vecLength
          vecY <- vecY/vecLength
          
          Xcorner1 <- -vecY*radius*0.25+centersX[i]
          Ycorner1 <- vecX*radius*0.25+centersY[i]
          
          Xcorner2 <- vecY*radius*0.25+centersX[i]
          Ycorner2 <- -vecX*radius*0.25+centersY[i]
          
          Xcorner3 <- vecY*radius*0.25+centersX[j]
          Ycorner3 <- -vecX*radius*0.25+centersY[j]
          
          Xcorner4 <- -vecY*radius*0.25+centersX[j]
          Ycorner4 <- vecX*radius*0.25+centersY[j]
          
          if(x$correspondence[i,j]>0){
            segments(x0=c(Xcorner1,Xcorner2),x1=c(Xcorner4,Xcorner3),
                     y0=c(Ycorner1,Ycorner2),y1=c(Ycorner4,Ycorner3))
          }else{
            segments(x0=c(Xcorner1,Xcorner2),x1=c(Xcorner4,Xcorner3),
                     y0=c(Ycorner1,Ycorner2),y1=c(Ycorner4,Ycorner3),lty=c(1,2))
          }
        }
      }
    }
  }
  
  #Draw the cycle representing trajectories
  PointyCircle(centersX,centersY,pointyCircles.colors=traj.colors,
               rad=radius,pointy=pointy,pointy.at=angles-pi/2,...)
  text(centersX,centersY,traj.names,col=traj.names.colors,font=2)
  
}

#' @rdname trajectoryConvergencePlot
#' @param x The positions of the pointy circles to be drawn on the x axis.
#' @param y The positions of the pointy circles to be drawn on the y axis.
#' @param rad The radius of the pointy circles to be drawn. Defaults to 1.
#' @param pointy Boolean. Should the circles be pointy. Defaults to FALSE.
#' @param pointyat The angle in radiant at which the point should be drawn if \code{pointy} is FALSE.
#' @param pointyCircles.colors The color filling the pointy circles. Defaults to "grey".
#' @param ... Other parameters for function \code{\link{polygon}}. 
#' @noRd
#' @keywords internal
PointyCircle <- function (x,
                          y,
                          rad = 1,
                          pointy = FALSE,
                          pointy.at = NULL,
                          pointyCircles.colors = "grey",
                          ...){
  if (length(y)!=length(x)){
    stop("x and y must have the same length")
  }
  if (length(x)>1){
    if (length(rad)==1){
      rad <- rep(rad,length(x))
    }else if (length(rad)!=length(x)){
      stop("rad must have the same length as x or length 1")
    }
    if (length(pointyCircles.colors)==1){
      pointyCircles.colors <- rep(pointyCircles.colors,length(x))
    }else if (length(pointyCircles.colors)!=length(x)){
      stop("pointyCircles.colors must have the same length as x or length 1")
    }
    if (pointy==T){
      if (length(pointy.at)==1){
        pointy.at <- rep(pointy.at,length(x))
      }else if (length(pointy.at)!=length(x)){
        stop("if pointy = T, pointy.at must have the same length as x or length 1")
      }
    }
  }
  
  for (i in 1:length(x)){
    if (pointy==T){
      pointy.at[i] <- pointy.at[i]%%(2*pi)
      
      xPoint <- cos(pointy.at[i])*sqrt(2)
      yPoint <- sin(pointy.at[i])*sqrt(2)
      
      baseX <- c(xPoint,cos(seq(pointy.at[i]+(pi/4),pointy.at[i]-(pi/4)+2*pi,0.1)))
      baseY <- c(yPoint,sin(seq(pointy.at[i]+(pi/4),pointy.at[i]-(pi/4)+2*pi,0.1)))
    }else{
      baseX <- cos(seq(-pi,pi,0.1))
      baseY <- sin(seq(-pi,pi,0.1))
    }
    polygon(x = baseX*rad[i]+x[i], y = baseY*rad[i]+y[i],col=pointyCircles.colors[i],...)
  }
}

#' @rdname trajectoryConvergencePlot
#' @param tau The tau value outputted by the Mann-Kendall test in function \code{\link{trajectoryConvergence}}.
#' @param conv.color The color used to mark convergent trajectories. Defaults to "red".
#' @param div.color The color used to mark divergent trajectories. Defaults to "blue".
#' @noRd
#' @keywords internal
ConvDivColor <- function (tau,
                          conv.color = "red",
                          div.color = "blue"){
  if (tau>0){
    colBase <- t(col2rgb(div.color))/255
    output <- rgb((1-colBase)*((1-abs(tau))^2)+colBase,maxColorValue = 1)
  }else{
    colBase <- t(col2rgb(conv.color))/255
    output <- rgb((1-colBase)*((1-abs(tau))^2)+colBase,maxColorValue = 1)
  }
}


