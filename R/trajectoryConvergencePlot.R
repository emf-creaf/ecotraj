#' Summary plot for trajectory convergence and divergence
#' 
#' Provides plots to represent trajectory convergence and divergence tests performed by the function \code{\link{trajectoryConvergence}}.
#' 
#' @encoding UTF-8
#' @name trajectoryCyclicalPlots
#' 
#' @details
#' The function \code{trajectoryConvergencePlot} provides ways to visualize pairwise convergence and divergence between trajectories using calls to function \code{\link{trajectoryConvergence}} which performs the tests.
#' In the plots, trajectories are represented by circles. The convergence or divergence between pairs of trajectories are represented by links. If convergence tests are symmetric, the links are simple. If the convergence tests are asymmetric, the links are displayed as half arrows pointing from the trajectory converging or diverging towards the trajectory being approached or diverged from.
#' The width and color hue of the links are proportional to the tau statistic of the Mann.Kendall test performed by the \code{\link{trajectoryConvergence}} function. 
#' The function \code{trajectoryConvergencePlot} also offers the possibility to plot both tests at the same time.
#' 
#' The function finally offers some possibilities more relevant in the context of cyclical ecological trajectory analysis (CETA) such as the possibility to make the circles representing trajectories "pointy" so that they communicate a cyclical organization. This is useful if fixed dates trajectories are being studied (see \code{\link{extractFixedDateTrajectories}}).
#' The function also allows to display "cyclical shifts" but values have to be provided by the user as they can be obtained in a variety of ways summarizing the outputs of \code{\link{cycleShifts}}. 
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' 
#' @seealso \code{\link{trajectoryConvergence}}
#' 
#' @examples
#' data("avoca")
#' avoca_D_man <- vegclust::vegdiststruct(avoca_strat, 
#'                                        method="manhattan", 
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
#' #More refined output with both type of tests and only plotting significant test results (p-value < 0.05)
#' trajectoryConvergencePlot(avoca_x,type="both",alpha.filter=0.05)
#' 
#' #Much more refined output with nicer colors bigger half arrows and personalized trajectory names
#' trajectoryConvergencePlot(avoca_x,type="both",alpha.filter=0.05,
#'                           #adjusting the attributes of the links between trajectories
#'                           half.arrows.size = 1.5,conv.color = "orangered",div.color="dodgerblue",
#'                           #controling the size of circles representing trajectories (links are enlarged to match it)
#'                           radius=1.2, 
#'                           #customizing the border of the circles representing trajectories
#'                           traj.colors = "black",border="white",lwd=2,
#'                           #giving trajectories new names and adjusting the names color
#'                           traj.names=LETTERS[1:8],traj.names.colors="white")
#'
#' 
#' @rdname trajectoryConvergencePlot
#' @param x An object of class \code{\link{trajectories}}.
#' @param type A string indicating the convergence test, either "pairwise.asymmetric", "pairwise.symmetric" or "both" (see \code{\link{trajectoryConvergence}}).
#' @param alpha.filter The minimum p-value for a link to be drawn (see \code{\link{trajectoryConvergence}}). Defaults to NULL (all links drawn).
#' @param traj.colors The colors for the trajectories (circles). Defaults to "grey".
#' @param traj.names The names of trajectories. Defaults to the names provided in \code{x}.
#' @param traj.names.colors The color of the names of trajectories on the circles. Defaults to "black".
#' @param ... Additional parameters passed to \code{\link{polygon}} to personalize the circles representing trajectories.
#' @param radius The radius of the circles representing trajectories. Defaults to 1.
#' @param conv.color The color used to mark convergent trajectories. Defaults to "red".
#' @param div.color The color used to mark divergent trajectories. Defaults to "blue".
#' @param half.arrows.size A multiplication coefficient for the size of the arrow heads when representing asymmetric tests results. Defaults to 1.
#' @param tau.links.transp The transparency of the links representing the tau statistic of the Mann.Kendall test (see \code{\link{trajectoryConvergence}}).
#' @param top A string indicating whether the top of the plotting area should contain a circle representing a trajectory ("top") or should be in between two circles ("between"). Defaults to "between".
#' @param pointy Boolean. Should the circles representing trajectories be made pointy? Useful in the context of CETA to represent fixed date trajectories (see \code{\link{trajectoryCyclical}}).
#' @param CS Experimental. Allows to add arrows representing cyclical shifts in the context of CETA.
#' @param CSinfConf Experimental. Lower confidence intervals for cyclical shifts.
#' @param CSsupConf Experimental. Upper confidence intervals for cyclical shifts.
#' @param coeffCS Experimental. A multiplication coefficient for the arrows representing cyclical shifts.
#' @param lwd.arrows.CS Experimental. The width of the arrows representing cyclical shifts.
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
                                       CS = NULL,
                                       CSinfConf = NULL,
                                       CSsupConf = NULL,
                                       coeffCS = 0.5,
                                       lwd.arrows.CS = 3){
  
  widthMult <- 1#this is a multiplication coefficient for the width of the "tau links" it will change to make room to display more links in case "both" is selected
  if (type=="pairwise.asymmetric"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.asymmetric")
  }else if (type=="pairwise.symmetric"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric")
  }else if (type=="both"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric")
    ConvTestAsym <- trajectoryConvergence(x,type="pairwise.asymmetric")
    widthMult <- 1/3
  }else{
    stop("invalid value for type")
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
  
  #Draw the cycle's dates
  PointyCircle(centersX,centersY,pointyCircles.colors=traj.colors,
               rad=radius,pointy=pointy,pointy.at=angles-pi/2,...)
  text(centersX,centersY,traj.names,col=traj.names.colors,font=2)
  
  #Add arrows for cyclical shifts if asked!
  if (!is.null(CS)){
    baseArrowsX <- centersX*(1+radius)
    baseArrowsY <- centersY*(1+radius)
    if ((!is.null(CSinfConf))&&(!is.null(CSsupConf))){
      arrows(x0=baseArrowsY*CSinfConf*coeffCS+baseArrowsX,
             y0=-baseArrowsX*CSinfConf*coeffCS+baseArrowsY,
             x1=baseArrowsY*CSsupConf*coeffCS+baseArrowsX,
             y1=-baseArrowsX*CSsupConf*coeffCS+baseArrowsY,
             angle=90,length=0.05,col="grey30",code=3,lwd=lwd.arrows.CS-1)
    }
    arrows(x0=baseArrowsX,y0=baseArrowsY,
           x1=baseArrowsY*CS*coeffCS+baseArrowsX,
           y1=-baseArrowsX*CS*coeffCS+baseArrowsY,
           lwd=lwd.arrows.CS,length=0.1)
    points(baseArrowsX,baseArrowsY,pch=16)
  }
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
    polygon(x = baseX*rad[i]+x[i], y = baseY*rad[i]+y[i],col=pointyCircles.colors,...)
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


