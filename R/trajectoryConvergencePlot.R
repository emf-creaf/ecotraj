trajectoryConvergencePlot <- function (x,
                                       conv.test="pairwise.asymmetric",
                                       col.Traj="grey",
                                       col.border.Traj="black",
                                       lwd.border.Traj=1,
                                       pointy=F,
                                       names.Traj=NULL,
                                       col.names.Traj="black",
                                       top="between",
                                       radius=1,
                                       col.conv="red",
                                       col.div="blue",
                                       tau.bands.transp=0.3,
                                       alpha.filter=NULL,
                                       CS=NULL,
                                       CSinfConf=NULL,
                                       CSsupConf=NULL,
                                       coeffCS=0.5,
                                       lwd.arrows.CS=3){
  
  widthMult <- 1#this is a multiplication coefficient for the width of the "tau bands" it will change to make room to display more links in case "both" is selected in 
  if (conv.test=="pairwise.asymmetric"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.asymmetric")
  }else if (conv.test=="pairwise.symmetric"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric")
  }else if (conv.test=="both"){
    ConvTest <- trajectoryConvergence(x,type="pairwise.symmetric")
    ConvTestAsym <- trajectoryConvergence(x,type="pairwise.asymmetric")
    widthMult <- 1/3
  }else{
    stop("conv.test invalid value for conv.test")
  }
  
  if (is.null(names.Traj)){
    names.Traj <- colnames(ConvTest$tau)
  }
  
  nTraj <- ncol(ConvTest$tau)
  #find the angles and centers of the shapes representing the trajectories
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
  if (conv.test=="both"){
    for (i in rownames(ConvTest$tau)){
      for (j in colnames(ConvTest$tau)){
        if (i!=j){
          if (is.null(alpha.filter)||(ConvTestAsym$p.value[i,j]<alpha.filter)){
            tau <- ConvTestAsym$tau[i,j]
            tauabs <- abs(tau)
            taucol <- ConvDivColor(tau,ConvCol=col.conv,DivCol=col.div)
            
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
            
            Xcorner3 <- vecY*radius*tauabs*widthMult+centersXj-(vecX*radius+2*vecX*radius*tauabs*widthMult)
            Ycorner3 <- -vecX*radius*tauabs*widthMult+centersYj-(vecY*radius+2*vecY*radius*tauabs*widthMult)
            
            Xcorner4 <- 1.5*vecY*radius*tauabs*widthMult+centersXj-(vecX*radius+2*vecX*radius*tauabs*widthMult)
            Ycorner4 <- -1.5*vecX*radius*tauabs*widthMult+centersYj-(vecY*radius+2*vecY*radius*tauabs*widthMult)
            
            Xcorner5 <- centersXj-vecX*radius
            Ycorner5 <- centersYj-vecY*radius
            
            #finally draw the half-arrow!
            polygon(x=c(Xcorner1,Xcorner2,Xcorner3,Xcorner4,Xcorner5),
                    y=c(Ycorner1,Ycorner2,Ycorner3,Ycorner4,Ycorner5),
                    col=rgb(t(col2rgb(taucol)/255),alpha=tau.bands.transp),border=taucol,lty=2)
          }
        }
      }
    }
  }
  if (conv.test!="pairwise.asymmetric"){
    count <- 1
    for (i in rownames(ConvTest$tau)[1:(nrow(ConvTest$tau)-1)]){
      count <- count+1
      for (j in colnames(ConvTest$tau)[count:ncol(ConvTest$tau)]){
        if (is.null(alpha.filter)||(ConvTest$p.value[i,j]<alpha.filter)){
          tau <- ConvTest$tau[i,j]
          tauabs <- abs(tau)
          taucol <- ConvDivColor(tau,ConvCol=col.conv,DivCol=col.div)
          
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
                  col=rgb(t(col2rgb(taucol)/255),alpha=tau.bands.transp),border=taucol)
        }
      }
    }
  }
  if (conv.test=="pairwise.asymmetric"){
    for (i in rownames(ConvTest$tau)){
      for (j in colnames(ConvTest$tau)){
        if (i!=j){
          if (is.null(alpha.filter)||(ConvTest$p.value[i,j]<alpha.filter)){
            tau <- ConvTest$tau[i,j]
            tauabs <- abs(tau)
            taucol <- ConvDivColor(tau,ConvCol=col.conv,DivCol=col.div)
            
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
            
            Xcorner3 <- vecY*radius*tauabs*widthMult+centersX[j]-(vecX*radius+2*vecX*radius*tauabs*widthMult)
            Ycorner3 <- -vecX*radius*tauabs*widthMult+centersY[j]-(vecY*radius+2*vecY*radius*tauabs*widthMult)
            
            Xcorner4 <- 1.5*vecY*radius*tauabs*widthMult+centersX[j]-(vecX*radius+2*vecX*radius*tauabs*widthMult)
            Ycorner4 <- -1.5*vecX*radius*tauabs*widthMult+centersY[j]-(vecY*radius+2*vecY*radius*tauabs*widthMult)
            
            Xcorner5 <- centersX[j]-vecX*radius
            Ycorner5 <- centersY[j]-vecY*radius
            
            #finally draw the half-arrow!
            polygon(x=c(Xcorner1,Xcorner2,Xcorner3,Xcorner4,Xcorner5),
                    y=c(Ycorner1,Ycorner2,Ycorner3,Ycorner4,Ycorner5),
                    col=rgb(t(col2rgb(taucol)/255),alpha=tau.bands.transp),border=taucol)
          }
        }
      }
    }
  }
  
  #Draw the cycle's dates
  PointyCircle(centersX,centersY,
               col.fill=col.Traj,col.border=col.border.Traj,rad=radius,lwd.border=lwd.border.Traj,
               pointy=pointy,pointy.at=angles-pi/2)
  text(centersX,centersY,names.Traj,col=col.names.Traj,font=2)
  
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


PointyCircle <- function (x,
                          y,
                          rad=1,
                          pointy=FALSE,
                          pointy.at=NULL,
                          col.fill="grey",
                          col.border="black",
                          lwd.border=1){
  if (length(y)!=length(x)){
    stop("x and y must have the same length")
  }
  if (length(x)>1){
    if (length(rad)==1){
      rad <- rep(rad,length(x))
    }else if (length(rad)!=length(x)){
      stop("rad must have the same length as x or length 1")
    }
    if (length(col.fill)==1){
      col.fill <- rep(col.fill,length(x))
    }else if (length(col.fill)!=length(x)){
      stop("col.fill must have the same length as x or length 1")
    }
    if (length(col.border)==1){
      col.border <- rep(col.border,length(x))
    }else if (length(col.border)!=length(x)){
      stop("col.border must have the same length as x or length 1")
    }
    if (length(lwd.border)==1){
      lwd.border <- rep(lwd.border,length(x))
    }else if (length(lwd.border)!=length(x)){
      stop("lwd.border must have the same length as x or length 1")
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
    polygon(x = baseX*rad[i]+x[i], y = baseY*rad[i]+y[i],
            col=col.fill[i],border=col.border[i],lwd=lwd.border[i])
  }
}


ConvDivColor <- function (x,
                          ConvCol="red",
                          DivCol="blue"){
  if (x>0){
    colBase <- t(col2rgb(DivCol))/255
    output <- rgb((1-colBase)*((1-abs(x))^2)+colBase,maxColorValue = 1)
  }else{
    colBase <- t(col2rgb(ConvCol))/255
    output <- rgb((1-colBase)*((1-abs(x))^2)+colBase,maxColorValue = 1)
  }
}


