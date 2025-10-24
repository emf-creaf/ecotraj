#' Heat map-like plots for Relative Trajectory Movement Assessment (RTMA)
#'
#' Function \code{trajectoryRMAPlot} provides heat map-like plots for Relative Trajectory Movement Assessment (RTMA) performed by function \code{\link{trajectoryRMA}}.
#' 
#' @param x An object of class \code{\link{RTMA}}.
#' @param mode The mode of trajectory relationship display (see details). Defaults to \code{"full"}, ignored if \code{relationships.colors} is specified.
#' @param relationships.colors Vector of user-chosen colors to represent trajectory relationships. Must have specific properties (see details). Overrides \code{mode}. 
#' @param traj.names The names of trajectories. Defaults to the names provided in \code{x}.
#' @param order.row A re-ordering (and potential selection) of the rows of the output. If provided without \code{order.col}, the re-ordering is also applied to columns.
#' @param order.col If \code{order.row} is provided, this parameter allows to set a different order or selection for columns of the final display (allowing rectangular plots).
#' @param vertical Flag to indicate if the top trajectory names should be rotated 90°, defaults to \code{FALSE}.
#' @param legend Flag to indicate if the legend should be plotted, defaults to \code{FALSE}.
#'
#' @details
#' Function \code{trajectoryRMAPlot} provides heat map-like plots for Relative Trajectory Movement Assessment (RTMA).
#' A key feature of the function is its different \code{mode} of representation, allowing to put more or less emphasis on some aspect of trajectories relative movements.
#' The 12 relative movement relationships recognized by RTMA may belong to three higher-order groups: the convergence group, the divergence group and the oriented group (see \code{\link{trajectoryRMA}} for more details).
#' The parameter \code{mode} allows to display targeted groups or combination of groups instead of the detailed relationships. Possible values for \code{mode} are:
#' \itemize{
#'     \item{\code{"full"}: Default value. Display the finest level relationships.}
#'     \item{\code{"convdiv"}: Displays and groups relationships belonging to the convergence and divergence groups.}
#'     \item{\code{"oriented"}: Displays and groups relationships belonging to the oriented group.}
#'     \item{\code{"crossed.groups"}: Displays and groups relationships belonging simultaneously to the oriented group and either the convergence or the divergence group.}
#'     \item{\code{"convdiv.complete"}: As \code{"convdiv"} but adding in details the relationships that are neither from the convergence or divergence group.}
#'     \item{\code{"oriented.complete"}: As \code{"oriented"} but adding in details the relationships that are not from the oriented group.}
#'     \item{\code{"crossed.groups.complete"}: As \code{"crossed.groups"} but adding in details the relationships that do not belong to both the oriented group and either the convergence or divergence group.}
#' }
#' Relationships belonging to the oriented group are asymmetric. Practically, this means that one trajectory is in front while the other is in the back.
#' In the \code{trajectoryRMAPlot} output, crossed cells indicate that the corresponding ROW trajectory is the trajectory in front.
#' 
#' COLORS: Each \code{mode} comes with default colors for the heat map. Nonetheless, the parameter \code{relationships.colors} allows user-defined colors instead.
#' The vectors of colors provided in \code{relationships.colors} must have length 21 with \code{\link{names}} corresponding to the names of the different relationships recognized by RTMA
#' (can be found in a \code{\link{RTMA}} object \code{x} with \code{x$dynamic_relationships_taxonomy$dynamic_relationship}).
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @references
#' Djeghri et al. (in preparation) Uncovering the relative movements of ecological trajectories.
#'
#' @seealso \code{\link{trajectoryRMA}}, \code{\link{trajectoryConvergencePlot}} 
#'
#' @examples
#' #Prepare data
#' data("avoca")
#' avoca_D_man <- vegclust::vegdiststruct(avoca_strat, 
#'                                        method ="manhattan", 
#'                                        transform = function(x){log(x+1)})
#' years <- c(1971, 1974, 1978, 1983, 1987, 1993, 1999, 2004, 2009)
#' avoca_times <- years[avoca_surveys]
#' avoca_x <- defineTrajectories(d = avoca_D_man,  
#'                               sites = avoca_sites, 
#'                               times = avoca_times)
#' #Perform RTMA
#' avoca_RTMA <- trajectoryRMA(avoca_x)
#' 
#' #Default (full) output
#' par(mar=c(2,2,2,12)) #New margin to accomodate the legend
#' trajectoryRMAPlot(avoca_RTMA,legend=TRUE)
#' 
#' #Play with different visualization modes of relationship groups
#' trajectoryRMAPlot(avoca_RTMA,mode="convdiv",legend=TRUE)
#' trajectoryRMAPlot(avoca_RTMA,mode="oriented",legend=TRUE)
#' trajectoryRMAPlot(avoca_RTMA,mode="crossed.groups",legend=TRUE)
#' 
#' 
#' @name trajectoryRMAPlot
#' @export
trajectoryRMAPlot <- function(x,
                              mode = "full",
                              relationships.colors = NULL,
                              traj.names = NULL,
                              order.row = NULL,
                              order.col = NULL,
                              vertical = FALSE,
                              legend = FALSE
                              ){
  if (!inherits(x,"RTMA")) stop("'x' must be of class RTMA.")
  if(is.null(relationships.colors)){
    if (mode%in%c("full",
                  "oriented.complete",
                  "convdiv.complete",
                  "crossed.groups.complete")){
      relationships.colors <- c("white","grey70","grey30",
                                "red3","red2",
                                "blue3","blue2",
                                "orange","orange",
                                "deepskyblue","deepskyblue",
                                "darkorange2",
                                "darkorange2",
                                "dodgerblue2",
                                "dodgerblue2",
                                "yellow2","yellow2",
                                "green3","green3",
                                "aquamarine2","aquamarine2")
      if (mode=="oriented.complete"){
        relationships.colors[!is.na(x$dynamic_relationships_taxonomy$oriented_group)] <- "green3"
      }else if (mode=="convdiv.complete"){
        relationships.colors[grep("convergence",x$dynamic_relationships_taxonomy$conv_div_group)] <- "red2"
        relationships.colors[grep("divergence",x$dynamic_relationships_taxonomy$conv_div_group)] <- "blue2"
      }else if (mode=="crossed.groups.complete"){
        conv <- grep("convergence",x$dynamic_relationships_taxonomy$conv_div_group)
        div <- grep("divergence",x$dynamic_relationships_taxonomy$conv_div_group)
        oriented <- grep("oriented",x$dynamic_relationships_taxonomy$oriented_group)
        relationships.colors[intersect(conv,oriented)] <- "yellow2"
        relationships.colors[intersect(div,oriented)] <- "aquamarine2"
      }
    }else if(mode=="oriented"){
      relationships.colors <- rep(NA,nrow(x$dynamic_relationships_taxonomy))
      relationships.colors[!is.na(x$dynamic_relationships_taxonomy$oriented_group)] <- "green3"
    }else if (mode=="convdiv"){
      relationships.colors <- rep(NA,nrow(x$dynamic_relationships_taxonomy))
      relationships.colors[grep("convergence",x$dynamic_relationships_taxonomy$conv_div_group)] <- "red2"
      relationships.colors[grep("divergence",x$dynamic_relationships_taxonomy$conv_div_group)] <- "blue2"
    }else if(mode=="crossed.groups"){
      relationships.colors <- rep(NA,nrow(x$dynamic_relationships_taxonomy))
      conv <- grep("convergence",x$dynamic_relationships_taxonomy$conv_div_group)
      div <- grep("divergence",x$dynamic_relationships_taxonomy$conv_div_group)
      oriented <- grep("oriented",x$dynamic_relationships_taxonomy$oriented_group)
      relationships.colors[intersect(conv,oriented)] <- "yellow2"
      relationships.colors[intersect(div,oriented)] <- "aquamarine2"
      
    }else stop("Invalid string for 'mode'.")
    names(relationships.colors) <- x$dynamic_relationships_taxonomy$dynamic_relationship
  }
  if(!is.null(traj.names)){
    colnames(x$dynamic_relationships) <- traj.names
    rownames(x$dynamic_relationships) <- traj.names
  }
  if(!is.null(order.row)){
    if (is.null(order.col)){
      order.col <- order.row
    }
    rels <- x$dynamic_relationships[order.row,order.col]
  }else{
    rels <- x$dynamic_relationships
  }
  
  ncols <- ncol(rels)
  nrows <- nrow(rels)
  
  #prepare square coordinates
  ybottom <- rep((nrows-1):0,ncols)
  ytop <- ybottom+1
  xleft <- rep(0:(ncols-1),each=nrows)
  xright <- xleft+1
  
  #prepare colors
  colors.squares <- relationships.colors[rels]
  colors.squares[is.na(colors.squares)] <- "black"
  
  #plot the squares
  plot(NA,NA,axes=FALSE,ylab="",xlab="",
       ylim=c(0,ncols),xlim=c(0,nrows),asp=1)
  
  rect(xleft=xleft,xright=xright,ybottom=ybottom,ytop=ytop,
       col=colors.squares,border=NA)
  
  
  #add information on oriented relationships
  if (mode!="convdiv"){
    oriented <- x$dynamic_relationships_taxonomy[rels,"oriented_group"]
    if(mode=="crossed.groups"){
      oriented[is.na(x$dynamic_relationships_taxonomy[rels,"conv_div_group"])] <- NA
    }else if(mode=="convdiv.complete"){
      oriented[-grep("pursuit",x$dynamic_relationships_taxonomy[rels,"dynamic_relationship"])] <- NA
    }
    front <- grep("(front)",oriented)
    
    segments(x0=xleft[front],y0=ybottom[front],
             x1=xright[front],y1=ytop[front],col="white",lwd=2,lend=1)
    segments(x0=xleft[front],y0=ybottom[front],
             x1=xright[front],y1=ytop[front],col="black",lwd=1,lend=1)
  }
  
  
  #put clean borders
  rect(xleft=xleft,xright=xright,ybottom=ybottom,ytop=ytop)
  
  #add trajectory names
  if (vertical){
    text(colnames(rels),y=nrows,x=1:ncols-0.5,xpd=NA,srt=90,adj=0)
  }else{
    text(colnames(rels),pos=3,y=nrows,x=1:ncols-0.5,xpd=NA)
  }
  text(rownames(rels),pos=2,x=0,y=nrows:1-0.5,xpd=NA)
  
  #add legend
  if (legend==T){
    names(relationships.colors) <- substr(names(relationships.colors),
                                          start=1,
                                          stop=(unlist(lapply(
                                            strsplit(names(relationships.colors), ''),
                                            function(x) which(x == "(")))-2))
    
    relationships.colors <- relationships.colors[unique(names(relationships.colors))]
    relationships.colors <- relationships.colors[c("neutral","parallel","antiparallel",
                                                   "convergence","weak convergence","approaching","weak approaching",
                                                   "catch-up","pursuit","escape",
                                                   "weak departing","departing","weak divergence","divergence")]
    relationships.colors[is.na(relationships.colors)] <- "black"
    
    leg <- legend(x=ncols,y=nrows,xpd=NA,bty="n",cex=0.7,title="Relative movement\nrelationships",title.font=2,text.font=1,
                  fill=relationships.colors[1:3],names(relationships.colors)[1:3])
    
    leg <- legend(x=ncols,y=leg$rect$top-leg$rect$h,title="Relative movement",title.font=2,text.font=1,title.col="white",
                  xpd=NA,bty="n",cex=0.7,
                  fill=relationships.colors[-(1:3)],names(relationships.colors)[-(1:3)])
    ydist <- (leg$text$y[1]-leg$text$y[2])/2
    rect(xleft=leg$rect$left+ncols/100,xright=leg$rect$left+leg$rect$w,
         ytop=c(leg$text$y[1]+ydist,leg$text$y[7]+ydist),
         ybottom=c(leg$text$y[5]-ydist,leg$text$y[11]-ydist),border=relationships.colors[c(4,14)],xpd=NA,lwd=2)
    text(c("Convergence group","Divergence group"),cex=0.8,col=relationships.colors[c(4,14)],
         pos=4,x=leg$rect$left+ncols/100,y=c(leg$text$y[1]+ydist*3,leg$text$y[11]-ydist*3),xpd=NA,font=4)
    rect(xleft=leg$rect$left+ncols/200,xright=leg$rect$left+leg$rect$w+ncols/200,
         ytop=leg$text$y[3]+ydist,ybottom=leg$text$y[9]-ydist,border=relationships.colors[9],lwd=2,xpd=NA)
    text("Oriented group",srt=90,cex=0.8,y=leg$text$y[6],x=leg$rect$left+leg$rect$w+ncols/200,xpd=NA,col=relationships.colors[9],font=4,adj=c(0.5,1))
  }
}
