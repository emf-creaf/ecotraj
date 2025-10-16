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
                                "red2","red3",
                                "blue2","blue3",
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
  
  if (legend==T){
    
    names(relationships.colors) <- substr(names(relationships.colors),
                                          start=1,
                                          stop=(unlist(lapply(
                                            strsplit(names(relationships.colors), ''),
                                            function(x) which(x == "(")))-2))
    
    relationships.colors <- relationships.colors[unique(names(relationships.colors))]
    legend(x=ncols,y=nrows,xpd=NA,bty="n",cex=0.7,
           fill=relationships.colors,names(relationships.colors))
  }
}
