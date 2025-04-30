#' Dynamic variation and variation decomposition
#'
#' \itemize{
#'   \item{Function \code{dynamicVariation} assesses the amount of dynamic variation observed across trajectories and the relative contribution of each of them.}
#'   \item{Function \code{variationDecomposition} performs a sum of squares decomposition of total variation in three components: (1) across trajectories (entities); (2) across time points;
#' (3) their interaction.}
#' }
#' 
#' @param x An object of class \code{trajectories} (or its children subclasses \code{fd.trajectories} or \code{cycles}).
#' @param ... Additional params to be passed to function \code{\link{trajectoryDistances}}.
#'
#' @details
#' Function \code{variationDecomposition} requires trajectories to be synchronous. The SS sum of \code{temporal} and \code{interaction} components correspond to the SS sum, across trajectories, of 
#' function \code{\link{trajectoryInternalVariation}}. 
#' 
#' @seealso \code{\link{defineTrajectories}}, \code{\link{is.synchronous}}, \code{\link{trajectoryDistances}}, \code{\link{trajectoryInternalVariation}}
#'
#' @returns 
#' \itemize{
#'  \item{Function \code{dynamicVariance} returns a list with three elements (dynamic sum of squares, dynamic variance and a vector of trajectory relative contributions)}
#'  \item{Function \code{variationDecomposition} returns a data frame with results (sum of squares, degrees of freedom and variance estimates) for each variance component and the total.}
#' }
#' 
#' @examples
#' #Description of entities and surveys
#' entities <- c("1","1","1","1","2","2","2","2","3","3","3","3")
#' surveys <- c(1,2,3,4,1,2,3,4,1,2,3,4)
#'   
#' #Raw data table
#' xy<-matrix(0, nrow=12, ncol=2)
#' xy[2,2]<-1
#' xy[3,2]<-2
#' xy[4,2]<-3
#' xy[5:6,2] <- xy[1:2,2]
#' xy[7,2]<-1.5
#' xy[8,2]<-2.0
#' xy[5:6,1] <- 0.25
#' xy[7,1]<-0.5
#' xy[8,1]<-1.0
#' xy[9:10,1] <- xy[5:6,1]+0.25
#' xy[11,1] <- 1.0
#' xy[12,1] <-1.5
#' xy[9:10,2] <- xy[5:6,2]
#' xy[11:12,2]<-c(1.25,1.0)
#' 
#' d <- dist(xy)
#' 
#' # Defines trajectories
#' x <- defineTrajectories(d, entities, surveys)
#' 
#' # Assessment of dynamic variation and individual trajectory contributions
#' dynamicVariation(x)
#' 
#' # Variation decomposition (entity, temporal and interaction) for synchronous 
#' # trajectories:
#' variationDecomposition(x)
#' 
#' # check the correspondence with internal variation
#' sum(variationDecomposition(x)[c("time", "interaction"),"ss"])
#' sum(trajectoryInternalVariation(x)$internal_ss)
#'
#' @name dynamicVariation
#' 
#' @export
dynamicVariation<-function(x, ...) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  D <- trajectoryDistances(x, ...)
  ss <- diag(.gowerCentered(D))
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  varDT <- sum(ss)/(nsite-1)
  RC <- ss/sum(ss)
  names(RC) <- siteIDs
  return(list(dynamic_ss = sum(ss), dynamic_variance = varDT, relative_contributions = RC))
}

#' @rdname dynamicVariation
#' @export
variationDecomposition<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  if(!is.synchronous(x)) stop("Trajectories need to be synchronous (all sites need to be surveyed the same times)")
  surveys <- x$metadata$surveys
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  
  d <- x$d
  Dmat <-as.matrix(d)
  n <- nrow(Dmat)
  I <- diag(n)
  
  # Anderson (2017). Permutational Multivariate Analysis of Variance (PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article ID: stat07841.
  G <- .gowerCentered(d)
  #model matrix
  M1 <- model.matrix(~a,data.frame(a = factor(sites)), contrasts = list(a = "contr.helmert"))
  M2 <- model.matrix(~a,data.frame(a = factor(surveys)), contrasts = list(a = "contr.helmert"))
  H1 <- M1%*%MASS::ginv(t(M1)%*%M1)%*%t(M1)
  H2 <- M2%*%MASS::ginv(t(M2)%*%M2)%*%t(M2)
  ss_total <- sum(diag(G))
  ss_entities <- sum(diag(t(H1)%*%G%*%H1))
  ss_time <- sum(diag(t(H2)%*%G%*%H2))
  ss_interaction <- sum(diag((I-t(H1) - t(H2))%*%G%*%(I-H1 - H2)))
  df <- data.frame(ss = rep(NA, 4), df = rep(NA, 4), variance = rep(NA,4))
  df$ss = c(ss_entities, ss_time, ss_interaction, ss_total)
  df$df = c(nsite-1, nsurveysite[1]-1, (nsite -1)*(nsurveysite[1]-1), nsite*nsurveysite[1]-1)
  df$variance = df$ss/df$df
  row.names(df) <-c("entities", "time", "interaction", "total")
  return(df)
}
