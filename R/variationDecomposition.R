#' Variation decomposition
#'
#' Performs a sum of squares decomposition of variation in three components: (1) across trajectory entities; (2) across time;
#' (3) interaction.
#' 
#' @param x An object of class \code{trajectories} (or its children subclasses \code{fd.trajectories} or \code{cycles}) where
#' trajectories are synchronous.
#'
#' @seealso \code{\link{defineTrajectories}}, \code{\link{is.synchronous}}
#' @returns A data frame with results (sum of squares, degrees of freedom and variance estimates) for each variance component and the total.
#' @export
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
#' d <- dist(xy)
#' 
#' # Defines trajectories
#' x <- defineTrajectories(d, sites, surveys)
#' 
#' # Variation decomposition
#' variationDecomposition(x)
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
  df <- data.frame(ss = rep(NA, 4), df = rep(NA, 4), var = rep(NA,4))
  df$ss = c(ss_entities, ss_time, ss_interaction, ss_total)
  df$df = c(nsite-1, nsurveysite[1]-1, (nsite -1)*(nsurveysite[1]-1), nsite*nsurveysite[1]-1)
  df$var = df$ss/df$df
  row.names(df) <-c("entities", "time", "interaction", "total")
  return(df)
}