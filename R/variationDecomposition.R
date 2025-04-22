#' Title
#'
#' @param x 
#'
#' @returns
#' @export
#'
#' @examples
variationDecomposition<-function(x) {
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
  surveys <- x$metadata$surveys
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
    exclude <- unique(c(exclude,which(x$metadata$internal==FALSE)))
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
    exclude <- unique(c(exclude,which(x$metadata$internal==FALSE)))
  } else {
    sites <- x$metadata$sites
  }
  siteIDs <- unique(sites)
  nsite <- length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) nsurveysite[i] = sum(sites==siteIDs[i])
  if(length(unique(nsurveysite))>1) stop("All sites need to be surveyed the same number of times")
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
  ss_space <- sum(diag(t(H1)%*%G%*%H1))
  ss_time <- sum(diag(t(H2)%*%G%*%H2))
  ss_spacetime <- sum(diag((I-t(H1) - t(H2))%*%G%*%(I-H1 - H2)))
  df <- data.frame(ss = rep(NA, 4), df = rep(NA, 4), var = rep(NA,4))
  df$ss = c(ss_space, ss_time, ss_spacetime, ss_total)
  df$df = c(nsite-1, nsurveysite[1]-1, (nsite -1)*(nsurveysite[1]-1), nsite*nsurveysite[1]-1)
  df$var = df$ss/df$df
  row.names(df) <-c("space", "time", "spacetime", "total")
  return(df)
}