#' Distances between segments
#'
#' Calculates the distance between pairs of directed segments of input trajectories.
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param distance.type The type of distance index to be calculated (Besse et al. 2016):
#'   \itemize{
#'     \item{\code{Hausdorff}: Hausdorff distance between two segments.}
#'     \item{\code{directed-segment}: Directed segment distance (default).}
#'     \item{\code{PPA}: Perpendicular-parallel-angle distance.}
#'   }
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#'
#' @details 
#' Details of calculations are given in De \enc{Cáceres}{Caceres} et al (2019). 
#' The input distance matrix \code{d} should ideally be metric. That is, all subsets of distance triplets should fulfill the triangle inequality (see utility function \code{\link{is.metric}}). 
#' All ETA functions that require metricity include a parameter '\code{add}', which by default is TRUE, meaning that whenever the triangle inequality is broken the minimum constant required to fulfill it is added to the three distances.
#' If such local (an hence, inconsistent across triplets) corrections are not desired, users should find another way modify \code{d} to achieve metricity, such as PCoA, metric MDS or non-metric MDS (see vignette 'Introduction to Ecological Trajectory Analysis'). 
#' If parameter '\code{add}' is set to FALSE and problems of triangle inequality exist, ETA functions may provide missing values in some cases where they should not.
#' 
#' @returns A list with the following elements:
#' \itemize{
#'   \item{\code{Dseg}: Distance matrix between segments.}
#'   \item{\code{Dini}: Distance matrix between initial points of segments.}
#'   \item{\code{Dfin}: Distance matrix between final points of segments.}
#'   \item{\code{Dinifin}: Distance matrix between initial points of one segment and the final point of the other.}
#'   \item{\code{Dfinini}: Distance matrix between final points of one segment and the initial point of the other.}
#' }
#' 
#' @export
#'
#' @references
#' Besse, P., Guillouet, B., Loubes, J.-M. & François, R. (2016). Review and perspective for distance based trajectory clustering. IEEE Trans. Intell. Transp. Syst., 17, 3306–3317.
#' 
#' De \enc{Cáceres}{Caceres} M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R & Hubbell S. (2019). 
#' Trajectory analysis in community ecology. Ecological Monographs 89, e01350.
#' 
segmentDistances<-function(x, distance.type ="directed-segment", add = TRUE) {
  distance.type <- match.arg(distance.type, c("directed-segment", "Hausdorff", "PPA"))
  if(!inherits(x, "trajectories")) stop("'x' should be of class `trajectories`")
  
  d <- x$d
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
  
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  nsurveysite<-numeric(nsite)
  for(i in 1:nsite) {
    nsurveysite[i] = sum(sites==siteIDs[i])
  }
  if(sum(nsurveysite==1)>0) stop("All sites need to be surveyed at least twice")
  dmat <- as.matrix(d)
  n <- nrow(dmat)
  nseg <- sum(nsurveysite)-nsite
  segnames <- character(nseg)
  cnt <- 1
  for(i in 1:nsite) {
    if(!is.null(surveys)) {
      surv <- surveys[sites==siteIDs[i]]
      surv <- sort(surv) #Surveys may not be in order
    }
    else surv <- 1:nsurveysite[i]
    for(j in 1:(nsurveysite[i]-1)) {
      segnames[cnt] <- paste0(siteIDs[i],"[",surv[j],"-",surv[j+1],"]")
      cnt <- cnt+1
    }
  }
  dsegmat <- matrix(0, nseg, nseg)
  rownames(dsegmat) <-segnames
  colnames(dsegmat) <-segnames
  dinisegmat <- dsegmat
  dfinsegmat <- dsegmat
  dinifinsegmat <- dsegmat
  
  os1 <- 1
  for(i1 in 1:nsite) {
    ind_surv1 <- which(sites==siteIDs[i1])
    #Surveys may not be in order
    if(!is.null(surveys)) ind_surv1 <- ind_surv1[order(surveys[sites==siteIDs[i1]])]
    for(s1 in 1:(nsurveysite[i1]-1)) {
      os2 <- 1
      for(i2 in 1:nsite) {
        ind_surv2 <- which(sites==siteIDs[i2])
        #Surveys may not be in order
        if(!is.null(surveys)) ind_surv2 <- ind_surv2[order(surveys[sites==siteIDs[i2]])]
        for(s2 in 1:(nsurveysite[i2]-1)) {
          #Select submatrix from dmat
          ind12 <- c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])
          # print(ind12)
          # print(c(os1, os2))
          dmat12 <- dmat[c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1]),
                         c(ind_surv1[s1],ind_surv1[s1+1],ind_surv2[s2],ind_surv2[s2+1])]
          dsegmat[os1,os2] <- .twoSegmentDistanceC(dmat12, type=distance.type, add)
          dsegmat[os2,os1] <- dsegmat[os1,os2]
          dinisegmat[os2,os1] <- dinisegmat[os1,os2]<-dmat[ind_surv1[s1],ind_surv2[s2]]
          dfinsegmat[os2,os1] <- dfinsegmat[os1,os2]<-dmat[ind_surv1[s1+1],ind_surv2[s2+1]]
          dinifinsegmat[os1,os2] <- dmat[ind_surv1[s1],ind_surv2[s2+1]]
          dinifinsegmat[os2,os1] <- dmat[ind_surv1[s1+1],ind_surv2[s2]]
          os2 <- os2+1
        }
      }
      os1 <- os1+1
    }
  }
  
  return(list(Dseg = as.dist(dsegmat), Dini=as.dist(dinisegmat), Dfin = as.dist(dfinsegmat),
              Dinifin=dinifinsegmat))
}