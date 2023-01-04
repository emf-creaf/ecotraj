.envelopeVar<-function(d, sites, surveys, ...){
  D_T <- trajectoryDistances(d = d, sites = sites, surveys = surveys, ...)
  r <- ncol(as.matrix(D_T))
  return(sum(as.vector(D_T)^2)/(r*(r-1)))
}
#' Reference envelopes
#' 
#' 
#' @encoding UTF-8
#' @name envelopes
#' @aliases segmentDistances trajectoryDistances trajectoryLengths trajectoryLengths2D trajectoryAngles trajectoryAngles2D
#'          trajectoryProjection trajectoryConvergence trajectoryDirectionality 
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states (see details).
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param surveys A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @param ... Additional parameters for function \code{\link{trajectoryDistances}}
#' 
#' @details 
#' 
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' @author Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la Baie de Saint-Brieuc
#' 
#' @references
#' 
#' @seealso \code{\link{trajectorymetrics}} 
#' 
#' @examples 
#'  
envelopeVariability<-function(d, sites, surveys = NULL, nboot.ci = NULL, alpha.ci = 0.05, ...){
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  if(is.null(nboot.ci)){
    return(.envelopeVar(d = d, sites = sites, surveys = surveys, ...))
  } else {
    if(is.integer(nboot.ci)) stop("'nboot.ci' must be an integer")
    if(nboot.ci < 1) stop("'nboot.ci' must be larger than 0")
    var_e <- .envelopeVar(d = d, sites = sites, surveys = surveys, ...)
    unique_sites <- unique(sites)
    nsites <- length(unique_sites)
    if(is.null(surveys)) {
      surveys <- numeric(0)
      for(i in 1:nsites) {
        surveys<-c(surveys, 1:sum(sites==sites[i]))
      }
    }    
    var_e_boot<-numeric(nboot.ci)    
    for(i in 1:nboot.ci){
      bsample <- sample(1:nsites, nsites, replace = TRUE)
      b_sites<-numeric(0)
      b_indices <- numeric(0)
      b_surveys<-numeric(0)
      for(j in 1:length(bsample)) {
        bsj <- bsample[j]
        nt <- sum(sites==bsj)
        b_indices <-c(b_indices, which(sites==bsj))
        b_sites <- c(b_sites, rep(j, nt))
        b_surveys <- c(b_surveys, surveys[sites==bsj])
      }
      b_d <- as.dist(as.matrix(d)[b_indices, b_indices])
      var_e_boot[i] <- .envelopeVar(d = b_d, sites = b_sites, surveys = b_surveys, ...)
    }
    return(c(Var = var_e, quantile(var_e_boot, probs = c((alpha.ci/2.0), 1.0 - (alpha.ci/2.0)))))
  }
}
compareAgainstEnvelope<-function(d, sites, surveys, envelope_sites, fuzzy_exponent, nboot = NULL, ...) {
  
}