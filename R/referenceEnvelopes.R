.trajectoryEnvelopeVar<-function(d, sites, surveys, ...){
  D_T <- trajectoryDistances(d = d, sites = sites, surveys = surveys, ...)
  r <- ncol(as.matrix(D_T))
  return(sum(as.vector(D_T)^2)/(r*(r-1)))
}

.stateEnvelopeVar <-function(d){
  r <- ncol(as.matrix(d))
  return(sum(as.vector(d)^2)/(r*(r-1)))
}
# Generates a bootstrap sample of the envelope variability
.bootstrapSamples<-function(ind, nboot) {
  bs <- matrix(nrow = nboot, ncol = length(ind))
  for(i in 1:nboot){
    bs[i, ] <- sample(ind, length(ind), replace = TRUE)
  }
  return(bs)
}

#' Reference envelopes
#' 
#' 
#' @encoding UTF-8
#' @name envelope
#' @aliases trajectoryEnvelopeVariability stateEnvelopeVariability compareToStateEnvelope compareToTrajectoryEnvelope
#' 
#' @param d A symmetric \code{\link{matrix}} or an object of class \code{\link{dist}} containing the distance values between pairs of ecosystem states (see details).
#' @param sites A vector indicating the site corresponding to each ecosystem state.
#' @param surveys A vector indicating the survey corresponding to each ecosystem state (only necessary when surveys are not in order).
#' @param nboot.ci Number of bootstrap samples for confidence intervals. If nboot.ci = NULL then confidence intervals are not estimated.
#' @param alpha.ci Error in confidence intervals.
#' @param envelope A vector indicating the set of sites that conform the reference envelope (other sites will be compared to the envelope)
#' @param m Fuzziness exponent for quality value assessment
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
trajectoryEnvelopeVariability<-function(d, sites, surveys = NULL, nboot.ci = NULL, alpha.ci = 0.05, ...){
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  if(is.null(nboot.ci)){
    return(.trajectoryEnvelopeVar(d = d, sites = sites, surveys = surveys, ...))
  } else {
    if(is.integer(nboot.ci)) stop("'nboot.ci' must be an integer")
    if(nboot.ci < 1) stop("'nboot.ci' must be larger than 0")
    var_e <- .trajectoryEnvelopeVar(d = d, sites = sites, surveys = surveys, ...)
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
      bsample <- sample(unique_sites, nsites, replace = TRUE)
      b_sites<-numeric(0)
      b_indices <- numeric(0)
      b_surveys<-numeric(0)
      for(j in 1:length(bsample)) {
        bsj <- bsample[j]
        nt <- sum(sites==bsj)
        b_indices <-c(b_indices, which(sites==bsj))
        b_sites <- c(b_sites, rep(bsj, nt))
        b_surveys <- c(b_surveys, surveys[sites==bsj])
      }
      b_d <- as.dist(as.matrix(d)[b_indices, b_indices])
      var_e_boot[i] <- .trajectoryEnvelopeVar(d = b_d, sites = b_sites, surveys = b_surveys, ...)
    }
    return(c(Var = var_e, quantile(var_e_boot, probs = c((alpha.ci/2.0), 1.0 - (alpha.ci/2.0)))))
  }
}

#' @rdname envelope
stateEnvelopeVariability<-function(d, nboot.ci = NULL, alpha.ci = 0.05){
  if(is.null(nboot.ci)){
    return(.stateEnvelopeVar(d = d))
  } else {
    if(is.integer(nboot.ci)) stop("'nboot.ci' must be an integer")
    if(nboot.ci < 1) stop("'nboot.ci' must be larger than 0")
    var_e <- .stateEnvelopeVar(d = d)
    var_e_boot<-numeric(nboot.ci)
    n <- length(rownames(as.matrix(d)))
    for(i in 1:nboot.ci){
      bsample <- sample(1:n, n, replace = TRUE)
      b_d <- as.dist(as.matrix(d)[bsample, bsample])
      var_e_boot[i] <- .stateEnvelopeVar(d = b_d)
    }
    return(c(Var = var_e, quantile(var_e_boot, probs = c((alpha.ci/2.0), 1.0 - (alpha.ci/2.0)))))
  }
}

#' @rdname envelope
compareToTrajectoryEnvelope<-function(d, sites, envelope, surveys = NULL, m = 1.5, 
                            nboot.ci = NULL, alpha.ci = 0.05,
                            ...) {
  if(length(sites)!=nrow(as.matrix(d))) stop("'sites' needs to be of length equal to the number of rows/columns in d")
  if(!is.null(surveys)) if(length(sites)!=length(surveys)) stop("'sites' and 'surveys' need to be of the same length")
  if(length(envelope)<2) stop("At least two sites must be part of the envelope")
  if(sum(envelope %in% sites)< length(envelope)) stop("Some elements in 'envelope' are not in 'sites'")
  unique_sites <- unique(sites) 
  non_envelope <- unique_sites[!(unique_sites %in% envelope)]
  if(length(non_envelope)==0) stop("At least one site must be outside the envelope")
  
  nsites <- length(unique_sites)
  if(is.null(surveys)) {
    surveys <- numeric(0)
    for(i in 1:nsites) {
      surveys<-c(surveys, 1:sum(sites==sites[i]))
    }
  }   
  D_T <- trajectoryDistances(d = d, sites = sites, surveys = surveys, ...)

  r <- length(envelope)
  sites_T <- colnames(as.matrix(D_T))
  sel_T_env <- sites_T %in% envelope
  D_T_env <- as.dist(as.matrix(D_T)[sel_T_env, sel_T_env])
  var_env <- sum(as.vector(D_T_env)^2)/(r*(r-1))
  
  D2E <- numeric(length(unique_sites))
  is_env <- unique_sites %in% envelope 
  for(i in 1:length(unique_sites)) {
    D_T_i <- as.vector(as.matrix(D_T)[sites_T == unique_sites[i], sel_T_env])
    D2E[i] <- sum(D_T_i^2)/r - 0.5*var_env
  }
  Q <- 1.0/(1.0 + (D2E/var_env)^(1/(m-1))) 
  return(data.frame(Site = unique_sites, Envelope = is_env, Var = var_env, SquaredDist = D2E, Q = Q))
}

#' @rdname envelope
compareToStateEnvelope<-function(d, envelope, m = 1.5, nboot.ci = NULL, alpha.ci = 0.05, ...) {
  if(length(envelope)<2) stop("At least two observations must be part of the envelope")
  obs <- rownames(as.matrix(d))
  if(sum(envelope %in% obs)< length(envelope)) stop("Some elements in 'envelope' are not in row names of 'd'")
  non_envelope <- obs[!(obs %in% envelope)]
  if(length(non_envelope)==0) stop("At least one site must be outside the envelope")
  
  sel_env <- obs %in% envelope
  d_env <- as.dist(as.matrix(d)[sel_env, sel_env])
  var_env <- .stateEnvelopeVar(d_env)

  if(is.null(nboot.ci)){
    D2E <- numeric(length(obs))
    for(i in 1:length(obs)) {
      d_i <- as.vector(as.matrix(d)[obs == obs[i], sel_env])
      D2E[i] <- sum(d_i^2)/length(d_i) - 0.5*var_env
    }
    Q <- 1.0/(1.0 + (D2E/var_env)^(1/(m-1))) 
    return(data.frame(Observation = obs, Envelope = sel_env, SquaredDist = D2E, Q = Q))
  } else {
    if(is.integer(nboot.ci)) stop("'nboot.ci' must be an integer")
    if(nboot.ci < 1) stop("'nboot.ci' must be larger than 0")
    bsm <- .bootstrapSamples(which(sel_env),nboot = nboot.ci)
    D2E <- numeric(length(obs))
    lci_d <- numeric(length(obs))
    uci_d <- numeric(length(obs))
    lci_q <- numeric(length(obs))
    uci_q <- numeric(length(obs))
    var_env_boot <- numeric(nboot.ci)
    D2E_i_boot <- numeric(nboot.ci)
    for(j in 1:nboot.ci) {
      var_env_boot[j] <- .stateEnvelopeVar(d = as.dist(as.matrix(d)[bsm[j,], bsm[j,]]))
    }
    for(i in 1:length(obs)) {
      d_i <- as.vector(as.matrix(d)[i, which(sel_env)])
      D2E[i] <- sum(d_i^2)/length(d_i) - 0.5*var_env
      for(j in 1:nboot.ci) {
        d_i_boot <-  as.vector(as.matrix(d)[i, bsm[j,]])
        D2E_i_boot[j] <- sum(d_i_boot^2)/length(d_i_boot) - 0.5*var_env_boot[j]
      }
      Q_i_boot <- 1.0/(1.0 + (D2E_i_boot/var_env_boot)^(1/(m-1))) 
      Q_i_boot[D2E_i_boot==0] <- 1.0
      lci_d[i] <- quantile(D2E_i_boot, probs = alpha.ci/2.0, na.rm=TRUE)
      uci_d[i] <- quantile(D2E_i_boot, probs = 1.0 - (alpha.ci/2.0), na.rm=TRUE)
      lci_q[i] <- quantile(Q_i_boot, probs = alpha.ci/2.0, na.rm=TRUE)
      uci_q[i] <- quantile(Q_i_boot, probs = 1.0 - (alpha.ci/2.0), na.rm=TRUE)
    }
    Q <- 1.0/(1.0 + (D2E/var_env)^(1/(m-1))) 
    return(data.frame(Observation = obs, Envelope = sel_env, SquaredDist = D2E, 
                      Lower_D = lci_d, Upper_D = uci_d, 
                      Q = Q, Lower_Q = lci_q, Upper_Q = uci_q))
  }
}