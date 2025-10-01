#' Relative Trajectory Movement Assessment (RTMA)
#'
#' RTMA
#' @param x An object of class \code{\link{trajectories}}.
#' @param nperm The number of permutations to be used in the dynamic correspondence test. Defaults to \code{999}.
#' @param full.output Flag to indicate that the full output of tests should be computed. Defaults to \code{TRUE}. Setting to FALSE will improve computation speed but yield incomplete outputs (see details).
#' @param add Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.

#' @name RTMA
#' @export
RTMA <- function(x,
                 alpha = 0.05,
                 nperm = 999,
                 full.output = TRUE,
                 add = TRUE
                 ){
  
  if(!inherits(x,"trajectories"))
    stop ("'x' should be of class 'trajectory'")
  if(inherits(x, "fd.trajectories")) {
    sites <- x$metadata$fdT
  } else if(inherits(x, "cycles")) {
    sites <- x$metadata$cycles
  } else if(inherits(x, "sections")) {
    sites <- x$metadata$sections
  } else {
    sites <- x$metadata$sites
  }
  trajs <- unique(sites)
  
  #preparing the output
  output <- list()
  
  #Empty scenario attribution matrix
  output$ScenarioAttribution <- matrix(NA,length(trajs),length(trajs))
  colnames(output$ScenarioAttribution) <- trajs
  rownames(output$ScenarioAttribution) <- trajs
  
  #convergence tests
  sym <- trajectoryConvergence(x,type="pairwise.symmetric",add = add)
  asym <- trajectoryConvergence(x,type="pairwise.asymmetric",add = add)
  
  #dynamic correspondence tests
  if (full.output == TRUE){
    Dcor <- trajectoryCorrespondence(x,nperm)
  }else if (full.output == FALSE){
    Dcor <- output$ScenarioAttribution
    x$metadata$sites <- sites#this is a little hack to handle correctly different class of trajectories
  }else{
    stop("full.output must be a logical flag.")
  }
  
  #Sidak correction procedure for alpha:
  alpha <- 1-(1-alpha)^(1/4)
  
  
  #ATTRIBUTING A SCENARIO
  for (j in trajs[1:(length(trajs)-1)]){
    for (k in trajs[(which(trajs==j)+1):length(trajs)]){
      
      if (sym$p.value[j,k]<=alpha){
        #this is the branch with significant symmetric convergence or divergence
        if (sym$tau[j,k]>0){
          #this is the branch with a divergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Divergence scenario (may be rare)
            output$ScenarioAttribution[k,j] <- "Weak Divergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]>0){
                #this is the Divergence scenario, all tests agree on divergence
                output$ScenarioAttribution[k,j] <- "Divergence"
              }else{
                #this a (very) unlikely scenario where both asymmetric test say convergence when the symmetric test says divergence
                output$ScenarioAttribution[k,j] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the Escape scenario with symmetric divergence, and disagreeing asymmetric tests
              output$ScenarioAttribution[k,j] <- "Escape"
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
              #this is the Departing scenario where there is one divergent asymmetric test and symmetric divergence
              output$ScenarioAttribution[k,j] <- "Departing"
            }else{
              #this is an unlikely scenario with one asymmetric test significant and opposed to the symmetric test
              output$ScenarioAttribution[k,j] <- "Other (opposed sym and one asym)"
            }
          }
          
          
          
        }else{
          #this is the branch with a convergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Convergence scenario (may be rare)
            output$ScenarioAttribution[k,j] <- "Weak Convergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]<0){
                #this is the Convergence scenario, all tests agree on convergence
                output$ScenarioAttribution[k,j] <- "Convergence"
              }else{
                #this a (very) unlikely scenario where both asymmetric test say divergence when the symmetric test says convergence
                output$ScenarioAttribution[k,j] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the CatchUp scenario with symmetric convergence, and disagreeing asymmetric tests
              output$ScenarioAttribution[k,j] <- "CatchUp"
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]<0){
              #this is the Approaching scenario where there is one convergent asymmetric test and symmetric convergence
              output$ScenarioAttribution[k,j] <- "Approaching"
            }else{
              #this is an unlikely scenario with one asymmetric test significant and opposed to the symmetric test
              output$ScenarioAttribution[k,j] <- "Other (opposed sym and one asym)"
            }
          }
        }
        
      }else{
        #this is the branch with non-significant symmetric test
        if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
          #this is the branch with no significant convergence tests
          
          #the first if checks if it is needed to do the dynamic correspondence test
          if (full.output==TRUE){
            Dcorbis <- Dcor[c(j,k),c(j,k)]
          }else{
            xbis <- subsetTrajectories(x,site_selection = c(j,k))
            Dcorbis <- trajectoryCorrespondence(xbis,nperm=nperm)
            Dcor[k,j] <- Dcorbis[2,1]
            Dcor[j,k] <- Dcorbis[1,2]
          }
          
          if (Dcorbis[2,1]<=alpha){
            #this is the branch where we have Parallelism
            if (Dcorbis[1,2]>0){
              #this is the Parallel scenario
              output$ScenarioAttribution[k,j] <- "Parallel"
            }else{
              #this is the Antiparallel scenario
              output$ScenarioAttribution[k,j] <- "Antiparallel"
            }
          }else{
            #this is the Neutral scenario (nothing significant)
            output$ScenarioAttribution[k,j] <- "Neutral"
          }
        }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
          #this is the branch where the asymmetric tests both are significant
          if (sign(asym$tau[j,k])==(sign(asym$tau[k,j]*(-1)))){
            #this is the Pursuit scenario (no significant symmetric test and opposed asymmetric tests)
            output$ScenarioAttribution[k,j] <- "Pursuit"
          }else{
            #this is an unlikely branch where the two asymmetric tests have agreement but the symmetric test is non-significant 
            if (asym$tau[j,k]<0){
              #this is the second version of the weak convergence scenario
              output$ScenarioAttribution[k,j] <- "Weak Convergence"
            }else{
              #this is the second version of the weak divergence scenario
              output$ScenarioAttribution[k,j] <- "Weak Divergence"
            }
          }
        }else{
          #this is the branch where only one asymmetric test is significant
          if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
            #this is the Weak Departing scenario where there is only one divergent asymmetric test significant
            output$ScenarioAttribution[k,j] <- "Weak Departing"
          }else{
            #this is the Weak Approaching scenario where there is only one convergent asymmetric test significant
            output$ScenarioAttribution[k,j] <- "Weak Approaching"
          }
        }
      }
    }
  }
  #finish building the output
  output$DetailedTests$SymmetricConvergence <- sym
  output$DetailedTests$AsymmetricConvergence <- asym
  output$DetailedTests$DynamicCorrespondence <- Dcor
  
  return(output)
}
