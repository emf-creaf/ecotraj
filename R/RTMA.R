#' Relative Trajectory Movement Assessment (RTMA)
#'
#' Relative Trajectory Movement Assessment (RTMA) is a method allowing testing and qualifying of the relative movements of ecological trajectories (e.g. "convergence", "parallel" etc., see details) as described in Djeghri et al. (in prep).
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param alpha The alpha level for the tests performed in RTMA. Defaults to \code{0.05}.
#' @param nperm Passed to function \code{\link{trajectoryCorrespondence}}. The number of permutations to be used in the dynamic correspondence test. Defaults to \code{999}.
#' @param full.output Flag to indicate that the full output of tests should be computed. Defaults to \code{TRUE}. Setting to FALSE will improve computation speed but yield incomplete outputs (see details).
#' @param add Passed to function \code{\link{trajectoryConvergence}}. Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#' 
#' @details
#' Function \code{RTMA} attributes a "scenario" of relative movement to pairs of ecological trajectories A and B by combining four tests:
#' \itemize{
#'     \item{Three convergence tests performed through internal callings of function \code{\link{trajectoryConvergence}}:
#'        \itemize{
#'             \item{The symmetric convergence test between trajectories A and B.}
#'             \item{The asymmetric convergence test assessing if trajectory A approaches trajectory B.}
#'             \item{The asymmetric convergence test assessing if trajectory B approaches trajectory A.}
#'        }
#'     }
#'     \item{One dynamic correspondence test performed through internal callings to function \code{\link{trajectoryCorrespondence}}.}
#'  } 
#' To account for multiple testing, \code{RTMA} performs internally a \enc{Šidák}{Sidak} (1967) correction on the alpha level provided in parameter \code{alpha}.
#' 
#' The results of the four tests (p-values and sign of statistic) are used to assign to each trajectory pair a "scenario" describing their relative movements. RTMA recognizes a total of 10 scenarios, some existing in "weak" variations:
#' \itemize{
#'     \item{\code{Convergence} scenario: The two trajectories converge, exists in a weak version.}
#'     \item{\code{Divergence} scenario: The two trajectories diverge, exists in a weak version.}
#'     \item{\code{Approaching} scenario: One trajectory approaches the other (likely relatively immobile), exists in a weak version, potentially indicative of an approach by the side.}
#'     \item{\code{Departing} scenario: One trajectory moves away from the other (likely relatively immobile), exists in a weak version, potentially indicative of a departure from the side.}
#'     \item{\code{Pursuit} scenario: The two trajectories follow each other.}
#'     \item{\code{Catchup} scenario: As \code{Pursuit} but the following trajectory moves faster.}
#'     \item{\code{Escape} scenario: As \code{Pursuit} but the leading trajectory is faster.}
#'     \item{\code{Parallel} scenario: The two trajectories travel side by side with broadly similar movements.}
#'     \item{\code{Antiparallel} scenario: As \code{Parallel} but the two trajectories travel in opposite directions.}
#'     \item{\code{Neutral} scenario:  The two trajectories have no particular movements relative to each other (effectively the Null Hypothesis).}
#'     }
#' In rare cases, unlikely scenarios (labelled \code{Other}) may occur. These involve contradictory patterns hard to interpret.
#' 
#' LIMITATIONS: RTMA has some limitations, in particular it uses trend tests not well suited to study trajectories pairs with changing relative movements (e.g. if two trajectories cross each other, they are first converging then diverging).
#' We advise users to not only rely on RTMA but to also visualize trajectories using function \code{\link{trajectoryPCoA}} for ecological interpretations. See Djeghri et al. (in prep) for more details.
#' Note also that, because RTMA incorporates a correction for multiple testing, it necessitates somewhat long trajectories to operate (minimum number of survey = 6 at alpha = 0.05).
#' 
#' COMPUTATION TIME: The dynamic correspondence tests performed in RTMA are computationally costly permutation tests only used when all three convergence tests are non-significant.
#' Function \code{RTMA} performs by default all tests but it is possible to only perform the tests useful for RTMA by setting \code{full.output = FALSE}.
#' This will reduce computation time but the details of the output of RTMA will not contain the information on all possible dynamic correspondence tests, only on relevant ones.
#' 
#' PLOTTING: Function \code{\link{trajectoryConvergencePlot}} provides options to plot the results of RTMA.
#' @returns
#' Function \code{RTMA} returns an object of class \code{\link{list}} containing:
#' \itemize{
#'     \item{\code{relationship}, a matrix containing the relative movements scenario attributed to each pair of trajectories.}
#'     \item{\code{symmetricConvergence}, a list containing the results of the symmetric convergence test.}
#'     \item{\code{asymmetricConvergence}, a list containing the results of the two asymmetric convergence tests.}
#'     \item{\code{dynamicCorrespondence}, a matrix containing the results of the the dynamic correspondence tests (partial if \code{full.out = FALSE}).}
#'     \item{\code{parameters}, a vector containing the parameters \code{alpha}, the \enc{Šidák}{Sidak} corrected \code{alpha}, and \code{nperm}}.
#'  }
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @references
#' Djeghri et al. (in preparation) Uncovering the relative movements of ecological trajectories.
#' 
#' \enc{Šidák}{Sidak}, Z. (1967) Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association 62:648-633.
#'
#' @seealso \code{\link{trajectoryConvergence}}, \code{\link{trajectoryCorrespondence}}, \code{\link{trajectoryConvergencePlot}} 
#'
#' @examples
#' #Obtain and format some trajectories
#' data("avoca")
#' avoca_D_man <- vegclust::vegdiststruct(avoca_strat, 
#'                                        method ="manhattan", 
#'                                        transform = function(x){log(x+1)})
#' years <- c(1971, 1974, 1978, 1983, 1987, 1993, 1999, 2004, 2009)
#' avoca_times <- years[avoca_surveys]
#' avoca_x <- defineTrajectories(d = avoca_D_man,  
#'                               sites = avoca_sites, 
#'                               times = avoca_times)
#'                               
#' #Visualize the trajectories
#' trajectoryPCoA(avoca_x,traj.colors = RColorBrewer::brewer.pal(8,"Accent"),length=0.1,lwd=2)
#' legend("bottomleft",bty="n",legend=1:8,col=RColorBrewer::brewer.pal(8,"Accent"),lwd=2,ncol=2)
#' 
#' #Perform RTMA
#' RTMA(avoca_x)
#' 
#' 
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
  #stop in case trajectories are not of same size
  if (any(table(sites)-mean(table(sites))!=0)) stop("RTMA only applies to trajectories with the same number of surveys")
  
  #Sidak correction procedure for alpha:
  alphaUncor <- alpha
  alpha <- 1-(1-alpha)^(1/4)
  
  #stop if trajectories are too short to yield significant results in the convergence tests given the chosen alpha
  if (as.numeric(Kendall::MannKendall(1:table(sites)[1])$sl)>alpha) stop("Trajectories are not long enough to yield significant tests given the alpha level")
  
  #check if trajectories are long enough to apply RTMA at this alpha level:
  
  trajs <- unique(sites)
  
  #preparing the output
  output <- list()
  
  #Empty scenario attribution matrix
  output$relationship <- matrix(NA,length(trajs),length(trajs))
  colnames(output$relationship) <- trajs
  rownames(output$relationship) <- trajs
  
  #convergence tests
  sym <- trajectoryConvergence(x,type="pairwise.symmetric",add = add)
  asym <- trajectoryConvergence(x,type="pairwise.asymmetric",add = add)
  
  #dynamic correspondence tests
  if (full.output == TRUE){
    Dcor <- trajectoryCorrespondence(x,nperm)
  }else if (full.output == FALSE){
    Dcor <- output$relationship
    x$metadata$sites <- sites#this is a little hack to handle correctly different class of trajectories
  }else{
    stop("full.output must be a logical flag.")
  }
  
  
  #ATTRIBUTING A SCENARIO
  for (j in trajs[1:(length(trajs)-1)]){
    for (k in trajs[(which(trajs==j)+1):length(trajs)]){
      
      if (sym$p.value[j,k]<=alpha){
        #this is the branch with significant symmetric convergence or divergence
        if (sym$tau[j,k]>0){
          #this is the branch with a divergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Divergence scenario (may be rare)
            output$relationship[k,j] <- "Weak Divergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]>0){
                #this is the Divergence scenario, all tests agree on divergence
                output$relationship[k,j] <- "Divergence"
              }else{
                #this a (very) unlikely scenario where both asymmetric test say convergence when the symmetric test says divergence
                output$relationship[k,j] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the Escape scenario with symmetric divergence, and disagreeing asymmetric tests
              output$relationship[k,j] <- "Escape"
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
              #this is the Departing scenario where there is one divergent asymmetric test and symmetric divergence
              output$relationship[k,j] <- "Departing"
            }else{
              #this is an unlikely scenario with one asymmetric test significant and opposed to the symmetric test
              output$relationship[k,j] <- "Other (opposed sym and one asym)"
            }
          }
          
          
          
        }else{
          #this is the branch with a convergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Convergence scenario (may be rare)
            output$relationship[k,j] <- "Weak Convergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]<0){
                #this is the Convergence scenario, all tests agree on convergence
                output$relationship[k,j] <- "Convergence"
              }else{
                #this a (very) unlikely scenario where both asymmetric test say divergence when the symmetric test says convergence
                output$relationship[k,j] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the CatchUp scenario with symmetric convergence, and disagreeing asymmetric tests
              output$relationship[k,j] <- "CatchUp"
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]<0){
              #this is the Approaching scenario where there is one convergent asymmetric test and symmetric convergence
              output$relationship[k,j] <- "Approaching"
            }else{
              #this is an unlikely scenario with one asymmetric test significant and opposed to the symmetric test
              output$relationship[k,j] <- "Other (opposed sym and one asym)"
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
              output$relationship[k,j] <- "Parallel"
            }else{
              #this is the Antiparallel scenario
              output$relationship[k,j] <- "Antiparallel"
            }
          }else{
            #this is the Neutral scenario (nothing significant)
            output$relationship[k,j] <- "Neutral"
          }
        }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
          #this is the branch where the asymmetric tests both are significant
          if (sign(asym$tau[j,k])==(sign(asym$tau[k,j]*(-1)))){
            #this is the Pursuit scenario (no significant symmetric test and opposed asymmetric tests)
            output$relationship[k,j] <- "Pursuit"
          }else{
            #this is an unlikely branch where the two asymmetric tests have agreement but the symmetric test is non-significant 
            if (asym$tau[j,k]<0){
              #this is the second version of the weak convergence scenario
              output$relationship[k,j] <- "Weak Convergence"
            }else{
              #this is the second version of the weak divergence scenario
              output$relationship[k,j] <- "Weak Divergence"
            }
          }
        }else{
          #this is the branch where only one asymmetric test is significant
          if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
            #this is the Weak Departing scenario where there is only one divergent asymmetric test significant
            output$relationship[k,j] <- "Weak/Side Departing"
          }else{
            #this is the Weak Approaching scenario where there is only one convergent asymmetric test significant
            output$relationship[k,j] <- "Weak/Side Approaching"
          }
        }
      }
    }
  }
  #finish building the output
  output$symmetricConvergence <- sym
  output$asymmetricConvergence <- asym
  output$dynamicCorrespondence <- Dcor
  output$parameters <- c(alphaUncor,alpha,nperm)
  names(output$parameters) <- c("alpha","alpha corrected","nperm")
  #define its class
  class(output) <- c("RTMA","list")
  
  return(output)
}
