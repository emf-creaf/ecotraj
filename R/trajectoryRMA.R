#' Relative Trajectory Movement Assessment (RTMA)
#'
#' Relative Trajectory Movement Assessment (RTMA) is a method allowing testing and qualifying of the relative movements of ecological trajectories (e.g. "convergence", "parallel" etc., see details) as described in Djeghri et al. (in prep). It is implemented in function \code{trajectoryRMA()}.
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param alpha The alpha level for the tests performed in RTMA. Defaults to \code{0.05}.
#' @param nperm Passed to function \code{\link{trajectoryCorrespondence}}. The number of permutations to be used in the dynamic correspondence test. Defaults to \code{999}.
#' @param full.output Flag to indicate that the full output of tests should be computed. Defaults to \code{TRUE}. Setting to FALSE will improve computation speed but yield incomplete outputs (see details).
#' @param add Passed to function \code{\link{trajectoryConvergence}}. Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#' 
#' @details
#' Function \code{trajectoryRMA} attributes a relationship to pairs of ecological trajectories A and B describing their relative movement. It does so by combining four tests:
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
#' To account for multiple testing, \code{trajectoryRMA} performs internally a \enc{Šidák}{Sidak} (1967) correction on the alpha level provided in parameter \code{alpha}.
#' 
#' The results of the four tests (p-values and sign of statistic) are used to assign to each trajectory pair a relationship describing their relative movements. RTMA recognizes a total of 12 relationships, some existing in "weak" variations:
#' \itemize{
#'     \item{\code{"Convergence"} relationship: The two trajectories converge. Exists in a weak version.}
#'     \item{\code{"Divergence"} relationship: The two trajectories diverge. Exists in a weak version.}
#'     \item{\code{"Approaching"} relationship: One trajectory approaches the other.}
#'     \item{\code{"Approaching-Stationary"} relationship: As \code{"Approaching"} but the trajectory approached is stationary relative to the approaching trajectory.}
#'     \item{\code{"Departing"} relationship: One trajectory moves away from the other.}
#'     \item{\code{"Departing-Stationary"} relationship: As \code{"Departing"} but the trajectory departed from is stationary relative to the departing trajectory.}
#'     \item{\code{"Pursuit"} relationship: The two trajectories follow each other.}
#'     \item{\code{"Catchup"} relationship: As \code{"Pursuit"} but the following trajectory moves faster.}
#'     \item{\code{"Escape"} relationship: As \code{"Pursuit"} but the leading trajectory is faster.}
#'     \item{\code{"Parallel"} relationship: The two trajectories travel side by side with broadly similar movements.}
#'     \item{\code{"Antiparallel"} relationship: As \code{"Parallel"} but the two trajectories travel in opposite directions.}
#'     \item{\code{"Neutral"} relationship: The two trajectories have no particular movements relative to each other (effectively the Null Hypothesis for RTMA).}
#'     }
#' Some relationships are asymmetric (e.g. in \code{"Pursuit"} there is a leading and a following trajectory).
#' In these asymmetric relationships the output of function \code{trajectoryRMA} gives the details (see Value section).
#' In rare cases, unlikely relationships (labelled \code{"Other"}) may occur. These involve contradictory patterns hard to interpret.
#' 
#' LIMITATIONS: RTMA has some limitations, in particular it uses trend tests not well suited to study trajectories pairs with changing relative movements (e.g. if two trajectories cross each other, they are first converging then diverging).
#' We advise users to not only rely on RTMA but to also visualize trajectories using function \code{\link{trajectoryPCoA}} for ecological interpretations. See Djeghri et al. (in prep) for more details.
#' Note also that, because RTMA incorporates a correction for multiple testing, it needs somewhat long trajectories to operate (minimum number of survey = 6 at alpha = 0.05).
#' 
#' COMPUTATION TIME: The dynamic correspondence tests performed in RTMA are computationally costly permutation tests only used when all three convergence tests are non-significant.
#' Function \code{trajectoryRMA} performs by default all tests but it is possible to only perform the tests useful for RTMA by setting \code{full.output = FALSE}.
#' This will reduce computation time but the details of the output of RTMA will not contain the information on all possible dynamic correspondence tests, only on relevant ones.
#' 
#' PLOTTING: Function \code{\link{trajectoryConvergencePlot}} provides options to plot the results of RTMA.
#' @returns
#' Function \code{trajectoryRMA} returns an object of class \code{\link{list}} containing:
#' \itemize{
#'     \item{\code{relationship}, a matrix containing the relative movements relationships attributed to each pair of trajectories.}
#'     \item{\code{symmetricConvergence}, a list containing the results of the symmetric convergence test.}
#'     \item{\code{asymmetricConvergence}, a list containing the results of the two asymmetric convergence tests.}
#'     \item{\code{dynamicCorrespondence}, a matrix containing the results of the the dynamic correspondence tests (partial if \code{full.out = FALSE}).}
#'     \item{\code{parameters}, a vector containing the parameters \code{alpha}, the \enc{Šidák}{Sidak} corrected \code{alpha}, and \code{nperm}}.
#'  }
#' In addition to the relationships recognized by RTMA, \code{relationship} provides details on asymmetric relationships (namely \code{"Approaching"}, \code{"Approaching-Stationary"}, \code{"Departing"}, \code{"Departing-Stationary"}, \code{"Pursuit"},  \code{"CatchUp"},  \code{"Escape"}).
#' In asymmetric relationships, the two trajectories have different behavior denoted in \code{relationship} by additional descriptive suffixes pasted on the relationship labels using \code{"_"} as a separator (e.g. \code{"Departing_Departer"}). In the matrix \code{"relationship"}, the suffixes apply to the ROW trajectory.
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
#' @name trajectoryRMA
#' @export
trajectoryRMA <- function(x,
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
  
  #Empty relationship attribution matrix
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
  
  
  #ATTRIBUTING A RELATIONSHIP
  for (j in trajs[1:(length(trajs)-1)]){
    for (k in trajs[(which(trajs==j)+1):length(trajs)]){
      
      if (sym$p.value[j,k]<=alpha){
        #this is the branch with significant symmetric convergence or divergence
        if (sym$tau[j,k]>0){
          #this is the branch with a divergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Divergence relationship (may be rare)
            output$relationship[k,j] <- "Weak Divergence"
            output$relationship[j,k] <- "Weak Divergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]>0){
                #this is the Divergence relationship, all tests agree on divergence
                output$relationship[k,j] <- "Divergence"
                output$relationship[j,k] <- "Divergence"
              }else{
                #this a (very) unlikely relationship where both asymmetric test say convergence when the symmetric test says divergence
                output$relationship[k,j] <- "Other (opposed sym and asym)"
                output$relationship[j,k] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the Escape relationship with symmetric divergence, and disagreeing asymmetric tests
              if (asym$tau[k,j]<0){
                output$relationship[k,j] <- "Escape_Follower"
                output$relationship[j,k] <- "Escape_Leader"
              }else{
                output$relationship[k,j] <- "Escape_Leader"
                output$relationship[j,k] <- "Escape_Follower"
              }
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
              #this is the Departing-Stationary relationship where there is one divergent asymmetric test and symmetric divergence
              if (asym$p.value[k,j]<=alpha){
                output$relationship[k,j] <- "Departing-Stationary_Departer"
                output$relationship[j,k] <- "Departing-Stationary_Departee"
              }else{
                output$relationship[k,j] <- "Departing-Stationary_Departee"
                output$relationship[j,k] <- "Departing-Stationary_Departer"
              }
            }else{
              #this is an unlikely relationship with one asymmetric test significant and opposed to the symmetric test
              output$relationship[k,j] <- "Other (opposed sym and one asym)"
              output$relationship[j,k] <- "Other (opposed sym and one asym)"
            }
          }
          
          
          
        }else{
          #this is the branch with a convergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Convergence relationship (may be rare)
            output$relationship[k,j] <- "Weak Convergence"
            output$relationship[j,k] <- "Weak Convergence"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]<0){
                #this is the Convergence relationship, all tests agree on convergence
                output$relationship[k,j] <- "Convergence"
                output$relationship[j,k] <- "Convergence"
              }else{
                #this a (very) unlikely relationship where both asymmetric test say divergence when the symmetric test says convergence
                output$relationship[k,j] <- "Other (opposed sym and asym)"
                output$relationship[j,k] <- "Other (opposed sym and asym)"
              }
            }else{
              #this is the CatchUp relationship with symmetric convergence, and disagreeing asymmetric tests
              output$relationship[k,j] <- "CatchUp"
              if (asym$tau[k,j]<0){
                output$relationship[k,j] <- "CatchUp_Follower"
                output$relationship[j,k] <- "CatchUp_Leader"
              }else{
                output$relationship[k,j] <- "CatchUp_Leader"
                output$relationship[j,k] <- "CatchUp_Follower"
              }
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]<0){
              #this is the Approaching-Stationary relationship where there is one convergent asymmetric test and symmetric convergence
              if (asym$p.value[k,j]<=alpha){
                output$relationship[k,j] <- "Approaching-Stationary_Approacher"
                output$relationship[j,k] <- "Approaching-Stationary_Approachee"
              }else{
                output$relationship[k,j] <- "Approaching-Stationary_Approachee"
                output$relationship[j,k] <- "Approaching-Stationary_Approacher"
              }
            }else{
              #this is an unlikely relationship with one asymmetric test significant and opposed to the symmetric test
              output$relationship[k,j] <- "Other (opposed sym and one asym)"
              output$relationship[j,k] <- "Other (opposed sym and one asym)"
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
              #this is the Parallel relationship
              output$relationship[k,j] <- "Parallel"
              output$relationship[j,k] <- "Parallel"
            }else{
              #this is the Antiparallel relationship
              output$relationship[k,j] <- "Antiparallel"
              output$relationship[j,k] <- "Antiparallel"
            }
          }else{
            #this is the Neutral relationship (nothing significant)
            output$relationship[k,j] <- "Neutral"
            output$relationship[j,k] <- "Neutral"
          }
        }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
          #this is the branch where the asymmetric tests both are significant
          if (sign(asym$tau[j,k])==(sign(asym$tau[k,j]*(-1)))){
            #this is the Pursuit relationship (no significant symmetric test and opposed asymmetric tests)
            if (asym$tau[k,j]<0){
              output$relationship[k,j] <- "Pursuit_Follower"
              output$relationship[j,k] <- "Pursuit_Leader"
            }else{
              output$relationship[k,j] <- "Pursuit_Leader"
              output$relationship[j,k] <- "Pursuit_Follower"
            }
          }else{
            #this is an unlikely branch where the two asymmetric tests have agreement but the symmetric test is non-significant 
            if (asym$tau[j,k]<0){
              #this is the second version of the weak convergence relationship
              output$relationship[k,j] <- "Weak Convergence"
              output$relationship[j,k] <- "Weak Convergence"
            }else{
              #this is the second version of the weak divergence relationship
              output$relationship[k,j] <- "Weak Divergence"
              output$relationship[j,k] <- "Weak Divergence"
            }
          }
        }else{
          #this is the branch where only one asymmetric test is significant
          if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
            #this is the Departing relationship where there is only one divergent asymmetric test significant
            if (asym$p.value[k,j]<=alpha){
              output$relationship[k,j] <- "Departing_Departer"
              output$relationship[j,k] <- "Departing_Departee"
            }else{
              output$relationship[k,j] <- "Departing_Departee"
              output$relationship[j,k] <- "Departing_Departer"
            }
          }else{
            #this is the Approaching relationship where there is only one convergent asymmetric test significant
            if (asym$p.value[k,j]<=alpha){
              output$relationship[k,j] <- "Approaching_Approacher"
              output$relationship[j,k] <- "Approaching_Approachee"
            }else{
              output$relationship[k,j] <- "Approaching_Approachee"
              output$relationship[j,k] <- "Approaching_Approacher"
            }
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
