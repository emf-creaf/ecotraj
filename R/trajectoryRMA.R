#' Relative Trajectory Movement Assessment (RTMA)
#'
#' Relative Trajectory Movement Assessment (RTMA) is a method for testing and qualifying of the relative movements of ecological trajectories (e.g. "convergence", "parallel" etc., see details) as described in Djeghri et al. (in prep). It is implemented in function \code{trajectoryRMA()}.
#' 
#' @param x An object of class \code{\link{trajectories}}.
#' @param alpha The alpha level for the tests performed in RTMA. Defaults to \code{0.05}.
#' @param nperm Passed to function \code{\link{trajectoryCorrespondence}}. The number of permutations to be used in the dynamic correspondence test. Defaults to \code{999}.
#' @param full.output Flag to indicate that the full output of tests should be computed. Defaults to \code{TRUE}. Setting to \code{FALSE} will improve computation speed but yield incomplete outputs (see details).
#' @param add Passed to function \code{\link{trajectoryConvergence}}. Flag to indicate that constant values should be added (local transformation) to correct triplets of distance values that do not fulfill the triangle inequality.
#'
#' @encoding UTF-8
#' @aliases RTMA trajectoryRMA
#' 
#' 
#' @details
#' Function \code{trajectoryRMA} attributes a dynamic relationship to pairs of ecological trajectories A and B describing their relative movement. It does so by combining four tests:
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
#' The results of the four tests (p-values and sign of statistic) are used to assign to each trajectory pair a relationship describing their relative movements. RTMA recognizes a total of 12 relationships, some existing in "weak" variations. 
#' The following five dynamic relationships are \emph{symmetric}, i.e. applying to the two trajectories without distinction of roles:
#' \itemize{
#'     \item{\code{"convergence"} - The two trajectories converge. Exists in a weak version.}
#'     \item{\code{"divergence"} - The two trajectories diverge. Exists in a weak version.}
#'     \item{\code{"parallel"} - The two trajectories travel side by side with broadly similar movements.}
#'     \item{\code{"antiparallel"} - As in \code{"parallel"} but the two trajectories travel in opposite directions.}
#'     \item{\code{"neutral"} - The two trajectories have no particular movements relative to each other (effectively the null hypothesis for RTMA).}
#' }
#' The following seven dynamic relationships are \emph{asymmetric} (e.g. in \code{"pursuit"} there is a leading and a following trajectory). In these asymmetric relationships the output of function \code{trajectoryRMA} gives the role of each trajectory (see Value section). A more general interpretation of asymmetry is to consider that the relationship is \emph{oriented} (see below, relationship groups).
#' \itemize{
#'     \item{\code{"approaching"} - One trajectory approaches the other.}
#'     \item{\code{"approaching-stationary"} - As in \code{"approaching"} but the trajectory approached is stationary relative to the approaching trajectory.}
#'     \item{\code{"departing"} - One trajectory moves away from the other.}
#'     \item{\code{"departing-stationary"} - As in \code{"departing"} but the trajectory departed from is stationary relative to the departing trajectory.}
#'     \item{\code{"pursuit"} - The two trajectories follow each other.}
#'     \item{\code{"catch-up"} - As in \code{"pursuit"} but the following trajectory moves faster.}
#'     \item{\code{"escape"} - As in \code{"pursuit"} but the leading trajectory is faster.}
#' }
#' 
#' In rare cases, unlikely relationships (labelled \code{"unknown"}, with a short description in brackets) may occur. These involve contradictory patterns hard to interpret.
#' 
#' RELATIONSHIP GROUPS: It is possible to further sort trajectory relationships in broad \emph{relationship groups} (not always mutually exclusive). Three such groups are recognized in RTMA:
#' \itemize{
#'    \item{The \code{"convergence group"}, includes relationships that display convergence in the broadest sense with a trend of diminishing distance between the two trajectories. Formally this group includes relationships of
#'    \code{"convergence"} and its weak version, \code{"approaching"},\code{"approaching-stationary"} and \code{"catch-up"}.}
#'    \item{The \code{"divergence group"}, includes relationships that display divergence in the broadest sense with a trend of increasing distance between the two trajectories. Formally this group includes relationships of
#'    \code{"divergence"} and its weak version, \code{"departing"},\code{"departing-stationary"} and \code{"escape"}.}
#'    \item{The \code{"oriented group"}, includes relationships that have, broadly speaking, a trajectory \emph{in front} and a trajectory \emph{in the back} implying an orientation to their relationship. This group includes all asymmetric relationships, formally:
#'    \code{"approaching"},\code{"approaching-stationary"}, \code{"departing"},\code{"departing-stationary"}, \code{"catch-up"}, \code{"escape"} and \code{"pursuit"}.}
#' }
#' Note that a given relationship may belong to two groups (either convergence or divergence group + oriented group) and that \code{"parallel"},\code{"antiparallel"} and \code{"neutral"} relationships stand on their own, not belonging to any groups.
#' In our experience, relationship groups have proven a useful conceptual tool to reveal large scale patterns particularly when adressing many trajectory relationships (see Djeghri et al. in prep).
#' 
#' LIMITATIONS: RTMA has some limitations, in particular it uses trend tests not well suited to study trajectories pairs with changing relative movements (e.g. if two trajectories cross each other, they are first converging then diverging).
#' We advise users to not only rely on RTMA but to also visualize trajectories using function \code{\link{trajectoryPCoA}} for ecological interpretations. See Djeghri et al. (in prep) for more details.
#' Note also that, because RTMA incorporates a correction for multiple testing, it needs somewhat long trajectories to operate (minimum number of survey = 6 at alpha = 0.05).
#' 
#' COMPUTATION TIME: The dynamic correspondence tests performed in RTMA are computationally costly permutation tests only used when all three convergence tests are non-significant.
#' Function \code{trajectoryRMA} performs by default all tests but it is possible to only perform the tests useful for RTMA by setting \code{full.output = FALSE}.
#' This will reduce computation time but the details of the output of RTMA will not contain the information on all possible dynamic correspondence tests, only on relevant ones.
#' 
#' PLOTTING: Functions \code{\link{trajectoryConvergencePlot}} and \code{\link{trajectoryRMAPlot}} provide options to plot the results of RTMA.
#' @returns
#' Function \code{trajectoryRMA} returns an object of classes \code{\link{list}} and \code{RTMA} containing:
#' \itemize{
#'     \item{\code{dynamic_relationships_taxonomy}: a data-frame containing the names of the relative movement relationships recognized by RTMA as well as corresponding relationship groups. This part of \code{trajectoryRMA} output is independent of the trajectories used as input and is primarily a bestiary (see details). It can be used to transform the \code{dynamic_relationships} matrix (see below) to focus on chosen relationship groups. }
#'     \item{\code{dynamic_relationships}: a matrix containing the relative movement relationships attributed to each pair of trajectories.}
#'     \item{\code{symmetric_convergence}: a list containing the results of the symmetric convergence test.}
#'     \item{\code{asymmetric_convergence}: a list containing the results of the two asymmetric convergence tests.}
#'     \item{\code{correspondence}: a matrix containing the results of the the dynamic correspondence tests (partial if \code{full.out = FALSE}).}
#'     \item{\code{parameters}: a vector containing the parameters \code{alpha}, the \enc{Šidák}{Sidak} \code{corrected_alpha}, and \code{nperm}}.
#'  }
#' In addition to the relationships recognized by RTMA, matrix \code{dynamic_relationships} provides the role of each trajectory in asymmetric relationships. 
#' The role is provided in parenthesis and applies to the trajectory of the ROW index. For example, \code{"approaching (approacher)"} means 
#' that the trajectory of the corresponding row is approaching the trajectory of the corresponding column, which will have \code{"approaching (target)"}.
#' In symmetric relationships, the wording \code{(symmetric)} is added to indicate that there is no distinction of roles.
#'
#' 
#' @author Nicolas Djeghri, UBO
#' @author Miquel De \enc{Cáceres}{Caceres}, CREAF
#' 
#' @references
#' Djeghri et al. (in preparation) Uncovering the relative movements of ecological trajectories.
#' 
#' \enc{Šidák}{Sidak}, Z. (1967) Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association 62:648-633.
#'
#' @seealso \code{\link{trajectoryConvergence}}, \code{\link{trajectoryCorrespondence}}, \code{\link{trajectoryConvergencePlot}}, \code{\link{trajectoryRMAPlot}} 
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
#' trajectoryRMA(avoca_x)
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

  #Empty relationship attribution matrix
  relationships <- matrix(NA,length(trajs),length(trajs))
  colnames(relationships) <- trajs
  rownames(relationships) <- trajs
  
  #convergence tests
  sym <- trajectoryConvergence(x,type="pairwise.symmetric",add = add)
  asym <- trajectoryConvergence(x,type="pairwise.asymmetric",add = add)
  
  #dynamic correspondence tests
  if (full.output == TRUE){
    Dcor <- trajectoryCorrespondence(x,nperm)
  } else if (full.output == FALSE){
    Dcor <- matrix(NA,length(trajs),length(trajs))
    colnames(Dcor) <- trajs
    rownames(Dcor) <- trajs
    x$metadata$sites <- sites #this is a little hack to handle correctly different class of trajectories
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
            relationships[k,j] <- "weak divergence (symmetric)"
            relationships[j,k] <- "weak divergence (symmetric)"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]>0){
                #this is the Divergence relationship, all tests agree on divergence
                relationships[k,j] <- "divergence (symmetric)"
                relationships[j,k] <- "divergence (symmetric)"
              }else{
                #this a (very) unlikely relationship where both asymmetric test say convergence when the symmetric test says divergence
                relationships[k,j] <- "unknown (opposed sym and asym)"
                relationships[j,k] <- "unknown (opposed sym and asym)"
              }
            }else{
              #this is the Escape relationship with symmetric divergence, and disagreeing asymmetric tests
              if (asym$tau[k,j]<0){
                relationships[k,j] <- "escape (follower)"
                relationships[j,k] <- "escape (leader)"
              }else{
                relationships[k,j] <- "escape (leader)"
                relationships[j,k] <- "escape (follower)"
              }
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
              #this is the Departing-Stationary relationship where there is one divergent asymmetric test and symmetric divergence
              if (asym$p.value[k,j]<=alpha){
                relationships[k,j] <- "departing-stationary (departer)"
                relationships[j,k] <- "departing-stationary (stationary origin)"
              }else{
                relationships[k,j] <- "departing-stationary (stationary origin)"
                relationships[j,k] <- "departing-stationary (departer)"
              }
            }else{
              #this is an unlikely relationship with one asymmetric test significant and opposed to the symmetric test
              relationships[k,j] <- "unknown (opposed sym and one asym)"
              relationships[j,k] <- "unknown (opposed sym and one asym)"
            }
          }
          
          
          
        }else{
          #this is the branch with a convergent symmetric test
          
          if ((asym$p.value[j,k]>alpha)&(asym$p.value[k,j]>alpha)){
            #this is the first version of the Weak Convergence relationship (may be rare)
            relationships[k,j] <- "weak convergence (symmetric)"
            relationships[j,k] <- "weak convergence (symmetric)"
          }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
            #this is the branch where the asymmetric tests both are significant
            if (sign(asym$tau[j,k])==sign(asym$tau[k,j])){
              if (asym$tau[j,k]<0){
                #this is the Convergence relationship, all tests agree on convergence
                relationships[k,j] <- "convergence (symmetric)"
                relationships[j,k] <- "convergence (symmetric)"
              }else{
                #this a (very) unlikely relationship where both asymmetric test say divergence when the symmetric test says convergence
                relationships[k,j] <- "unknown (opposed sym and asym)"
                relationships[j,k] <- "unknown (opposed sym and asym)"
              }
            } else {
              #this is the Catch-up relationship with symmetric convergence, and disagreeing asymmetric tests
              if (asym$tau[k,j]<0){
                relationships[k,j] <- "catch-up (follower)"
                relationships[j,k] <- "catch-up (leader)"
              }else{
                relationships[k,j] <- "catch-up (leader)"
                relationships[j,k] <- "catch-up (follower)"
              }
            }
          }else{
            #this is the branch where only one asymmetric test is significant
            if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]<0){
              #this is the Approaching-Stationary relationship where there is one convergent asymmetric test and symmetric convergence
              if (asym$p.value[k,j]<=alpha){
                relationships[k,j] <- "approaching-stationary (approacher)"
                relationships[j,k] <- "approaching-stationary (stationary target)"
              }else{
                relationships[k,j] <- "approaching-stationary (stationary target)"
                relationships[j,k] <- "approaching-stationary (approacher)"
              }
            }else{
              #this is an unlikely relationship with one asymmetric test significant and opposed to the symmetric test
              relationships[k,j] <- "unknown (opposed sym and one asym)"
              relationships[j,k] <- "unknown (opposed sym and one asym)"
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
              relationships[k,j] <- "parallel (symmetric)"
              relationships[j,k] <- "parallel (symmetric)"
            }else{
              #this is the Antiparallel relationship
              relationships[k,j] <- "antiparallel (symmetric)"
              relationships[j,k] <- "antiparallel (symmetric)"
            }
          }else{
            #this is the Neutral relationship (nothing significant)
            relationships[k,j] <- "neutral (symmetric)"
            relationships[j,k] <- "neutral (symmetric)"
          }
        }else if ((asym$p.value[j,k]<=alpha)&(asym$p.value[k,j]<=alpha)){
          #this is the branch where the asymmetric tests both are significant
          if (sign(asym$tau[j,k])==(sign(asym$tau[k,j]*(-1)))){
            #this is the Pursuit relationship (no significant symmetric test and opposed asymmetric tests)
            if (asym$tau[k,j]<0){
              relationships[k,j] <- "pursuit (follower)"
              relationships[j,k] <- "pursuit (leader)"
            }else{
              relationships[k,j] <- "pursuit (leader)"
              relationships[j,k] <- "pursuit (follower)"
            }
          }else{
            #this is an unlikely branch where the two asymmetric tests have agreement but the symmetric test is non-significant 
            if (asym$tau[j,k]<0){
              #this is the second version of the weak convergence relationship
              relationships[k,j] <- "weak convergence (symmetric)"
              relationships[j,k] <- "weak convergence (symmetric)"
            }else{
              #this is the second version of the weak divergence relationship
              relationships[k,j] <- "weak divergence (symmetric)"
              relationships[j,k] <- "weak divergence (symmetric)"
            }
          }
        }else{
          #this is the branch where only one asymmetric test is significant
          if (asym$tau[c(j,k),c(j,k)][which(asym$p.value[c(j,k),c(j,k)]<=alpha)]>0){
            #this is the Departing relationship where there is only one divergent asymmetric test significant
            if (asym$p.value[k,j]<=alpha){
              relationships[k,j] <- "departing (departer)"
              relationships[j,k] <- "departing (origin)"
            }else{
              relationships[k,j] <- "departing (origin)"
              relationships[j,k] <- "departing (departer)"
            }
          }else{
            #this is the Approaching relationship where there is only one convergent asymmetric test significant
            if (asym$p.value[k,j]<=alpha){
              relationships[k,j] <- "approaching (approacher)"
              relationships[j,k] <- "approaching (target)"
            }else{
              relationships[k,j] <- "approaching (target)"
              relationships[j,k] <- "approaching (approacher)"
            }
          }
        }
      }
    }
  }
  
  #preparing the output
  output <- list()
  output$dynamic_relationships <- relationships
  output$symmetric_convergence <- sym
  output$asymmetric_convergence <- asym
  output$correspondence <- Dcor
  output$parameters <- c(alphaUncor,alpha,nperm)
  names(output$parameters) <- c("alpha","corrected_alpha","nperm")
  
  output[["dynamic_relationships_taxonomy"]] <- data.frame(dynamic_relationship=c("neutral (symmetric)","parallel (symmetric)","antiparallel (symmetric)",
                                                                                "convergence (symmetric)","weak convergence (symmetric)",
                                                                                "divergence (symmetric)","weak divergence (symmetric)",
                                                                                "approaching (approacher)","approaching (target)",
                                                                                "departing (departer)","departing (origin)",
                                                                                "approaching-stationary (approacher)",
                                                                                "approaching-stationary (stationary target)",
                                                                                "departing-stationary (departer)",
                                                                                "departing-stationary (stationary origin)",
                                                                                "catch-up (leader)","catch-up (follower)",
                                                                                "pursuit (leader)","pursuit (follower)",
                                                                                "escape (leader)","escape (follower)"),
                                                         conv_div_group=c(NA,NA,NA,
                                                                          "convergence group","convergence group",
                                                                          "divergence group","divergence group",
                                                                          "convergence group","convergence group",
                                                                          "divergence group","divergence group",
                                                                          "convergence group",
                                                                          "convergence group",
                                                                          "divergence group",
                                                                          "divergence group",
                                                                          "convergence group","convergence group",
                                                                          NA,NA,
                                                                          "divergence group","divergence group"),
                                                         oriented_group=c(NA,NA,NA,
                                                                          NA,NA,
                                                                          NA,NA,
                                                                          "oriented group (back)","oriented group (front)",
                                                                          "oriented group (front)","oriented group (back)",
                                                                          "oriented group (back)",
                                                                          "oriented group (front)",
                                                                          "oriented group (front)",
                                                                          "oriented group (back)",
                                                                          "oriented group (front)","oriented group (back)",
                                                                          "oriented group (front)","oriented group (back)",
                                                                          "oriented group (front)","oriented group (back)"))
  rownames(output$dynamic_relationships_taxonomy) <- output$dynamic_relationships_taxonomy$dynamic_relationship
  
  #define its class
  class(output) <- c("RTMA","list")
  
  return(output)
}
