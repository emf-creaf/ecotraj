#' Glomel vegetation dataset
#'
#' Vegetation data set to illustrate Ecological Quality Assessment (EQA) 
#' 
#' The nature reserve of Landes et Marais de Glomel (Brittany, France) is
#' composed of temperate Atlantic wet heaths whose reference state is commonly
#' considered dominated by plant communities associated to acid, nutrient poor soils
#' that are at least seasonally water logged and dominated by \emph{Erica tetralix} and 
#' \emph{E. ciliaris}. The data set consists of 23 rows and 46 columns. 
#' The first five stations (rows) were used to define the reference
#' envelope, and the next 18 stations (rows) where those for which the conservation status 
#' was to be assessed.
#' 
#' @encoding UTF-8
#' @format Glomel is an object of class data.frame composed of 23 observations and 46 variables.
#'
#' \describe{
#'  \item{ID}{Station ID.}
#'  \item{Ref}{Logical flag to indicate stations used to define the reference envelope.}
#'  \item{Complementary}{Comments regarding the quality of the ecosystem.}
#'  \item{...}{Percent cover values (derived from Braun-Blanquet ordinal scale) for 43 species of vascular plants.}
#' }
#'
#' @name glomel
#' @aliases glomel 
#' @docType data
#' @author XXX
#' @keywords data
#' 
#' @seealso \code{\link{envelope}}
NULL


#' Glenan dataset
#'
#' Maerl bed data set to illustrate Ecological Quality Assessment (EQA) 
#' 
#' Experimental data set built by Tauran et al. (2020) to study the impact of fishing dredges 
#' and varying fishing pressures on maerl beds, in the bay of Brest (Brittany, France). 
#' 
#' @encoding UTF-8
#' @format Glenan is an object of class data.frame composed of 32 observations and 252 variables.
#'
#' \describe{
#'  \item{Abundance.\emph{x}}{Abundance (number of individuals) of each taxon \emph{x}}
#'  \item{Surveys}{Indicates different Maerl bed surveys.}
#'  \item{Treatment}{Combinations of fishing dredges and pressure levels. 'CTRL' stands for control. Fishing dredges are: 
#'    \itemize{
#'       \item{(1) a clam dredge (CD), 70 to 90 kg, 1.5 m wide, 40 teeth of 11 cm each;}
#'       \item{(2) a queen scallop dredge (QSD), 120 kg,1.8 m wide, with a blade;}
#'       \item{(3) a king scallop dredge (KSD), 190 kg, 1.8 m wide, 18 teeth of 10 cm each every 9 cm.}
#'    }
#'  }
#' }
#'
#' @name glenan
#' @aliases glenan 
#' @docType data
#' @keywords data
#' 
#' @references 
#' Tauran, A., Dubreuil, J., Guyonnet, B., Grall, J., 2020. Impact of fishing gears and fishing intensities on maerl beds: An experimental approach. Journal of Experimental Marine Biology and Ecology 533, 151472. https://doi.org/10.1016/j.jembe.2020.151472
#' 
#' @seealso \code{\link{envelope}}
NULL
