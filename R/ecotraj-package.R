#' ecotraj: Ecological Trajectory Analysis
#'
#' Assists ecologists in the analysis of temporal changes of ecosystems, defined as trajectories on a chosen multivariate space
#'
#' @name ecotraj-package
#' @aliases ecotraj ecotraj-package
#' @docType package
#' @author \strong{Maintainer}: Miquel De Cáceres
#' \email{miquelcaceres@@gmail.com}
#' [ORCID](https://orcid.org/0000-0001-7132-2080)
#'
#' Authors: \itemize{
#' \item{ Anthony Sturbois}
#' }
#' Contributors: \itemize{
#' \item{ Javier De la Casa}
#' }
#'  
#' @seealso Useful links: \itemize{ \item{
#' \url{https://emf-creaf.github.io/ecotraj/index.html}} }
#'
#' @references De Caceres et al., 2019 (\doi{10.1002/ecm.1350}), Sturbois et al., 2021 (\doi{10.1016/j.ecolmodel.2020.109400})
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Kendall MannKendall
#' @importFrom MASS ginv
#' @importFrom graphics arrows text
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom stats as.dist cmdscale dist model.matrix quantile
#' @importFrom utils combn setTxtProgressBar txtProgressBar
#' @useDynLib ecotraj, .registration = TRUE
## usethis namespace: end
NULL