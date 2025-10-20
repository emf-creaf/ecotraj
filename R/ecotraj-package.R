#' ecotraj: Ecological Trajectory Analysis
#'
#' Analysis of temporal changes (i.e. dynamics) of ecological entities, defined as trajectories on a chosen multivariate space
#'
#' @name ecotraj-package
#' @aliases ecotraj ecotraj-package
#' @docType package
#' @author \strong{Maintainer}: Miquel De Cáceres
#' \email{miquelcaceres@@gmail.com}
#' [ORCID](https://orcid.org/0000-0001-7132-2080)
#'
#' Authors: \itemize{
#' \item{ Nicolas Djeghri
#' [ORCID](https://orcid.org/0000-0001-5740-3386)}
#' \item{ Anthony Sturbois
#' [ORCID](https://orcid.org/0000-0002-9219-4468)}
#' }
#' Contributors: \itemize{
#' \item{ Javier De la Casa}
#' }
#'  
#' @seealso Useful links: \itemize{ \item{
#' \url{https://emf-creaf.github.io/ecotraj/index.html}} }
#'
#' @references De Caceres et al., 2019 (\doi{10.1002/ecm.1350}), Sturbois et al., 2021 (\doi{10.1016/j.ecolmodel.2020.109400}), Sturbois et al., 2023 (\doi{10.1002/ecs2.4726}).
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Kendall MannKendall
#' @importFrom MASS ginv
#' @importFrom graphics arrows text rect segments points par polygon
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices col2rgb rgb
#' @importFrom stats as.dist cmdscale dist model.matrix quantile
#' @importFrom utils combn
#' @useDynLib ecotraj, .registration = TRUE
## usethis namespace: end
NULL