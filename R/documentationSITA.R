#' furseals dataset
#'
#' This is a subset of a data sets from Kernaléguen et al. (2015). 
#' 
#' Briefly, fur seals the Antarctic fur seal Arctocephalus gazella and subantarctic fur seal A. tropicalis whisker SI values yield unique long-term information
#' on individual behaviour which integrates the spatial, trophic and temporal dimensions of the ecological niche. The foraging strategies of this two species of sympatric 
#' fur seals were examined in the winter 2001/2002 at Crozet, Amsterdam and Kerguelen Islands (Southern Ocean) using the stable isotope values of serially sampled whiskers. 
#' The subset of the initial data set is composed of consecutive whisker sections (3 mm-long) starting from the proximal (facial) end, with the most recently synthesized tissue remaining 
#' under the skin. Only individuals (n = 47) with whiskers totalizing at least 30 sections were selected in the initail data, and only those 30 sections were selected.
#' 
#' @encoding UTF-8
#' @format furseals is an object of class data.frame composed of 1414 observations and 8 variables.
#'
#' \describe{
#'  \item{ID_SITA}{Fur seal ID used by Sturbois et al. (under review), from 1 to 47}
#'  \item{ID}{Fur seal ID used by Kernaléguen et al. (2015) in the initial data set.}
#'  \item{Species}{Fur seal species: the Antarctic fur seal Arctocephalus gazella or the subantarctic fur seal A. tropicalis.}
#'  \item{Sexe}{Fur seal gender, either 'Male' or 'Female'.}
#'  \item{Time}{Number of the whisker sections from 1 to 30.}
#'  \item{Place}{Breeding place: Crozet, Amsterdam or Kerguelen}
#'  \item{d13C}{δ 13C value}
#'  \item{d15N}{δ 15N value}
#' }
#'
#' @name furseals
#' @aliases furseals 
#' @docType data
#' @author Kernaléguen, L., Arnould, J.P.Y., Guinet, C., Cherel, Y.
#' @keywords data
#' @references 
#' \enc{Kernaléguen}{Kernaleguen}, L., Arnould, J.P.Y., Guinet, C., Cherel, Y., 2015. Determinants of individual foraging specialization inlarge marine vertebrates, the Antarctic and subantarctic fur seals. Journal of Animal Ecology 1081–1091.
NULL
#' pike dataset
#'
#' This data sets comes from Cucherousset et al. (2013). 
#' 
#' Briefly, Cucherousset et al. (2013) released 192 individually tagged, hatchery-raised, juvenile pike (Esox lucius L.) with variable initial trophic position (fin δ 13C/δ 15N values). 
#' Based on δ values, individuals were classified into zooplanktivorous (δ 15N < 10 ‰) and piscivorous (δ 15N > 10 ‰) as cannibalism is commonly observed in this species. 
#' Individuals were released in a temporarily flooded grassland where pike eggs usually hatch of the Brière marsh (France) to identify the determinants of juvenile natal departure. 
#' The release site was connected through a unique point to an adjacent pond used as a nursery habitat. Fish were continuously recaptured when migrating from flooded grassland 
#' to adjacent pond. Recaptured individuals (n = 29) were anaesthetized, checked for tags, measured for fork length, fin-clipped to quantify changes in δ 13C and δ 15N values, 
#' and released.
#'
#' @encoding UTF-8
#' @format pike is an object of class dataframe composed of 58 observations of 10 variables.
#' 
#' \describe{
#'  \item{trophic_status_initial}{Initial trophic status at release}
#'  \item{ID}{ID used for each individual by Cucherousset et al. (2013)}
#'  \item{Time}{Time of the stable isotope measurement: 1 (Release) or 2 (Departure)}
#'  \item{Time_L}{Time of the stable isotope measurement as string, either 'Release' or 'Departure'}
#'  \item{Date}{Date of release (common for all individuals) or recapture (variable dependind of the date of departure)}
#'  \item{Size_mm}{Size (length) of juvenile pike, in mm}
#'  \item{d13C}{δ 13C values}
#'  \item{d15N}{δ 15N values}
#'  \item{Residence_time}{Number of days between the release and the recapture} 
#'  \item{Trophic_status_final}{Trophic status at the end of the study}
#' }
#'
#' @name pike
#' @aliases pike 
#' @docType data
#' @author Cucherousset, J., Paillisson, J.-M., Roussel, J.-M.
#' @keywords data
#' @references 
#' Cucherousset, J., Paillisson, J.-M., Roussel, J.-M., 2013. Natal departure timing from spatially varying environments is dependent of individual ontogenetic status. Naturwissenschaften 100, 761–768.
NULL

#' isoscape dataset
#'
#' This data sets is a subset from Espinasse et al. (2020). 
#' 
#' Briefly, Espinasse et al. (2020) tested the application of isoscapes modelled from satellite data to the description of secondary production in the Northeast pacific. 
#' The output model fits in a 0.25° x 0.25° spatial grid covering the region spanning from 46 to 62°N and from 195 to 235°E and supporting δ 13C and δ 15N isoscapes 
#' from 1998 to 2017. The subset is composed of modelled  δ 13C and δ 15N values of a 1° x 1° spatial grid from the original modelled dataset for 2013 and 2015. 
#'
#' @encoding UTF-8
#' @format isoscape is an object of class dataframe composed of 978 observations of 6 variables.
#' \describe{
#'  \item{Latitude}{Latitude coordinate of the station, in degrees}
#'  \item{Longitude}{Longitude coordinate of the station, in degrees}
#'  \item{d13C}{δ 13C modelled value}
#'  \item{d15N}{δ 15N modelled value}
#'  \item{station}{station ID}
#'  \item{Year}{Year corresponding to modelled stable isotope values}
#' }
#'
#' @name isoscape
#' @aliases isoscape
#' @docType data
#' @author Espinasse, B., Hunt, B.P.V., Batten, S.D., Pakhomov, E.A.
#' @keywords data
#' @seealso heatmapdata 
#' @references 
#' Espinasse, B., Hunt, B.P.V., Batten, S.D., Pakhomov, E.A., 2020. Defining isoscapes in the Northeast Pacific as an index of ocean productivity. Global Ecol Biogeogr 29, 246–261.
NULL

#' heatmapdata dataset
#'
#' Espinasse et al. (2020) tested the application of isoscapes modelled from satellite data to the description of secondary production in the Northeast pacific. 
#' The output model fits in a 0.25° x 0.25° spatial grid covering the region spanning from 46 to 62°N and from 195 to 235°E and supporting δ 13C and δ 15N isoscapes 
#' from 1998 to 2017.
#'
#' This data sets is composed of trajectory metrics calculated by Sturbois et al. (2021) for all stations within all inter-annual consecutive periods between 1998 and 2017 
#' calculated from the whole data set of Espinasse et al. (2020) for a 1° x 1° spatial grid.
#'
#' @encoding UTF-8
#' @format heatmapdata is an object of class dataframe composed of 9206 observations of 9 variables.
#' 
#' \describe{
#'  \item{Latitude}{Latitude coordinate of the station, in degrees}
#'  \item{Longitude}{Longitude coordinate of the station, in degrees}
#'  \item{d13C}{δ13C modelled value}
#'  \item{d15N}{δ15N modelled value}
#'  \item{station}{Station ID}
#'  \item{Years}{Period corresponding to the calculation of trajectory metrics}
#'  \item{Angles}{Angle α (i.e direction) in the stable isotope space}
#'  \item{Lengths}{Net change values (i.e direction) in the stable isotope space}
#'  \item{Angles2}{Angle alpha values (i.e direction) in the stable isotope space transformed for a potential use with function \code{geom_spoke}}
#' }
#'
#' @name heatmapdata
#' @aliases heatmapdata
#' @docType data
#' @author Espinasse, B., Hunt, B.P.V., Batten, S.D., Pakhomov, E.A.
#' @keywords data
#' @seealso isoscape 
#' @references 
#' Espinasse, B., Hunt, B.P.V., Batten, S.D., Pakhomov, E.A., 2020. Defining isoscapes in the Northeast Pacific as an index of ocean productivity. Global Ecol Biogeogr 29, 246–261.
NULL




