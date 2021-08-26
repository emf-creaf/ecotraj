% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documentationSITA.R
\docType{data}
\name{pike}
\alias{pike}
\title{pike dataset}
\format{
pike is an object of class dataframe composed of 58 observations of 10 variables.

\describe{
 \item{trophic_status_initial}{Initial trophic status at release}
 \item{ID}{ID used for each individual by Cucherousset et al. (2013)}
 \item{Time}{Time of the stable isotope measurement: 1 (Release) or 2 (Departure)}
 \item{Time_L}{Time of the stable isotope measurement as string, either 'Release' or 'Departure'}
 \item{Date}{Date of release (common for all individuals) or recapture (variable dependind of the date of departure)}
 \item{Size_mm}{Size (length) of juvenile pike, in mm}
 \item{d13C}{δ13C values}
 \item{d15N}{δ15N values}
 \item{Residence_time}{Number of days between the release and the recapture} 
 \item{Trophic_status_final}{Trophic status at the end of the study}
}
}
\description{
This data sets comes from Cucherousset et al. (2013).
}
\details{
Briefly, Cucherousset et al. (2013) released 192 individually tagged, hatchery-raised, juvenile pike (Esox lucius L.) with variable initial trophic position (fin δ13C/δ15N values). 
Based on δ15N values, individuals were classified into zooplanktivorous (δ15N < 10 ‰) and piscivorous (δ15N > 10 ‰) as cannibalism is commonly observed in this species. 
Individuals were released in a temporarily flooded grassland where pike eggs usually hatch of the Brière marsh (France) to identify the determinants of juvenile natal departure. 
The release site was connected through a unique point to an adjacent pond used as a nursery habitat. Fish were continuously recaptured when migrating from flooded grassland 
to adjacent pond. Recaptured individuals (n = 29) were anaesthetized, checked for tags, measured for fork length, fin-clipped to quantify changes in δ13C and δ15N values, 
and released.
}
\references{
Cucherousset, J., Paillisson, J.-M., Roussel, J.-M., 2013. Natal departure timing from spatially varying environments is dependent of individual ontogenetic status. Naturwissenschaften 100, 761–768.
}
\author{
Cucherousset, J., Paillisson, J.-M., Roussel, J.-M.
}
\keyword{data}