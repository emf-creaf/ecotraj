####
# Utilities for circular statistics
# Inspired from package 'circular', to avoid dependencies from this package, orphaned by CRAN
####
.meancircular <- function(x, na.rm=FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x)==0) {
    warning("No observations (after removing missing values)")
    return(NA)
  }
  sinr = 0.0
  cosr = 0.0
  circmean = NA
  for(i in 1:length(x)) {
    sinr = sinr+ sin(x[i]*(pi/180))
    cosr = cosr+ cos(x[i]*(pi/180))
  }
  if (sqrt(sinr^2 + cosr^2)/length(x) > 0.0000001) circmean = (180/pi)*atan2(sinr, cosr)
  return(circmean)
}
.rhocircular <- function(x, na.rm=FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x)==0) {
    warning("No observations (after removing missing values)")
    return(NA)
  }
  sinr <- sum(sin(x*(pi/180)))
  cosr <- sum(cos(x*(pi/180)))
  return(sqrt(sinr^2 + cosr^2)/length(x))
}
.sdcircular <- function(x, na.rm=FALSE) {
  rbar <- .rhocircular(x, na.rm = na.rm)
  circsd <- sqrt(-2*log(rbar))
  return(circsd)
}
