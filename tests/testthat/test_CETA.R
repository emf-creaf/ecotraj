
#Let's define our toy sampling times:
timesToy <- 0:30 #The sampling times of the time series
cycleDurationToy <- 10 #The duration of the cycles (i.e. the periodicity of the time series)
datesToy <- timesToy%%cycleDurationToy #The dates associated to each times

#And state where the sampling occurred, for now let's only use one site "A"
sitesToy <- rep(c("A"),length(timesToy))

#Then prepare some toy data:
#Prepare a noise and trend term to make the data more interesting
noise <- 0.05
trend <- 0.05

#Make cyclical data (note that we apply the trend only to x:
x <- sin((timesToy*2*pi)/cycleDurationToy)+rnorm(length(timesToy),mean=0,sd=noise)+trend*timesToy
y <- cos((timesToy*2*pi)/cycleDurationToy)+rnorm(length(timesToy),mean=0,sd=noise)
matToy <- cbind(x,y)

#And express it as a distance matrix (ETA is based on distances, increasing its generality)
dToy <- dist(matToy)

# Define cyclical trajectory
xToy <- defineTrajectories(dToy, sites = sitesToy, times = timesToy)

test_that("Fixed date trajectories can be build and analyzed",{
  fdtrajToy <- extractFixedDateTrajectories(xToy,
                                            cycleDuration = cycleDurationToy)
  expect_s3_class(fdtrajToy, "fd.trajectories")
  expect_type(trajectoryDirectionality(fdtrajToy), "double")
  expect_type(trajectoryVariability(fdtrajToy), "double")
  expect_s3_class(trajectoryDistances(fdtrajToy), "dist")
})

test_that("Cycles can be build and analyzed",{
  cycleToy <- extractCycles(xToy,
                            cycleDuration = cycleDurationToy)
  expect_s3_class(cycleToy, "cycles")
  expect_s3_class(trajectoryLengths(cycleToy), "data.frame")
  expect_type(trajectoryDirectionality(cycleToy), "double")
  expect_type(trajectoryVariability(cycleToy), "double")
  expect_s3_class(trajectoryDistances(cycleToy), "dist")

  expect_type(cycleConvexity(xToy, cycleDuration = cycleDurationToy), "double")
  expect_s3_class(cycleShifts(xToy, cycleDuration = cycleDurationToy), "data.frame")
  
  cycleToy_cent <- centerTrajectories(cycleToy,
                                      exclude = which(cycleToy$metadata$internal==FALSE))
  expect_s3_class(cycleToy_cent, "cycles")
  expect_s3_class(trajectoryDistances(cycleToy_cent), "dist")
})


