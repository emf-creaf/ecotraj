
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

test_that("Fixed date trajectories can be build, subset and analyzed",{
  fdtrajToy <- extractFixedDateTrajectories(xToy,
                                            cycleDuration = cycleDurationToy)
  expect_s3_class(fdtrajToy, "fd.trajectories")
  expect_s3_class(subsetTrajectories(fdtrajToy, site_selection = "A"), "fd.trajectories")
  expect_s3_class(subsetTrajectories(fdtrajToy, subtrajectory_selection = "A_fdT_3"), "fd.trajectories")
  expect_s3_class(centerTrajectories(fdtrajToy), "fd.trajectories")
  expect_s3_class(smoothTrajectories(fdtrajToy), "fd.trajectories")
  expect_s3_class(trajectoryLengths(fdtrajToy), "data.frame")
  expect_s3_class(trajectorySpeeds(fdtrajToy), "data.frame")
  expect_type(trajectoryDirectionality(fdtrajToy), "double")
  expect_type(trajectoryVariability(fdtrajToy), "double")
  expect_s3_class(trajectoryMetrics(fdtrajToy), "data.frame")
  expect_type(segmentDistances(fdtrajToy), "list")
  expect_s3_class(trajectoryDistances(fdtrajToy), "dist")
})

test_that("Cycles can be build, subset and analyzed",{
  cycleToy <- extractCycles(xToy,
                            cycleDuration = cycleDurationToy)
  expect_s3_class(cycleToy, "cycles")
  expect_s3_class(subsetTrajectories(cycleToy, site_selection = "A"), "cycles")
  expect_s3_class(subsetTrajectories(cycleToy, subtrajectory_selection = "A_C3"), "cycles")
  expect_s3_class(centerTrajectories(cycleToy), "cycles")
  expect_s3_class(smoothTrajectories(cycleToy), "cycles")
  expect_s3_class(trajectoryLengths(cycleToy), "data.frame")
  expect_s3_class(trajectorySpeeds(cycleToy), "data.frame")
  expect_warning(trajectoryMetrics(cycleToy))
  expect_s3_class(cycleMetrics(xToy,
                               cycleDuration = cycleDurationToy), "data.frame")
  expect_type(trajectoryDirectionality(cycleToy), "double")
  expect_type(trajectoryVariability(cycleToy), "double")
  expect_type(segmentDistances(cycleToy), "list")
  expect_s3_class(trajectoryDistances(cycleToy), "dist")

  expect_type(cycleConvexity(xToy, cycleDuration = cycleDurationToy), "double")
  expect_s3_class(cycleShifts(xToy, cycleDuration = cycleDurationToy), "data.frame")
  
  cycleToy_cent <- centerTrajectories(cycleToy,
                                      exclude = which(cycleToy$metadata$internal==FALSE))
  expect_s3_class(cycleToy_cent, "cycles")
  expect_s3_class(trajectoryDistances(cycleToy_cent), "dist")
})


