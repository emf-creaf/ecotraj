sites <- c("1","1","1","1","2","2","2","2","3","3","3","3")
surveys <- c(1,2,3,4,1,2,3,4,1,2,3,4)
times <- c(1.0,2.2,3.1,4.2,1.0,1.5,2.8,3.9,1.6,2.8,3.9,4.3)
xy<-matrix(0, nrow=12, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4,2]<-3
xy[5:6,2] <- xy[1:2,2]
xy[7,2]<-1.5
xy[8,2]<-2.0
xy[5:6,1] <- 0.25
xy[7,1]<-0.5
xy[8,1]<-1.0
xy[9:10,1] <- xy[5:6,1]+0.25
xy[11,1] <- 1.0
xy[12,1] <-1.5
xy[9:10,2] <- xy[5:6,2]
xy[11:12,2]<-c(1.25,1.0)
d <- dist(xy)

x <- defineTrajectories(d, sites, surveys, times)


test_that("Trajectories can be analyzed",{
  expect_s3_class(trajectoryLengths(x), "data.frame")
  expect_s3_class(trajectoryLengths(x, relativeToInitial = TRUE), "data.frame")
  expect_s3_class(trajectoryLengths(x, all = TRUE), "data.frame")
  expect_s3_class(trajectoryLengths(x, relativeToInitial = TRUE, all = TRUE), "data.frame")
  expect_equal(trajectoryLengths(x), trajectoryLengths2D(xy, sites, surveys))
  expect_s3_class(trajectorySpeeds(x), "data.frame")
  expect_equal(trajectorySpeeds(x), trajectorySpeeds2D(xy, sites, surveys, times))
  expect_s3_class(trajectoryAngles(x), "data.frame")
  expect_s3_class(trajectoryAngles(x, relativeToInitial = TRUE), "data.frame")
  expect_s3_class(trajectoryAngles(x, all = TRUE), "data.frame")
  expect_s3_class(trajectoryAngles(x, relativeToInitial = TRUE, all = TRUE), "data.frame")
  expect_type(trajectoryDirectionality(x), "double")
  expect_s3_class(trajectoryVariability(x), "data.frame")
  expect_type(segmentDistances(x), "list")
  expect_s3_class(trajectoryDistances(x, distance.type = "DSPD"), "dist")
  expect_s3_class(trajectoryDistances(x, distance.type = "SPD"), "dist")
  expect_s3_class(trajectoryDistances(x, distance.type = "TSPD"), "dist")
  expect_type(trajectoryConvergence(x, type = "pairwise.asymmetric"), "list")
  expect_type(trajectoryConvergence(x, type = "pairwise.symmetric"), "list")
  expect_error(trajectoryConvergence(x, type = "multiple"))
  expect_s3_class(trajectoryMetrics(x), "data.frame")
  expect_s3_class(trajectoryShifts(x), "data.frame")
  expect_s3_class(trajectoryWindowMetrics(x, 1), "data.frame")
})

# Shuffle surveys and check if the results are the same
xy2 <- xy
surveys2 <- surveys
times2  <- times
sites2 <- sites
temp = xy2[5,]
xy2[5,] = xy2[6,]
xy2[6,] = temp
surveys2[5] = 2
surveys2[6] = 1
temp = times2[5]
times2[5] = times2[6]
times2[6] = temp
x2 <- defineTrajectories(dist(xy2), sites2, surveys2, times2)

test_that("Trajectory analysis gives the same result after shuffling surveys",{
  expect_equal(trajectoryLengths(x), trajectoryLengths(x2))
  expect_equal(trajectoryLengths2D(xy, sites, surveys), trajectoryLengths2D(xy2, sites2, surveys2))
  expect_equal(trajectorySpeeds(x), trajectorySpeeds(x2))
  expect_equal(trajectorySpeeds2D(xy, sites, surveys, times), trajectorySpeeds2D(xy2, sites2, surveys2, times2))
  expect_equal(trajectoryAngles(x), trajectoryAngles(x2))
  expect_equal(trajectoryDirectionality(x), trajectoryDirectionality(x2))
  expect_equal(trajectoryVariability(x), trajectoryVariability(x2))
  expect_equal(segmentDistances(x), segmentDistances(x2))
  expect_equal(trajectoryDistances(x, distance.type = "DSPD"), trajectoryDistances(x2, distance.type = "DSPD"))
  expect_equal(trajectoryDistances(x, distance.type = "SPD"), trajectoryDistances(x2, distance.type = "SPD"))
  expect_equal(trajectoryDistances(x, distance.type = "TSPD"), trajectoryDistances(x2, distance.type = "TSPD"))
  expect_equal(trajectoryConvergence(x), trajectoryConvergence(x2))
  expect_equal(trajectoryShifts(x), trajectoryShifts(x2))
})
