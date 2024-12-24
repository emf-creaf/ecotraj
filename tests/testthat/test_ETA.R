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
  expect_s3_class(trajectorySpeeds(x), "data.frame")
  expect_s3_class(trajectoryAngles(x), "data.frame")
  expect_s3_class(trajectoryAngles(x, relativeToInitial = TRUE), "data.frame")
  expect_s3_class(trajectoryAngles(x, all = TRUE), "data.frame")
  expect_s3_class(trajectoryAngles(x, relativeToInitial = TRUE, all = TRUE), "data.frame")
  expect_type(trajectoryDirectionality(x), "double")
  expect_type(trajectoryVariability(x), "double")
  expect_type(segmentDistances(x), "list")
  expect_s3_class(trajectoryDistances(x), "dist")
  expect_type(trajectoryConvergence(x), "list")
})