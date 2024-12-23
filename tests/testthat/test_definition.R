library(ecotraj)

#Description of sites and surveys
sites = c("1","1","1","2","2","2")
surveys=c(1,2,3,1,2,3)

#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1

d <- dist(xy)


test_that("Trajectories are well defined",{
  x <- defineTrajectories(d, sites, surveys)
  expect_s3_class(x, "trajectories")
  expect_s3_class(subsetTrajectories(x, "2"), "trajectories")
})
