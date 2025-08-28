library(ecotraj)

#Description of sites, surveys and times
sites <- c("1","1","1","2","2","2")
surveys<-c(1,2,3,1,2,3)
times <- c(1.2, 2.1, 3.3, 1.3, 2.0, 3.1)
times_perm <- c(1.2, 2.1, 3.3, 1.3, 3.1, 2.0) # Permutes 5 and 6
surveys_perm<-c(1,2,3,1,3,2)

#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1
xy_perm <- rbind(xy[1:4,], xy[6,], xy[5,])

d <- dist(xy)
d_perm <- dist(xy_perm)

test_that("Trajectories are well defined",{
  expect_s3_class(defineTrajectories(d, sites), "trajectories")
  expect_s3_class(defineTrajectories(d, sites, surveys), "trajectories")
  expect_s3_class(defineTrajectories(d, sites, surveys, times), "trajectories")
  expect_s3_class(defineTrajectories(d, sites, times = times), "trajectories")
  expect_s3_class(subsetTrajectories(defineTrajectories(d, sites, surveys, times), "2"), "trajectories")
  expect_s3_class(subsetTrajectories(defineTrajectories(d, sites, surveys, times), 
                                     site_selection = c("1", "2"), 
                                     survey_selection = c(1,3)), "trajectories")
  expect_s3_class(defineTrajectories(as.matrix(d), sites), "trajectories")
  expect_s3_class(defineTrajectories(as.matrix(d), sites, surveys), "trajectories")
  expect_s3_class(defineTrajectories(as.matrix(d), sites, surveys, times), "trajectories")
  expect_s3_class(defineTrajectories(as.matrix(d), sites, times = times), "trajectories")
  t <- defineTrajectories(as.matrix(d_perm), sites, times = times_perm)
  expect_equal(t$metadata$surveys, surveys_perm)
})


test_that("Trajectories are wrongly defined",{
  expect_error(defineTrajectories(d, sites[1:5]))
  expect_error(defineTrajectories(d, sites, surveys[1:5]))
  expect_error(defineTrajectories(d, sites, surveys, times[1:5]))
  expect_error(defineTrajectories(d, sites, as.character(surveys)))
  expect_error(defineTrajectories(d, sites, surveys, as.character(times)))
})

test_that("Trajectories can be transformed", {
  x <- defineTrajectories(d, sites, surveys)
  expect_s3_class(centerTrajectories(x),"trajectories")
  expect_s3_class(smoothTrajectories(x),"trajectories")
  expect_s3_class(interpolateTrajectories(x, c(1.5, 2.5)),"trajectories")
})

