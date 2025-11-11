# Synchronicity in trajectory observations

Checks whether trajectories are synchronous, meaning that observation
times are equal

## Usage

``` r
is.synchronous(x)
```

## Arguments

- x:

  An object of class `trajectories` (or its children subclasses
  `fd.trajectories` or `cycles`)

## Value

A boolean indicating whether trajectories are synchronous

## See also

[`defineTrajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)

## Examples

``` r
#Description of sites and surveys
sites <- c("1","1","1","2","2","2")
surveys <- c(1,2,3,1,2,3)
  
#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1

#Synchronous trajectories
x1 <- defineTrajectories(dist(xy), sites, surveys)
is.synchronous(x1)
#> [1] TRUE

# Non synchronous trajectories
x2 <- defineTrajectories(dist(xy[1:5,]), sites[1:5], surveys[1:5])
is.synchronous(x2)
#> [1] FALSE
```
