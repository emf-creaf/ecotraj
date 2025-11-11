# Trajectory subsetting

Subsets data structures for trajectory analysis

## Usage

``` r
subsetTrajectories(
  x,
  site_selection = NULL,
  subtrajectory_selection = NULL,
  survey_selection = NULL,
  window_selection = NULL
)
```

## Arguments

- x:

  An object of class `trajectories` (or its children subclasses
  `fd.trajectories` or `cycles`)

- site_selection:

  A character vector indicating the subset of entity (site) trajectories
  to be selected (if NULL, all sites are included).

- subtrajectory_selection:

  A character vector indicating the subset of cycles or fixed date
  trajectories to be selected (only used when `x` is of class
  `fd.trajectories` or `cycles`).

- survey_selection:

  An integer vector indicating the subset of surveys to be included (if
  NULL, all surveys are included).

- window_selection:

  An ordered pair of time values (e.g. `c(lower, upper)`) to subset the
  observations to a time window.

## Value

An object (list) of class `trajectories` (or its children subclasses
`fd.trajectories` or `cycles`), depending on the input.

## Details

When using function `subsetTrajectories` on cycles or fixed-date
trajectories then the parameter `site_selection` applies to sites (hence
allows selecting multiple cycles or fixed-date trajectories). Specific
cycles or fixed-date trajectories can be selected using
`trajectory_selection`.

## See also

[`defineTrajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md),
[`trajectoryCyclical`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)

## Examples

``` r
#Description of entities surveys and times
entities <- c("1","1","1","2","2","2")
surveys <- c(1,2,3,1,2,3)
times <- c(10, 20, 35, 10, 20, 35)
  
#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1

d <- dist(xy)

# Defines trajectories
x <- defineTrajectories(d, entities, surveys, times = times)
x
#> $d
#>          1        2        3        4        5
#> 2 1.000000                                    
#> 3 2.000000 1.000000                           
#> 4 0.500000 1.118034 2.061553                  
#> 5 1.118034 0.500000 1.118034 1.000000         
#> 6 2.236068 1.414214 1.000000 2.061553 1.118034
#> 
#> $metadata
#>   sites surveys times
#> 1     1       1    10
#> 2     1       2    20
#> 3     1       3    35
#> 4     2       1    10
#> 5     2       2    20
#> 6     2       3    35
#> 
#> attr(,"class")
#> [1] "trajectories" "list"        

# Extracts (subset) second trajectory
x_2 <- subsetTrajectories(x, "2")
x_2
#> $d
#>          1        2
#> 2 1.000000         
#> 3 2.061553 1.118034
#> 
#> $metadata
#>   sites surveys times
#> 1     2       1    10
#> 2     2       2    20
#> 3     2       3    35
#> 
#> attr(,"class")
#> [1] "trajectories" "list"        

# Extracts window corresponding to observation times 20, 35
x_3 <- subsetTrajectories(x, window_selection = c(20, 35))
x_3
#> $d
#>          1        2        3
#> 2 1.000000                  
#> 3 0.500000 1.118034         
#> 4 1.414214 1.000000 1.118034
#> 
#> $metadata
#>   sites surveys times
#> 1     1       2    20
#> 2     1       3    35
#> 3     2       2    20
#> 4     2       3    35
#> 
#> attr(,"class")
#> [1] "trajectories" "list"        
```
