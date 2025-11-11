# Trajectory definition

Defines data structures for trajectory analysis

## Usage

``` r
defineTrajectories(d, sites, surveys = NULL, times = NULL)
```

## Arguments

- d:

  A symmetric [`matrix`](https://rdrr.io/r/base/matrix.html) or an
  object of class [`dist`](https://rdrr.io/r/stats/dist.html) containing
  the distance values between pairs of ecological states..

- sites:

  A character vector indicating the ecological entity (site, individual,
  community) corresponding to each ecological state (other types are
  converted to character).

- surveys:

  An integer vector indicating the survey corresponding to each
  ecological state (only necessary when surveys are not in order and
  `times` is not provided).

- times:

  A numeric vector indicating survey times.

## Value

An object (list) of class `trajectories` with the following elements:

- `d`: An object of class [`dist`](https://rdrr.io/r/stats/dist.html)
  containing relationships between ecological states

- `metadata`: A data frame describing trajectory states, with the
  following columns:

  - `sites`: A character vector indicating the ecological entity
    corresponding to each ecological state.

  - `surveys`: An integer vector indicating the survey corresponding to
    each ecological state.

  - `times`: A numeric vector indicating survey times.

## Details

If `surveys` is not provided, but `times` is available, surveys will be
taken as the order of times. Otherwise, `surveys` will be assumed to be
in order for all the occurrences of the same value of `sites`. If
`times` is not provided, then it is made equal to `surveys`.

## See also

[`subsetTrajectories`](https://emf-creaf.github.io/ecotraj/reference/subsetTrajectories.md)

## Examples

``` r
#Description of entities (sites) and surveys
entities <- c("1","1","1","2","2","2")
surveys <- c(1,2,3,1,2,3)
  
#Raw data table
xy<-matrix(0, nrow=6, ncol=2)
xy[2,2]<-1
xy[3,2]<-2
xy[4:6,1] <- 0.5
xy[4:6,2] <- xy[1:3,2]
xy[6,1]<-1

d <- dist(xy)

# Defines trajectories
x <- defineTrajectories(d, entities, surveys)
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
#> 1     1       1     1
#> 2     1       2     2
#> 3     1       3     3
#> 4     2       1     1
#> 5     2       2     2
#> 6     2       3     3
#> 
#> attr(,"class")
#> [1] "trajectories" "list"        
```
