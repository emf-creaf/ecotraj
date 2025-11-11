# Dynamic variation and variation decomposition

- Function `dynamicVariation` assesses the amount of dynamic variation
  observed across trajectories and the relative contribution of each of
  them.

- Function `variationDecomposition` performs a sum of squares
  decomposition of total variation in three components: (1) across
  trajectories (entities); (2) across time points; (3) their
  interaction.

## Usage

``` r
dynamicVariation(x, ...)

variationDecomposition(x)
```

## Arguments

- x:

  An object of class `trajectories` (or its children subclasses
  `fd.trajectories` or `cycles`).

- ...:

  Additional params to be passed to function
  [`trajectoryDistances`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md).

## Value

- Function `dynamicVariance` returns a list with three elements (dynamic
  sum of squares, dynamic variance and a vector of trajectory relative
  contributions)

- Function `variationDecomposition` returns a data frame with results
  (sum of squares, degrees of freedom and variance estimates) for each
  variance component and the total.

## Details

Function `variationDecomposition` requires trajectories to be
synchronous. The SS sum of `temporal` and `interaction` components
correspond to the SS sum, across trajectories, of function
[`trajectoryInternalVariation`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md).

## See also

[`defineTrajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md),
[`is.synchronous`](https://emf-creaf.github.io/ecotraj/reference/is.synchronous.md),
[`trajectoryDistances`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md),
[`trajectoryInternalVariation`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)

## Examples

``` r
#Description of entities and surveys
entities <- c("1","1","1","1","2","2","2","2","3","3","3","3")
surveys <- c(1,2,3,4,1,2,3,4,1,2,3,4)
  
#Raw data table
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

# Defines trajectories
x <- defineTrajectories(d, entities, surveys)

# Assessment of dynamic variation and individual trajectory contributions
dynamicVariation(x)
#> $dynamic_ss
#> [1] 0.7114251
#> 
#> $dynamic_variance
#> [1] 0.3557125
#> 
#> $relative_contributions
#>          1          2          3 
#> 0.51366204 0.06351261 0.42282535 
#> 

# Variation decomposition (entity, temporal and interaction) for synchronous 
# trajectories:
variationDecomposition(x)
#>                    ss df  variance
#> entities     2.489583  2 1.2447917
#> time         7.453125  3 2.4843750
#> interaction  1.718750  6 0.2864583
#> total       11.661458 11 1.0601326

# check the correspondence with internal variation
sum(variationDecomposition(x)[c("time", "interaction"),"ss"])
#> [1] 9.171875
sum(trajectoryInternalVariation(x)$internal_ss)
#> [1] 9.171875
```
