# Functions for Cyclical Ecological Trajectory Analysis

The Cyclical extension of Ecological Trajectory Analysis (CETA) aims at
allowing ETA to describe ecological trajectories presenting cyclical
dynamics such as seasonal or day/night cycles. We call such trajectories
"cyclical". CETA operates by subdividing cyclical trajectories into two
types of sub-trajectories of interest: cycles and fixed-date
trajectories.

- Cycles are sub-trajectories joining the ecological states belonging to
  the same cycle.

- Fixed-date trajectories are sub-trajectories joining the ecological
  states of the same date in different cycles (e.g. in a multi-annual
  cyclical trajectory with seasonality, a fixed-date trajectory might
  join all the ecological states associated with the January months of
  the different years).

We recommend reading the vignette on CETA prior to use it.The CETA
functions provided here achieve one of two goals:

1.  Reformatting data to analyze either cycles or fixed-date
    trajectories. The reformatted data can then be fed into existing ETA
    functions to obtain desired metrics (although special care need to
    be taken with cycles, see details).

2.  Providing new metrics relevant to cycles complementing other ETA
    functions.

## Usage

``` r
extractCycles(
  x,
  cycleDuration,
  dates = NULL,
  startdate = NA,
  externalBoundary = "end",
  minEcolStates = 3
)

extractFixedDateTrajectories(
  x,
  cycleDuration,
  dates = NULL,
  fixedDate = NULL,
  namesFixedDate = NULL,
  minEcolStates = 2
)

cycleConvexity(
  x,
  cycleDuration,
  dates = NULL,
  startdate = NA,
  externalBoundary = "end",
  minEcolStates = 3,
  add = TRUE
)

cycleShifts(
  x,
  cycleDuration,
  dates = NULL,
  datesCS = NULL,
  centering = TRUE,
  minEcolStates = 3,
  add = TRUE
)

cycleMetrics(
  x,
  cycleDuration,
  dates = NULL,
  startdate = NA,
  externalBoundary = "end",
  minEcolStates = 3,
  add = TRUE
)
```

## Arguments

- x:

  An object of class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  describing a cyclical trajectory.

- cycleDuration:

  A value indicating the duration of a cycle. Must be in the same units
  as times.

- dates:

  An optional vector indicating the dates (\< `cycleDuration`)
  corresponding to each ecosystem state. Must be in the same units as
  times. Defaults to times modulo cycleDuration (see details).

- startdate:

  An optional value indicating at which date the cycles must begin. Must
  be in the same units as times. Defaults to `min(dates)`.

- externalBoundary:

  An optional string, either `"end"` or `"start"`, indicating whether
  the start or end of the cycles must be considered "external". Defaults
  to `"end"`.

- minEcolStates:

  An optional integer indicating the minimum number of ecological states
  to return a fixed-date trajectory. Fixed-date trajectories comprising
  less ecological states than minEcolStates are discarded and do not
  appear in the output of the function. Defaults to 2.

- fixedDate:

  An optional vector of dates for which fixed-date trajectories must be
  computed. Defaults to `unique(dates)`, resulting in returning all
  possible fixed-date trajectories.

- namesFixedDate:

  An optional vector of names associated to each `fixedDate`. Defaults
  to `round(fixedDate,2)`.

- add:

  Flag to indicate that constant values should be added (local
  transformation) to correct triplets of distance values that do not
  fulfill the triangle inequality.

- datesCS:

  An optional vector indicating the dates for which a cyclical shift
  must be computed. Default to `unique(dates)` resulting in the
  computation of all possible cyclical shifts.

- centering:

  An optional boolean. Should the cycles be centered before computing
  cyclical shifts? Defaults to `TRUE`.

## Value

Function `extractCycles` returns the base information needed to describe
cycles. Its outputs are meant to be used as input for other ETA
functions. Importantly, within cycles, ecological states can be
considered "internal" or "external". Some operations and metrics within
ETA use all ecological states whereas others use only "internal" ones
(see details). Function `extractCycles` returns an object of class
`cycles` containing:

- `d`: an object of class [`dist`](https://rdrr.io/r/stats/dist.html),
  the new distance matrix describing the cycles. To take in account
  ecological states that are both the end of a cycle and the start of
  another,`d` contains duplications. As compared to the input matrix,
  `d` may present deletions of ecological states that do not belong to
  any cycles (e.g. due to `minEcolStates`))

- `metadata`: an object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) describing the
  ecological states in `d` with columns:

  - `sites`: the sites associated to each ecological states.

  - `Cycles`: the names of the cycle each ecological states belongs to.
    The cycle name is built by combining the site name with C1, C2,
    C3... in chronological order.

  - `surveys`: renumbering of the surveys to describe individual Cycles.

  - `times`: the times associated to each ecological states.

  - `internal`: a boolean vector with `TRUE` indicating "internal"
    ecological states whereas `FALSE` indicates "external" ecological
    states. This has implications for how the outputs of `extractCycles`
    are treated by other ETA functions (see details).

  - `dates`: the dates associated to each ecological states.

- `interpolationInfo`: an output that only appear if ecological states
  have been interpolated. It is used internally by plotting functions
  (see
  [`cyclePCoA`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclicalPlots.md))
  but is not intended to be of interest to the end user.

Function `extractFixedDateTrajectories` returns the base information
needed to describe fixed-date trajectories. Its outputs are meant to be
used as inputs for other ETA functions in order to obtain desired
metrics. Function `extractFixedDateTrajectories` returns an object of
class `fd.trajectories` containing:

- `d`: an object of class [`dist`](https://rdrr.io/r/stats/dist.html),
  the new distance matrix describing the fixed-date trajectories. As
  compared to the input matrix, `d` may present deletions of ecological
  states that do not belong to any fixed-date trajectories (e.g. due to
  `minEcolStates`))

- `metadata`: an object of class
  [`data.frame`](https://rdrr.io/r/base/data.frame.html) describing the
  ecological states in `d` with columns:

  - `sites`: the sites to each ecological states.

  - `fdT`: the names of the fixed-date trajectory each ecological states
    belongs to. The fixed-date trajectory name is built by combining the
    site name with "fdT" and the name of the fixed date (from
    `namesFixedDate`).

  - `surveys`: renumbering of the surveys to describe individual fixed
    date trajectories.

  - `times`: the times associated to each ecological states.

  - `dates`: the dates associated to each ecological states.

Function `cycleConvexity` returns the a vector containing values between
0 and 1 describing the convexity of cycles. Importantly, outputs of
`extractCycles` should not be used as inputs for `cycleConvexity` (see
details).

Function `cycleShifts` returns an object of class
[`data.frame`](https://rdrr.io/r/base/data.frame.html) describing
cyclical shifts (i.e. advances and delays). Importantly, outputs of
`extractCycles` should not be used as inputs for `cycleShifts` (see
details). The columns of the
[`data.frame`](https://rdrr.io/r/base/data.frame.html) are:

- `site`: the site for which each cycle shift has been computed.

- `dateCS`: the date for which a cycle shift has been computed.

- `timeCS`: the time of the ecological state for which a cycle shift has
  been computed (i.e. the time associated to the projected ecological
  state).

- `timeRef`: the time associated to the reference ecological state.

- `timeScale`: the time difference between the reference and the
  projected ecological state.

- `cyclicalShift`: the cyclical shift computed (an advance if positive,
  a delay if negative) in the same units as the times input.

Function `cycleMetrics` returns a data frame where rows are cycles and
columns are different cycle metrics.

## Details

CETA functions:

- Function `extractCycles` reformats an object of class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  describing one or more cyclical trajectories into a new object of
  class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  designed for the analysis cycles.

- Function `extractFixedDateTrajectories` reformats an object of class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  describing one or more cyclical trajectories into a new object of
  class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  designed for the analysis fixed-date trajectories.

- Function `cycleConvexity` computes the "convexity" of the cycles
  embedded in one or more cyclical trajectories.

- Function `cycleShifts` computes the cyclical shifts (i.e. advances and
  delays) that can be obtain from one or more cyclical trajectories.

CETA is a little more time-explicit than the rest of ETA. Hence the
parameter `times` is needed to initiate the CETA approach (classical ETA
functions can work from `surveys` which is only ordinal). CETA also
distinguishes between times and dates. Times represent linear time
whereas dates represent circular time (e.g. the month of year). Dates
are circular variables, coming back to zero when reaching their maximum
value `cycleDuration` corresponding to the duration of a cycle. In CETA,
dates are by default assumed to be `times` modulo `cycleDuration`. This
should fit many applications but if this is not the case (i.e. if there
is an offset between times and dates), dates can be specified. `dates`
however need to remain compatible with `times` and `cycleDuration` (i.e.
(times modulo cycleDuration) - (dates modulo cycleDuration) needs to be
a constant).

IMPORTANT: Cycles within CETA comprises both "internal" and "external"
ecological states (see the output of function `extractCycles`). This
distinction is a solution to what we call the "December-to-January
segment problem". Taking the example of a monthly resolved multi-annual
time series, a way to make cycles would be to take the set of ecological
states representing months from January to December of each year.
However, this omits the segment linking December of year Y to January of
year Y+1. However, including this segments means having two January
months in the same cycle. The proposed solution in CETA (in the case of
this specific example) is to set the January month of year Y+1 as
"external". "external" ecological states need a specific handling for
some operation in ETA, namely:

- Centering where external ecological states must be excluded from
  computation but included nonetheless in the procedure. This is handled
  automatically by the function
  [`centerTrajectories`](https://emf-creaf.github.io/ecotraj/reference/transformTrajectories.md).

- Trajectory internal variability, where external ecological states must
  be excluded. This handled directly by the
  [`trajectoryInternalVariation`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  function.

- Visualization through principal coordinate analysis of the cycles. The
  dedicated function
  [`cyclePCoA`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclicalPlots.md)
  must be preferred over
  [`trajectoryPCoA`](https://emf-creaf.github.io/ecotraj/reference/trajectoryPlot.md).

As a general rule the outputs of `extractCycles` should be used as
inputs in other, non-CETA function (e.g. `trajectoryDistances`). There
is three important exceptions to that rule: the functions
`cycleConvexity`, `cycleShifts` and `cycleMetrics`. Instead, the inputs
of these three functions should parallel the inputs of `extractCycles`
in a given analysis. For `cycleConvexity`, this is because convexity
uses angles obtained from the whole cyclical trajectory, and not only
the cycles. For `cycleShifts`, this is because cyclical shifts are not
obtained with respect to a particular set of cycles. For `cycleMetrics`,
this is because it calls `cycleConvexity`. The function instead compute
the most adapted set of cycles to obtain the metric.

Note: Function `cycleShifts` is computation intensive for large data
sets, it may not execute immediately.

Further information and detailed examples of the use of CETA functions
can be found in the associated vignette.

## References

Djeghri et al. (under review) Going round in cycles, but going
somewhere: Ecological Trajectory Analysis as a tool to decipher
seasonality and other cyclical dynamics.

## See also

[`trajectoryCyclicalPlots`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclicalPlots.md),
[`cycleShiftArrows`](https://emf-creaf.github.io/ecotraj/reference/cycleShiftArrows.md),
[`trajectoryMetrics`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md),
[`trajectoryComparison`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md),
[`trajectoryConvergencePlot`](https://emf-creaf.github.io/ecotraj/reference/trajectoryConvergencePlot.md)

## Author

Nicolas Djeghri, UBO

Miquel De Cáceres, CREAF

## Examples

``` r
#First build a toy dataset with:
#The sampling times of the time series
timesToy <- 0:30 

#The duration of the cycles (i.e. the periodicity of the time series)
cycleDurationToy <- 10 

#The sites sampled (only one named "A")
sitesToy <- rep(c("A"),length(timesToy)) 

#And prepare a trend term
trend <- 0.05

#Build cyclical data (note that we apply the trend only to x):
x <- sin((timesToy*2*pi)/cycleDurationToy)+trend*timesToy
y <- cos((timesToy*2*pi)/cycleDurationToy)
matToy <- cbind(x,y)

#And express it as distances:
dToy <- dist(matToy)

#Make it an object of class trajectory:
cyclicalTrajToy <- defineTrajectories(d = dToy,
                                      sites = sitesToy,
                                      times = timesToy)

#At this stage, cycles and / or fixed date trajectories are not isolated.
#This done with the two CETA "extract" functions:
cyclesToy <- extractCycles(x = cyclicalTrajToy,
                           cycleDuration = cycleDurationToy)
fdTrajToy <- extractFixedDateTrajectories(x = cyclicalTrajToy,
                                          cycleDuration = cycleDurationToy)

#The output of these functions can be used as input
#for other ETA functions to get metrics of interest
#such as trajectory length:
trajectoryLengths(x = cyclesToy)
#>             S1        S2        S3       S4        S5        S6       S7
#> A_C1 0.6657661 0.6486857 0.6200532 0.590033 0.5706904 0.5706904 0.590033
#> A_C2 0.6657661 0.6486857 0.6200532 0.590033 0.5706904 0.5706904 0.590033
#> A_C3 0.6657661 0.6486857 0.6200532 0.590033 0.5706904 0.5706904 0.590033
#>             S8        S9       S10     Path
#> A_C1 0.6200532 0.6486857 0.6657661 6.190457
#> A_C2 0.6200532 0.6486857 0.6657661 6.190457
#> A_C3 0.6200532 0.6486857 0.6657661 6.190457
trajectoryLengths(x = fdTrajToy)
#>          S1  S2  S3 Path
#> A_fdT_0 0.5 0.5 0.5  1.5
#> A_fdT_1 0.5 0.5  NA  1.0
#> A_fdT_2 0.5 0.5  NA  1.0
#> A_fdT_3 0.5 0.5  NA  1.0
#> A_fdT_4 0.5 0.5  NA  1.0
#> A_fdT_5 0.5 0.5  NA  1.0
#> A_fdT_6 0.5 0.5  NA  1.0
#> A_fdT_7 0.5 0.5  NA  1.0
#> A_fdT_8 0.5 0.5  NA  1.0
#> A_fdT_9 0.5 0.5  NA  1.0

#or distances between trajectories:
trajectoryDistances(x = cyclesToy)
#>           A_C1      A_C2
#> A_C2 0.4602096          
#> A_C3 0.8276663 0.4602096
trajectoryDistances(x = fdTrajToy)
#>           A_fdT_0   A_fdT_1   A_fdT_2   A_fdT_3   A_fdT_4   A_fdT_5   A_fdT_6
#> A_fdT_1 0.3072093                                                            
#> A_fdT_2 0.8665898 0.5780759                                                  
#> A_fdT_3 1.4323048 1.1644260 0.6200532                                        
#> A_fdT_4 1.8510122 1.6249720 1.1457171 0.5618814                              
#> A_fdT_5 2.0155644 1.8313038 1.4389310 0.9356465 0.3826877                    
#> A_fdT_6 1.8510122 1.7685615 1.5710039 1.2479634 0.8255705 0.3826877          
#> A_fdT_7 1.4323048 1.5044289 1.5356710 1.4521130 1.2479634 0.9356465 0.5618814
#> A_fdT_8 0.8665898 1.0704423 1.3521130 1.5356710 1.5710039 1.4389310 1.1457171
#> A_fdT_9 0.3072093 0.5255705 1.0704423 1.5044289 1.7685615 1.8313038 1.6249720
#>           A_fdT_7   A_fdT_8
#> A_fdT_1                    
#> A_fdT_2                    
#> A_fdT_3                    
#> A_fdT_4                    
#> A_fdT_5                    
#> A_fdT_6                    
#> A_fdT_7                    
#> A_fdT_8 0.6200532          
#> A_fdT_9 1.1644260 0.5780759

#In addition CETA adds two additional specific metrics.
#that require the same inputs as function extractCycles():
cycleConvexity(x = cyclicalTrajToy,
               cycleDuration = cycleDurationToy)
#> A_C1 A_C2 A_C3 
#>   NA    1    1 
#The NA with the first cycle, is expected:
#Cycle convexity cannot be computed right at the boundary of the time series
cycleShifts(x = cyclicalTrajToy,
            cycleDuration = cycleDurationToy)
#>    sites dateCS timeCS timeRef timeScale cyclicalShift
#> 1      A      0     20      10        10  0.000000e+00
#> 2      A      1     21      11        10  0.000000e+00
#> 3      A      2     22      12        10 -1.776357e-15
#> 4      A      3     23      13        10  0.000000e+00
#> 5      A      4     24      14        10  0.000000e+00
#> 6      A      5     15       5        10 -1.776357e-15
#> 7      A      5     25       5        20  0.000000e+00
#> 8      A      5     25      15        10  0.000000e+00
#> 9      A      6     16       6        10  0.000000e+00
#> 10     A      7     17       7        10 -8.881784e-16
#> 11     A      8     18       8        10  0.000000e+00
#> 12     A      9     19       9        10  0.000000e+00
#Note that because our cycles are perfectly regular here, the cyclicalShift
#computed are all 0 (or close because of R's computing approximations)

#Subsetting cycles and fixed date trajectories:
subsetTrajectories(cyclesToy,
                   subtrajectory_selection = "A_C1") 
#> $d
#>            1         2         3         4         5         6         7
#> 2  0.6657661                                                            
#> 3  1.2578463 0.6486857                                                  
#> 4  1.7105119 1.2102150 0.6200532                                        
#> 5  1.9731062 1.6249720 1.1486130 0.5900330                              
#> 6  2.0155644 1.8501135 1.5346716 1.0962457 0.5706904                    
#> 7  1.8317650 1.8640587 1.7442756 1.4761035 1.0755705 0.5706904          
#> 8  1.4404147 1.6687507 1.7639284 1.7021130 1.4761035 1.0962457 0.5900330
#> 9  0.8838104 1.2897072 1.6021130 1.7639284 1.7442756 1.5346716 1.1486130
#> 10 0.2354979 0.7755705 1.2897072 1.6687507 1.8640587 1.8501135 1.6249720
#> 11 0.5000000 0.2354979 0.8838104 1.4404147 1.8317650 2.0155644 1.9731062
#>            8         9        10
#> 2                               
#> 3                               
#> 4                               
#> 5                               
#> 6                               
#> 7                               
#> 8                               
#> 9  0.6200532                    
#> 10 1.2102150 0.6486857          
#> 11 1.7105119 1.2578463 0.6657661
#> 
#> $metadata
#>    sites cycles surveys times dates internal
#> 1      A   A_C1       1     0     0     TRUE
#> 2      A   A_C1       2     1     1     TRUE
#> 3      A   A_C1       3     2     2     TRUE
#> 4      A   A_C1       4     3     3     TRUE
#> 5      A   A_C1       5     4     4     TRUE
#> 6      A   A_C1       6     5     5     TRUE
#> 7      A   A_C1       7     6     6     TRUE
#> 8      A   A_C1       8     7     7     TRUE
#> 9      A   A_C1       9     8     8     TRUE
#> 10     A   A_C1      10     9     9     TRUE
#> 11     A   A_C1      11    10     0    FALSE
#> 
#> attr(,"class")
#> [1] "cycles"       "trajectories" "list"        
subsetTrajectories(fdTrajToy,
                   subtrajectory_selection = c("A_fdT_2","A_fdT_4"))
#> $d
#>          1        2        3        4        5
#> 2 1.148613                                    
#> 3 0.500000 1.353729                           
#> 4 1.142821 0.500000 1.148613                  
#> 5 1.000000 1.686966 0.500000 1.353729         
#> 6 1.338943 1.000000 1.142821 0.500000 1.148613
#> 
#> $metadata
#>   sites     fdT surveys times dates
#> 1     A A_fdT_2       1     2     2
#> 2     A A_fdT_4       1     4     4
#> 3     A A_fdT_2       2    12     2
#> 4     A A_fdT_4       2    14     4
#> 5     A A_fdT_2       3    22     2
#> 6     A A_fdT_4       3    24     4
#> 
#> attr(,"class")
#> [1] "fd.trajectories" "trajectories"    "list"           
                
#General metrics describing the geometry of cycles:
cycleMetrics(x = cyclicalTrajToy,
             cycleDuration = cycleDurationToy)
#>   cycle site  n t_start t_end   length mean_speed mean_angle convexity
#> 1  A_C1    A 10       0    10 6.190457  0.6190457         NA        NA
#> 2  A_C2    A 10      10    20 6.190457  0.6190457         36         1
#> 3  A_C3    A 10      20    30 6.190457  0.6190457         36         1
#>   internal_ss internal_variance
#> 1    8.667408         0.9630454
#> 2    8.667408         0.9630454
#> 3    8.667408         0.9630454
             
```
