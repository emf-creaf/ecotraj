# Transform trajectories

The following functions are provided to transform trajectories:

- Function `smoothTrajectories` performs multivariate smoothing on
  trajectory data using a Gaussian kernel.

- Function `centerTrajectories` shifts all trajectories to the center of
  the multivariate space and returns a modified distance matrix.

- Function `averageTrajectories` creates an "average" trajectory where
  the position of each observation is the average of the position of the
  corresponding observations of the input trajectories.

- Function `interpolateTrajectories` relocates trajectory ecological
  states to those corresponding to input times, via interpolation.

## Usage

``` r
smoothTrajectories(
  x,
  survey_times = NULL,
  kernel_scale = 1,
  fixed_endpoints = TRUE
)

centerTrajectories(x, exclude = integer(0))

averageTrajectories(
  x,
  group = NULL,
  keep_members = FALSE,
  output_name = "average"
)

interpolateTrajectories(x, times)
```

## Arguments

- x:

  An object of class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  (or of a sub-class such as
  [`cycles`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)).
  Function `averageTrajectories` requires synchronous input
  trajectories.

- survey_times:

  A vector indicating the survey time for all surveys (if `NULL`, time
  between consecutive surveys is considered to be one)

- kernel_scale:

  Scale of the Gaussian kernel, related to survey times

- fixed_endpoints:

  A logical flag to force keeping the location of trajectory endpoints
  unmodified

- exclude:

  An integer vector indicating observations that are excluded from
  trajectory centroid computation. Note: for objects of class
  [`cycles`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md),
  `external` are excluded by default.

- group:

  Character vector of the sites (trajectories) to be averaged. If `NULL`
  all trajectories contribute to the average

- keep_members:

  Boolean flag to keep the group member trajectories in the result

- output_name:

  A string with the name for the average trajectory

- times:

  A numeric vector indicating new observation times for trajectories.
  Values should be comprised between time limits of the original
  trajectories.

## Value

A modified object of class
[`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md),
where distance matrix has been transformed. When calling
`interpolateTrajectories` and `averageTrajectories`, also the number of
observations and metadata is likely to be affected.

## Details

We recommend reading the article "Transforming trajectories" on the
package website prior to use these functions.

Details of calculations for trajectory centering are given in De Cáceres
et al (2019). Functions `centerTrajectories` and `averageTrajectories`
perform centering/averaging of trajectories using matrix algebra as
explained in Anderson (2017).

When using transformation functions on objects of class
[`cycles`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
or
[`fd.trajectories`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md),
the corresponding transformations are applied to trajectory subsections
(e.g. cycles) instead of being applied to the whole trajectory.

## References

De Cáceres M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit
R & Hubbell S. (2019). Trajectory analysis in community ecology.
Ecological Monographs 89, e01350.

Anderson (2017). Permutational Multivariate Analysis of Variance
(PERMANOVA). Wiley StatsRef: Statistics Reference Online. 1-15. Article
ID: stat07841.

## See also

[`trajectoryPlot`](https://emf-creaf.github.io/ecotraj/reference/trajectoryPlot.md)
[`trajectoryMetrics`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
[`trajectoryComparison`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)

## Author

Miquel De Cáceres, CREAF

Nicolas Djeghri, UBO
