# Transform trajectories

The following functions are provided to transform trajectories:

- Function `smoothTrajectories` performs multivariate smoothing on
  trajectory data using a Gaussian kernel.

- Function `centerTrajectories` shifts all trajectories to the center of
  the multivariate space and returns a modified distance matrix.

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

interpolateTrajectories(x, times)
```

## Arguments

- x:

  An object of class
  [`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md).

- survey_times:

  A vector indicating the survey time for all surveys (if `NULL`, time
  between consecutive surveys is considered to be one)

- kernel_scale:

  Scale of the Gaussian kernel, related to survey times

- fixed_endpoints:

  A logical flag to force keeping the location of trajectory endpoints
  unmodified

- exclude:

  An integer vector indicating sites that are excluded from trajectory
  centroid computation. Note: for objects of class
  [`cycles`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md),
  `external` are excluded by default.

- times:

  A numeric vector indicating new observation times for trajectories.
  Values should be comprised between time limits of the original
  trajectories.

## Value

A modified object of class
[`trajectories`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md),
where distance matrix has been transformed. When calling
`interpolateTrajectories`, also the number of observations and metadata
is likely to be affected.

## Details

Details of calculations are given in De Cáceres et al (2019). Function
`centerTrajectories` performs centering of trajectories using matrix
algebra as explained in Anderson (2017).

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

## Author

Miquel De Cáceres, CREAF

Nicolas Djeghri, UBO
