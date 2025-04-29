
# Ecological trajectory analysis <a href="https://emf-creaf.github.io/ecotraj/"><img src="man/figures/logo.png" align="right" height="139" alt="ecotraj website" /></a>

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/ecotraj)](https://cran.r-project.org/package=ecotraj)
[![](https://cranlogs.r-pkg.org/badges/ecotraj)](https://cran.rstudio.com/web/packages/ecotraj/index.html)
[![R-CMD-check](https://github.com/emf-creaf/ecotraj/workflows/R-CMD-check/badge.svg)](https://github.com/emf-creaf/ecotraj/actions)

## Introduction

Package `ecotraj` is a package providing Ecological Trajectory Analysis
(ETA), a framework to assist ecologists in the analysis of temporal
changes in ecosystems defined on a chosen multivariate space.

ETA is related to the following publications:

- Initial presentation of the **Ecological Trajectory Analysis (ETA)**
  framework (De Cáceres et al. 2019).

- ETA extension of **metrics** and **visualisation modes** (Sturbois et
  al. 2021a).

- ETA extension to **cyclical** trajectory data (Djeghri et al. in
  prep.).

- Application to a multivariate space defined by **stable isotope** data
  (Sturbois et al. 2021b).

- Application of ETA framework for the **assessment of ecological
  quality** (Sturbois et al. 2023).

## Package installation and documentation

Package `ecotraj` can be found at
[CRAN](https://cran.r-project.org/package=ecotraj). In addition, the
latest stable `ecotraj` R package can also be installed from GitHub as
follows:

``` r
remotes::install_github("emf-creaf/ecotraj")
```

Although the package comes with function documentation, articles
explaining how to use the package can be found at
<https://emf-creaf.github.io/ecotraj/>.

## Package functions

The following lists the main package functions, grouped by analysis
framework:

**Functions for Ecological Trajectory Analysis (ETA)**

| Name | Short description |
|----|----|
| `defineTrajectories()` | Defines trajectory objects to be analyzed. |
| `subsetTrajectories()` | Subsets trajectories corresponding to selected set of sites or surveys from a dataset containing multiple trajectories. |
| `smoothTrajectories()` | Performs multivariate smoothing on trajectory data using a Gaussian kernel. |
| `centerTrajectories()` | Shifts all trajectories to the center of the multivariate space and returns a modified distance matrix. |
| `trajectoryLengths()` | Calculates lengths of directed segments and total path length of trajectories. |
| `trajectoryLengths2D()` | Calculates lengths of directed segments and total path lengths of trajectories from 2D coordinates given as input. |
| `trajectorySpeeds()` | Calculates speeds of directed segments and total path speed of trajectories. |
| `trajectorySpeeds2D()` | Calculates speeds of directed segments and total path speed of trajectories from 2D coordinates given as input. |
| `trajectoryAngles()` | Calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points. |
| `trajectoryAngles2D()` | Calculates the angle between consecutive pairs of directed segments or between segments of ordered triplets of points from 2D coordinates. |
| `trajectoryConvergence()` | Performs the Mann-Kendall trend test on (1) the distances between trajectories; (2) the distance between points of one trajectory to the other; or (3) the variance of states among trajectories. |
| `trajectoryDirectionality()` | Assesses the level directionality of each trajectory. |
| `trajectoryVariability()` | Assesses the internal variability within each trajectory. |
| `trajectoryDistances()` | Estimates distances between pairs of ecosystem trajectories. |
| `trajectoryShifts()` | Calculates trajectory shifts (i.e. advances and delays) between trajectories assumed to follow a similar path but with different speeds or time lags. |
| `trajectoryPlot()` | Draws trajectories in a 2D-plot corresponding to the input coordinates. |
| `trajectoryPCoA()` | Performs principal coordinates analysis and draws trajectories in the ordination scatterplot. |
| `trajectoryMetrics()` | Calculates several (whole path) metrics on trajectories. |
| `trajectoryWindowMetrics()` | Calculates several metrics on subtrajectories defined using moving windows. |

**Functions specific for Cyclical Ecological Trajectory Analysis
(CETA)**

| Name | Short description |
|----|----|
| `extractCycles()` | Extracts (i.e. subsets) cycles from cyclical trajectories. |
| `extractFixedDateTrajectories()` | Extracts fixed-date trajectories from cyclical trajectories. |
| `cycleShifts()` | Estimates cyclical shifts (i.e. advances or delays). |
| `cycleConvexity()` | Estimates the degree of convexity of cycles. |
| `cycleMetrics()` | Calculates several (whole path) metrics on cycles. |
| `cyclePCoA()` | Performs Principal Coordinates Analysis and draws cycles in the ordination scatter plot. |
| `fixedDateTrajectoryPCoA()` | Performs Principal Coordinates Analysis and draws fixed date trajectories in the ordination scatter plot. |

**Functions for Ecological Quality Assessment (EQA)**

| Name | Short description |
|----|----|
| `stateEnvelopeVariability()` | Assess the variability of a set of *states* taken as reference envelope. |
| `trajectoryEnvelopeVariability()` | Assess the variability of a set of *trajectories* taken as reference envelope. |
| `compareToStateEnvelope()` | Evaluate the ecological quality of stations/observations with respect to a set of *states* taken as reference envelope. |
| `compareToTrajectoryEnvelope()` | Evaluate the ecological quality of stations/observations with respect to a set of *trajectories* taken as reference envelope. |

## Related package

The R package [**ecoregime**](https://mspinillos.github.io/ecoregime/),
by Martina Sánchez-Pinillos, implements the Ecological Dynamic Regimes
(EDR) framework to characterize and compare groups of ecological
trajectories in multidimensional spaces defined by ecosystem state
variables.

## References

- De Cáceres, M., Coll, L., Legendre, P., Allen, R.B., Wiser, S.K.,
  Fortin, M.J., Condit, R. & Hubbell, S. (2019). Trajectory analysis in
  community ecology. Ecological Monographs 89, e01350
  (<https://doi.org/10.1002/ecm.1350>).

- Djeghri et al (in preparation) Going round in cycles, but going
  somewhere: Ecological Trajectory Analysis as a tool to decipher
  seasonality and other cyclical dynamics.

- Sturbois, A., De Cáceres, M., Sánchez-Pinillos, M., Schaal, G.,
  Gauthier, O., Le Mao, P., Ponsero, A., & Desroy, N. (2021a). Extending
  community trajectory analysis : New metrics and representation.
  Ecological Modelling, 440, 109400
  (<https://doi.org/10.1016/j.ecolmodel.2020.109400>).

- Sturbois, A., Cucherousset, J., De Cáceres, M., Desroy, N., Riera, P.,
  Carpentier, A., Quillien, N., Grall, J., Espinasse, B., Cherel, Y.,
  Schaal, G. (2021). Stable Isotope Trajectory Analysis (SITA) : A new
  approach to quantify and visualize dynamics in stable isotope studies.
  Ecological Monographs, 92, e1501 (<https://doi.org/10.1002/ecm.1501>).

- Sturbois, A., De Cáceres, M., Bifolchi, A., Bioret, F., Boyé, A.,
  Gauthier, O., Grall, J., Grémare, A., Labrune, C., Robert, A., Schaal,
  G., Desroy, N. (2023). Ecological Quality Assessment: a general
  multivariate framework to report the quality of ecosystems and their
  dynamics with respect to reference conditions. Ecosphere 14, e4726
  (<https://doi.org/10.1002/ecs2.4726>).
