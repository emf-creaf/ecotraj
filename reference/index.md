# Package index

## Datasets

Example datasets

- [`avoca`](https://emf-creaf.github.io/ecotraj/reference/avoca.md)
  [`avoca_sites`](https://emf-creaf.github.io/ecotraj/reference/avoca.md)
  [`avoca_strat`](https://emf-creaf.github.io/ecotraj/reference/avoca.md)
  [`avoca_surveys`](https://emf-creaf.github.io/ecotraj/reference/avoca.md)
  : Avoca permanent plot dataset
- [`furseals`](https://emf-creaf.github.io/ecotraj/reference/furseals.md)
  : furseals dataset
- [`pike`](https://emf-creaf.github.io/ecotraj/reference/pike.md) : pike
  dataset
- [`heatmapdata`](https://emf-creaf.github.io/ecotraj/reference/heatmapdata.md)
  : heatmapdata dataset
- [`isoscape`](https://emf-creaf.github.io/ecotraj/reference/isoscape.md)
  : isoscape dataset
- [`glomel`](https://emf-creaf.github.io/ecotraj/reference/glomel.md) :
  Glomel vegetation dataset
- [`glenan`](https://emf-creaf.github.io/ecotraj/reference/glenan.md) :
  Glenan dataset
- [`northseaZoo`](https://emf-creaf.github.io/ecotraj/reference/northseaZoo.md)
  : North Sea zooplankton dataset

## Trajectory data

Set of functions to define and subset trajectories

- [`defineTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/defineTrajectories.md)
  : Trajectory definition
- [`subsetTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/subsetTrajectories.md)
  : Trajectory subsetting

## Trajectory metrics

Set of functions to calculate trajectory metrics

- [`trajectoryLengths()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryLengths2D()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectorySpeeds()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectorySpeeds2D()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryAngles()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryAngles2D()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryDirectionality()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryInternalVariation()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryMetrics()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  [`trajectoryWindowMetrics()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md)
  : Trajectory metrics

## Trajectory comparison

Set of functions to compare trajectories

- [`trajectoryProjection()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryProjection.md)
  : Trajectory projection
- [`segmentDistances()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)
  [`trajectoryDistances()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)
  [`trajectoryConvergence()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)
  [`trajectoryCorrespondence()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)
  [`trajectoryShifts()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)
  : Trajectory comparison
- [`trajectoryRMA()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryRMA.md)
  : Relative Trajectory Movement Assessment (RTMA)
- [`trajectoryRMAPlot()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryRMAPlot.md)
  : Heat map-like plots for Relative Trajectory Movement Assessment
  (RTMA)
- [`dynamicVariation()`](https://emf-creaf.github.io/ecotraj/reference/dynamicVariation.md)
  [`variationDecomposition()`](https://emf-creaf.github.io/ecotraj/reference/dynamicVariation.md)
  : Dynamic variation and variation decomposition

## Cyclical trajectory analysis

Set of functions for cyclical trajectory analysis

- [`extractCycles()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
  [`extractFixedDateTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
  [`cycleConvexity()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
  [`cycleShifts()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
  [`cycleMetrics()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)
  : Functions for Cyclical Ecological Trajectory Analysis

## Trajectory plots

Set of plotting functions

- [`trajectoryPCoA()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryPlot.md)
  [`trajectoryPlot()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryPlot.md)
  : Trajectory plots
- [`cyclePCoA()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclicalPlots.md)
  [`fixedDateTrajectoryPCoA()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclicalPlots.md)
  : Cyclical trajectory plots
- [`trajectoryConvergencePlot()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryConvergencePlot.md)
  : Summary plot for trajectory convergence and divergence
- [`cycleShiftArrows()`](https://emf-creaf.github.io/ecotraj/reference/cycleShiftArrows.md)
  : Displaying cycle shifts

## Ecological Quality Assessment

Set of functions to determine the membership to state/trajectory
envelopes

- [`trajectoryEnvelopeVariability()`](https://emf-creaf.github.io/ecotraj/reference/referenceEnvelopes.md)
  [`stateEnvelopeVariability()`](https://emf-creaf.github.io/ecotraj/reference/referenceEnvelopes.md)
  [`compareToTrajectoryEnvelope()`](https://emf-creaf.github.io/ecotraj/reference/referenceEnvelopes.md)
  [`compareToStateEnvelope()`](https://emf-creaf.github.io/ecotraj/reference/referenceEnvelopes.md)
  : Ecological quality assessment

## Transform functions

Transforming trajectories

- [`smoothTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/transformTrajectories.md)
  [`centerTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/transformTrajectories.md)
  [`interpolateTrajectories()`](https://emf-creaf.github.io/ecotraj/reference/transformTrajectories.md)
  : Transform trajectories

## Miscellaneous

Miscellaneous functions

- [`is.metric()`](https://emf-creaf.github.io/ecotraj/reference/is.metric.md)
  : Metricity
- [`is.synchronous()`](https://emf-creaf.github.io/ecotraj/reference/is.synchronous.md)
  : Synchronicity in trajectory observations
