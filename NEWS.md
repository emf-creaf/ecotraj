# ecotraj 1.1.2
* New functions trajectoryConvergencePlot and cycleShiftArrows by N. Djeghri

# ecotraj 1.1.1
* Trajectory subsetting according to time window limits
* Bug correction in trajectory definition (missing surveys but existing times)
* Bug correction in trajectory asymmetric tests (where done as symmetric previously)

# ecotraj 1.1.0
* New functions is.metric() and is.synchronous().
* Update of function trajectoryInternalVariation() (former trajectoryVariability()).
* New function dynamicVariation().
* New function variationDecomposition() for synchronous trajectories.
* New function interpolateTrajectories() for interpolation to given times.
* Function trajectoryConvergence() modified to allow a global convergence test (type = "multiple").

# ecotraj 1.0.0
* Trajectory data structures have been introduced.
* New functions for cyclical trajectory analysis (Djeghri et al, in prep.).
* New ETA metrics trajectoryVariability() and trajectorySpeeds().
* Permutational test added for trajectoryDirectionality().
* New function trajectoryMetrics() for evaluation of multiple metrics.
* New function trajectoryWindowMetrics() for metrics evaluated over moving windows.
* New function trajectoryShifts() to compare trajectories that are similar in shape but differ in speed.
* New dissimilarity metric for trajectories, taking into account differences in time.

# ecotraj 0.2.1
* Update of function 'centerTrajectories' to add flexibility in trajectory centroid definition.
* New vignette to illustrate centering.
* Improvement of trajectoryProjection() to return relative position within selected segments, also when orthogonal projection does not exist.

# ecotraj 0.2.0
* New function 'smoothTrajectories' added.
* Update of vignette 'Introduction to to Ecological Trajectory Analysis (ETA)'.
* Correction of squared distances to centroid (GitHub issue #3).
* EQA between individual states and dynamic envelopes (experimental).
