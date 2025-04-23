# ecotraj 1.1.0
* New functions is.metric() and is.synchronous().
* Update of function trajectoryVariability().
* New function variationDecomposition() for synchronous trajectories.

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
