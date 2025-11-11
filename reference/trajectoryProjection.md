# Trajectory projection

Performs an projection of a set of target points onto a specified
trajectory and returns the distance to the trajectory (i.e. rejection)
and the relative position of the projection point within the trajectory.

## Usage

``` r
trajectoryProjection(
  d,
  target,
  trajectory,
  tol = 1e-06,
  add = TRUE,
  force = TRUE
)
```

## Arguments

- d:

  A symmetric [`matrix`](https://rdrr.io/r/base/matrix.html) or an
  object of class [`dist`](https://rdrr.io/r/stats/dist.html) containing
  the distance values between pairs of ecological states (see details).

- target:

  An integer vector of the ecological states to be projected.

- trajectory:

  An integer vector of the ecological states conforming the trajectory
  onto which target states are to be projected.

- tol:

  Numerical tolerance value to determine that projection of a point lies
  within the trajectory.

- add:

  Flag to indicate that constant values should be added (local
  transformation) to correct triplets of distance values that do not
  fulfill the triangle inequality.

- force:

  Flag to indicate that when projection falls out of the reference
  trajectory for a given, the closest point in the trajectory will be
  used.

## Value

A data frame with the following columns:

- `distanceToTrajectory`: Distances to the trajectory, i.e. rejection.
  If there is no orthogonal projection the distance corresponds to the
  minimum distance to the trajectory.

- `segment`: Segment that includes the projected point or the closest
  state.

- `relativeSegmentPosition`: Relative position of the projected point
  within the segment, i.e. values from 0 to 1 with 0 representing the
  start of the segment and 1 representing its end.

- `relativeTrajectoryPosition`: Relative position of the projected point
  within the trajectory, i.e. values from 0 to 1 with 0 representing the
  start of the trajectory and 1 representing its end.

## Author

Miquel De Cáceres, CREAF
