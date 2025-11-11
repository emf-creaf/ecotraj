# Metricity

Checks whether the input dissimilarity matrix is metric (i.e. all
triplets fulfill the triangle inequality).

## Usage

``` r
is.metric(x, tol = 1e-04)
```

## Arguments

- x:

  Either an object of class `trajectories`, a symmetric
  [`matrix`](https://rdrr.io/r/base/matrix.html) or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the distance
  values between pairs of ecological states.

- tol:

  Tolerance value for metricity

## Value

A boolean indicating metric property

## Author

Miquel De Cáceres, CREAF
