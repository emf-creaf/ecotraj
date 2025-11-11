# Ecological quality assessment

Functions to assess the variability of ecological reference envelopes
and to assess the ecological quality of target stations/observations
with respect to reference envelopes (Sturbois et al., under review).

## Usage

``` r
trajectoryEnvelopeVariability(
  d,
  sites,
  surveys = NULL,
  envelope = NULL,
  nboot.ci = NULL,
  alpha.ci = 0.05,
  ...
)

stateEnvelopeVariability(d, envelope = NULL, nboot.ci = NULL, alpha.ci = 0.05)

compareToTrajectoryEnvelope(
  d,
  sites,
  envelope,
  surveys = NULL,
  m = 1.5,
  comparison_target = "trajectories",
  distances_to_envelope = FALSE,
  distance_percentiles = FALSE,
  ...
)

compareToStateEnvelope(
  d,
  envelope,
  m = 1.5,
  nboot.ci = NULL,
  alpha.ci = 0.05,
  distances_to_envelope = FALSE,
  distance_percentiles = FALSE,
  ...
)
```

## Arguments

- d:

  A symmetric [`matrix`](https://rdrr.io/r/base/matrix.html) or an
  object of class [`dist`](https://rdrr.io/r/stats/dist.html) containing
  the distance values between pairs of ecological states (see details).

- sites:

  A vector indicating the site corresponding to each ecological state.

- surveys:

  A vector indicating the survey corresponding to each ecological state
  (only necessary when surveys are not in order).

- envelope:

  A vector indicating the set of sites that conform the reference
  envelope (other sites will be compared to the envelope)

- nboot.ci:

  Number of bootstrap samples for confidence intervals. If nboot.ci =
  NULL then confidence intervals are not estimated.

- alpha.ci:

  Error in confidence intervals.

- ...:

  Additional parameters for function
  [`trajectoryDistances`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.md)

- m:

  Fuzziness exponent for quality value assessment

- comparison_target:

  String indicating the component to be compared to the reference
  envelope. Either 'trajectories' (to compare complete trajectories) or
  'states' (to compare individual trajectory states).

- distances_to_envelope:

  Flag to indicate that distances to envelope should be included in the
  result

- distance_percentiles:

  Flag to include the percentage of distances to the envelope (among
  sites corresponding to the reference) that are smaller than that of
  the site.

## Value

- Functions `stateEnvelopeVariability` and
  `trajectoryEnvelopeVariability` are used to assess the variability of
  reference envelopes.

- Functions `compareToStateEnvelope` and `compareToTrajectoryEnvelope`
  return data frame with columns identifying the envelope and the Q
  statistic for the ecological quality with respect to the envelope. If
  `nboot.ci != NULL` extra columns are added to indicate the boundaries
  of a confidence interval for Q, built using bootstrap samples of the
  reference envelope.

## Details

Functions `stateEnvelopeVariability` and `trajectoryEnvelopeVariability`
are used to assess the variability of reference envelopes. Functions
`compareToStateEnvelope` and `compareToTrajectoryEnvelope` are used to
evaluate the ecological quality of stations/observations with respect to
a predefined reference envelope.

## References

Sturbois, A., De Cáceres, M., Bifolchi, A., Bioret, F., Boyé, A.,
Gauthier, O., Grall, J., Grémare, A., Labrune, C., Robert, A., Schaal,
G., Desroy, N. (2023). Ecological Quality Assessment: a general
multivariate framework to report the quality of ecosystems and their
dynamics with respect to reference conditions. Ecosphere.

## See also

[`trajectoryMetrics`](https://emf-creaf.github.io/ecotraj/reference/trajectoryMetrics.md),
[`glomel`](https://emf-creaf.github.io/ecotraj/reference/glomel.md)

## Author

Miquel De Cáceres, CREAF

Anthony Sturbois, Vivarmor nature, Réserve Naturelle nationale de la
Baie de Saint-Brieuc

## Examples

``` r
data(glomel)
 
# Extract compositional data matrix
glomel_comp <- as.matrix(glomel[,!(names(glomel) %in% c("ID", "Ref", "Complementary"))])
rownames(glomel_comp) <- glomel$ID
 
# Calculate Bray-Curtis distance matrix 
glomel_bc <- vegan::vegdist(glomel_comp, method = "bray")
 
# Define reference envelope (5 stations) by observation ID
glomel_env <- glomel$ID[glomel$Ref]
 
# Assess quality with respect to reference envelope
compareToStateEnvelope(glomel_bc, glomel_env)
#>    Observation Envelope          Q
#> 1           1A     TRUE 1.00000000
#> 2           1B     TRUE 0.81060634
#> 3           1C     TRUE 0.56214189
#> 4           5A     TRUE 1.00000000
#> 5           13     TRUE 0.88208650
#> 6            2    FALSE 0.57224382
#> 7           3A    FALSE 0.39919080
#> 8           3B    FALSE 1.00000000
#> 9           3C    FALSE 0.42318855
#> 10          5B    FALSE 0.28846018
#> 11          6A    FALSE 0.04261028
#> 12          6C    FALSE 0.34551240
#> 13          6D    FALSE 0.02853182
#> 14          7A    FALSE 0.04810891
#> 15          8A    FALSE 0.15160202
#> 16          8B    FALSE 0.11181164
#> 17          8C    FALSE 0.02057370
#> 18           9    FALSE 0.09946578
#> 19          10    FALSE 0.81557507
#> 20          11    FALSE 0.56752876
#> 21          12    FALSE 1.00000000
#> 22         14A    FALSE 0.66577245
#> 23         14B    FALSE 0.09988729
```
