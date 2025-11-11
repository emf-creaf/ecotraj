# North Sea zooplankton dataset

A multi-annual (1958-2021), monthly resolved dataset of zooplankton
community composition in the Northern and Southern North Sea used to
illustrate Cyclical Ecological Trajectory Analysis (CETA)

## Format

northseaZoo is an object of class
[`list`](https://rdrr.io/r/base/list.html) composed of 3 objects:

- Hellinger:

  a [`data.frame`](https://rdrr.io/r/base/data.frame.html) containing
  Hellinger-transformed zooplankton taxa abundances.

- times:

  a vector indicating the date (in year) associated to each line in
  `Hellinger`.

- sites:

  a vector indicating the site (`"NNS"` = Northern North Sea, `"SNS"` =
  Southern North Sea) associated to each line in `Hellinger`.

## Details

The data describes the zooplankton community in the North Sea sampled by
the [Continuous Plankton Recorder (CPR)](https://www.cprsurvey.org/)
survey. The CPR survey operates through towing of CPR samplers across
commercial routes of merchant ships (plankton silk mesh = 270 microm,
sampling depth = 5-10 m). When brought back to the laboratory, plankton
is counted and identified taxonomically following standardized
protocols. The raw data provided by the survey
([doi:10.17031/66f12be296d70](https://doi.org/10.17031/66f12be296d70) ).
was reformated into two monthly-resolved time series of the commonest
zooplankton taxa in the Northern North Sea (`"NNS"`) and the Southern
North Sea (`"SNS"`). During data processing, a smoothing was performed
by taking a rolling average (for each month, 5 values were averaged: a 3
months window + the corresponding month of the previous and next years).
The abundances were finally Hellinger-transformed, making them amenable
to ecological diversity study.

## See also

[`trajectoryCyclical`](https://emf-creaf.github.io/ecotraj/reference/trajectoryCyclical.md)

## Author

Nicolas Djeghri, Université de Bretagne Occidentale, France

Pierre Hélaouët and CPR survey staff, Marine Biological Association,
United Kingdom
