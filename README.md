ecotraj
================

## Introduction

Package `ecotraj` is a package providing a framework to assist
ecologists in the analysis of temporal changes in ecosystems defined on
a chosen multivariate space. It is related to the following
publications:

-   First presentation of the community trajectory analysis framework:
    De Cáceres et al. (2019) (<https://doi.org/10.1002/ecm.1350>).

-   Extension of community trajectory analysis: Sturbois et al. (2021)
    (<https://doi.org/10.1016/j.ecolmodel.2020.109400>).

## Package installation

Package `ecotraj` can be found at [CRAN](https://cran.r-project.org/).
In addition, the latest stable `ecotraj` R package can also be installed
from GitHub as follows:

``` r
devtools::install_github("emf-creaf/ecotraj")
```

Additionally, users can have help to run package functions directly as
package vignettes, by forcing their inclusion in installation:

``` r
devtools::install_github("emf-creaf/ecotraj", 
                         build_opts = c("--no-resave-data", "--no-manual"),
                         build_vignettes = TRUE)
```

## References

-   De Cáceres, M., Coll, L., Legendre, P., Allen, R.B., Wiser, S.K.,
    Fortin, M.J., Condit, R. & Hubbell, S.. (2019). Trajectory analysis
    in community ecology. Ecological Monographs.

-   Sturbois, A., De Cáceres, M., Sánchez-Pinillos, M., Schaal, G.,
    Gauthier, O., Le Mao, P., Ponsero, A., & Desroy, N. (2021).
    Extending community trajectory analysis : New metrics and
    representation. Ecological Modelling, 440, 109400.
