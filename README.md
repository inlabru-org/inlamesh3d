
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inlamesh3d

<!-- badges: start -->
<!-- badges: end -->

The goal of `inlamesh3d` is to add 3D mesh functionality to fmesher and
INLA. When the methods have been tested, they will be incorporated into
the main packages instead.

For online documentation, see
<https://inlabru-org.github.io/inlamesh3d/>

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("inlabru-org/inlamesh3d", build_vignettes = TRUE)
```

## Example

You must load `INLA` before loading `inlamesh3d` to get the correct
overrides.

``` r
suppressPackageStartupMessages(library(INLA))
library(inlamesh3d)
#> Loading required package: fmesher
#> 
#> Attaching package: 'inlamesh3d'
#> The following objects are masked from 'package:INLA':
#> 
#>     inla.spde2.matern, inla.spde2.pcmatern
loc <- matrix(rnorm(12), 4, 3)
tv <- geometry::delaunayn(loc)
mesh <- inla.mesh3d(loc, tv)
str(mesh)
#> List of 4
#>  $ manifold: chr "R3"
#>  $ n       : int 4
#>  $ loc     : num [1:4, 1:3] -0.7103 -0.7264 -1.0812 1.5949 -0.0406 ...
#>  $ graph   :List of 1
#>   ..$ tv: int [1, 1:4] 1 2 4 3
#>  - attr(*, "class")= chr [1:2] "inla_mesh_3d" "inla.mesh"
```

New version of `inla.spde2.matern`, as well as `fm_basis()`,
`fm_evaluator()`, and `fm_fem()` methods are now available to use in
much the same way as for 2D meshes.

For a full example, see `vignette("spde3d")`.
