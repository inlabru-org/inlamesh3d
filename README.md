
<!-- README.md is generated from README.Rmd. Please edit that file -->

# inlamesh3d

<!-- badges: start -->

<!-- badges: end -->

The goal of inlamesh3d is to add 3D mesh functionality to INLA. When the
methods have been tested, they will be incorporated into the main
package instead.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("finnlindgren/inlamesh3d", build_vignettes = TRUE)
```

## Example

You must load `INLA` before loading `inlamesh3d` to get the correct
overrides.

``` r
library(INLA)
#> Loading required package: Matrix
#> Loading required package: sp
#> Loading required package: parallel
#> Loading required package: foreach
#> This is INLA_20.07.27 built 2020-07-30 15:35:07 UTC.
#>  - See www.r-inla.org/contact-us for how to get help.
#>  - To enable PARDISO sparse library; see inla.pardiso()
#>  - Save 196.6Mb of storage running 'inla.prune()'
library(inlamesh3d)
#> 
#> Attaching package: 'inlamesh3d'
#> The following objects are masked from 'package:INLA':
#> 
#>     inla.mesh.fem, inla.spde.make.A, inla.spde2.matern
loc <- matrix(rnorm(12), 4, 3)
tv <- geometry::delaunayn(loc)
mesh <- inla.mesh3d(loc, tv)
str(mesh)
#> List of 4
#>  $ manifold: chr "R3"
#>  $ n       : int 4
#>  $ loc     : num [1:4, 1:3] -0.0805 -0.4462 0.7394 1.1667 0.0331 ...
#>  $ graph   :List of 1
#>   ..$ tv: int [1, 1:4] 1 3 4 2
#>  - attr(*, "class")= chr "inla_mesh_3d"
```

New versions of `inla.spde2.matern`, `inla.spde.make.A` and
`inla.mesh.fem` are now available to use in much the same way as for 2D
meshes.

For a full example, see `vignette("spde3d")`.
