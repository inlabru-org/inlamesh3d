
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
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
#> This is INLA_99.99.9999 built 2023-06-23 15:51:37 UTC.
#>  - See www.r-inla.org/contact-us for how to get help.
library(inlamesh3d)
#> 
#> Attaching package: 'inlamesh3d'
#> The following objects are masked from 'package:INLA':
#> 
#>     inla.mesh.fem, inla.spde.make.A, inla.spde2.matern,
#>     inla.spde2.pcmatern
loc <- matrix(rnorm(12), 4, 3)
tv <- geometry::delaunayn(loc)
mesh <- inla.mesh3d(loc, tv)
str(mesh)
#> List of 4
#>  $ manifold: chr "R3"
#>  $ n       : int 4
#>  $ loc     : num [1:4, 1:3] -0.217 1.68 -0.126 -0.822 -1.091 ...
#>  $ graph   :List of 1
#>   ..$ tv: int [1, 1:4] 3 1 2 4
#>  - attr(*, "class")= chr "inla_mesh_3d"
```

New versions of `inla.spde2.matern`, `inla.spde.make.A` and
`inla.mesh.fem` are now available to use in much the same way as for 2D
meshes.

For a full example, see `vignette("spde3d")`.
