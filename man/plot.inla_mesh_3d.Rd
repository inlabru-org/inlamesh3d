% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mesh3d.R
\name{plot.inla_mesh_3d}
\alias{plot.inla_mesh_3d}
\title{Draw 3D mesh}
\usage{
\method{plot}{inla_mesh_3d}(
  x,
  include = c(FALSE, TRUE, TRUE),
  alpha = c(0.9, 0.3, 0.1),
  col = c("black", "blue", "red"),
  t_sub = NULL,
  size = 5,
  lwd = 2,
  fill_col = NULL,
  add = FALSE,
  ...
)
}
\arguments{
\item{x}{A \code{inla_mesh_3d} object}

\item{include}{Whether to draw vertices, edges, and triangles, Default: c(FALSE, TRUE, TRUE)}

\item{alpha}{Opaqueness of vertices, edges, and triangles, Default: c(0.9, 0.3, 0.1)}

\item{col}{Colour of vertices, edges, and triangles, Default: \code{c("black", "blue", "red")}}

\item{t_sub}{A subset of tetrahedron indices to include, Default: NULL
(draw all tetrahedra)}

\item{size}{Vertex size, in pixels}

\item{lwd}{Edge width, in pixels}

\item{add}{If TRUE, adds to an existing \code{rgl} device, Default: FALSE}

\item{...}{Further parameters passed through to the rgl plotting functions}
}
\description{
Draw a 3D mesh with rgl
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
