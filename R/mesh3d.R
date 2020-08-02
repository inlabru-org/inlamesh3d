#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param loc PARAM_DESCRIPTION
#' @param tv PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.mesh3d
#' @aliases inla_mesh_3d
#' @export

inla.mesh3d <- function(loc, tv) {
  stopifnot(ncol(loc) == 3)
  stopifnot(ncol(tv) == 4)
  stopifnot(max(as.vector(tv)) <= nrow(loc))
  mesh <- list(
    manifold = "R3",
    n = nrow(loc),
    loc = loc,
    graph = list(tv = tv)
  )
  class(mesh) <- "inla_mesh_3d"
  mesh
}


# FEM ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.mesh.fem
#' @export

inla.mesh.fem <- function(mesh,
                          ...) {
  UseMethod("inla.mesh.fem", mesh)
}

#' @export
#' @rdname inla.mesh.fem

inla.mesh.fem.inla_mesh_3d <- function(mesh, ...) {
  inla.mesh3d.fem(mesh, ...)
}

#' @export
#' @rdname inla.mesh.fem

inla.mesh.fem.default <- function(mesh, ...) {
  INLA::inla.mesh.fem(mesh, ...)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param e1 PARAM_DESCRIPTION
#' @param e2 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname row_cross_product
#' @export

row_cross_product <- function(e1, e2) {
  cbind(
    e1[, 2] * e2[, 3] - e1[, 3] * e2[, 2],
    e1[, 3] * e2[, 1] - e1[, 1] * e2[, 3],
    e1[, 1] * e2[, 2] - e1[, 2] * e2[, 1]
  )
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param e1 PARAM_DESCRIPTION
#' @param e2 PARAM_DESCRIPTION
#' @param e3 PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname row_volume_product
#' @export

row_volume_product <- function(e1, e2, e3) {
  rowSums(row_cross_product(e1, e2) * e3)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.mesh.fem
#' @export
#' @importFrom Matrix sparseMatrix diag
inla.mesh3d.fem <- function(mesh, ...) {
  v1 <- mesh$loc[mesh$graph$tv[, 1], , drop = FALSE]
  v2 <- mesh$loc[mesh$graph$tv[, 2], , drop = FALSE]
  v3 <- mesh$loc[mesh$graph$tv[, 3], , drop = FALSE]
  v4 <- mesh$loc[mesh$graph$tv[, 4], , drop = FALSE]
  e1 <- v2 - v1
  e2 <- v3 - v2
  e3 <- v4 - v3
  e4 <- v1 - v4
  vols_t <- abs(row_volume_product(e1, e2, e3)) / 6

  c0 <- Matrix::sparseMatrix(
    i = as.vector(mesh$graph$tv),
    j = as.vector(mesh$graph$tv),
    x = rep(vols_t / 4, times = 4),
    dims = c(mesh$n, mesh$n)
  )
  vols_v <- Matrix::diag(c0)

  # Sign changes for b2 and b4 for consistent in/out vector orientation
  b1 <- row_cross_product(e2, e3)
  b2 <- -row_cross_product(e3, e4)
  b3 <- row_cross_product(e4, e1)
  b4 <- -row_cross_product(e1, e2)

  g_i <- g_j <- g_x <- c()
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    GG <- rbind(
      b1[tt, , drop = FALSE],
      b2[tt, , drop = FALSE],
      b3[tt, , drop = FALSE],
      b4[tt, , drop = FALSE]
    )
    g_i <- c(g_i, rep(mesh$graph$tv[tt, ], each = 4))
    g_j <- c(g_j, rep(mesh$graph$tv[tt, ], times = 4))
    g_x <- c(g_x, as.vector((GG %*% t(GG)) / vols_t[tt] / 36))
  }
  g1 <- Matrix::sparseMatrix(
    i = g_i,
    j = g_j,
    x = g_x,
    dims = c(mesh$n, mesh$n)
  )

  list(c0 = c0, g1 = g1, g2 = g1 %*% (g1 / vols_v), va = vols_v, vt = vols_t)
}


# Barycentric coordinates ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.spde.make.A
#' @export

inla.spde.make.A <- function(mesh,
                              ...) {
  UseMethod("inla.spde.make.A", mesh)
}

#' @export
#' @rdname inla.spde.make.A

inla.spde.make.A.inla_mesh_3d <- function(mesh, ...) {
  inla.mesh3d.make.A(mesh, ...)
}

#' @export
#' @rdname inla.spde.make.A

inla.spde.make.A.default <- function(mesh, ...) {
  INLA::inla.spde.make.A(mesh, ...)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @param loc PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.mesh3d.bary
#' @export

inla.mesh3d.bary <- function(mesh, loc) {
  stopifnot(inherits(mesh, "inla_mesh_3d"))
  loc_1 <- rbind(t(loc), 1)
  bary <- matrix(-Inf, nrow(loc), 4)
  vt <- integer(nrow(loc))
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    loc_t <- mesh$loc[mesh$graph$tv[tt, ], , drop = FALSE]
    # Barycentric coordinates fulfil
    # 1) rbind(t(loc), 1) = rbind(t(loc_t), 1) * w
    # 2) w >= 0
    w <- solve(rbind(t(loc_t), 1), loc_1)
    idx <- which(apply(w, 2, min) > apply(bary, 1, min))
    if (length(idx) > 0) {
      bary[idx, ] <- t(w[, idx, drop = FALSE])
      vt[idx] <- tt
    }
  }
  list(bary = bary, vt = vt)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @param loc PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[Matrix]{sparseMatrix}}
#' @rdname inla.mesh3d.make.A
#' @export
#' @importFrom Matrix sparseMatrix
inla.mesh3d.make.A <- function(mesh, loc, ...) {
  stopifnot(inherits(mesh, "inla_mesh_3d"))
  bary <- inla.mesh3d.bary(mesh, loc)
  # Accept 0.1% relative distance outside the closest tetrahedron
  ok <- apply(bary$bary, 1, min) > -1e-3
  ii <- which(ok)
  if (length(ii) < nrow(loc)) {
    warning(paste0(
      "Some locations were not found inside the mesh: \n(",
      paste0(which(!ok), collapse = ", "),
      ")"
    ))
  }
  A <- Matrix::sparseMatrix(
    dims = c(nrow(loc), mesh$n),
    i = rep(ii, 4),
    j = as.vector(mesh$graph$tv[bary$vt[ii], ]),
    x = as.vector(bary$bary[ii, ])
  )
  A
}


# SPDE object ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname inla.spde2.matern
#' @export

inla.spde2.matern <- function(mesh,
                          ...) {
  UseMethod("inla.spde2.matern", mesh)
}

#' @export
#' @rdname inla.spde2.matern

inla.spde2.matern.inla_mesh_3d <- function(mesh, ...) {
  inla.spde2.matern3d(mesh, ...)
}

#' @export
#' @rdname inla.spde2.matern

inla.spde2.matern.default <- function(mesh, ...) {
  INLA::inla.spde2.matern(mesh, ...)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @param param PARAM_DESCRIPTION, Default: NULL
#' @param constr PARAM_DESCRIPTION, Default: FALSE
#' @param extraconstr.int PARAM_DESCRIPTION, Default: NULL
#' @param extraconstr PARAM_DESCRIPTION, Default: NULL
#' @param fractional.method PARAM_DESCRIPTION, Default: c("parsimonious", "null")
#' @param n.iid.group PARAM_DESCRIPTION, Default: 1
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link{param2.matern3d}},\code{\link{param2.matern_calc}}
#' @rdname inla.spde2.matern
#' @keywords Internal
#' @importFrom Matrix Diagonal
inla.spde2.matern3d <-
  function(mesh,
           param = NULL,
           constr = FALSE,
           extraconstr.int = NULL,
           extraconstr = NULL,
           fractional.method = c("parsimonious", "null"),
           n.iid.group = 1,
           ...) {
    INLA:::inla.require.inherits(mesh, "inla_mesh_3d", "'mesh'")
    fractional.method <- match.arg(fractional.method)

    if (is.null(param)) {
      stop("'param' must be provided. See ?param2.matern3d")
    } else {
      deprecated <-
        c(
          "B.tau", "B.kappa",
          "prior.variance.nominal",
          "prior.range.nominal",
          "prior.tau", "prior.kappa",
          "theta.prior.mean", "theta.prior.prec"
        )
      deprecated <- deprecated[deprecated %in% names(list(...))]
      if (length(deprecated) > 0) {
        warning(paste("'param' specified;  ",
          "Ignoring deprecated parameter(s) ",
          paste(deprecated, collapse = ", "), ".",
          sep = ""
        ))
      }
    }

    alpha <- param$alpha
    n.spde <- mesh$n
    fem <- inla.mesh3d.fem(mesh)

    if (alpha == 2) {
      B.phi0 <- param$B.tau
      B.phi1 <- 2 * param$B.kappa
      M0 <- fem$c0
      M1 <- fem$g1
      M2 <- fem$g2
    } else if (alpha == 1) {
      B.phi0 <- param$B.tau
      B.phi1 <- param$B.kappa
      M0 <- fem$c0
      M1 <- fem$g1 * 0
      M2 <- fem$g1
    } else if (!param$is.stationary) {
      stop("Non-stationary Matern with fractional alpha is not implemented.")
    } else if ((alpha < 2) && (alpha > 1)) {
      if (fractional.method == "parsimonious") {
        lambda <- alpha - floor(alpha)
        b <- matrix(c(1, 0, 0, 1, 1, 0, 1, 2, 1), 3, 3) %*%
          solve(
            matrix(1 / (c(4:2, 3:1, 2:0) + lambda), 3, 3),
            1 / (c(4:2) + lambda - alpha)
          )
      } else if (fractional.method == "null") {
        b <- c(1, alpha, alpha * (alpha - 1) / 2)
      } else {
        stop(paste("Unknown fractional.method '", fractional.method,
          "'.",
          sep = ""
        ))
      }
      B.phi0 <- param$B.tau + (alpha - 2) * param$B.kappa
      B.phi1 <- 2 * param$B.kappa
      M0 <- fem$c0 * b[1]
      M1 <- fem$g1 * b[2] / 2
      M2 <- fem$g2 * b[3]
    } else if ((alpha < 1) && (alpha > 0)) {
      if (fractional.method == "parsimonious") {
        lambda <- alpha - floor(alpha)
        b <- matrix(c(1, 0, 1, 1), 2, 2) %*%
          solve(
            matrix(1 / (c(2:1, 1:0) + lambda), 2, 2),
            1 / (c(2:1) + lambda - alpha)
          )
      } else if (fractional.method == "null") {
        b <- c(1, alpha)
      } else {
        stop(paste("Unknown fractional.method '", fractional.method,
          "'.",
          sep = ""
        ))
      }
      B.phi0 <- param$B.tau + (alpha - 1) * param$B.kappa
      B.phi1 <- param$B.kappa
      M0 <- fem$c0 * b[1]
      M1 <- fem$g1 * 0
      M2 <- fem$g1 * b[2]
    } else {
      stop(paste("Unsupported alpha value (", alpha,
        "). Supported values are 0 < alpha <= 2",
        sep = ""
      ))
    }

    if (n.iid.group == 1) {
      spde <-
        INLA::inla.spde2.generic(
          M0 = M0, M1 = M1, M2 = M2,
          B0 = B.phi0, B1 = B.phi1, B2 = 1,
          theta.mu = param$theta.prior.mean,
          theta.Q = param$theta.prior.prec,
          transform = "identity",
          BLC = param$BLC
        )
    } else {
      if (nrow(B.phi0) > 1) {
        B.phi0 <- kronecker(matrix(1, n.iid.group, 1), B.phi0)
      }
      if (nrow(B.phi1) > 1) {
        B.phi1 <- kronecker(matrix(1, n.iid.group, 1), B.phi1)
      }
      spde <-
        INLA::inla.spde2.generic(
          M0 = kronecker(Matrix::Diagonal(n.iid.group), M0),
          M1 = kronecker(Matrix::Diagonal(n.iid.group), M1),
          M2 = kronecker(Matrix::Diagonal(n.iid.group), M2),
          B0 = B.phi0, B1 = B.phi1, B2 = 1,
          theta.mu = param$theta.prior.mean,
          theta.Q = param$theta.prior.prec,
          transform = "identity",
          BLC = param$BLC
        )
    }
    spde$model <- "matern"
    spde$BLC <- param$BLC

    if (constr || !is.null(extraconstr.int) || !is.null(extraconstr)) {
      A.constr <- matrix(numeric(0), 0, n.spde * n.iid.group)
      e.constr <- matrix(numeric(0), 0, 1)
      if (constr) {
        A.constr <- rbind(
          A.constr,
          matrix(
            fem$va / n.iid.group,
            1, n.spde * n.iid.group
          )
        )
        e.constr <- rbind(e.constr, 0)
      }
      if (!is.null(extraconstr.int)) {
        if (ncol(extraconstr.int$A) == n.spde) {
          A.constr <-
            rbind(
              A.constr,
              kronecker(
                matrix(1 / n.iid.group, 1, n.iid.group),
                as.matrix(extraconstr.int$A %*% fem$c0)
              )
            )
        } else {
          A.constr <-
            rbind(
              A.constr,
              as.matrix(extraconstr.int$A %*%
                kronecker(
                  Matrix::Diagonal(n.iid.group),
                  fem$c0
                ))
            )
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr.int$e))
      }
      if (!is.null(extraconstr)) {
        if (ncol(extraconstr$A) == n.spde) {
          A.constr <-
            rbind(
              A.constr,
              kronecker(
                matrix(1 / n.iid.group, 1, n.iid.group),
                as.matrix(extraconstr$A)
              )
            )
        } else {
          A.constr <- rbind(A.constr, as.matrix(extraconstr$A))
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr$e))
      }

      spde$f$constr <- FALSE
      spde$f$extraconstr <- list(A = A.constr, e = e.constr)
    }

    ## Attach the mesh, so downstream code can have access
    spde$mesh <- mesh

    return(invisible(spde))
  }


# Parameter construction ####

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param mesh PARAM_DESCRIPTION
#' @inheritDotParams param2.matern_calc -dim -dof
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname param2.matern
#' @export

param2.matern <- function(mesh,
                          ...) {
  UseMethod("param2.matern", mesh)
}

#' @export
#' @rdname param2.matern

param2.matern.inla_mesh_3d <- function(mesh, ...) {
  param2.matern_calc(dim = 3, dof = mesh$n, ...)
}

#' @export
#' @rdname param2.matern

param2.matern.inla.mesh <- function(mesh, ...) {
  param2.matern_calc(dim = 2, dof = mesh$n, ...)
}

#' @export
#' @rdname param2.matern

param2.matern.inla.mesh.1d <- function(mesh, ...) {
  param2.matern_calc(dim = 1, dof = mesh$m, ...)
}

B_range_sigma_to_tau_kappa <- function(B_range, B_sigma, d, nu) {
  alpha <- nu + d/2
  log_kappa0 <- log(8 * nu) / 2 - B_range[, 1]
  log_tau0 <- (lgamma(nu) - lgamma(alpha) - d / 2 * log(4 * pi)) / 2 -
    B_sigma[, 1] - nu * log_kappa0
  B_tau <- cbind(log_tau0,
                 -B_sigma[, -1, drop = FALSE] +
                   nu * B_range[, -1, drop = FALSE])
  B_kappa <- cbind(log_kappa0,
                   -B_range[, -1, drop = FALSE])
  list(B_tau = B_tau, B_kappa = B_kappa)
}

construct_prior <- function(B_range, B_sigma,
                            prior_range, prior_sigma,
                            d, nu,
                            prior_theta = NULL) {
  n_theta <- ncol(B_range) - 1L
  if (!is.null(prior_theta)) {
    stopifnot(!is.null(prior_theta[["mean"]]))
    stopifnot(!is.null(prior_theta[["prec"]]))
  } else {
    if (is.null(prior_range)) {
      stop("'prior_range' or 'prior_theta' must be provided.")
    }
    if (is.null(prior_sigma)) {
      stop("'prior_sigma' or 'prior_theta' must be provided.")
    }

    if (n_theta > 0) {
      prior_theta <- list()
      # B_range[, -1, drop = FALSE] %*% mean on average mean_lrange - B_range[, 1]
      # B_sigma[, -1, drop = FALSE] %*% mean on average mean_lsigma - B_sigma[, 1]
      prior_theta$mean <-
        qr.solve(rbind(B_range[, -1, drop = FALSE],
                       B_sigma[, -1, drop = FALSE]),
                 c(log(prior_range[[1]]) - B_range[, 1],
                   log(prior_sigma[[1]]) - B_sigma[, 1]))
      # B_range[, -1, drop = FALSE]^2 %*% variance on average sd_lrange^2
      # B_sigma[, -1, drop = FALSE]^2 %*% variance on average sd_lsigma^2
      # prec = 1 / variance
      sd_lrange <- log(prior_range[[2]]) / qnorm(0.99)
      sd_lsigma <- log(prior_sigma[[2]]) / qnorm(0.99)
      prior_theta$prec <-
        1 / qr.solve(rbind(B_range[, -1, drop = FALSE]^2,
                           B_sigma[, -1, drop = FALSE]^2),
                     c(sd_lrange^2 - rep(0, nrow(B_range)),
                       sd_lsigma^2 - rep(0, nrow(B_sigma))))
    } else {
      prior_theta$mean <- rep(0, n_theta) ## Empty vector
      prior_theta$prec <- matrix(0, n_theta, n_theta) ## Empty matrix
    }
  }

  # Ensure 'prec' is a matrix
  prior_theta$prec <- as.matrix(prior_theta$prec)
  if (ncol(prior_theta$prec) == 1) {
    prior_theta$prec <-
      diag(as.vector(prior_theta$prec), n_theta, n_theta)
  } else if ((nrow(prior_theta$prec) != n_theta) ||
             (ncol(prior_theta$prec) != n_theta)) {
    stop(paste("Size of prior_theta$prec is (",
               paste(dim(prior_theta$prec), collapse = ",", sep = ""),
               ") but should be (",
               paste(c(n_theta, n_theta), collapse = ",", sep = ""),
               ")."))
  }

  prior_theta
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dim PARAM_DESCRIPTION, Default: NULL
#' @param dof PARAM_DESCRIPTION, Default: NULL
#' @param alpha PARAM_DESCRIPTION, Default: 2
#' @param prior_range `list(median, factor)`,
#'   Default: `factor=5`
#' @param prior_sigma `list(median, factor)`,
#'   Default: `factor=5`
#' @param B_range PARAM_DESCRIPTION, Default: `matrix(c(0, 1, 0), 1, 3)`
#' @param B_sigma PARAM_DESCRIPTION, Default: `matrix(c(0, 0, 1), 1, 3)`
#' @param prior_theta If not `NULL`, `prior_range` and `prior_sigma` will only
#' be used to set the offset columns of `B_range` and `B_sigma`, and not to set
#' the prior distribution for the `theta` parameters.
#'   A list with elements `mean` and `prec`.
#'   Default: `NULL`
#' @param \dots PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details The priors for `range` (the spatial correlation range parameter)
#' and `sigma` (the field standard deviation parameter) are constructed so that
#' the prior median is given by `median`, and
#' \deqn{P(median/factor < param < median*factor) = 0.98}
#' which means that the std.dev. of `log(param)` is `log(factor)/qnorm(0.99)`.
#'
#' Remark: In the old `param2.matern.orig` function, the internal parameter
#' scale used a default prior precision \eqn{0.1}. This corresponds to
#' `factor=exp(qnorm((1+0.98)/2)/sqrt(0.1))`, that equals \eqn{1566.435}.
#' This large value has lead to numerical and other problems for many models.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname param2.matern
#' @export
#'

param2.matern_calc <-
  function(dim = NULL,
           dof = NULL,
           alpha = 2,
           prior_range = NULL,
           prior_sigma = NULL,
           B_range = matrix(c(0, 1, 0), 1, 3),
           B_sigma = matrix(c(0, 0, 1), 1, 3),
           prior_theta = NULL,
           ...) {

    if (is.null(dim)) {
      stop("'dim' must not be NULL.")
    }
    if (is.null(dof)) {
      stop("'dof' must not be NULL.")
    }
    if (is.null(B_range)) {
      stop("B_range must not be NULL.")
    }
    if (is.null(B_sigma)) {
      stop("B_sigma must not be NULL.")
    }
    stopifnot(!is.null(prior_range))
    stopifnot(!is.null(prior_sigma))

    if (!is.list(prior_range)) {
      prior_range <- list(prior_range)
    }
    if (length(prior_range) < 2) {
      prior_range <- c(prior_range, list(5))
    }
    if (!is.list(prior_sigma)) {
      prior_sigma <- list(prior_sigma)
    }
    if (length(prior_sigma) < 2) {
      prior_sigma <- c(prior_sigma, list(5))
    }
    is.stationary <- (nrow(B_range) == 1) && (nrow(B_sigma) == 1) &&
      (length(prior_range[[1]]) == 1) && (length(prior_sigma[[1]]) == 1)

    d <- dim
    nu <- alpha - d / 2
    stopifnot(nu > 0)

    if (is.null(dof)) {
      dof <- max(nrow(B_range), nrow(B_sigma),
                 length(prior_range[[1]]), length(prior_sigma[[1]]))
    }
    n_theta <- ncol(B_range) - 1L

    B_range <- INLA:::inla.spde.homogenise_B_matrix(B_range, dof, n_theta)
    B_sigma <- INLA:::inla.spde.homogenise_B_matrix(B_sigma, dof, n_theta)
    B_range[, 1] <- log(prior_range[[1]])
    B_sigma[, 1] <- log(prior_sigma[[1]])
    B_tau_kappa <- B_range_sigma_to_tau_kappa(B_range = B_range,
                                              B_sigma = B_sigma,
                                              d = d,
                                              nu = nu)
    B_tau <- B_tau_kappa$B_tau
    B_kappa <- B_tau_kappa$B_kappa

    if (is.stationary) {
      B_tau <- B_tau[1, , drop = FALSE]
      B_kappa <- B_kappa[1, , drop = FALSE]
    }

    if (n_theta > 0) {
      B_theta <- cbind(0, diag(1, n_theta, n_theta))
      rownames(B_theta) <- rownames(B_theta, do.NULL = FALSE, prefix = "theta.")
    } else {
      B_theta <- NULL
    }

    rownames(B_tau) <- rownames(B_tau, do.NULL = FALSE, prefix = "tau.")
    rownames(B_kappa) <- rownames(B_kappa, do.NULL = FALSE, prefix = "kappa.")
    B_variance <- 2 * B_sigma
    rownames(B_variance) <-
      rownames(B_variance, do.NULL = FALSE, prefix = "variance.nominal.")
    rownames(B_sigma) <-
      rownames(B_sigma, do.NULL = FALSE, prefix = "sigma.nominal.")
    rownames(B_range) <-
      rownames(B_range, do.NULL = FALSE, prefix = "range.nominal.")
    BLC <- rbind(B_theta, B_range, B_sigma, B_tau, B_kappa, B_variance)

    ## Construct prior.
    prior_theta <- construct_prior(B_range, B_sigma,
                                   prior_range, prior_sigma,
                                   d = d, nu = nu,
                                   prior_theta)

    param <-
      list(alpha = alpha,
           is.stationary = is.stationary,
           B_range = B_range, B_sigma = B_sigma,
           B.tau = B_tau, B.kappa = B_kappa, BLC = BLC,
           theta.prior.mean = prior_theta$mean,
           theta.prior.prec = prior_theta$prec)
    return(param)
  }
