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
    loc = as.matrix(loc),
    graph = list(tv = as.matrix(tv))
  )
  class(mesh) <- "inla_mesh_3d"
  mesh
}


# FEM ----

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

  g_i <- g_j <- g_x <- numeric(nrow(mesh$graph$tv) * 16)
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    GG <- rbind(
      b1[tt, , drop = FALSE],
      b2[tt, , drop = FALSE],
      b3[tt, , drop = FALSE],
      b4[tt, , drop = FALSE]
    )
    ii <- (tt - 1) * 16 + seq_len(16)
    g_i[ii] <- rep(mesh$graph$tv[tt, ], each = 4)
    g_j[ii] <- rep(mesh$graph$tv[tt, ], times = 4)
    g_x[ii] <- as.vector((GG %*% t(GG)) / vols_t[tt] / 36)
  }
  g1 <- Matrix::sparseMatrix(
    i = g_i,
    j = g_j,
    x = g_x,
    dims = c(mesh$n, mesh$n)
  )

  list(c0 = c0, g1 = g1, g2 = g1 %*% (g1 / vols_v), va = vols_v, vt = vols_t)
}


#' @rdname inla.mesh.fem
#' @export
#' @importFrom Matrix sparseMatrix diag
inla.mesh3d.volumes <- function(mesh, ...) {
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

  list(t = vols_t, v = vols_v)
}


# Barycentric coordinates ----

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

inla.mesh3d.bary <- function(mesh, loc, divide_along = 1, t_subset = NULL) {
  stopifnot(inherits(mesh, "inla_mesh_3d"))
  loc_1 <- rbind(t(loc), 1)
  bary <- matrix(-Inf, nrow(loc), 4)
  vt <- integer(nrow(loc))
  if (is.null(t_subset)) {
    t_subset <- seq_len(nrow(mesh$graph$tv))
  }
  if (length(t_subset) > 1000) {
    t_order <- order(mesh$loc[mesh$graph$tv[t_subset, 1], divide_along, drop = FALSE])

    t_left <- t_subset[t_order[seq_len(ceiling(length(t_subset) / 2))]]
    max_left <- max(mesh$loc[as.vector(mesh$graph$tv[t_left, ]), divide_along])
    loc_left <- which(loc[, divide_along] <= max_left)
    message(paste0("Splitting ", divide_along,
                   " into L: #T = ", length(t_left),
                   " #loc = ", length(loc_left)))
    bary_left <- inla.mesh3d.bary(
      mesh,
      loc[loc_left, , drop = FALSE],
      divide_along = (divide_along %% 3) + 1,
      t_subset = t_left
    )

    idx <- which(apply(bary_left$bary, 1, min) >
                   apply(bary[loc_left, , drop = FALSE], 1, min))
    if (length(idx) > 0) {
      bary[loc_left[idx], ] <- bary_left$bary[idx, , drop = FALSE]
      vt[loc_left[idx]] <- bary_left$vt[idx]
    }

    t_right <- setdiff(t_subset, t_left)
    if (length(t_right) > 0) {
      min_right <- min(mesh$loc[as.vector(mesh$graph$tv[t_right, ]), divide_along])
      loc_right <- which(loc[, divide_along] >= min_right)
      message(paste0("Splitting ", divide_along,
                     " into R: #T = ", length(t_right),
                     " #loc = ", length(loc_right)))
      bary_right <- inla.mesh3d.bary(
        mesh,
        loc[loc_right, , drop = FALSE],
        divide_along = (divide_along %% 3) + 1,
        t_subset = t_right
      )


      idx <- which(apply(bary_right$bary, 1, min) >
                     apply(bary[loc_right, , drop = FALSE], 1, min))
      if (length(idx) > 0) {
        bary[loc_right[idx], ] <- bary_right$bary[idx, , drop = FALSE]
        vt[loc_right[idx]] <- bary_right$vt[idx]
      }
    }

  } else {
    message(paste0("Handling: #T = ", length(t_subset),
                   ", #loc = ", nrow(loc),
                   ", #loc/#T = ", nrow(loc) / length(t_subset)))
    time0 <- proc.time()
    for (tt in t_subset) {
      loc_t <- mesh$loc[mesh$graph$tv[tt, ], , drop = FALSE]
      # Barycentric coordinates fulfil
      # 1) rbind(t(loc), 1) = rbind(t(loc_t), 1) * w
      # 2) w >= 0
      try({
        w <- solve(rbind(t(loc_t), 1), loc_1)
        idx <- which(apply(w, 2, min) > apply(bary, 1, min))
        if (length(idx) > 0) {
          bary[idx, ] <- t(w[, idx, drop = FALSE])
          vt[idx] <- tt
        }
      },
      silent = TRUE)
    }
    time1 <- proc.time()
    time <- time1 - time0
    message(paste0(
      "Handled: Time/(10^3 Tetra) = ",
      signif(time[1] / (length(t_subset) / 1e3), 3)
    ))
  }
  list(bary = bary, vt = vt)
}

inla.mesh3d.bary.old <- function(mesh, loc) {
  stopifnot(inherits(mesh, "inla_mesh_3d"))
  loc_1 <- rbind(t(loc), 1)
  bary <- matrix(-Inf, nrow(loc), 4)
  vt <- integer(nrow(loc))
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    loc_t <- mesh$loc[mesh$graph$tv[tt, ], , drop = FALSE]
    # Barycentric coordinates fulfil
    # 1) rbind(t(loc), 1) = rbind(t(loc_t), 1) * w
    # 2) w >= 0
    print(tt)
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
  mesh$loc <- as.matrix(mesh$loc)
  loc <- as.matrix(loc)
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


# SPDE object ----

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





# SPDE with pcmatern priors ----

#' @title Matern SPDE model object with PC prior for INLA
#'
#' @description
#' Create an `inla.spde2` model object for a Matern model, using a PC
#' prior for the parameters.
#'
#' This method constructs a Matern SPDE model, with spatial range \eqn{\rho}
#' and standard deviation parameter \eqn{\sigma}.  In the parameterisation
#'
#' \deqn{(\kappa^2-\Delta)^{\alpha/2}(\tau }{(kappa^2-Delta)^(alpha/2) (tau
#' x(u)) = W(u)}\deqn{ x(u))=W(u)}{(kappa^2-Delta)^(alpha/2) (tau x(u)) = W(u)}
#'
#' the spatial scale parameter \eqn{\kappa=\sqrt{8\nu}/\rho}, where
#' \eqn{\nu=\alpha-d/2}, and \eqn{\tau} is proportional to \eqn{1/\sigma}.
#'
#' Stationary models are supported for \eqn{0 < \alpha \leq 2}{0 < alpha <= 2},
#' with spectral approximation methods used for non-integer \eqn{\alpha}, with
#' approximation method determined by `fractional.method`.
#'
#' Integration and other general linear constraints are supported via the
#' `constr`, `extraconstr.int`, and `extraconstr` parameters,
#' which also interact with `n.iid.group`.
#'
#' @details
#' The joint PC prior density for the spatial range, \eqn{\rho}, and the
#' marginal standard deviation, \eqn{\sigma}, and is \deqn{ }{p(rho, sigma) =
#' (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S sigma) }\deqn{ \pi(\rho,
#' \sigma) = }{p(rho, sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S
#' sigma) }\deqn{ \frac{d \lambda_\rho}{2} \rho^{-1-d/2} \exp(-\lambda_\rho
#' \rho^{-d/2}) }{p(rho, sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S
#' exp(-S sigma) }\deqn{ \lambda_\sigma\exp(-\lambda_\sigma \sigma) }{p(rho,
#' sigma) = (d R)/2 rho^(-1-d/2) exp(-R rho^(-d/2)) S exp(-S sigma) } where
#' \eqn{\lambda_\rho}{R} and \eqn{\lambda_\sigma}{S} are hyperparameters that
#' must be determined by the analyst. The practical approach for this in INLA
#' is to require the user to indirectly specify these hyperparameters through
#' \deqn{P(\rho < \rho_0) = p_\rho} and \deqn{P(\sigma > \sigma_0) = p_\sigma}
#' where the user specifies the lower tail quantile and probability for the
#' range (\eqn{\rho_0} and \eqn{p_\rho}) and the upper tail quantile and
#' probability for the standard deviation (\eqn{\sigma_0} and
#' \eqn{\alpha_\sigma}).
#'
#' This allows the user to control the priors of the parameters by supplying
#' knowledge of the scale of the problem. What is a reasonable upper magnitude
#' for the spatial effect and what is a reasonable lower scale at which the
#' spatial effect can operate? The shape of the prior was derived through a
#' construction that shrinks the spatial effect towards a base model of no
#' spatial effect in the sense of distance measured by Kullback-Leibler
#' divergence.
#'
#' The prior is constructed in two steps, under the idea that having a spatial
#' field is an extension of not having a spatial field. First, a spatially
#' constant random effect (\eqn{\rho = \infty}) with finite variance is more
#' complex than not having a random effect (\eqn{\sigma = 0}). Second, a
#' spatial field with spatial variation (\eqn{\rho < \infty}) is more complex
#' than the random effect with no spatial variation. Each of these extensions
#' are shrunk towards the simpler model and, as a result, we shrink the spatial
#' field towards the base model of no spatial variation and zero variance
#' (\eqn{\rho = \infty} and \eqn{\sigma = 0}).
#'
#' The details behind the construction of the prior is presented in Fuglstad,
#' et al. (2016) and is based on the PC prior framework (Simpson, et al.,
#' 2015).
#'
#' @param mesh The mesh to build the model on, as an [inla.mesh()] or
#' [inla.mesh.1d()] object.
#' @param alpha Fractional operator order, \eqn{0<\alpha\leq 2}{0 < alpha <= 2}
#' supported, for \eqn{\nu=\alpha-d/2>0}.
#' @param param Further model parameters. Not currently used.
#' @param constr If `TRUE`, apply an integrate-to-zero constraint.
#' Default `FALSE`.
#' @param extraconstr.int Field integral constraints.
#' @param extraconstr Direct linear combination constraints on the basis
#' weights.
#' @param fractional.method Specifies the approximation method to use for
#' fractional (non-integer) `alpha` values. `'parsimonious'` gives an
#' overall approximate minimal covariance error, `'null'` uses
#' approximates low-order properties.
#' @param n.iid.group If greater than 1, build an explicitly iid replicated
#' model, to support constraints applied to the combined replicates, for
#' example in a time-replicated spatial model. Constraints can either be
#' specified for a single mesh, in which case it's applied to the average of
#' the replicates (`ncol(A)` should be `mesh$n` for 2D meshes,
#' `mesh$m` for 1D), or as general constraints on the collection of
#' replicates (`ncol(A)` should be `mesh$n * n.iid.group` for 2D
#' meshes, `mesh$m * n.iid.group` for 1D).
#' @param prior.range A length 2 vector, with `(range0,Prange)` specifying
#' that \eqn{P(\rho < \rho_0)=p_\rho}, where \eqn{\rho} is the spatial range of
#' the random field. If `Prange` is `NA`, then `range0` is used
#' as a fixed range value.
#' @param prior.sigma A length 2 vector, with `(sigma0,Psigma)` specifying
#' that \eqn{P(\sigma > \sigma_0)=p_\sigma}, where \eqn{\sigma} is the marginal
#' standard deviation of the field.  If `Psigma` is `NA`, then
#' `sigma0` is used as a fixed range value.
#' @return An `inla.spde2` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.2d()], [inla.mesh.create()],
#' [inla.mesh.1d()], [inla.mesh.basis()],
#' [inla.spde2.matern()], [inla.spde2.generic()]
#' @references Fuglstad, G.-A., Simpson, D., Lindgren, F., and Rue, H. (2016)
#' Constructing Priors that Penalize the Complexity of Gaussian Random Fields.
#' arXiv:1503.00256
#'
#' Simpson, D., Rue, H., Martins, T., Riebler, A., and Sørbye, S. (2015)
#' Penalising model component complexity: A principled, practical approach to
#' constructing priors. arXiv:1403.4630
#'
#' @export inla.spde2.pcmatern
inla.spde2.pcmatern <- function(mesh,
                              ...) {
  UseMethod("inla.spde2.pcmatern", mesh)
}

#' @export
#' @method inla.spde2.pcmatern inla_mesh_3d
#' @rdname inla.spde2.pcmatern

inla.spde2.pcmatern.inla_mesh_3d <- function(mesh, ...) {
  inla.spde2.pcmatern3d(mesh, ...)
}

#' @export
#' @rdname inla.spde2.pcmatern
#' @method inla.spde2.pcmatern default

inla.spde2.pcmatern.default <- function(mesh, ...) {
  INLA::inla.spde2.pcmatern(mesh, ...)
}

#' @export
#' @rdname inla.spde2.pcmatern
inla.spde2.pcmatern3d <-
  function(mesh,
           alpha = 2,
           param = NULL,
           constr = FALSE,
           extraconstr.int = NULL,
           extraconstr = NULL,
           fractional.method = c("parsimonious", "null"),
           n.iid.group = 1,
           prior.range = NULL,
           prior.sigma = NULL)
  {
    ## Implementation of PC prior for standard deviation and range
    ##    - Sets the parametrization to range and standard deviation
    ##    - Sets prior according to hyperparameters for range   : prior.range
    ##                                              and std.dev.: prior.sigma
    ## Calls inla.spde2.matern to construct the object, then changes the prior
    if (inherits(mesh, "inla_mesh_3d")) {
      d <- 3
    } else if (inherits(mesh, "inla.mesh")) {
      d <- 2
    } else if (inherits(mesh, "inla.mesh.1d")) {
      d <- 1
    } else {
      stop(paste("Unknown mesh class '",
                 paste(class(mesh), collapse=",", sep=""),
                 "'.", sep=""))
    }

    if (missing(prior.range) || is.null(prior.range) ||
        !is.vector(prior.range) || (length(prior.range) != 2)) {
      stop("'prior.range' should be a length 2 vector 'c(range0,tailprob)' or a fixed range specified with 'c(range,NA)'.")
    }
    if (missing(prior.sigma) || is.null(prior.sigma) ||
        !is.vector(prior.sigma) || (length(prior.sigma) != 2)) {
      stop("'prior.sigma' should be a length 2 vector 'c(sigma0,tailprob)' or a fixed sigma specified with 'c(sigma,NA)'.")
    }
    if (prior.range[1] <= 0){
      stop("'prior.range[1]' must be a number greater than 0 specifying a spatial range")
    }
    if (prior.sigma[1] <= 0){
      stop("'prior.sigma[1]' must be a number greater than 0 specifying a standard deviation")
    }
    if (!is.na(prior.range[2]) &&
        ((prior.range[2] <= 0) || (prior.range[2] >= 1))) {
      stop("'prior.range[2]' must be a probaility strictly between 0 and 1 (or NA to specify a fixed range)")
    }
    if (!is.na(prior.sigma[2]) &&
        ((prior.sigma[2] <= 0) || (prior.sigma[2] >= 1))) {
      stop("'prior.sigma[2]' must be a probaility strictly between 0 and 1 (or NA to specify a fixed sigma)")
    }

    nu <- alpha-d/2
    if (nu <= 0) {
      stop(paste("Smoothness nu = alpha-dim/2 = ", nu,
                 ", but must be > 0.", sep=""))
    }

    spde   <- inla.spde2.matern(mesh = mesh,
                                param =
                                  param2.matern(
                                    mesh,
                                    alpha = alpha,
                                    prior_range = 1,
                                    prior_sigma = 1,
                                    ),
                                constr = constr,
                                extraconstr.int = extraconstr.int,
                                extraconstr = extraconstr,
                                fractional.method = fractional.method,
                                n.iid.group = n.iid.group)

    ## Calculate hyperparameters
    is.fixed.range <- is.na(prior.range[2])
    if (is.fixed.range) {
      lam1 <- 0
      initial.range <- log(prior.range[1])
    } else {
      lam1 <- -log(prior.range[2])*prior.range[1]^(d/2)
      initial.range <- log(prior.range[1]) + 1
    }

    is.fixed.sigma <- is.na(prior.sigma[2])
    if (is.fixed.sigma){
      lam2 <- 0
      initial.sigma <- log(prior.sigma[1])
    } else{
      lam2 <- -log(prior.sigma[2])/prior.sigma[1]
      initial.sigma <- log(prior.sigma[1]) - 1
    }

    pcmatern.param = c(lam1, lam2, d)

    ## Change prior information
    spde$f$hyper.default <-
      list(theta1=list(prior="pcmatern",
                       param=pcmatern.param,
                       initial=initial.range,
                       fixed=is.fixed.range),
           theta2=list(initial=initial.sigma,
                       fixed=is.fixed.sigma))

    ## Change the model descriptor
    spde$model = "pcmatern"

    invisible(spde)
  }




# Parameter construction ----

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
      sd_lrange <- log(prior_range[[2]]) / stats::qnorm(0.99)
      sd_lsigma <- log(prior_sigma[[2]]) / stats::qnorm(0.99)
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
#' which means that the std.dev. of `log(param)` is `log(factor)/stats::qnorm(0.99)`.
#'
#' Remark: In the old `param2.matern.orig` function, the internal parameter
#' scale used a default prior precision \eqn{0.1}. This corresponds to
#' `factor=exp(stats::qnorm((1+0.98)/2)/sqrt(0.1))`, that equals \eqn{1566.435}.
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



# Plotting

#' @title Draw 3D mesh
#' @description Draw a 3D mesh with rgl
#' @param x A `inla_mesh_3d` object
#' @param include Whether to draw vertices, edges, and triangles, Default: c(FALSE, TRUE, TRUE)
#' @param alpha Opaqueness of vertices, edges, and triangles, Default: c(0.9, 0.3, 0.1)
#' @param col Colour of vertices, edges, and triangles, Default: `c("black", "blue", "red")`
#' @param t_sub A subset of tetrahedron indices to include, Default: NULL
#' (draw all tetrahedra)
#' @param size Vertex size, in pixels
#' @param lwd Edge width, in pixels
#' @param add If TRUE, adds to an existing `rgl` device, Default: FALSE
#' @param ... Further parameters passed through to the rgl plotting functions
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @rdname plot.inla_mesh_3d

plot.inla_mesh_3d <- function(x,
                              include = c(FALSE, TRUE, TRUE),
                              alpha = c(0.9, 0.3, 0.1),
                              col = c("black", "blue", "red"),
                              t_sub =  NULL,
                              size = 5,
                              lwd = 2,
                              fill_col = NULL,
                              add = FALSE, ...) {
  stopifnot(requireNamespace("rgl", quietly = TRUE))
  if (!add) {
    dev <- rgl::open3d()
    rgl::view3d(0, 0, fov = 0)
  }
  else {
    dev <- NULL
  }
  if (is.null(t_sub)) {
    tetrav <- x$graph$tv
  } else {
    tetrav <- x$graph$tv[t_sub, , drop = FALSE]
  }
  # Triangles
  triv <- rbind(tetrav[, -1, drop = FALSE],
                tetrav[, -2, drop = FALSE],
                tetrav[, -3, drop = FALSE],
                tetrav[, -4, drop = FALSE])
  triv <- unique(t(apply(triv, 1, sort)))
  # Edges
  edgev <- rbind(triv[, -1, drop = FALSE],
                 triv[, -2, drop = FALSE],
                 triv[, -3, drop = FALSE])
  edgev <- unique(t(apply(edgev, 1, sort)))
  # Plot
  triv <- as.vector(t(triv))
  edgev <- as.vector(t(edgev))
  if (include[3]) {
    if (is.null(fill_col)) {
      rgl::triangles3d(x$loc[triv, , drop = FALSE],
                       lwd = lwd, color = col[3], alpha = alpha[3], ...)
    } else if (length(fill_col) == nrow(x$loc)) {
      rgl::triangles3d(x$loc[triv, , drop = FALSE],
                       lwd = lwd, color = fill_col[triv],
                       alpha = alpha[3], ...)
    } else {
      stop("Only per-vertex fill colours implemented.")
    }
  }
  if (include[2]) {
    rgl::lines3d(x$loc[edgev, , drop = FALSE],
                 lwd = lwd, color = col[2], alpha = alpha[2], ...)
  }
  if (include[1]) {
    idx <- unique(as.vector(tetrav))
    rgl::points3d(x$loc[idx, , drop = FALSE],
                  size = size, lwd = lwd, color = col[1], alpha = alpha[1],
                  ...)
  }
}

