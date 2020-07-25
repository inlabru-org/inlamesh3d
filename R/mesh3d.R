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
  class(mesh) <- "inla.mesh3d"
  mesh
}

row_cross_product <- function(e1, e2) {
  cbind(e1[, 2] * e2[, 3] - e1[, 3] * e2[, 2],
        e1[, 3] * e2[, 1] - e1[, 1] * e2[, 3],
        e1[, 1] * e2[, 2] - e1[, 2] * e2[, 1])
}

row_volume_product <- function(e1, e2, e3) {
  rowSums(row_cross_product(e1, e2) * e3)
}

inla.mesh3d.fem <- function(mesh) {
  v1 <- mesh$loc[mesh$graph$tv[, 1], , drop = FALSE]
  v2 <- mesh$loc[mesh$graph$tv[, 2], , drop = FALSE]
  v3 <- mesh$loc[mesh$graph$tv[, 3], , drop = FALSE]
  v4 <- mesh$loc[mesh$graph$tv[, 4], , drop = FALSE]
  e1 <- v2 - v1
  e2 <- v3 - v2
  e3 <- v4 - v3
  e4 <- v1 - v4
  vols_t <- abs(row_volume_product(e1, e2, e3)) / 6

  c0 <- Matrix::sparseMatrix(i = as.vector(mesh$graph$tv),
                             j = as.vector(mesh$graph$tv),
                             x = rep(vols_t / 4, times = 4),
                             dims = c(mesh$n, mesh$n))
  vols_v <- Matrix::diag(c0)

  # Sign changes for b2 and b4 for consistent in/out vector orientation
  b1 <- row_cross_product(e2, e3)
  b2 <- -row_cross_product(e3, e4)
  b3 <- row_cross_product(e4, e1)
  b4 <- -row_cross_product(e1, e2)

  g_i <- g_j <- g_x <- c()
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    GG <- rbind(b1[tt, , drop = FALSE],
               b2[tt, , drop = FALSE],
               b3[tt, , drop = FALSE],
               b4[tt, , drop = FALSE])
    g_i <- c(g_i, rep(mesh$graph$tv[tt, ], each = 4))
    g_j <- c(g_j, rep(mesh$graph$tv[tt, ], times = 4))
    g_x <- c(g_x, as.vector((GG %*% t(GG)) / vols_t[tt] / 36))
  }
  g1 <- Matrix::sparseMatrix(i = g_i,
                             j = g_j,
                             x = g_x,
                             dims = c(mesh$n, mesh$n))

  list(c0 = c0, g1 = g1, g2 = g1 %*% (g1 / vols_v), va = vols_v, vt = vols_t)
}

inla.mesh3d.bary <- function(mesh, loc) {
  stopifnot(inherits(mesh, "inla.mesh3d"))
  loc_1 <- rbind(t(loc), 1)
  bary <- matrix(0, nrow(loc), 4)
  vt <- integer(nrow(loc))
  found <- rep(FALSE, nrow(loc))
  for (tt in seq_len(nrow(mesh$graph$tv))) {
    loc_t <- loc[mesh$graph$tv[tt, ], , drop = FALSE]
    # Barycentric coordinates fulfil
    # 1) rbind(t(loc), 1) = rbind(t(loc_t), 1) * w
    # 2) w >= 0
    w <- solve(rbind(t(loc_t), 1), loc_1)
    idx <- which(colSums(w >= 0) == 4)
    if (length(idx) > 0) {
      idx_v <- which(!found)[idx]
      bary[idx_v, ] <- t(w[, idx, drop = FALSE])
      vt[idx_v] <- tt
      found[which(!found)[idx]] <- TRUE
      loc_1 <- loc_1[, !found, drop = FALSE]
    }
  }
  list(bary = bary, vt = vt)
}

inla.mesh3d.make.A <- function(mesh, loc) {
  stopifnot(inherits(mesh, "inla.mesh3d"))
  bary <- inla.mesh3d.bary(mesh, loc)
  ok <- bary$vt > 0
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

param2.matern3d <-
  function(mesh,
           alpha=2,
           B.tau = matrix(c(0, 1, 0), 1, 3),
           B.kappa = matrix(c(0, 0, 1), 1, 3),
           prior.variance.nominal = 1,
           prior.range.nominal = NULL,
           prior.tau = NULL,
           prior.kappa = NULL,
           theta.prior.mean = NULL,
           theta.prior.prec = 0.1) {
    INLA::inla.require.inherits(mesh, "inla.mesh3d", "'mesh'")
    if (is.null(B.tau))
      stop("B.tau must not be NULL.")
    if (is.null(B.kappa))
      stop("B.kappa must not be NULL.")
    is.stationary = (nrow(B.kappa)==1) && (nrow(B.tau)==1)

    d <- 3
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = mesh$n
    n.theta = ncol(B.kappa) - 1L

    B.kappa = INLA::inla.spde.homogenise_B_matrix(B.kappa, n.spde, n.theta)
    B.tau = INLA::inla.spde.homogenise_B_matrix(B.tau, n.spde, n.theta)

    B.prec =
      cbind(+ lgamma(alpha.nominal) - lgamma(nu.nominal)
            + (d/2)*log(4*pi) + 2*nu.nominal*B.kappa[,1]
            + 2*B.tau[,1],
            2*nu.nominal*B.kappa[,-1,drop=FALSE]
            + 2*B.tau[,-1,drop=FALSE] )
    if (is.stationary) {
      B.tau = B.tau[1,,drop=FALSE]
      B.kappa = B.kappa[1,,drop=FALSE]
      B.prec = B.prec[1,,drop=FALSE]
    }
    B.variance = -B.prec
    B.range =
      cbind(0.5*log(8*nu.nominal)-B.kappa[,1],
            -B.kappa[,-1,drop=FALSE])

    if (n.theta>0) {
      B.theta = cbind(0,diag(1, n.theta))
      rownames(B.theta) <- rownames(B.theta, do.NULL=FALSE, prefix="theta.")
    } else {
      B.theta = NULL
    }
    rownames(B.tau) <- rownames(B.tau, do.NULL=FALSE, prefix="tau.")
    rownames(B.kappa) <- rownames(B.kappa, do.NULL=FALSE, prefix="kappa.")
    rownames(B.variance) <-
      rownames(B.variance, do.NULL=FALSE, prefix="variance.nominal.")
    rownames(B.range) <-
      rownames(B.range, do.NULL=FALSE, prefix="range.nominal.")
    BLC = rbind(B.theta, B.tau, B.kappa, B.variance, B.range)


    ## Construct priors.
    if (is.null(theta.prior.prec)) {
      theta.prior.prec = diag(0.1, n.theta, n.theta)
    } else {
      theta.prior.prec = as.matrix(theta.prior.prec)
      if (ncol(theta.prior.prec) == 1) {
        theta.prior.prec =
          diag(as.vector(theta.prior.prec), n.theta, n.theta)
      }
      if ((nrow(theta.prior.prec) != n.theta) ||
          (ncol(theta.prior.prec) != n.theta)) {
        stop(paste("Size of theta.prior.prec is (",
                   paste(dim(theta.prior.prec), collapse=",", sep=""),
                   ") but should be (",
                   paste(c(n.theta, n.theta), collapse=",", sep=""),
                   ")."))
      }
    }

    if (is.null(theta.prior.mean)) {
      if (is.null(prior.range.nominal)) {
        stop("'prior.range.nominal' must be provided.")
      } else {
        if (!is.numeric(prior.range.nominal) ||
            (length(prior.range.nominal) != 1)) {
          stop(paste0("'prior.range.nominal' must be NULL or a single scalar value.\n",
                      "Did you intend to supply 'prior.range' to inla.spde2.pcmatern instead?"))
        }
      }

      if (is.null(prior.kappa)) {
        prior.kappa = sqrt(8*nu.nominal)/prior.range.nominal
      }
      if (is.null(prior.tau)) {
        prior.tau =
          sqrt(gamma(nu.nominal)/gamma(alpha.nominal)/
                 (4*pi*prior.kappa^(2*nu.nominal)*prior.variance.nominal))
      }

      if (n.theta>0) {
        theta.prior.mean =
          qr.solve(rbind(B.tau[,-1,drop=FALSE], B.kappa[,-1,drop=FALSE]),
                   c(log(prior.tau) - B.tau[,1],
                     log(prior.kappa) - B.kappa[,1]))
      } else {
        theta.prior.mean = rep(0, n.theta) ## Empty vector
      }
    }

    param =
      list(is.stationary=is.stationary,
           B.tau=B.tau, B.kappa=B.kappa, BLC=BLC,
           theta.prior.mean=theta.prior.mean,
           theta.prior.prec=theta.prior.prec)
    return(param)
  }

inla.spde2.matern3d <-
  function(mesh,
           param = NULL,
           constr = FALSE,
           extraconstr.int = NULL,
           extraconstr = NULL,
           fractional.method = c("parsimonious", "null"),
           n.iid.group = 1,
           ...)
  {
    INLA::inla.require.inherits(mesh, "inla.mesh.3d", "'mesh'")
    fractional.method = match.arg(fractional.method)

    if (is.null(param)) {
      stop("'param' must be provided. See ?param2.matern3d")
    } else {
      deprecated =
        !c(missing(B.tau), missing(B.kappa),
           missing(prior.variance.nominal),
           missing(prior.range.nominal),
           missing(prior.tau), missing(prior.kappa),
           missing(theta.prior.mean), missing(theta.prior.prec))
      deprecated =
        c("B.tau", "B.kappa",
          "prior.variance.nominal",
          "prior.range.nominal",
          "prior.tau", "prior.kappa",
          "theta.prior.mean", "theta.prior.prec")[deprecated]
      if (length(deprecated) > 0) {
        warning(paste("'param' specified;  ",
                      "Ignoring deprecated parameter(s) ",
                      paste(deprecated, collapse=", "), ".", sep=""))
      }
    }

    d <- 3
    nu = alpha-d/2
    nu.nominal = max(0.5, nu)
    alpha.nominal = max(nu.nominal+d/2, alpha)

    n.spde = meh$n
    n.theta = ncol(B.kappa)-1L

    # d == 3
    fem <- inla.mesh3d.fem(mesh)

    if (alpha==2) {
      B.phi0 = param$B.tau
      B.phi1 = 2*param$B.kappa
      M0 = fem$c0
      M1 = fem$g1
      M2 = fem$g2
    } else if (alpha==1) {
      B.phi0 = param$B.tau
      B.phi1 = param$B.kappa
      M0 = fem$c0
      M1 = fem$g1*0
      M2 = fem$g1
    } else if (!param$is.stationary) {
      stop("Non-stationary Matern with fractional alpha is not implemented.")
    } else if ((alpha<2) && (alpha>1)) {
      if (fractional.method == "parsimonious") {
        lambda = alpha-floor(alpha)
        b = matrix(c(1,0,0, 1,1,0, 1,2,1),3,3) %*%
          solve(matrix(1/(c(4:2, 3:1, 2:0)+lambda), 3, 3),
                1/(c(4:2)+lambda-alpha))
      } else if (fractional.method == "null") {
        b = c(1,alpha,alpha*(alpha-1)/2)
      } else {
        stop(paste("Unknown fractional.method '", fractional.method,
                   "'.", sep=""))
      }
      B.phi0 = param$B.tau + (alpha-2)*param$B.kappa
      B.phi1 = 2*param$B.kappa
      M0 = fem$c0*b[1]
      M1 = fem$g1*b[2]/2
      M2 = fem$g2*b[3]
    } else if ((alpha<1) && (alpha>0)) {
      if (fractional.method == "parsimonious") {
        lambda = alpha-floor(alpha)
        b = matrix(c(1,0,1,1),2,2) %*%
          solve(matrix(1/(c(2:1, 1:0)+lambda), 2, 2),
                1/(c(2:1)+lambda-alpha))
      } else if (fractional.method == "null") {
        b = c(1,alpha)
      } else {
        stop(paste("Unknown fractional.method '", fractional.method,
                   "'.", sep=""))
      }
      B.phi0 = param$B.tau + (alpha-1)*param$B.kappa
      B.phi1 = param$B.kappa
      M0 = fem$c0*b[1]
      M1 = fem$g1*0
      M2 = fem$g1*b[2]
    } else {
      stop(paste("Unsupported alpha value (", alpha,
                 "). Supported values are 0 < alpha <= 2", sep=""))
    }

    if (n.iid.group == 1) {
      spde =
        INLA::inla.spde2.generic(M0=M0, M1=M1, M2=M2,
                                 B0=B.phi0, B1=B.phi1, B2=1,
                                 theta.mu = param$theta.prior.mean,
                                 theta.Q = param$theta.prior.prec,
                                 transform = "identity",
                                 BLC = param$BLC)
    } else {
      if (nrow(B.phi0) > 1) {
        B.phi0 <- kronecker(matrix(1, n.iid.group, 1), B.phi0)
      }
      if (nrow(B.phi1) > 1) {
        B.phi1 <- kronecker(matrix(1, n.iid.group, 1), B.phi1)
      }
      spde =
        inla.spde2.generic(M0=kronecker(Matrix::Diagonal(n.iid.group), M0),
                           M1=kronecker(Matrix::Diagonal(n.iid.group), M1),
                           M2=kronecker(Matrix::Diagonal(n.iid.group), M2),
                           B0=B.phi0, B1=B.phi1, B2=1,
                           theta.mu = param$theta.prior.mean,
                           theta.Q = param$theta.prior.prec,
                           transform = "identity",
                           BLC = param$BLC)
    }
    spde$model = "matern"
    spde$BLC = param$BLC

    if (constr || !is.null(extraconstr.int) || !is.null(extraconstr)) {
      A.constr = matrix(numeric(0), 0, n.spde*n.iid.group)
      e.constr = matrix(numeric(0), 0, 1)
      if (constr) {
        A.constr <- rbind(A.constr,
                          matrix(colSums(fem$c1)/n.iid.group,
                                 1, n.spde*n.iid.group))
        e.constr <- rbind(e.constr, 0)
      }
      if (!is.null(extraconstr.int)) {
        if (ncol(extraconstr.int$A) == n.spde) {
          A.constr <-
            rbind(A.constr,
                  kronecker(matrix(1/n.iid.group, 1, n.iid.group),
                            as.matrix(extraconstr.int$A %*% fem$c1)))
        } else {
          A.constr <-
            rbind(A.constr,
                  as.matrix(extraconstr.int$A %*%
                              kronecker(Matrix::Diagonal(n.iid.group),
                                        fem$c1)))
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr.int$e))
      }
      if (!is.null(extraconstr)) {
        if (ncol(extraconstr$A) == n.spde) {
          A.constr <-
            rbind(A.constr,
                  kronecker(matrix(1/n.iid.group, 1, n.iid.group),
                            as.matrix(extraconstr$A)))
        } else {
          A.constr <- rbind(A.constr, as.matrix(extraconstr$A))
        }
        e.constr <- rbind(e.constr, as.matrix(extraconstr$e))
      }

      spde$f$constr = FALSE
      spde$f$extraconstr = list(A=A.constr, e=e.constr)
    }

    ## Attach the mesh, so downstream code can have access
    spde$mesh <- mesh

    return(invisible(spde))
  }
