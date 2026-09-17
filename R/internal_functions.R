#' Start a pool of mirai daemons with pmxcopula loaded, skipping the
#' restart if a matching pool is already running
#'
#' @param cores Number of daemons to run.
#' @noRd
ensure_daemons <- function(cores) {
  if (mirai::status()$connections != cores) {
    mirai::daemons(cores)
    mirai::everywhere(library(pmxcopula))
  }
}

#' Calculate Correlations
#'
#' @param data A data.frame with rows corresponding to observations and columns corresponding to variables.
#' @param pairs_matrix Matrix with 2 column and each row containing a pair of
#' variable names. If set to NULL, every possible variable
#' pair is included.
#'
#' @return A data.frame containing the correlation calculated for the dataset.
#'
#' @importFrom combinat combn
#' @noRd
calc_correlation <- function(data, pairs_matrix = NULL){

  if (is.null(pairs_matrix)) {
    pairs_matrix <- t(combinat::combn(colnames(data), 2))
  }

  cor <- matrix(data = NA, nrow = nrow(pairs_matrix), ncol = 5)
  colnames(cor) <- c("statistic", "Var1", "Var2","var_pair", "cor")
  for (i in 1:nrow(pairs_matrix)) {
    cor_pearson <- cor(x = data[, pairs_matrix[i,1]],
                       y = data[, pairs_matrix[i,2]],
                       use = "pairwise.complete.obs",
                       method = "pearson")

    cor[i,] <- c("correlation", pairs_matrix[i,1], pairs_matrix[i,2],
                 paste(pairs_matrix[i,1], pairs_matrix[i,2], sep = "-"),
                 cor_pearson)
  }

  cor <- cor |> as.data.frame() |> dplyr::mutate(cor = as.numeric(cor))
  return(cor)

}

#' Kernel functional estimate for an isotropic bivariate normal kernel
#'
#' R driver around the Rcpp core cpp_kfe_isotropic_batch(). Matches
#' ks:::kfe()'s output (validated against it in data-raw/validate_fast_hpi.R).
#'
#' @noRd
fast_kfe2d <- function(x, g, r) {
  ind_mat <- ks:::dmvnorm.deriv(x = rep(0, 2), deriv.order = r, only.index = TRUE, deriv.vec = TRUE)
  uniq_vals <- cpp_kfe_isotropic_batch(x, g, r, TRUE)
  # cpp_kfe_isotropic_batch()'s output is ordered (r,0), (r-1,1), ..., (0,r)
  key <- paste(ind_mat[, 1], ind_mat[, 2])
  ukey <- paste(r - (0:r), 0:r)
  uniq_vals[match(key, ukey)]
}

#' Kernel functional estimate for an isotropic bivariate normal kernel,
#' binned
#'
#' Same as fast_kfe2d(), but for n above ks:::default.bflag()'s ~500
#' threshold, where ks:::kfe() itself switches to a binned (linear-binning +
#' FFT-convolution) estimate instead of the exact O(n^2) sum. The binning
#' itself (ks:::binning()) is already fast and reused unchanged; only the
#' kernel-grid evaluation and convolution are done via the Rcpp core
#' cpp_kfe_isotropic_binned_batch(). Matches ks:::kfe(..., binned = TRUE)'s
#' output (validated against it in data-raw/validate_fast_hpi.R).
#'
#' @noRd
fast_kfe2d_binned <- function(x, g, r) {
  bp <- ks:::binning(x, H = g^2 * diag(2))
  # ks:::kdde.binned.nd() uses sum(bin.par$counts), not nrow(x) - equal for
  # unweighted data (always the case here) but kept exact regardless.
  n <- sum(bp$counts)
  delta <- vapply(bp$eval.points, function(e) (max(e) - min(e)) / (length(e) - 1), numeric(1))
  ind_mat <- ks:::dmvnorm.deriv(x = rep(0, 2), deriv.order = r, only.index = TRUE, deriv.vec = TRUE)
  uniq_vals <- cpp_kfe_isotropic_binned_batch(bp$counts, delta, g, r, n)
  key <- paste(ind_mat[, 1], ind_mat[, 2])
  ukey <- paste(r - (0:r), 0:r)
  uniq_vals[match(key, ukey)]
}

#' Fast bivariate plug-in bandwidth selector
#'
#' A faster ks::Hpi(x, nstage = 2, pilot = "samse", pre = "sphere") - its
#' own default for 2D data. Reuses ks's own pilot-estimation and
#' optimization code (ks:::gsamse, ks:::invvec, ks:::nur, ...) and only
#' replaces the kernel functional estimation step - via fast_kfe2d() (exact,
#' n <= ks:::default.bflag()'s ~500 threshold) or fast_kfe2d_binned()
#' (binned, above it), matching ks::Hpi()'s own choice there exactly.
#' Validated against ks::Hpi() in data-raw/validate_fast_hpi.R. Falls back
#' to plain ks::Hpi() for anything not 2D or too few points.
#'
#' @param x A matrix or data.frame with 2 columns.
#'
#' @return A 2x2 bandwidth matrix, as ks::Hpi() would return.
#' @noRd
fast_Hpi2d <- function(x) {
  tryCatch({
    x <- as.matrix(x)
    if (ncol(x) != 2 || nrow(x) < 10) stop("fast path only validated for 2D data")

    n <- nrow(x); d <- 2
    binned <- ks:::default.bflag(d = 2, n = n)
    x.star <- ks:::pre.sphere(x)
    S12 <- ks:::matrix.sqrt(var(x))

    S.star <- var(x.star)
    g6.star <- ks:::gsamse(S.star, n = n, modr = 6)
    psihat6.star <- if (binned) fast_kfe2d_binned(x.star, g6.star, 6) else fast_kfe2d(x.star, g6.star, 6)
    g.star <- ks:::gsamse(S.star, n = n, modr = 4, nstage = 2, psihat = psihat6.star)
    psihat.star <- if (binned) fast_kfe2d_binned(x.star, g.star, 4) else fast_kfe2d(x.star, g.star, 4)

    psi2r4.mat <- ks:::invvec(psihat.star)
    Hstart <- ks::Hns(x = x.star, deriv.order = 0)

    Idr <- diag(1)
    pi.temp <- function(vechH) {
      H <- ks:::invvech(vechH) %*% ks:::invvech(vechH)
      Hinv <- chol2inv(chol(H))
      IdrvH <- Idr %x% ks:::vec(H)
      int.var <- 1 / (det(H)^(1/2) * n) * ks:::nur(r = 0, A = Hinv, mu = rep(0, d), Sigma = diag(d)) *
        2^(-d - 0) * pi^(-d/2)
      pi.val <- int.var + 1/4 * sum(diag(t(IdrvH) %*% psi2r4.mat %*% IdrvH))
      drop(pi.val)
    }

    Hstart <- ks:::matrix.sqrt(Hstart)
    result <- stats::optim(ks:::vech(Hstart), pi.temp, method = "BFGS", control = list(trace = 0, REPORT = 1))
    H <- ks:::invvech(result$par) %*% ks:::invvech(result$par)
    S12 %*% H %*% S12
  }, error = function(e) ks::Hpi(x))
}

#' Exact contour levels for a bivariate kde
#'
#' A faster ks::contourLevels.kde(fhat, cont = 1:99, approx = FALSE): same
#' quantile-of-fitted-density calculation, with the O(n^2) density
#' evaluation done via cpp_dmvnorm_mixture_eval() instead of ks::kde().
#' Validated in data-raw/validate_fast_hpi.R.
#'
#' @param x A matrix with 2 columns (the data the kde was fitted on).
#' @param H The same bandwidth matrix the kde was fitted with.
#' @param w Weights, as passed to ks::kde(); NULL means unweighted.
#'
#' @return A named numeric vector, as ks::contourLevels() returns.
#' @noRd
fast_contour_levels <- function(x, H, w = NULL) {
  x <- as.matrix(x)
  if (is.null(w)) w <- rep(1, nrow(x))
  dobs <- cpp_dmvnorm_mixture_eval(x, x, H, w)
  stats::quantile(dobs, prob = (100 - (1:99)) / 100)
}

#' Unbinned kde grid evaluation for a bivariate kde
#'
#' A faster ks::kde(x, H, compute.cont = FALSE) in the unbinned case (n
#' below ks:::default.bflag()'s ~500 threshold for 2D - see fast_kde()).
#' Builds the same evenly-spaced 151x151 grid ks:::make.grid.ks() would,
#' then evaluates the density on it via cpp_kde_grid_truncated(), which
#' replicates ks:::kde.grid.2d()'s own support-window-truncated computation
#' (see ks:::make.supp(), ks:::find.gridpts()), just compiled instead of
#' going through repeated mvtnorm::dmvnorm() calls. Validated in
#' data-raw/validate_fast_hpi.R.
#'
#' Only returns the fields this package's call sites actually use (x,
#' eval.points, estimate, H, gridded), not a full "kde"-class object.
#'
#' @param x A matrix or data.frame with 2 columns.
#' @param H Bandwidth matrix.
#' @param w Weights; NULL means unweighted.
#'
#' @noRd
fast_kde_grid <- function(x, H, w = NULL) {
  x <- as.matrix(x)
  if (is.null(w)) w <- rep(1, nrow(x))

  Hsqrt <- ks:::matrix.sqrt(H)
  tol <- 3.7 * diag(Hsqrt)
  xmin <- apply(x, 2, min) - tol
  xmax <- apply(x, 2, max) + tol
  gx <- seq(xmin[1], xmax[1], length.out = 151)
  gy <- seq(xmin[2], xmax[2], length.out = 151)

  estimate <- cpp_kde_grid_truncated(gx, gy, x, H, w, tol[1], tol[2])

  list(x = x, eval.points = list(gx, gy), estimate = estimate, H = H,
       gridtype = c("linear", "linear"), gridded = TRUE, w = w, names = colnames(x))
}

#' Fast 2D kernel density estimate for contour extraction
#'
#' Same as ks::kde(), including its own choice of binned vs unbinned grid
#' evaluation (ks:::default.bflag(), not overridden here), except:
#'
#' - if H is not supplied, it's computed via fast_Hpi2d() instead of
#'   ks::kde()'s own internal ks::Hpi() call.
#' - for compute.cont = TRUE, approx.cont = FALSE and n <= 500 (the case
#'   ks::kde() itself would use its slow unbinned path for, and this
#'   package's hot path always uses), the grid evaluation and exact contour
#'   levels are computed via fast_kde_grid()/fast_contour_levels() instead
#'   of ks::kde()'s own (much slower) versions of the same thing.
#'
#' Every other combination falls through to plain ks::kde(), untouched.
#'
#' @param data A data.frame or matrix with 2 columns.
#' @param ... Passed on to ks::kde(), e.g. H, compute.cont, approx.cont, w.
#'
#' @return An object of class "kde" (or, in the fast-path case, a plain
#' list with the same fields this package's call sites actually use).
#' @noRd
fast_kde <- function(data, ...) {
  dots <- list(...)
  if (is.null(dots$H) && NCOL(data) == 2) dots$H <- fast_Hpi2d(data)

  fast_path <- isTRUE(dots$compute.cont) && identical(dots$approx.cont, FALSE) &&
    NCOL(data) == 2 && NROW(data) <= 500
  if (fast_path) {
    kd <- fast_kde_grid(data, dots$H, dots$w)
    kd$cont <- fast_contour_levels(data, dots$H, dots$w)
    kd
  } else {
    do.call(ks::kde, c(list(data), dots))
  }
}
