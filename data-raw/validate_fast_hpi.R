# Validates the Rcpp-accelerated bivariate kde/Hpi replacements in
# R/internal_functions.R and src/deriv_sum.cpp against ks's own (much
# slower) reference implementations, using the package's bundled example
# data. Not run as part of R CMD check (this package's tests/ directory is
# itself excluded from the build, see .Rbuildignore) - this is a manual
# regression check to re-run after modifying deriv_sum.cpp/internal_functions.R,
# or after upgrading the ks package, since fast_Hpi2d()/fast_kde() call
# several ks:::-prefixed internals that could change between ks versions.
# Last validated against ks 1.15.1 (see DESCRIPTION's ks (>= 1.15.1) floor -
# bump both together after re-running this against a newer ks release).
#
# All reported diffs should stay at floating-point-noise level (~1e-10 or
# smaller); anything larger signals a real discrepancy, not rounding.

devtools::load_all(".")
data(pediatric_sim)
pairs_matrix <- matrix(c("CREA", "CREA", "AGE", "BW"), 2, 2)

cat("=== fast_Hpi2d() vs ks::Hpi(), all bundled replicates ===\n")
max_rel <- 0
for (b in 1:100) {
  for (p in 1:nrow(pairs_matrix)) {
    sim_data_b <- na.omit(pediatric_sim[pediatric_sim$simulation_nr == b, pairs_matrix[p, ]])
    H_true <- ks::Hpi(sim_data_b)
    H_fast <- fast_Hpi2d(sim_data_b)
    max_rel <- max(max_rel, max(abs((H_fast - H_true) / H_true)))
  }
}
cat("max relative diff in H:", max_rel, "\n\n")

cat("=== fast_kde() (grid + exact contour levels) vs ks::kde(), all bundled replicates ===\n")
max_abs_est <- 0
max_rel_cont <- 0
for (b in 1:100) {
  for (p in 1:nrow(pairs_matrix)) {
    sim_data_b <- na.omit(pediatric_sim[pediatric_sim$simulation_nr == b, pairs_matrix[p, ]])
    H <- fast_Hpi2d(sim_data_b)
    kd_true <- ks::kde(sim_data_b, H = H, compute.cont = TRUE, approx.cont = FALSE)
    kd_fast <- fast_kde(sim_data_b, H = H, compute.cont = TRUE, approx.cont = FALSE)
    max_abs_est <- max(max_abs_est, max(abs(kd_fast$estimate - kd_true$estimate)))
    max_rel_cont <- max(max_rel_cont, max(abs((kd_fast$cont - kd_true$cont) / kd_true$cont), na.rm = TRUE))
  }
}
cat("max abs diff in grid density estimate:", max_abs_est, "\n")
cat("max relative diff in exact contour levels:", max_rel_cont, "\n\n")

# The bundled pediatric_sim replicates alone didn't previously catch a real
# off-by-one in cpp_kde_grid_truncated()'s support-window truncation (ceil vs
# floor on the lower bound - see git history), so also check other bundled
# datasets, weighted kde, and some edge cases (n near the fast-path
# boundaries, near-singular covariance, ties, skew) that are more likely to
# exercise grid cells right at the truncation boundary.
cat("=== fast_Hpi2d() vs ks::Hpi(), other bundled datasets ===\n")
data(pediatric_3cov); data(NHANES_12cov); data(MIMIC_5cov); data(mix_data)

check_dataset <- function(df, name, group_col = NULL) {
  num_cols <- names(df)[sapply(df, is.numeric)]
  if (length(num_cols) < 2) return(invisible())
  pairs <- combn(num_cols, 2, simplify = FALSE)
  max_rel <- 0
  groups <- if (!is.null(group_col) && group_col %in% names(df)) unique(df[[group_col]]) else list(NULL)
  for (g in groups[seq_len(min(length(groups), 15))]) {
    sub <- if (is.null(g)) df else df[df[[group_col]] == g, ]
    for (p in pairs) {
      d <- unique(na.omit(sub[, p]))
      if (nrow(d) < 15) next
      H_true <- tryCatch(ks::Hpi(as.matrix(d)), error = function(e) NULL)
      if (is.null(H_true)) next
      max_rel <- max(max_rel, max(abs((fast_Hpi2d(d) - H_true) / H_true)), na.rm = TRUE)
    }
  }
  cat("  ", name, ": max relative diff in H:", max_rel, "\n")
}
check_dataset(pediatric_3cov, "pediatric_3cov")
check_dataset(NHANES_12cov, "NHANES_12cov")
check_dataset(MIMIC_5cov, "MIMIC_5cov")
check_dataset(mix_data, "mix_data")

cat("\n=== fast_kde() vs ks::kde(), weighted + edge cases ===\n")
check_kde <- function(d, H, w, label) {
  kd_true <- ks::kde(d, H = H, w = w, compute.cont = TRUE, approx.cont = FALSE)
  kd_fast <- fast_kde(d, H = H, w = w, compute.cont = TRUE, approx.cont = FALSE)
  cat("  ", label, ": grid abs diff:", max(abs(kd_fast$estimate - kd_true$estimate)),
      " contour rel diff:", max(abs((kd_fast$cont - kd_true$cont) / kd_true$cont), na.rm = TRUE), "\n")
}

set.seed(1)
n <- 200
d <- data.frame(x = rnorm(n), y = 0.6 * rnorm(n) + rnorm(n))
w <- runif(n, 0.5, 2); w <- w / mean(w)
check_kde(d, fast_Hpi2d(d), w, "weighted, n=200")

set.seed(3)
d_500 <- data.frame(x = rnorm(500), y = rnorm(500) + 0.3 * rnorm(500))
check_kde(d_500, fast_Hpi2d(d_500), rep(1, 500), "n=500 (fast-path boundary)")

set.seed(5)
n <- 150; base <- rnorm(n)
d_corr <- data.frame(x = base, y = base + rnorm(n, sd = 0.01))
check_kde(d_corr, fast_Hpi2d(d_corr), rep(1, n), "near-singular covariance")

set.seed(6)
d_ties <- data.frame(x = round(rnorm(200), 1), y = round(rnorm(200) * 2, 0))
check_kde(d_ties, fast_Hpi2d(d_ties), rep(1, 200), "many tied values")

set.seed(7)
d_skew <- data.frame(x = rgamma(200, shape = 2), y = rexp(200) + rnorm(200, sd = 0.3))
check_kde(d_skew, fast_Hpi2d(d_skew), rep(1, 200), "skewed/heavy-tailed")

# n > ks:::default.bflag()'s ~500 threshold: ks::Hpi()'s own default
# silently switches to a binned pilot estimate there (a different algorithm,
# not just a faster path) - fast_Hpi2d() reimplements that too (linear
# binning via ks:::binning(), unchanged, + FFT convolution via
# cpp_kfe_isotropic_binned_batch()), matching ks::Hpi()'s choice exactly
# rather than falling back. Caught two real bugs getting here (see git
# history): the n > 500 case not being handled at all, and once handled, a
# missing /n^2 normalization + a placement-offset parity mismatch from R's
# round-half-to-even round() in the FFT convolution.
#
# Note: because psi2r4.mat feeds a BFGS optimizer, even inputs matching to
# ~1e-14 can occasionally produce an H differing by ~1e-6 to ~1e-7 on
# numerically ill-conditioned (highly anisotropic) datasets - reproduced
# with ks's own code alone by perturbing its own psi4 by ~1e-14 and
# rerunning its own optim() call, so this is inherent BFGS sensitivity, not
# a sign of a bug. Diffs at that level are expected on some inputs; only
# diffs much larger than that (like the two bugs above) need investigating.
cat("\n=== fast_Hpi2d() vs ks::Hpi(), n > 500 (binned path) ===\n")
for (n in c(501, 600, 1000, 2000, 15000)) {
  set.seed(n)
  d_big <- data.frame(x = rnorm(n), y = 0.4 * rnorm(n) + rnorm(n))
  H_true <- ks::Hpi(d_big)
  H_fast <- fast_Hpi2d(d_big)
  cat("  n =", n, ": max relative diff in H:", max(abs((H_fast - H_true) / H_true)), "\n")
}
