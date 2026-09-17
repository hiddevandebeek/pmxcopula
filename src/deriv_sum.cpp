#include <Rcpp.h>
using namespace Rcpp;

// Probabilist's Hermite polynomial He_n(z), n in 0..6.
static inline double he(int n, double z) {
  switch (n) {
    case 0: return 1.0;
    case 1: return z;
    case 2: return z*z - 1.0;
    case 3: return z*z*z - 3.0*z;
    case 4: return z*z*z*z - 6.0*z*z + 3.0;
    case 5: return z*z*z*z*z - 10.0*z*z*z + 15.0*z;
    case 6: return z*z*z*z*z*z - 15.0*z*z*z*z + 45.0*z*z - 15.0;
  }
  return NA_REAL;
}

//' Kernel functional estimate for an isotropic bivariate normal kernel
//'
//' Computes the same O(n^2) pairwise-difference sum as
//' ks:::dmvnorm.deriv.sum(x, Sigma = g^2*diag(2), deriv.order = r1 + r2,
//' kfe = kfe, binned = FALSE) for one derivative-index combination
//' (r1, r2), using the closed-form isotropic Hermite-polynomial derivative
//' (D^(r1,r2) of an isotropic bivariate normal collapses to a product of
//' 1D Hermite polynomials) instead of ks's general anisotropic tensor
//' recursion. This is what ks::kfe() ends up computing whenever it's
//' called with an isotropic pilot bandwidth, which is the only case
//' ks::Hpi() uses for its "samse" pilot (the default for 2D data).
//' Validated in data-raw/validate_fast_hpi.R.
//'
//' @noRd
// [[Rcpp::export]]
double cpp_kfe_isotropic(NumericMatrix x, double g, int r1, int r2, bool kfe) {
  int n = x.nrow();
  int r = r1 + r2;
  double invg = 1.0 / g;
  double sign = (r % 2 == 0) ? 1.0 : -1.0;
  double C = sign * std::pow(g, -r) / (2.0 * M_PI * g * g);

  const double* x0 = &x(0, 0);
  const double* x1 = &x(0, 1);

  double s = 0.0;
  for (int i = 0; i < n; i++) {
    double xi0 = x0[i], xi1 = x1[i];
    for (int j = i + 1; j < n; j++) {
      double u = (xi0 - x0[j]) * invg;
      double v = (xi1 - x1[j]) * invg;
      double expo = std::exp(-0.5 * (u*u + v*v));
      s += he(r1, u) * he(r2, v) * expo;
    }
  }
  double hm0 = he(r1, 0.0) * he(r2, 0.0); // He_n(0) times exp(0) = He_n(0)
  double total = C * (2.0 * s + (double)n * hm0);
  if (kfe) total /= ((double)n * (double)n);
  return total;
}

//' Same as cpp_kfe_isotropic(), but for all derivative-index combinations
//' of a given order r at once (r=6 -> 7 combos, r=4 -> 5 combos), sharing
//' the pairwise-distance and exp() work across combinations instead of
//' redoing it per combination - only the Hermite polynomial factor differs
//' between them.
//'
//' @noRd
// [[Rcpp::export]]
NumericVector cpp_kfe_isotropic_batch(NumericMatrix x, double g, int r, bool kfe) {
  int n = x.nrow();
  int ncombo = r + 1;
  double invg = 1.0 / g;
  double sign = (r % 2 == 0) ? 1.0 : -1.0;
  double C = sign * std::pow(g, -r) / (2.0 * M_PI * g * g);

  const double* x0 = &x(0, 0);
  const double* x1 = &x(0, 1);

  std::vector<double> s(ncombo, 0.0);

  for (int i = 0; i < n; i++) {
    double xi0 = x0[i], xi1 = x1[i];
    for (int j = i + 1; j < n; j++) {
      double u = (xi0 - x0[j]) * invg;
      double v = (xi1 - x1[j]) * invg;
      double expo = std::exp(-0.5 * (u*u + v*v));
      for (int k = 0; k < ncombo; k++) {
        s[k] += he(r - k, u) * he(k, v) * expo;
      }
    }
  }

  NumericVector total(ncombo);
  for (int k = 0; k < ncombo; k++) {
    double hm0 = he(r - k, 0.0) * he(k, 0.0);
    double tot = C * (2.0 * s[k] + (double)n * hm0);
    if (kfe) tot /= ((double)n * (double)n);
    total[k] = tot;
  }
  return total;
}

//' Weighted bivariate normal mixture density, evaluated at a set of points
//'
//' Computes the same thing as ks:::kde.points(x, H, eval.points, w)$estimate
//' (equivalently ks:::dmvnorm.mixt(x = eval.points, mus = x, Sigmas = H
//' repeated n times, props = w/n)$estimate) - the O(n^2) evaluation
//' ks::contourLevels.kde() falls back to for exact (approx = FALSE)
//' contour levels. H is a general (not necessarily isotropic) 2x2 matrix
//' here, unlike the isotropic-only kfe functions above.
//'
//' @noRd
// [[Rcpp::export]]
NumericVector cpp_dmvnorm_mixture_eval(NumericMatrix eval_points, NumericMatrix x, NumericMatrix H, NumericVector w) {
  int n = x.nrow();
  int ne = eval_points.nrow();
  double Ha = H(0, 0), Hb = H(0, 1), Hd = H(1, 1);
  double det = Ha * Hd - Hb * Hb;
  double norm_const = 1.0 / (2.0 * M_PI * std::sqrt(det)) / (double)n;
  double inv_det = 1.0 / det;

  const double* x0 = &x(0, 0);
  const double* x1 = &x(0, 1);
  const double* e0 = &eval_points(0, 0);
  const double* e1 = &eval_points(0, 1);
  const double* wp = &w[0];

  NumericVector out(ne);
  for (int k = 0; k < ne; k++) {
    double yk0 = e0[k], yk1 = e1[k];
    double s = 0.0;
    for (int i = 0; i < n; i++) {
      double dx = yk0 - x0[i];
      double dy = yk1 - x1[i];
      double maha = (dx * dx * Hd - 2.0 * dx * dy * Hb + dy * dy * Ha) * inv_det;
      s += wp[i] * std::exp(-0.5 * maha);
    }
    out[k] = norm_const * s;
  }
  return out;
}

//' Unbinned bivariate kde grid evaluation, with the same support-window
//' truncation ks::kde() itself uses
//'
//' Computes the same thing as ks:::kde.grid.2d(x, H, gridx, ...)$estimate.
//' ks's version, for every data point, calls mvtnorm::dmvnorm() only on
//' the grid points within a box of half-width supp * diag(matrix.sqrt(H))
//' around that point (see ks:::make.supp(), ks:::find.gridpts()), skipping
//' cells far enough away for the Gaussian kernel to be negligible there.
//' This replicates that same box (tolx, toly, computed in R exactly as ks
//' does) using direct index arithmetic instead of a linear index search
//' (valid since the grid is evenly spaced), and the same closed-form
//' density for points inside the box, without going through
//' mvtnorm::dmvnorm() or R-level loops.
//'
//' @param gridx,gridy Grid coordinates, evenly spaced (as ks::kde() uses by
//' default).
//' @param x,H,w As elsewhere - the fitted data, bandwidth and weights.
//' @param tolx,toly Support window half-widths, supp * diag(matrix.sqrt(H)).
//'
//' @noRd
// [[Rcpp::export]]
NumericMatrix cpp_kde_grid_truncated(NumericVector gridx, NumericVector gridy,
                                      NumericMatrix x, NumericMatrix H, NumericVector w,
                                      double tolx, double toly) {
  int n = x.nrow();
  int nx = gridx.size();
  int ny = gridy.size();
  double Ha = H(0, 0), Hb = H(0, 1), Hd = H(1, 1);
  double det = Ha * Hd - Hb * Hb;
  double norm_const = 1.0 / (2.0 * M_PI * std::sqrt(det)) / (double)n;
  double inv_det = 1.0 / det;

  double gx0 = gridx[0], gy0 = gridy[0];
  double dgx = (gridx[nx - 1] - gx0) / (double)(nx - 1);
  double dgy = (gridy[ny - 1] - gy0) / (double)(ny - 1);

  const double* x0 = &x(0, 0);
  const double* x1 = &x(0, 1);
  const double* wp = &w[0];

  NumericMatrix out(nx, ny);

  for (int i = 0; i < n; i++) {
    double xi0 = x0[i], xi1 = x1[i];
    double wi = wp[i];

    // ks:::find.gridpts() takes the last grid point <= each support bound on
    // both sides (floor, not ceil, on the lower bound too) - see
    // data-raw/validate_fast_hpi.R for how this was caught.
    int gi_min = (int)std::floor((xi0 - tolx - gx0) / dgx);
    int gi_max = (int)std::floor((xi0 + tolx - gx0) / dgx);
    int gj_min = (int)std::floor((xi1 - toly - gy0) / dgy);
    int gj_max = (int)std::floor((xi1 + toly - gy0) / dgy);
    if (gi_min < 0) gi_min = 0;
    if (gj_min < 0) gj_min = 0;
    if (gi_max > nx - 1) gi_max = nx - 1;
    if (gj_max > ny - 1) gj_max = ny - 1;

    for (int gi = gi_min; gi <= gi_max; gi++) {
      double dx = gx0 + gi * dgx - xi0;
      for (int gj = gj_min; gj <= gj_max; gj++) {
        double dy = gy0 + gj * dgy - xi1;
        double maha = (dx * dx * Hd - 2.0 * dx * dy * Hb + dy * dy * Ha) * inv_det;
        out(gi, gj) += wi * norm_const * std::exp(-0.5 * maha);
      }
    }
  }
  return out;
}
