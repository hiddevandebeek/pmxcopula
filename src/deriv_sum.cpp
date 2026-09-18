#include <Rcpp.h>
#include <complex>
#include <algorithm>
using namespace Rcpp;
typedef std::complex<double> cd;

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

// Fills out[0..r] with He_0(z)..He_r(z) via the standard 3-term recurrence
// He_{n+1}(z) = z*He_n(z) - n*He_{n-1}(z), instead of r+1 separate calls to
// he() that each recompute powers of z from scratch. Exact, not an
// approximation - same values as he(), just fewer flops when several
// orders are needed for the same z (see cpp_kfe_isotropic_batch()).
static inline void he_array(int r, double z, double* out) {
  out[0] = 1.0;
  if (r == 0) return;
  out[1] = z;
  for (int m = 1; m < r; m++) out[m + 1] = z * out[m] - m * out[m - 1];
}

//' Kernel functional estimate for an isotropic bivariate normal kernel, for
//' all derivative-index combinations of a given order r at once (r=6 -> 7
//' combos, r=4 -> 5 combos)
//'
//' Computes the same O(n^2) pairwise-difference sum as
//' ks:::dmvnorm.deriv.sum(x, Sigma = g^2*diag(2), deriv.order = r, kfe =
//' kfe, binned = FALSE), for every derivative-index combination of order r
//' at once, sharing the pairwise-distance and exp() work across
//' combinations instead of redoing it per combination - only the Hermite
//' polynomial factor differs between them. Uses the closed-form isotropic
//' Hermite-polynomial derivative (D^(r1,r2) of an isotropic bivariate
//' normal collapses to a product of 1D Hermite polynomials) instead of
//' ks's general anisotropic tensor recursion - this is what ks::kfe() ends
//' up computing whenever it's called with an isotropic pilot bandwidth,
//' which is the only case ks::Hpi() uses for its "samse" pilot (the
//' default for 2D data). Validated in data-raw/validate_fast_hpi.R.
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
  std::vector<double> heu(ncombo), hev(ncombo);

  for (int i = 0; i < n; i++) {
    double xi0 = x0[i], xi1 = x1[i];
    for (int j = i + 1; j < n; j++) {
      double u = (xi0 - x0[j]) * invg;
      double v = (xi1 - x1[j]) * invg;
      double expo = std::exp(-0.5 * (u*u + v*v));
      he_array(r, u, heu.data());
      he_array(r, v, hev.data());
      for (int k = 0; k < ncombo; k++) {
        s[k] += heu[r - k] * hev[k] * expo;
      }
    }
  }

  std::vector<double> he0(ncombo);
  he_array(r, 0.0, he0.data());
  NumericVector total(ncombo);
  for (int k = 0; k < ncombo; k++) {
    double hm0 = he0[r - k] * he0[k];
    double tot = C * (2.0 * s[k] + (double)n * hm0);
    if (kfe) tot /= ((double)n * (double)n);
    total[k] = tot;
  }
  return total;
}

// In-place iterative radix-2 Cooley-Tukey FFT, n must be a power of 2.
// Matches R's fft() to floating-point precision (validated in
// data-raw/validate_fast_hpi.R) - used instead of calling fft() from R since
// R's own implementation is the actual bottleneck in ks:::kdde.binned.nd()
// for the grid sizes ks:::binning() produces.
static void fft1d(cd* a, int n, bool invert) {
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    double ang = 2 * M_PI / len * (invert ? 1 : -1);
    cd wlen(std::cos(ang), std::sin(ang));
    for (int i = 0; i < n; i += len) {
      cd w(1);
      for (int j2 = 0; j2 < len / 2; j2++) {
        cd u = a[i + j2], v = a[i + j2 + len / 2] * w;
        a[i + j2] = u + v;
        a[i + j2 + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) for (int i = 0; i < n; i++) a[i] /= n;
}

// data is a flat, row-major nr x nc buffer (data[i*nc+j] == element (i,j)) -
// a single contiguous allocation instead of nr separate row vectors. Rows
// are contiguous in this layout, so the row pass runs fft1d() in place with
// no copy; only the column pass (strided, not contiguous) needs a
// temporary buffer.
static void fft2d(cd* data, int nr, int nc, bool invert) {
  for (int i = 0; i < nr; i++) fft1d(data + i * nc, nc, invert);

  std::vector<cd> col(nr);
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) col[i] = data[i * nc + j];
    fft1d(col.data(), nr, invert);
    for (int i = 0; i < nr; i++) data[i * nc + j] = col[i];
  }
}

// R's round() is round-half-to-even (IEC 60559), unlike C's round() - needed
// for ks:::symconv.nd()'s placement offset below to match exactly.
static int r_round(double x) {
  double f = std::floor(x);
  double diff = x - f;
  if (diff < 0.5) return (int)f;
  if (diff > 0.5) return (int)f + 1;
  int fi = (int)f;
  return (fi % 2 == 0) ? fi : fi + 1;
}

//' Isotropic bivariate kernel functional estimate via linear binning + FFT
//' convolution
//'
//' Computes the same thing as ks:::kfe(x, G = g^2*diag(2), binned = TRUE,
//' deriv.order = r, deriv.vec = TRUE, add.index = FALSE) for r in {4, 6} -
//' the pilot-bandwidth kernel functional ks::Hpi() itself uses once n
//' exceeds ks:::default.bflag()'s ~500 threshold for 2D data. Mirrors
//' ks:::kdde.binned.nd() + ks:::symconv.nd() + the binned branch of
//' ks:::dmvnorm.deriv.sum() exactly: the kernel-grid evaluation (isotropic
//' closed form, same as cpp_kfe_isotropic_batch()) and the FFT convolution
//' are done here instead of via repeated R-level dmvnorm.deriv()/fft()
//' calls, which is what's actually slow in ks's own binned path (the linear
//' binning itself, via ks:::binning(), is already fast and is called from R
//' unchanged - see fast_kfe2d_binned() in R/internal_functions.R).
//' Validated in data-raw/validate_fast_hpi.R.
//'
//' @param counts Binned mass distribution (ks:::binning(x, H)$counts).
//' @param delta Bin grid spacing (dx, dy).
//' @param g,r Isotropic pilot bandwidth and derivative order.
//' @param n Sample size (sum of bin counts/weights).
//'
//' @noRd
// [[Rcpp::export]]
NumericVector cpp_kfe_isotropic_binned_batch(NumericMatrix counts, NumericVector delta,
                                              double g, int r, double n) {
  int M1 = counts.nrow(), M2 = counts.ncol();
  double lambda = g; // isotropic: sqrt(eigenvalues of g^2*I) = g
  int L1 = std::min((int)std::ceil((4 + r) * lambda / delta[0]), M1 - 1);
  int L2 = std::min((int)std::ceil((4 + r) * lambda / delta[1]), M2 - 1);
  int N1 = 2 * L1 - 1, N2 = 2 * L2 - 1;

  int P1 = 1; while (P1 < M1 + N1) P1 <<= 1;
  int P2 = 1; while (P2 < M2 + N2) P2 <<= 1;

  // ks:::symconv.nd() recomputes its own placement offset as
  // round(N/2)+1 - not the same as the L used to build the kernel offset
  // grid below, and can differ from it by parity.
  int Lc1 = r_round(N1 / 2.0) + 1;
  int Lc2 = r_round(N2 / 2.0) + 1;

  std::vector<cd> gcounts_pad(P1 * P2, cd(0, 0));
  for (int i = 0; i < M1; i++)
    for (int j = 0; j < M2; j++)
      gcounts_pad[(Lc1 - 1 + i) * P2 + (Lc2 - 1 + j)] = cd(counts(i, j), 0);
  fft2d(gcounts_pad.data(), P1, P2, false);

  double invg = 1.0 / g;
  double sign = (r % 2 == 0) ? 1.0 : -1.0;
  double C = sign * std::pow(g, -r) / (2.0 * M_PI * g * g);

  int ncombo = r + 1;
  NumericVector total(ncombo);
  int cell = P1 * P2;

  for (int k = 0; k < ncombo; k++) {
    int r1 = r - k, r2 = k;
    std::vector<cd> keval_pad(cell, cd(0, 0));

    for (int i = 0; i < N1; i++) {
      double u = delta[0] * (i - (L1 - 1)) * invg;
      double heu = he(r1, u);
      for (int j = 0; j < N2; j++) {
        double v = delta[1] * (j - (L2 - 1)) * invg;
        double val = C * heu * he(r2, v) * std::exp(-0.5 * (u * u + v * v)) / n;
        keval_pad[i * P2 + j] = cd(val, 0);
      }
    }
    fft2d(keval_pad.data(), P1, P2, false);

    std::vector<cd> prod(cell);
    for (int idx = 0; idx < cell; idx++)
      prod[idx] = keval_pad[idx] * gcounts_pad[idx];
    fft2d(prod.data(), P1, P2, true);

    // valid region: rows/cols [N1-1 .. N1-1+M1-1] x [N2-1 .. N2-1+M2-1],
    // then sum(counts*n*est) and an overall /n^2 (dmvnorm.deriv.sum()'s
    // kfe=TRUE normalization), folded together here.
    double s = 0.0;
    for (int i = 0; i < M1; i++)
      for (int j = 0; j < M2; j++)
        s += counts(i, j) * n * prod[(N1 - 1 + i) * P2 + (N2 - 1 + j)].real();
    total[k] = s / (n * n);
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
