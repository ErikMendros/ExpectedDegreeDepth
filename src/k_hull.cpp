
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

using namespace Rcpp;

namespace {

inline double ch(int n, int k) {
  if (k < 0 || k > n || n < 0) return 0.0;
  return R::choose(static_cast<double>(n), static_cast<double>(k));
}

std::vector<arma::uword> order_ratio(const arma::vec& ratio) {
  std::vector<arma::uword> ord(ratio.n_elem);
  std::iota(ord.begin(), ord.end(), 0);

  // R's order() places NaN/NA values last by default. stable_sort also
  // preserves the original order of ties.
  std::stable_sort(ord.begin(), ord.end(), [&](arma::uword a, arma::uword b) {
    const double xa = ratio[a];
    const double xb = ratio[b];
    const bool na = std::isnan(xa);
    const bool nb = std::isnan(xb);
    if (na != nb) return !na;  // finite/Inf before NaN
    if (na && nb) return false;
    return xa < xb;
  });

  return ord;
}

arma::mat shift_data(const arma::mat& X, const arma::vec& x) {
  arma::mat Z = X;
  Z.each_row() -= x.t();
  return Z;
}


double univariate_SDK_impl(const arma::mat& X,
                           const arma::vec& x,
                           int k) {
  const int n = static_cast<int>(X.n_rows);

  int L = 0;

  for (int i = 0; i < n; ++i) {
    if (X(i, 0) < x[0]) {
      ++L;
    }
  }

  const double denom = ch(n, k);

  return 1.0
  - ch(L, k) / denom
  - ch(n - L, k) / denom;
}


double bivariate_SDK_impl(const arma::mat& X, const arma::vec& x, int k) {
  arma::mat Z = shift_data(X, x);
  const int n = static_cast<int>(Z.n_rows);
  if (n == 0) return R_NaN;

  std::vector<int> cond(n);
  arma::vec ratio(n);

  for (int i = 0; i < n; ++i) {
    cond[i] = (Z(i, 1) > 0.0) ? 1 : 0;
    ratio[i] = -Z(i, 0) / Z(i, 1);
  }

  const std::vector<arma::uword> ord = order_ratio(ratio);
  std::vector<int> s(n);
  for (int i = 0; i < n; ++i) s[i] = cond[ord[i]];

  int init = 0;
  for (int i = 1; i < n; ++i) init += s[i];

  double SUM = ch((s[0] == 1) ? init : (n - 1 - init), k - 1);

  for (int j = 1; j < n; ++j) {
    if (s[j - 1] == 0) ++init;
    if (s[j] == 1) --init;

    const int side = (s[j] == 1) ? init : (n - 1 - init);
    SUM += ch(side, k - 1);
  }

  return 1.0 - SUM / ch(n, k);
}

double trivariate_SDK_impl(const arma::mat& X, const arma::vec& x, int k) {
  arma::mat Z = shift_data(X, x);
  const int n = static_cast<int>(Z.n_rows);
  if (n < 2) return R_NaN;

  arma::mat Y(n, 3);
  int n_up = 0;

  for (int i = 0; i < n; ++i) {
    Y(i, 0) = Z(i, 0) / Z(i, 2);
    Y(i, 1) = Z(i, 1) / Z(i, 2);
    Y(i, 2) = (Z(i, 2) > 0.0) ? 1.0 : ((Z(i, 2) < 0.0) ? -1.0 : 0.0);
    if (Z(i, 2) > 0.0) ++n_up;
  }

  const int n_down = n - n_up;
  double SDk = ch(n, k) - ch(n_up, k) - ch(n_down, k);
  arma::vec Counts(n - 1, arma::fill::zeros);

  for (int idx = 0; idx < n; ++idx) {
    const int m = n - 1;
    arma::vec ratio(m);
    std::vector<int> labels(m);
    std::vector<int> color(m);

    int p = 0;
    for (int r = 0; r < n; ++r) {
      if (r == idx) continue;

      const double X1 = Y(r, 0) - Y(idx, 0);
      const double X2 = Y(r, 1) - Y(idx, 1);
      const double X3 = Y(r, 2);

      color[p] = (Y(r, 2) != Y(idx, 2)) ? 1 : 0;
      labels[p] = (X2 * X3 > 0.0) ? 1 : 0;
      ratio[p] = -X1 / X2;
      ++p;
    }

    const std::vector<arma::uword> ord = order_ratio(ratio);
    std::vector<int> seq(m), clr(m);
    for (int i = 0; i < m; ++i) {
      seq[i] = labels[ord[i]];
      clr[i] = color[ord[i]];
    }

    int EL = 0;
    for (int i = 1; i < m; ++i) EL += seq[i];
    Counts[EL] += clr[0];

    for (int j = 0; j < n - 2; ++j) {
      if (seq[j] == 0) ++EL;
      if (seq[j + 1] == 1) --EL;
      Counts[EL] += clr[j + 1];
    }
  }

  double SUM = 0.0;
  for (int l = 1; l <= n - 1; ++l) {
    SUM += Counts[l - 1] * (ch(l - 1, k - 2) + ch(n - l - 1, k - 2));
  }

  SDk -= 0.25 * SUM;
  return SDk / ch(n, k);
}

arma::vec check_intersections_impl(const arma::mat& P, const arma::mat& Y) {
  const int d = static_cast<int>(P.n_cols);
  const int k = d - 2;

  if (static_cast<int>(P.n_rows) != k) {
    stop("P must have d-2 rows and d columns.");
  }
  if (static_cast<int>(Y.n_cols) != d) {
    stop("P and Y must have the same number of columns.");
  }

  arma::vec result(Y.n_rows, arma::fill::zeros);
  if (k == 0) {
    result.ones();
    return result;
  }

  arma::mat A = P.cols(0, k - 1).t();
  arma::mat A_inv;
  if (!arma::inv(A_inv, A)) {
    stop("Singular matrix encountered in check_intersections().");
  }

  for (arma::uword j = 0; j < Y.n_rows; ++j) {
    arma::vec u = Y.row(j).cols(0, k - 1).t();
    arma::vec tt = A_inv * (-u);

    bool ok = true;
    for (arma::uword r = 0; r < tt.n_elem; ++r) {
      if (tt[r] < -1e-10) {
        ok = false;
        break;
      }
    }
    result[j] = ok ? 1.0 : 0.0;
  }

  return result;
}

arma::vec HCounts_impl(const arma::mat& X, const arma::vec& x) {
  arma::mat Z = shift_data(X, x);
  const int n1 = static_cast<int>(Z.n_rows);
  const int d = static_cast<int>(Z.n_cols);
  const int r = d - 2;
  const int m = n1 - d + 2;

  if (d < 3) stop("HCounts() requires d >= 3.");
  if (r < 0 || r > n1 || m <= 0) {
    stop("Not enough points for HCounts() after filtering.");
  }

  arma::vec counts(m, arma::fill::zeros);
  std::vector<int> current(r);

  auto process_current = [&]() {
    arma::uvec cur(r);
    std::vector<char> selected(n1, 0);
    for (int i = 0; i < r; ++i) {
      cur[i] = static_cast<arma::uword>(current[i]);
      selected[current[i]] = 1;
    }

    arma::uvec rem(m);
    int pos = 0;
    for (int i = 0; i < n1; ++i) {
      if (!selected[i]) rem[pos++] = static_cast<arma::uword>(i);
    }

    arma::mat P = Z.rows(cur);
    arma::mat Yrem = Z.rows(rem);
    arma::vec vec01 = check_intersections_impl(P, Yrem);

    arma::mat Q, R;
    if (!arma::qr(Q, R, P.t())) {
      stop("QR decomposition failed in HCounts().");
    }
    if (Q.n_cols < static_cast<arma::uword>(d)) {
      stop("QR decomposition did not return a complete Q matrix.");
    }

    // R: Q[, (d-1):d] -> C++ zero-based columns d-2 and d-1.
    arma::mat A = Q.cols(d - 2, d - 1);
    arma::mat Y2 = Yrem * A;

    std::vector<int> cond(m);
    arma::vec ratio(m);
    for (int i = 0; i < m; ++i) {
      cond[i] = (Y2(i, 1) > 0.0) ? 1 : 0;
      ratio[i] = -Y2(i, 0) / Y2(i, 1);
    }

    const std::vector<arma::uword> ord = order_ratio(ratio);
    std::vector<int> s(m);
    std::vector<double> ord01(m);
    for (int i = 0; i < m; ++i) {
      s[i] = cond[ord[i]];
      ord01[i] = vec01[ord[i]];
    }

    int init = 0;
    for (int i = 1; i < m; ++i) init += s[i];
    counts[init] += ord01[0];

    for (int j = 1; j < m; ++j) {
      if (s[j - 1] == 0) ++init;
      if (s[j] == 1) --init;
      counts[init] += ord01[j];
    }
  };

  std::function<void(int, int)> generate = [&](int start, int depth) {
    if (depth == r) {
      process_current();
      return;
    }

    const int need = r - depth;
    for (int i = start; i <= n1 - need; ++i) {
      current[depth] = i;
      generate(i + 1, depth + 1);
    }
  };

  generate(0, 0);
  return counts;
}

double SDK_impl(const arma::mat& X, const arma::vec& x, int k) {
  const int n = static_cast<int>(X.n_rows);
  const int d = static_cast<int>(X.n_cols);

  if (static_cast<int>(x.n_elem) != d) {
    stop("length(x) must equal ncol(X).");
  }
  if (d < 1) stop("X must have at least one column.");

  if (d == 1) return univariate_SDK_impl(X, x, k);
  if (d == 2) return bivariate_SDK_impl(X, x, k);
  if (d == 3) return trivariate_SDK_impl(X, x, k);

  arma::mat Xproj = X.cols(0, d - 3);
  arma::vec xproj = x.subvec(0, d - 3);
  const double sdk_proj = SDK_impl(Xproj, xproj, k);

  arma::vec HC = HCounts_impl(X, x);
  const int expected = n - d + 2;
  if (static_cast<int>(HC.n_elem) != expected) {
    stop("After filtering shifted rows with first coordinate 0, HCounts has a different length than n-d+2. This mirrors a degeneracy in the original R code; check first-coordinate ties with x.");
  }

  double SUM = 0.0;
  for (int j = 1; j <= expected; ++j) {
    SUM += HC[j - 1] * (ch(j - 1, k - d + 1) + ch(n - d + 2 - j, k - d + 1));
  }

  return sdk_proj - SUM / (2.0 * (d - 1) * ch(n, k));
}

} // namespace



// [[Rcpp::export]]
double KHD_cpp(const arma::mat& X,
               const arma::vec& x,
               int k) {
  return SDK_impl(X, x, k);
}
