#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

// ------------------------------------------------------------
// Original specialized 3D implementation. Intentionally kept
// unchanged; called directly by R and by odd-dimensional recursion.
// ------------------------------------------------------------

// ============================================================
// ClosedForm
// ============================================================

inline double ClosedForm_cpp(double a, double n)
{
  return 2.0/(6 + 11*a + 6*a*a + a*a*a)
  + (2 + 2*a)/(a - n)
  - 10.0/(n - 1)
  + 10.0/n
  + 2*n/(6 + 11*a + 6*a*a + a*a*a)
  + (a + 2)/(n - 1 - a)
  + a/(n - a + 1);
}


// ============================================================
// ED_3d_cpp
// ============================================================

// [[Rcpp::export]]
double ED_3d_cpp(NumericMatrix Y, NumericVector x)
{
  int N = Y.nrow();
  int d = Y.ncol();

  // Shift data: X = sweep(Y, 2, x)
  NumericMatrix X(N, d);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < d; j++) {
      X(i,j) = Y(i,j) - x[j];
    }
  }

  // Remove rows where X[,1] == 0
  std::vector<int> keep;
  keep.reserve(N);

  for (int i = 0; i < N; i++) {
    if (X(i,0) != 0)
      keep.push_back(i);
  }

  int n = keep.size();

  NumericMatrix X2(n, d);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      X2(i,j) = X(keep[i], j);
    }
  }

  // Y <- cbind(X[,1]/X[,3], X[,2]/X[,3], sign(X[,3]))
  NumericMatrix Yp(n, 3);

  int n1 = 0;

  for (int i = 0; i < n; i++) {

    double z = X2(i,2);

    Yp(i,0) = X2(i,0) / z;
    Yp(i,1) = X2(i,1) / z;
    Yp(i,2) = (z > 0 ? 1 : -1);

    if (z > 0)
      n1++;
  }

  int n0 = n - n1;

  // Compute ED correction term
  double ed = n - 3
  - ((n+1.0)/(n+1.0-n1)
       - (1 + (double)n1/n
       + R::choose(n1,2)/R::choose(n,2)
       + R::choose(n1,3)/R::choose(n,3)))
       - ((n+1.0)/(n+1.0-n0)
       - (1 + (double)n0/n
       + R::choose(n0,2)/R::choose(n,2)
       + R::choose(n0,3)/R::choose(n,3)));

       // Counts vector
       std::vector<double> Counts(n-1, 0.0);

       // Main loop
       for (int idx = 0; idx < n; idx++) {

         std::vector<double> X1, X2v, X3;
         std::vector<int> CLR, Labels;

         X1.reserve(n-1);
         X2v.reserve(n-1);
         X3.reserve(n-1);
         CLR.reserve(n-1);
         Labels.reserve(n-1);

         for (int j = 0; j < n; j++) {

           if (j == idx)
             continue;

           double dx1 = Yp(j,0) - Yp(idx,0);
           double dx2 = Yp(j,1) - Yp(idx,1);
           double dx3 = Yp(j,2);

           X1.push_back(dx1);
           X2v.push_back(dx2);
           X3.push_back(dx3);

           CLR.push_back(dx3 != Yp(idx,2));
           Labels.push_back(dx2 * dx3 > 0);
         }

         int m = n - 1;

         // Precompute ratio to avoid repeated division in sort comparator
         std::vector<double> ratio(m);

         for (int i = 0; i < m; i++) {
           ratio[i] = -X1[i] / X2v[i];
         }

         // Order = order(ratio)
         std::vector<int> ord(m);

         for (int i = 0; i < m; i++)
           ord[i] = i;

         std::sort(
           ord.begin(),
           ord.end(),
           [&](int a, int b) {
             return ratio[a] < ratio[b];
           }
         );

         // Apply ordering
         std::vector<int> Seq(m), Color(m);

         for (int i = 0; i < m; i++) {
           Seq[i] = Labels[ord[i]];
           Color[i] = CLR[ord[i]];
         }

         int EL = 0;

         for (int i = 1; i < m; i++)
           EL += Seq[i];

         Counts[EL] += Color[0];

         for (int j = 0; j < m-1; j++) {

           if (Seq[j] == 0)
             EL++;

           if (Seq[j+1] == 1)
             EL--;

           Counts[EL] += Color[j+1];
         }
       }

       // SUM = Σ Counts[l] * ClosedForm(l-1, n)
       double SUM = 0.0;

       for (int l = 1; l < n; l++) {
         SUM += Counts[l-1] * ClosedForm_cpp(l-1, n);
       }

       ed -= 0.25 * SUM;

       return 1.0 / (1 + n - ed);
}



// ============================================================
// Expected degree depth is evaluated using the sample recursion.
// For d >= 4, we do not aggregate the k-hull depths.
// ============================================================

// ============================================================
// Binomial coefficient
// ============================================================

inline double choose_cpp(int n, int k)
{
  if (k < 0 || k > n || n < 0)
    return 0.0;

  if (k == 0 || k == n)
    return 1.0;

  k = std::min(k, n-k);

  double ans = 1.0;

  for (int i = 1; i <= k; ++i)
    ans *= (double)(n-k+i) / (double)i;

  return ans;
}


// ============================================================
// Solve linear system A x = b by Gaussian elimination
// with partial pivoting.
// ============================================================

bool solve_linear_system(
    std::vector<double> A,
    std::vector<double> b,
    int n,
    std::vector<double>& x
)
{
  x.assign(n, 0.0);

  // Augmented matrix
  std::vector<double> M(n*(n+1));

  for (int i = 0; i < n; ++i) {

    for (int j = 0; j < n; ++j)
      M[i*(n+1)+j] = A[i*n+j];

    M[i*(n+1)+n] = b[i];
  }

  // Elimination
  for (int col = 0; col < n; ++col) {

    int pivot = col;

    double best =
      std::abs(M[col*(n+1)+col]);

    for (int i = col+1; i < n; ++i) {

      double val =
        std::abs(M[i*(n+1)+col]);

      if (val > best) {
        best = val;
        pivot = i;
      }
    }

    if (best < 1e-14)
      return false;

    if (pivot != col) {

      for (int j = col; j <= n; ++j) {

        std::swap(
          M[col*(n+1)+j],
           M[pivot*(n+1)+j]
        );
      }
    }

    for (int i = col+1; i < n; ++i) {

      double factor =
        M[i*(n+1)+col] /
          M[col*(n+1)+col];

      if (factor == 0.0)
        continue;

      for (int j = col; j <= n; ++j) {

        M[i*(n+1)+j] -=
          factor * M[col*(n+1)+j];
      }
    }
  }

  // Back substitution
  for (int i = n-1; i >= 0; --i) {

    double rhs =
      M[i*(n+1)+n];

    for (int j = i+1; j < n; ++j) {

      rhs -=
        M[i*(n+1)+j] * x[j];
    }

    x[i] =
      rhs / M[i*(n+1)+i];
  }

  return true;
}


// ============================================================
// Complete Householder QR
// ============================================================

NumericMatrix complete_QR_Q(
    const NumericMatrix& M
)
{
  int d = M.nrow();
  int r = M.ncol();

  NumericMatrix R(M);
  NumericMatrix Q(d,d);

  for (int i = 0; i < d; ++i) {

    for (int j = 0; j < d; ++j)

      Q(i,j) =
        (i == j ? 1.0 : 0.0);
  }

  for (int k = 0; k < r; ++k) {

    double normx = 0.0;

    for (int i = k; i < d; ++i)

      normx =
        std::hypot(normx, R(i,k));

    if (normx < 1e-15)
      continue;

    double alpha =
      (R(k,k) >= 0.0)
      ? -normx
    : normx;

    std::vector<double> v(d-k);

    for (int i = k; i < d; ++i)

      v[i-k] =
        R(i,k);

    v[0] -= alpha;

    double vnorm2 = 0.0;

    for (double z : v)
      vnorm2 += z*z;

    if (vnorm2 < 1e-30)
      continue;

    for (int j = k; j < r; ++j) {

      double dot = 0.0;

      for (int i = k; i < d; ++i)

        dot +=
          v[i-k] * R(i,j);

      double factor =
        2.0 * dot / vnorm2;

      for (int i = k; i < d; ++i)

        R(i,j) -=
          factor * v[i-k];
    }

    for (int i = 0; i < d; ++i) {

      double dot = 0.0;

      for (int j = k; j < d; ++j)

        dot +=
          Q(i,j) * v[j-k];

      double factor =
        2.0 * dot / vnorm2;

      for (int j = k; j < d; ++j)

        Q(i,j) -=
          factor * v[j-k];
    }
  }

  return Q;
}


// ============================================================
// BIVARIATE SDK (CORRECTED)
// ============================================================

NumericVector bivariate_SDK_cpp(
    const NumericMatrix& X,
    const IntegerVector& K
)
{
  int n = X.nrow();
  int nk = K.size();

  NumericVector result(nk);

  std::vector<double> ratio(n);
  std::vector<int> cond(n);

  for (int i = 0; i < n; ++i) {

    cond[i] =
      (X(i,1) > 0.0) ? 1 : 0;

    ratio[i] =
      -X(i,0) / X(i,1);
  }

  std::vector<int> ord(n);

  for (int i = 0; i < n; ++i)
    ord[i] = i;

  std::sort(
    ord.begin(),
    ord.end(),
    [&](int a, int b) {
      return ratio[a] < ratio[b];
    }
  );

  std::vector<int> s(n);

  for (int i = 0; i < n; ++i)
    s[i] = cond[ord[i]];

  int init = 0;

  for (int i = 1; i < n; ++i)
    init += s[i];

  std::vector<double> counts(n, 0.0);

  counts[0] =
    (s[0] == 1)
    ? init
  : (n - 1 - init);

  for (int j = 1; j < n; ++j) {

    if (s[j-1] == 0)
      ++init;

    if (s[j] == 1)
      --init;

    counts[j] =
    (s[j] == 1)
      ? init
    : (n - 1 - init);
  }

  for (int q = 0; q < nk; ++q) {

    int k = K[q];

    double SUM = 0.0;

    for (int ell = 0; ell < n; ++ell)

      SUM +=
        choose_cpp(
          (int)counts[ell],
                     k-1
        );

    result[q] =
      1.0 -
      SUM / choose_cpp(n,k);
  }

  return result;
}


// ============================================================
// TRIVARIATE SDK
// ============================================================

NumericVector trivariate_SDK_cpp(
    const NumericMatrix& X,
    const IntegerVector& K
)
{
  int n = X.nrow();
  int nk = K.size();

  NumericVector result(nk);
  NumericMatrix Y(n,3);

  int n_up = 0;

  for (int i = 0; i < n; ++i) {

    double z = X(i,2);

    Y(i,0) =
      X(i,0) / z;

    Y(i,1) =
      X(i,1) / z;

    Y(i,2) =
      (z > 0.0 ? 1.0 : -1.0);

    if (z > 0.0)
      ++n_up;
  }

  int n_down = n-n_up;

  std::vector<double> Counts(n-1, 0.0);

  for (int idx = 0; idx < n; ++idx) {

    int m = n-1;

    std::vector<double> X1(m);
    std::vector<double> X2(m);
    std::vector<double> X3(m);

    std::vector<int> CLR(m);
    std::vector<int> Labels(m);

    int pos = 0;

    for (int j = 0; j < n; ++j) {

      if (j == idx)
        continue;

      X1[pos] =
        Y(j,0) - Y(idx,0);

      X2[pos] =
        Y(j,1) - Y(idx,1);

      X3[pos] =
        Y(j,2);

      CLR[pos] =
        (Y(j,2) != Y(idx,2));

      Labels[pos] =
        (X2[pos] * X3[pos] > 0.0);

      ++pos;
    }

    std::vector<double> ratio(m);

    for (int i = 0; i < m; ++i)

      ratio[i] =
        -X1[i] / X2[i];

        std::vector<int> ord(m);

        for (int i = 0; i < m; ++i)
          ord[i] = i;

        std::sort(
          ord.begin(),
          ord.end(),
          [&](int a, int b) {
            return ratio[a] < ratio[b];
          }
        );

        std::vector<int> Seq(m);
        std::vector<int> Color(m);

        for (int i = 0; i < m; ++i) {

          Seq[i] =
            Labels[ord[i]];

          Color[i] =
            CLR[ord[i]];
        }

        int EL = 0;

        for (int i = 1; i < m; ++i)
          EL += Seq[i];

        Counts[EL] += Color[0];

        for (int j = 0; j < m-1; ++j) {

          if (Seq[j] == 0)
            ++EL;

          if (Seq[j+1] == 1)
            --EL;

          Counts[EL] += Color[j+1];
        }
  }

  for (int q = 0; q < nk; ++q) {

    int k = K[q];

    double base =
      choose_cpp(n,k)
      - choose_cpp(n_up,k)
      - choose_cpp(n_down,k);

      double SUM = 0.0;

      for (int ell = 0; ell <= n-2; ++ell) {

        double term =
          choose_cpp(ell,k-2)
        + choose_cpp(n-2-ell,k-2);

        SUM +=
        Counts[ell] * term;
      }

      result[q] =
        (base - 0.25 * SUM)
        / choose_cpp(n,k);
  }

  return result;
}


// ============================================================
// check_intersections
// ============================================================

std::vector<int> check_intersections_cpp(
    const NumericMatrix& P,
    const NumericMatrix& Y
)
{
  // The (d-2)x(d-2) matrix is constant for all remaining points.
  // Factor it once instead of performing a fresh Gaussian elimination
  // (and allocating 3+ vectors) for each observation.
  const int k = P.nrow();
  const int m = Y.nrow();
  std::vector<int> result(m, 0);
  std::vector<double> lu(k*k);
  std::vector<int> piv(k);
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j)
      lu[i*k+j] = P(j,i);

  for (int col = 0; col < k; ++col) {
    int row = col;
    double maxabs = std::abs(lu[col*k+col]);
    for (int i = col+1; i < k; ++i) {
      const double cur = std::abs(lu[i*k+col]);
      if (cur > maxabs) { maxabs = cur; row = i; }
    }
    if (maxabs < 1e-14) return result;
    piv[col] = row;
    if (row != col)
      for (int j = 0; j < k; ++j)
        std::swap(lu[col*k+j], lu[row*k+j]);
    for (int i = col+1; i < k; ++i) {
      lu[i*k+col] /= lu[col*k+col];
      const double z = lu[i*k+col];
      for (int j = col+1; j < k; ++j)
        lu[i*k+j] -= z*lu[col*k+j];
    }
  }
  std::vector<double> b(k);
  for (int p = 0; p < m; ++p) {
    for (int i = 0; i < k; ++i) b[i] = -Y(p,i);
    for (int col = 0; col < k; ++col)
      std::swap(b[col], b[piv[col]]);
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < i; ++j)
        b[i] -= lu[i*k+j] * b[j];
    bool inside = true;
    for (int i = k-1; i >= 0; --i) {
      for (int j = i+1; j < k; ++j)
        b[i] -= lu[i*k+j] * b[j];
      b[i] /= lu[i*k+i];
      if (b[i] < -1e-10) inside = false;
    }
    result[p] = inside ? 1 : 0;
  }
  return result;
}


// ============================================================
// HCounts
// ============================================================

void HCounts_recursive(
    const NumericMatrix& X,
    int d,
    int target,
    int start,
    int depth,
    std::vector<int>& current,
    std::vector<double>& counts
)
{
  int n = X.nrow();

  if (depth == target) {

    NumericMatrix P(target,d);
    std::vector<char> selected(n,0);

    for (int i = 0; i < target; ++i) {

      selected[current[i]] = 1;

      for (int j = 0; j < d; ++j)

        P(i,j) =
          X(current[i],j);
    }

    int m = n-target;

    NumericMatrix Yremaining(m,d);

    int pos = 0;

    for (int i = 0; i < n; ++i) {

      if (selected[i])
        continue;

      for (int j = 0; j < d; ++j)

        Yremaining(pos,j) =
          X(i,j);

      ++pos;
    }

    std::vector<int> vec01 =
      check_intersections_cpp(
        P,
        Yremaining
      );

    NumericMatrix PT(d,target);

    for (int i = 0; i < target; ++i)

      for (int j = 0; j < d; ++j)

        PT(j,i) =
          P(i,j);

    NumericMatrix Q =
      complete_QR_Q(PT);

    NumericMatrix A(d,2);

    for (int i = 0; i < d; ++i) {

      A(i,0) =
        Q(i,d-2);

      A(i,1) =
        Q(i,d-1);
    }

    NumericMatrix Y(m,2);

    for (int i = 0; i < m; ++i) {

      Y(i,0) = 0.0;
      Y(i,1) = 0.0;

      for (int j = 0; j < d; ++j) {

        Y(i,0) +=
          Yremaining(i,j) * A(j,0);

        Y(i,1) +=
          Yremaining(i,j) * A(j,1);
      }
    }

    std::vector<int> cond(m);
    std::vector<double> ratio(m);

    for (int i = 0; i < m; ++i) {

      cond[i] =
        (Y(i,1) > 0.0) ? 1 : 0;

      ratio[i] =
        -Y(i,0) / Y(i,1);
    }

    std::vector<int> ord(m);

    for (int i = 0; i < m; ++i)
      ord[i] = i;

    std::sort(
      ord.begin(),
      ord.end(),
      [&](int a, int b) {
        return ratio[a] < ratio[b];
      }
    );

    std::vector<int> s(m);
    std::vector<int> ord01(m);

    for (int i = 0; i < m; ++i) {

      s[i] =
        cond[ord[i]];

      ord01[i] =
        vec01[ord[i]];
    }

    int init = 0;

    for (int i = 1; i < m; ++i)

      init += s[i];

    counts[init] +=
      ord01[0];

    for (int j = 1; j < m; ++j) {

      if (s[j-1] == 0)
        ++init;

      if (s[j] == 1)
        --init;

      counts[init] +=
      ord01[j];
    }

    return;
  }

  int max_i =
    n - (target-depth);

  for (int i = start; i <= max_i; ++i) {

    current[depth] = i;

    HCounts_recursive(
      X,
      d,
      target,
      i+1,
      depth+1,
      current,
      counts
    );
  }
}


std::vector<double> HCounts_cpp(
    const NumericMatrix& X
)
{
  int n = X.nrow();
  int d = X.ncol();

  int target = d-2;

  std::vector<double> counts(
      n-d+2,
      0.0
  );

  std::vector<int> current(target);

  HCounts_recursive(
    X,
    d,
    target,
    0,
    0,
    current,
    counts
  );

  return counts;
}


// ============================================================
// GENERAL SDK
// ============================================================

NumericVector SDK_cpp(
    const NumericMatrix& X,
    const IntegerVector& K
)
{
  int n = X.nrow();
  int d = X.ncol();
  int nk = K.size();

  NumericVector result(nk);

  if (d == 2)
    return bivariate_SDK_cpp(X,K);

  if (d == 3)
    return trivariate_SDK_cpp(X,K);

  NumericMatrix Xproj(n,d-2);

  for (int i = 0; i < n; ++i)

    for (int j = 0; j < d-2; ++j)

      Xproj(i,j) =
        X(i,j);

  NumericVector sdk_proj =
    SDK_cpp(Xproj,K);

  std::vector<double> HC =
    HCounts_cpp(X);

  int m = n-d+2;

  for (int q = 0; q < nk; ++q) {

    int k = K[q];

    int r = k-d+1;

    double SUM = 0.0;

    for (int j = 1; j <= m; ++j) {

      double term =
        choose_cpp(j-1,r)
      + choose_cpp(m-j,r);

      SUM +=
      HC[j-1] * term;
    }

    result[q] =
      sdk_proj[q]
    - SUM /
      (2.0 * (d-1)
         * choose_cpp(n,k));
  }

  return result;
}


// ============================================================
// Reciprocal binomial coefficient. Avoid forming binomial values,
// which can exceed the floating-point range for moderate n.
// ============================================================

inline double inverse_choose_ed(int m, int k) {
  if (k < 0 || k > m) return 0.0;
  k = std::min(k, m-k);
  double result = 1.0;
  for (int j = 1; j <= k; ++j)
    result *= static_cast<double>(j) / (m-k+j);
  return result;
}

// ============================================================
// Uncorrected reciprocal ED, T_d = sum_{k=0}^n (1 - SD^(k)).
// The final finite-sample correction is applied by the R wrapper.
// X is already centered at the query point.
// ============================================================

static double ed_reciprocal_recursive(const NumericMatrix& X) {
  const int n = X.nrow();
  const int d = X.ncol();

  // In general position, no subset of <=d observations can
  // contain the origin. This also handles an empty sample.
  if (n <= d) return static_cast<double>(n) + 1.0;

  if (d == 1) {
    int nplus = 0;
    for (int i = 0; i < n; ++i)
      nplus += (X(i,0) > 0.0);
    const int nminus = n - nplus;
    return (n+1.0)/(nplus+1.0)
      + (n+1.0)/(nminus+1.0) - 1.0;
  }

  if (d == 3) {
    // Stop here in odd dimension: reuse the ORIGINAL specialized 3D
    // implementation rather than projecting once more to dimension 1.
    // ED_3d_cpp returns 1/T_3, so invert to obtain T_3.
    NumericVector origin(3);
    return 1.0 / ED_3d_cpp(X, origin);
  }

  if (d == 2) {
    // For each directed ray 0 -> X_i, count the number of rays
    // strictly to its left using a circular angular sweep.
    const double pi = std::acos(-1.0);
    std::vector<double> angles(n);
    for (int i = 0; i < n; ++i) {
      double a = std::atan2(X(i,1), X(i,0));
      if (a < 0.0) a += 2.0*pi;
      angles[i] = a;
    }
    std::sort(angles.begin(), angles.end());
    std::vector<double> twice(2*n);
    for (int i = 0; i < n; ++i) {
      twice[i] = angles[i];
      twice[i+n] = angles[i] + 2.0*pi;
    }
    double total = 1.0;
    int j = 1;
    for (int i = 0; i < n; ++i) {
      if (j < i+1) j = i+1;
      while (j < i+n && twice[j] - twice[i] < pi) ++j;
      const int ell = j-i-1;
      const double a = n-ell;
      total += (n+1.0)/(a*(a+1.0));
    }
    return total;
  }

  // For d >= 4, project onto the orthogonal complement of
  // psi = span(e_{d-1}, e_d), then apply the sample ED recursion.
  NumericMatrix projected(n, d-2);
  for (int i = 0; i < n; ++i)
    for (int a = 0; a < d-2; ++a)
      projected(i,a) = X(i,a);

  double total = ed_reciprocal_recursive(projected);
  const std::vector<double> counts = HCounts_cpp(X);

  // HCounts_cpp counts each (d-1)-subset exactly (d-1) times,
  // once for each choice of the (d-2)-element anchor subset.
  // counts[ell] is the multiplicity of the halfspace count ell.
  const int maxell = n-d+1;
  double weighted = 0.0;
  for (int ell = 0; ell <= maxell; ++ell) {
    const double inv1 = inverse_choose_ed(n-ell+1, d);
    const double inv2 = inverse_choose_ed(d+ell, d);
    weighted += counts[ell]*(inv1 + inv2);
  }
  return total + (n+1.0)/(2.0*d*(d-1.0))*weighted;
}

// [[Rcpp::export]]
double ED_all_d_cpp(NumericMatrix X, NumericVector x) {
  const int nrows = X.nrow();
  const int d = X.ncol();
  if (d < 1 || x.size() != d)
    Rcpp::stop("Dimensions of query point and data do not match.");

  // Translate to the origin. Remove observations exactly equal to
  // the query point, matching the existing leave-one-out R wrapper.
  std::vector<double> kept;
  kept.reserve(static_cast<size_t>(nrows)*d);
  for (int i = 0; i < nrows; ++i) {
    bool is_query = true;
    for (int a = 0; a < d; ++a) {
      if (X(i,a) != x[a]) { is_query = false; break; }
    }
    if (is_query) continue;
    for (int a = 0; a < d; ++a)
      kept.push_back(X(i,a) - x[a]);
  }
  const int n = static_cast<int>(kept.size()/d);
  NumericMatrix centered(n,d);
  for (int i = 0; i < n; ++i)
    for (int a = 0; a < d; ++a)
      centered(i,a) = kept[i*d+a];

  // For d >= 4, an orthogonal change of coordinates makes the fixed
  // reference planes used in the recursion non-special. In particular,
  // the first coordinate of an observation may be zero without making
  // the projection degenerate. The ED is orthogonally invariant.
  // Normalize rows first: radial scaling of each observation by a
  // positive constant also leaves all hull-containment events intact.
  if (d >= 4 && n > d) {
    for (int i = 0; i < n; ++i) {
      double length = 0.0;
      for (int a = 0; a < d; ++a)
        length = std::hypot(length, centered(i,a));
      if (length == 0.0)
        Rcpp::stop("Zero observation remained after centering.");
      for (int a = 0; a < d; ++a)
        centered(i,a) /= length;
    }
    // A fixed sequence of Givens rotations. This is deterministic and
    // has negligible cost relative to halfspace counting for fixed d.
    for (int a = 0; a < d; ++a) {
      for (int b = a+1; b < d; ++b) {
        const double phi = 0.413 + 0.117*(a+1) + 0.071*(b+1);
        const double cs = std::cos(phi), sn = std::sin(phi);
        for (int i = 0; i < n; ++i) {
          const double u = centered(i,a), v = centered(i,b);
          centered(i,a) = cs*u - sn*v;
          centered(i,b) = sn*u + cs*v;
        }
      }
    }
  }

  const double T = ed_reciprocal_recursive(centered);
  return 1.0/T;
}
