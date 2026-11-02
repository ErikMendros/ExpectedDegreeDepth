# Reference implementation by exhaustive enumeration (small n only).
contains_origin_hybrid <- function(A, x) {
  d <- ncol(A)
  k <- nrow(A)
  if (k < d + 1L) return(FALSE)
  for (ids in combn(k, d + 1L, simplify = FALSE)) {
    B <- sweep(A[ids, , drop = FALSE], 2L, x, "-")
    Z <- rbind(t(B), rep(1, d + 1L))
    w <- tryCatch(solve(Z, c(rep(0, d), 1)), error = function(e) NULL)
    if (!is.null(w) && all(w >= -1e-9)) return(TRUE)
  }
  FALSE
}

brute_ed_hybrid <- function(X, x) {
  d <- ncol(X)
  if (nrow(X)) {
    same <- rowSums(sweep(X, 2L, x, "==")) == d
    X <- X[!same, , drop = FALSE]
  }
  n <- nrow(X)
  T <- n + 1.0
  if (n >= d + 1L) {
    for (k in (d + 1L):n) {
      covering <- vapply(combn(n, k, simplify = FALSE), function(ids) {
        contains_origin_hybrid(X[ids, , drop = FALSE], x)
      }, logical(1))
      T <- T - mean(covering)
    }
  }
  1 / T - 1 / (n + 1)
}

test_that("ED agrees with exhaustive enumeration in dimensions 1-7", {
  set.seed(431)
  for (d in 1:7) {
    n <- d + 4L
    A <- matrix(rnorm(d*d), ncol = d)
    # The original 3D routine has coordinate-axis restrictions;
    # keep data generic when using d=3 directly.
    if (d %in% c(4L, 5L)) A[1, 1] <- 0
    X <- rbind(A, -colSums(A), matrix(rnorm((n-d-1)*d), ncol = d))
    x <- rep(0, d)
    expect_equal(unname(ED(x, X)), brute_ed_hybrid(X, x), tolerance = 1e-10)
  }
})

test_that("3D wrapper uses the original ED_3d_cpp", {
  set.seed(11)
  X <- matrix(rnorm(120), ncol = 3)
  origin <- c(0, 0, 0)
  expected <- max(0, ED_3d_cpp(X, origin) - 1/(nrow(X)+1))
  expect_equal(unname(ED(origin, X)), round(expected, 12), tolerance = 1e-11)
})

test_that("multiple queries and leave-one-out agree with separate calls", {
  set.seed(17)
  X <- matrix(rnorm(24), ncol = 3)
  Q <- rbind(rep(0, 3), X[1, ], rep(1.7, 3))
  expect_equal(unname(ED(Q, X)),
               vapply(seq_len(nrow(Q)), function(i) ED(Q[i, ], X), numeric(1)))
})

test_that("not enough observations gives zero depth", {
  for (d in 1:7) {
    X <- matrix(rnorm(d*d), ncol = d)
    expect_equal(unname(ED(rep(0, d), X)), 0)
  }
})
