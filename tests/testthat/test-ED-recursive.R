# Independent reference: for general-position data, any covering set
# contains a covering simplex of d+1 observations (Caratheodory).
contains_origin <- function(A, x) {
  d <- ncol(A)
  k <- nrow(A)
  if (k < d + 1L) return(FALSE)
  for (idx in combn(k, d + 1L, simplify = FALSE)) {
    Z <- rbind(t(sweep(A[idx, , drop = FALSE], 2L, x, FUN = "-")),
               rep(1, d + 1L))
    w <- tryCatch(solve(Z, c(rep(0, d), 1)), error = function(e) NULL)
    if (!is.null(w) && all(w >= -1e-9)) return(TRUE)
  }
  FALSE
}

brute_force_ed <- function(X, x) {
  d <- ncol(X)
  if (nrow(X)) {
    is_query <- rowSums(sweep(X, 2L, x, FUN = "==")) == d
    X <- X[!is_query, , drop = FALSE]
  }
  n <- nrow(X)
  T <- n + 1.0
  if (n >= d + 1L) {
    for (k in (d + 1L):n) {
      subsets <- combn(n, k, simplify = FALSE)
      covering <- vapply(subsets, function(ids)
        contains_origin(X[ids, , drop = FALSE], x), logical(1))
      T <- T - mean(covering)
    }
  }
  1 / T - 1 / (n + 1)
}

test_that("recursive ED agrees with exhaustive enumeration", {
  set.seed(431)
  for (d in 1:5) {
    n <- d + 3L
    A <- matrix(rnorm(d * d), ncol = d)
    if (d >= 2L) A[1, 1] <- 0  # axis-aligned but nonzero observation
    X <- rbind(A, -colSums(A), matrix(rnorm((n-d-1)*d), ncol=d))
    x <- rep(0, d)
    expect_equal(unname(ED(x, X)), brute_force_ed(X, x), tolerance = 1e-10)
    expect_equal(unname(ED(X[1, ], X)), brute_force_ed(X, X[1, ]),
                 tolerance = 1e-10)
  }
})

test_that("zero and outside-hull depths are zero", {
  for (d in 1:5) {
    X <- matrix(seq_len(d * d), ncol = d)  # n=d
    expect_equal(unname(ED(rep(0, d), X)), 0)
    X <- matrix(runif((d+2)*d, min = 1, max = 2), ncol = d)
    expect_equal(unname(ED(rep(0, d), X)), 0)
  }
})

test_that("multiple query points and leave-one-out agree with separate calls", {
  set.seed(41)
  X <- matrix(rnorm(24), ncol = 3)
  Q <- rbind(rep(0, 3), X[1, ], rep(2, 3))
  expect_equal(unname(ED(Q, X)),
               vapply(seq_len(nrow(Q)), function(i) ED(Q[i, ], X), numeric(1)))
})
