#' K-Hull Depth
#'
#' Computes the k-hull depth of one or more query points with respect
#' to a multivariate data set.
#'
#' For a d-dimensional data set, setting `k = d + 1` gives the
#' simplicial depth.
#'
#' When a query point is also contained in `data`, that observation is
#' excluded from the reference data set before its depth is computed.
#'
#' @param x A numeric vector or matrix containing the query point(s).
#'   If `x` is a matrix, each row represents one query point.
#' @param data A numeric matrix containing the data points, with one
#'   observation per row.
#' @param k Integer giving the number of points used to form the
#'   random convex hull. For nondegenerate d-dimensional data,
#'   `k` should satisfy `d + 1 <= k <= n`.
#'
#' @return A numeric vector containing the k-hull depth of each
#'   query point in `x`.
#'
#' @details
#' The k-hull depth is the probability that a query point is contained
#' in the convex hull of `k` observations drawn from the reference
#' data set.
#'
#' For `k = d + 1`, where `d` is the dimension of the data,
#' this is the simplicial depth.
#'
#' @examples
#' set.seed(123)
#'
#' n <- 30
#' d <- 3
#'
#' data <- matrix(rnorm(n * d), ncol = d)
#' x <- c(0, 0, 0)
#'
#' # Simplicial depth
#' KHD(x, data, k = d + 1)
#'
#' # 6-hull depth
#' KHD(x, data, k = 6)
#'
#' @export
KHD <- function(x, data, k) {

  d <- ncol(data)

  if (is.null(d)) {
    stop("`data` must be a matrix.")
  }

  if (!is.numeric(data)) {
    stop("`data` must be numeric.")
  }

  # Same convention as ED():
  # vector -> one-row matrix,
  # matrix -> each row is one query point.
  if (is.null(dim(x))) {

    if (length(x) != d) {
      stop("Length of `x` must equal ncol(data).")
    }

    x <- matrix(x, nrow = 1)

  } else {

    x <- as.matrix(x)

    if (ncol(x) != d) {
      stop("`x` and `data` must have the same number of columns.")
    }
  }

  if (!is.numeric(x)) {
    stop("`x` must be numeric.")
  }

  if (
    length(k) != 1 ||
    !is.finite(k) ||
    k != floor(k)
  ) {
    stop("`k` must be a single integer.")
  }

  k <- as.integer(k)

  if (k < d + 1L) {
    stop("`k` must be at least d + 1.")
  }

  m <- nrow(x)
  depth <- numeric(m)

  for (i in seq_len(m)) {

    # Same leave-one-out rule as ED().
    # Remove exact matches of the current query point.
    idx <- which(
      rowSums(
        t(t(data) == x[i, ])
      ) == d
    )

    if (length(idx) > 0) {
      data_i <- data[-idx, , drop = FALSE]
    } else {
      data_i <- data
    }

    n_i <- nrow(data_i)

    if (k > n_i) {
      stop(
        sprintf(
          "`k` is %d but query point %d has only %d reference observations.",
          k,
          i,
          n_i
        )
      )
    }

    depth[i] <- KHD_cpp(
      data_i,
      x[i, ],
      k
    )
  }

  round(depth, 12)
}
