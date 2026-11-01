#' Expected Degree Depth
#'
#' Computes the expected degree depth of one or more query points with
#' respect to a multivariate data set, using the recursive finite-sample
#' identity rather than summing the individual k-hull depths.
#'
#' @param x A numeric vector or matrix containing the query point(s).
#'   If `x` is a matrix, each row represents one query point.
#' @param data A numeric matrix containing the data points, with one
#'   observation per row.
#'
#' @return A numeric vector containing the expected degree depth of
#'   each query point in `x`.
#'
#' @examples
#' n <- 100
#' d <- 2
#' data <- matrix(rnorm(n * d), ncol = d)
#' x <- matrix(0, nrow = 1, ncol = d)
#' ED(x, data)
#'
#' @export
ED <- function(x, data) {
  if (!is.matrix(data) || !is.numeric(data) || ncol(data) < 1L)
    stop("'data' must be a numeric matrix with at least one column.")

  d <- ncol(data)
  if (!is.numeric(x) || length(x) %% d != 0L)
    stop("'x' must contain a multiple of ncol(data) numeric values.")
  if (!is.null(dim(x)) && (!is.matrix(x) || ncol(x) != d))
    stop("Each row of 'x' must have ncol(data) coordinates.")

  # Preserve the previous behavior for numeric vectors and query matrices.
  x <- matrix(x, ncol = d)
  storage.mode(x) <- "double"
  storage.mode(data) <- "double"
  if (any(!is.finite(x)) || any(!is.finite(data)))
    stop("'x' and 'data' must contain only finite values.")

  depth <- numeric(nrow(x))
  for (i in seq_len(nrow(x))) {
    # The original package uses leave-one-out evaluation when the
    # query is exactly equal to a data observation.
    if (nrow(data) > 0L) {
      is_query <- rowSums(sweep(data, 2L, x[i, ], FUN = "==")) == d
      data_i <- data[!is_query, , drop = FALSE]
    } else {
      data_i <- data
    }
    n <- nrow(data_i)

    # C++ returns 1/T_d, with T_d = sum_k (1 - SD^(k)).
    # Apply the finite-sample correction exactly once.
    uncorrected <- ED_all_d_cpp(data_i, x[i, ])
    depth[i] <- max(0, uncorrected - 1 / (n + 1))
  }
  round(depth, 12)
}

# Internal compatibility helpers: they return uncorrected depth,
# just as the original ED_1d and ED_2d functions did.
ED_1d <- function(X, x) ED_all_d_cpp(as.matrix(X), as.numeric(x))
ED_2d <- function(X, x) ED_all_d_cpp(as.matrix(X), as.numeric(x))
