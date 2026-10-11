#' Expected Degree Depth
#'
#' Computes the expected degree depth of one or more query points with
#' respect to a multivariate data set.
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
#'
#' ED(x, data)
#'
#' @export
ED <- function(x, data){

  d <- ncol(data)
  x <- matrix(x, ncol=d)
  m <- nrow(x)
  depth <- numeric(m)

  for(i in seq_len(m)){

    # Remove x[i,] from data if present.
    # Uses exact row matching.
    idx <- which(rowSums(t(t(data) == x[i,])) == d)

    if(length(idx) > 0){
      data_i <- data[-idx, , drop = FALSE]
    } else {
      data_i <- data
    }

    n <- nrow(data_i)

    # Compute depth depending on dimension.
    if(d == 1){
      depth[i] <- ED_1d(data_i, x[i,]) - 1/(n+1)
    } else if(d == 2){
      depth[i] <- ED_2d(data_i, x[i,]) - 1/(n+1)
    } else if(d == 3){
      depth[i] <- ED_3d_cpp(data_i, x[i,]) - 1/(n+1)
    } else {
      depth[i] <- ED_all_d_cpp(data_i, x[i,]) - 1/(n+1)
    }
  }

  round(depth, 12)
}


ED_1d <- function(X, x){

  n <- nrow(X)
  ell <- sum(X > x)

  return(
    1 / (
      2 +
        (ell * (ell - 1)) /
        (n * (n - ell + 1)) +
        ((n - ell) * ((n - ell) - 1)) /
        (n * (ell + 1))
    )
  )
}


ED_2d <- function(X, x) {

  X <- sweep(X, 2, x)
  X <- X[X[, 1] != 0, , drop = FALSE]

  cond <- X[, 2] > 0
  ratio <- -X[, 1] / X[, 2]
  s <- cond[order(ratio)]

  n <- length(s)

  if (n < 3L)
    stop("At least three observations are required.")

  # Precompute counts
  delta <- (s[-n] == 0) - (s[-1] == 1)

  ones_vec <- numeric(n)
  ones_vec[1] <- sum(s[-1])

  ones_vec[2:n] <- ones_vec[1] + cumsum(delta)

  ell <- ifelse(s == 0, (n - 1) - ones_vec, ones_vec)

  # Use double precision to avoid integer overflow
  N <- as.double(n)

  # ED: always computed
  ed <- 1 / (
    3 + sum(
      (3*N + 1 - 2*ell) * (ell - 1) * ell /
        ((N - 1) * N * (N - ell) * (N + 1 - ell))
    )
  )
  return(ed)
}
