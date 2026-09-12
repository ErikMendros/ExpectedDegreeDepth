test_that("KHD accepts a vector query point", {

  set.seed(123)

  d <- 3
  n <- 12
  k <- 4

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- c(0, 0, 0)

  result <- KHD(x, data, k)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_true(is.finite(result))
})


test_that("KHD accepts multiple query points", {

  set.seed(123)

  d <- 3
  n <- 12
  m <- 4
  k <- 4

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- matrix(
    rnorm(m * d),
    nrow = m,
    ncol = d
  )

  result <- KHD(x, data, k)

  expect_length(result, m)
  expect_true(all(is.finite(result)))
})


test_that("KHD matrix evaluation equals separate evaluations", {

  set.seed(123)

  d <- 3
  n <- 12
  m <- 4
  k <- 4

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- matrix(
    rnorm(m * d),
    nrow = m,
    ncol = d
  )

  together <- KHD(x, data, k)

  separately <- vapply(
    seq_len(m),
    function(i) KHD(x[i, ], data, k),
    numeric(1)
  )

  expect_equal(
    together,
    separately,
    tolerance = 1e-10
  )
})


test_that("KHD uses leave-one-out when query belongs to data", {

  set.seed(123)

  d <- 3
  n <- 10
  k <- 4

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- data[1, ]

  result1 <- KHD(
    x,
    data,
    k
  )

  result2 <- KHD(
    x,
    data[-1, , drop = FALSE],
    k
  )

  expect_equal(
    result1,
    result2,
    tolerance = 1e-10
  )
})


test_that("KHD(data, data) gives leave-one-out depths", {

  set.seed(123)

  d <- 2
  n <- 8
  k <- 3

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  together <- KHD(
    data,
    data,
    k
  )

  separately <- vapply(
    seq_len(n),
    function(i) {
      KHD(
        data[i, ],
        data[-i, , drop = FALSE],
        k
      )
    },
    numeric(1)
  )

  expect_equal(
    together,
    separately,
    tolerance = 1e-10
  )
})


test_that("KHD returns values between zero and one", {

  set.seed(123)

  for (d in 1:5) {

    n <- max(15, d + 3)
    k <- d + 1

    data <- matrix(
      rnorm(n * d),
      ncol = d
    )

    x <- rep(0, d)

    result <- KHD(
      x,
      data,
      k
    )

    expect_true(is.finite(result))

    # Small numerical tolerance
    expect_true(result >= -1e-10)
    expect_true(result <= 1 + 1e-10)
  }
})
