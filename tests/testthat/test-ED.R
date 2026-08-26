test_that("ED accepts a vector as a query point", {

  set.seed(123)

  d <- 3
  n <- 100

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- c(0, 0, 0)

  result <- ED(x, data)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_true(is.finite(result))
})


test_that("ED accepts a matrix of multiple query points", {

  set.seed(123)

  d <- 3
  n <- 100
  m <- 5

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- matrix(
    rnorm(m * d),
    nrow = m,
    ncol = d
  )

  result <- ED(x, data)

  expect_length(result, m)
  expect_type(result, "double")
  expect_true(all(is.finite(result)))
})


test_that("ED gives the same result for a vector and a one-row matrix", {

  set.seed(123)

  d <- 4
  n <- 30

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- c(0, 0, 0, 0)

  result_vector <- ED(x, data)

  result_matrix <- ED(
    matrix(x, nrow = 1),
    data
  )

  expect_equal(
    result_vector,
    result_matrix,
    tolerance = 1e-12
  )
})


test_that("ED handles a query point that is contained in the data", {

  set.seed(123)

  d <- 3
  n <- 100

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- data[1, ]

  result <- ED(x, data)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_true(is.finite(result))
})


test_that("ED handles a query point that is not contained in the data", {

  set.seed(123)

  d <- 4
  n <- 30

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- rep(100, d)

  result <- ED(x, data)

  expect_length(result, 1)
  expect_type(result, "double")
  expect_true(is.finite(result))
})


test_that("ED returns one value per query point", {

  set.seed(123)

  d <- 4
  n <- 30
  m <- 10

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- matrix(
    rnorm(m * d),
    nrow = m,
    ncol = d
  )

  result <- ED(x, data)

  expect_length(result, m)
  expect_true(all(is.finite(result)))
})


test_that("ED returns values between zero and one", {

  set.seed(123)

  d <- 5
  n <- 30
  m <- 10

  data <- matrix(
    rnorm(n * d),
    ncol = d
  )

  x <- matrix(
    rnorm(m * d),
    nrow = m,
    ncol = d
  )

  result <- ED(x, data)

  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})


test_that("ED works in dimensions one through six", {

  set.seed(123)

  for (d in 1:6) {

    n <- 20

    data <- matrix(
      rnorm(n * d),
      ncol = d
    )

    x <- matrix(
      rep(0, d),
      nrow = 1
    )

    result <- ED(x, data)

    expect_length(result, 1)
    expect_true(is.finite(result))
    expect_true(result >= 0)
    expect_true(result <= 1)
  }
})
