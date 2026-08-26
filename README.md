# ExpectedDegreeDepth

An R package for computing the Expected Degree Depth (ED) of multivariate data points in any dimension.

## Installation

You can install the development version directly from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ErikMendros/ExpectedDegreeDepth")
```

## Usage

The main function is `ED()`.

It computes the expected degree depth of one or more query points with respect to a given data set. It is assumed that the data set is in general position. When a query point is also part of the data set, it is excluded from the reference data set when computing its depth. Thus, for `ED(data, data)`, the depth of each observation \(X_i\) is computed with respect to the data set excluding \(X_i\).

### Example

```r
library(ExpectedDegreeDepth)

set.seed(123)

# Generate a four-dimensional data set
n <- 50
d <- 4

data <- matrix(
  rnorm(n * d),
  ncol = d
)

# Query point
x <- c(0, 0, 0, 0)

# Compute expected degree depth
ED(x, data)
```

The query point can also be supplied as a matrix:

```r
x <- matrix(c(0, 0, 0, 0), nrow = 1)

ED(x, data)
```

Multiple query points can be evaluated simultaneously:

```r
x <- matrix(
  rnorm(10 * d),
  nrow = 10,
  ncol = d
)

ED(x, data)
```

The function returns one expected degree depth value for each query point.


## Dimensions

The implementation supports arbitrary dimensions.

Specialized implementations are used for dimensions 1, 2, and 3, while a general-dimensional implementation is used for dimensions greater than 3.


## Implementation

The package combines R and C++ implementations.

The computationally intensive parts are implemented in C++ using Rcpp.

The main R-level interface is:

```r
ED(x, data)
```

## Requirements

* R
* Rcpp

## Development

The package uses:

* Rcpp for the C++ implementation
* roxygen2 for documentation
* testthat for automated testing

To run the package tests:

```r
devtools::test()
```

To run the complete package check:

```r
devtools::check()
```

To regenerate documentation:

```r
devtools::document()
```

## Author

Erik Mendroš
