# gsignal Internal functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# isScalar()

test_that("parameters to isScalar() are correct", {
  expect_error(isScalar())
  expect_error(isScalar(1, 2))
})

test_that("isScalar() returns TRUE or FALSE", {
  expect_equal(isScalar(1), TRUE)
  expect_equal(isScalar(c(1,2)), FALSE)
  expect_equal(isScalar(NULL), FALSE)
  expect_equal(isScalar(matrix(c(1,2,3,4),2,2)), FALSE)
  expect_equal(isScalar("t"), TRUE)
  expect_equal(isScalar("test"), FALSE)
  expect_equal(isScalar(c("test", "ing")), FALSE)
  expect_equal(isScalar(complex(real = 1, imaginary = 1)), TRUE)
})

# -----------------------------------------------------------------------
# isPosscal()

test_that("parameters to isPosscal() are correct", {
  expect_error(isPosscal())
  expect_error(isPosscal(1, 2))
})

test_that("isPosscal returns TRUE or FALSE", {
  expect_equal(isPosscal(1), TRUE)
  expect_equal(isPosscal(-1), FALSE)
  expect_equal(isPosscal(c(1,2)), FALSE)
  expect_equal(isPosscal(NULL), FALSE)
  expect_equal(isPosscal(matrix(c(1,2,3,4),2,2)), FALSE)
  expect_equal(isPosscal("t"), FALSE)
  expect_equal(isPosscal(complex(real = 1, imaginary = 1)), FALSE)
  expect_equal(isPosscal(Re(complex(real = 1, imaginary = 1))), TRUE)
})

# -----------------------------------------------------------------------
# isWhole()

test_that("parameters to isWhole() are correct", {
  expect_error(isWhole())
  expect_error(isWhole(1, 2, 3))
})

test_that("isWhole() returns TRUE or FALSE", {
  expect_equal(isWhole(1), TRUE)
  expect_equal(isWhole(-1), TRUE)
  expect_equal(isWhole(-12.5), FALSE)
  expect_equal(isWhole(c(1,2)), TRUE)
  expect_equal(isWhole(NULL), FALSE)
  expect_equal(isWhole(matrix(c(1,2,3,4),2,2)), TRUE)
  expect_equal(isWhole("t"), FALSE)
  expect_equal(isWhole(complex(real = 1, imaginary = 1)), TRUE)
  expect_equal(isWhole(complex(real = 1.1, imaginary = 1)), FALSE)
  expect_equal(isWhole(complex(real = 1, imaginary = 1.1)), FALSE)
  expect_equal(isWhole(Re(complex(real = 1, imaginary = 1))), TRUE)
})

# -----------------------------------------------------------------------
# unfactor()

test_that("parameters to unfactor() are correct", {
  expect_error(unfactor())
  expect_error(unfactor(1, 2))
  expect_warning(unfactor(as.factor("test")))
})

test_that("unfactor returns integer levels of factor", {
  expect_equal(unfactor(as.factor(c(1,2,3))), c(1,2,3))
  expect_equal(unfactor(1), NULL)
  expect_equal(unfactor("test"), NULL)
})

# -----------------------------------------------------------------------
# sinc()

test_that("parameters to sinc() are correct", {
  expect_error(sinc())
  expect_error(sinc(1, 2))
})

test_that("sinc() returns correct values", {
  expect_equal(sinc(0), 1L)
  expect_equal(round(sum(sinc(1:1e6)), 4), 0)
})
