# gsignal Standard Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# detrend()

test_that("parameters to detrend() are correct", {
  expect_error(detrend())
  expect_error(detrend('invalid'))
  expect_error(detrend(1:10, -1))
  expect_error(detrend(1:10, 'invalid'))
  expect_error(detrend(1:10, 0, 1))
})

test_that("detrend() tests are correct", {
  N <- 32
  x <- seq(0, N - 1, 1) / N + 2
  y <- detrend (x)
  expect_true(all(abs(y) < 1e-10))
  
  N <- 32
  t <- seq(0, N - 1, 1) / N
  x <- t * t + 2
  y <- detrend (x, 2)
  expect_true(all(abs(y) < 1e-10))
  
  N <- 32
  t <- seq(0, N - 1, 1) / N
  x <- cbind(t, 4 * t - 3)
  y <- detrend (x)
  expect_true(all(abs(y) < 1e-10))
})

