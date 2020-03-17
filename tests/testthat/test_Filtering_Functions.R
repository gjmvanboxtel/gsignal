# gsignal filteringn functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# filtfilt()

test_that("parameters to filtfilt() are correct", {
  expect_error(filtfilt())
  expect_error(filtfilt(1, 2, 3))
  expect_error(filtfilt(0, 0, 1:10))
  expect_error(filtfilt(1, 2, array(1:8, c(2, 2, 2))))
  expect_error(filtfilt(1, 1, c('invalid', 'invalid')))
})

test_that("filtfilt() tests are correct", {
  expect_that(filtfilt(1, 1, 1:2), equals(c(1,2)))
  expect_that(filtfilt(1, 2, 1:2), equals(c(0.25, 0.50)))
  expect_that(filtfilt(2, 1, 1:2), equals(c(4, 8)))
  x <- runif(100)
  y <- filtfilt(1, 1, x)
  expect_that(length(y), equals(length(x)))
  x <- matrix(runif(200), 100, 2)
  y <- filtfilt(1, 1, x)
  expect_that(ncol(y), equals(ncol(x)))
  expect_that(nrow(y), equals(nrow(x)))
})
  