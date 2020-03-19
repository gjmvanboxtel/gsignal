# gsignal filtering functions
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

# -----------------------------------------------------------------------
# filtic()

test_that("parameters to filtic() are correct", {
  expect_error(filtic())
  expect_error(filtic(1))
  expect_error(filtic(1, 2))
  expect_error(filtic(1, 2, 3, 4, 5))
  expect_error(filtic(0, 0, 'invalid'))
})

test_that("filtic() tests are correct", {

  # Simple low pass filter
  b <- c(0.25, 0.25)
  a <- c(1.0, -0.5)
  expect_that(filtic(b, a, 1, 1), equals(0.75))
  
  # Simple high pass filter
  b <- c(0.25, -0.25)
  a <- c(1.0, 0.5)
  expect_that(filtic(b, a, 0, 1), equals(-0.25))

  # Second order cases
  # bs <- butter(2, 0.4)
  b <- c(0.2065721, 0.4131442, 0.2065721)
  a <- c(1.0000000, -0.3695274,  0.1958157)
  x <- y <- c(1, 1)
  expect_that(filtic(b, a, y, x), equals(c(0.7934280, 0.0107564), tolerance = 1e-7))
  N <- 1000;
  xx <- cos(2 * pi * seq(0, N-1, length.out = N)/8)
  yy <- filter(b, a, xx)
  x <- xx[seq(N, N - 1, -1)]
  y <- yy[seq(N, N - 1, -1)]
  zf <- filtic(b, a, y, x)
  expect_that(filtic(b, a, y, x), equals(c( 0.4039015, 0.1625113), tolerance = 1e-7))
  
})
  