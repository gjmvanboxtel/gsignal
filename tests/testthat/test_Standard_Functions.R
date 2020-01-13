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

# -----------------------------------------------------------------------
# ifft() and imvfft()

test_that("parameters to ifft() are correct", {
  expect_error(ifft())
  expect_error(ifft('invalid'))
  expect_error(ifft(1, 2))
})

test_that("ifft() tests are correct", {
  expect_equal(ifft(stats::fft(1:5)), 1:5)
  expect_equal(ifft(stats::fft(c(1+5i, 2+3i, 3+2i, 4+6i, 5+2i))), c(1+5i, 2+3i, 3+2i, 4+6i, 5+2i))
  expect_equal(imvfft(stats::mvfft(matrix(1:20, 4, 5))), matrix(1:20, 4, 5))
})

