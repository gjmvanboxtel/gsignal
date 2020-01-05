# gsignal Signal Measurement
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# findpeaks()

test_that("parameters to findpeaks() are correct", {
  expect_error(findpeaks())
  expect_error(findpeaks(1))
  expect_error(findpeaks(c(1, 1)))
  expect_error(findpeaks(complex(3, 1, 1)))
  expect_error(findpeaks(c(1, 1, 1), MinPeakHeight = -1))
  expect_error(findpeaks(c(1, 1, 1), MinPeakDistance = -1))
  expect_error(findpeaks(c(1, 1, 1), MinPeakWidth = -1))
  expect_error(findpeaks(c(1, 1, 1), DoubleSided = -1))
})

test_that("findpeaks() tests are correct", {

  expect_null(findpeaks (c(1, 1, 1)))
  expect_null(findpeaks (t(c(1, 1, 1))))

  # Test Matlab/Octave compatibility
  p <- findpeaks(c(34, 134, 353, 64, 134, 14, 56, 67, 234, 143, 64, 575, 8657))
  expect_equal(p$pks, c(353, 134, 234))

  
  ## Test for bug #45056
  ## Test input vector is an oversampled sinusoid with clipped peaks
  x <- pmin (3, cos (2*pi*c(0:8000) / 600) + 2.01)
  expect_equal(findpeaks(x)$pks, rep(3L, 27))
  
})

