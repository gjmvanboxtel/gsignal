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

# -----------------------------------------------------------------------
# peak2peak()

test_that("parameters to peak2peak() are correct", {
  expect_error(peak2peak())
  expect_error(peak2peak('invalid'))
  expect_error(peak2peak(1, 2, 3))
  expect_error(peak2peak(1, 1.5))
  expect_error(peak2peak(1, -1))
})

test_that("peak2peak() tests are correct", {
  
  x <- c(1:5)
  expect_equal(peak2peak(x), 4)
  
  x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
  expect_equal(peak2peak(x), c(2, 100, 1000))
  expect_equal(peak2peak(x, 1), c(999, 1498, 1997))
  
  x <- array(c(1, 1.5, 2, 100, 150, 200, 1000, 1500, 2000, 10000, 15000, 20000), c(2,3,2))
  expect_equal(peak2peak(x, 1), c(14999.0, 19998.5))
  expect_equal(peak2peak(x, 2), c(1499, 9998, 19850))
  expect_equal(peak2peak(x, 3), c(199, 19000))

  x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
  expect_equal(peak2peak(x), 4+8i)

})
