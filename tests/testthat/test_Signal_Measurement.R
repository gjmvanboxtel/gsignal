# gsignal Signal Measurement
library(gsignal)
library(testthat)

tol <- 1e-6

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
  x <- pmin(3, cos (2*pi*c(0:8000) / 600) + 2.01)
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

# -----------------------------------------------------------------------
# peak2rms()

test_that("parameters to peak2rms() are correct", {
  expect_error(peak2rms())
  expect_error(peak2rms('invalid'))
  expect_error(peak2rms(1, 2, 3))
  expect_error(peak2rms(1, 1.5))
  expect_error(peak2rms(1, -1))
})

test_that("peak2rms() tests are correct", {
  
  expect_equal(peak2rms(1), 1L)
  expect_equal(peak2rms(-5), 1L)
  
  x <- c(1:5)
  expect_equal(peak2rms(x), 5 / sqrt(11), tolerance = tol)
  
  x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
  expect_equal(peak2rms(x), c(3/sqrt(14/3), 200/sqrt(72500/3), 2000/sqrt(7250000/3)),
               tolerance = tol)
  expect_equal(peak2rms(x, 1), c(1000/sqrt(1010001/3), 1500/sqrt(2272504/3), 2000/sqrt(4040009/3)),
               tolerance = tol)
  
  x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
  expect_equal(peak2rms(x), 1.552125, tolerance = tol)
  
})

# -----------------------------------------------------------------------
# rms()

test_that("parameters to rms() are correct", {
  expect_error(rms())
  expect_error(rms('invalid'))
  expect_error(rms(1, 2, 3))
  expect_error(rms(1, 1.5))
  expect_error(rms(1, -1))
})

test_that("rms() tests are correct", {
  
  expect_equal(rms(0), 0L)
  expect_equal(rms(1), 1L)
  expect_equal(rms(c(1, 2, -1)), sqrt(2), tolerance = tol)
  
  x <- c(1:5)
  expect_equal(rms(x), sqrt(11), tolerance = tol)
  
  x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
  expect_equal(rms(x), c(sqrt(14/3), sqrt(72500/3), sqrt(7250000/3)),
               tolerance = tol)
  expect_equal(rms(x, 1), c(sqrt(336667), sqrt(2272504/3), sqrt(1346670)),
               tolerance = tol)
  
  x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
  expect_equal(rms(x), 6.63325, tolerance = tol)
  
})


# -----------------------------------------------------------------------
# rssq()

test_that("parameters to rssq() are correct", {
  expect_error(rssq())
  expect_error(rssq('invalid'))
  expect_error(rssq(1, 2, 3))
  expect_error(rssq(1, 1.5))
  expect_error(rssq(1, -1))
})

test_that("rssq() tests are correct", {
  
  expect_equal(rssq(1), 1L)
  expect_equal(rssq(-5), 5L)
  expect_equal(rssq(c(1, 2, -1)), sqrt(6), tolerance = tol)
  
  x <- c(1:5)
  expect_equal(rssq(x), sqrt(55), tolerance = tol)
  
  x <- matrix(c(1,2,3, 100, 150, 200, 1000, 1500, 2000), 3, 3)
  expect_equal(rssq(x), c(sqrt(14), sqrt(72500), sqrt(7250000)),
               tolerance = tol)
  expect_equal(rssq(x, 1), c(sqrt(1010001), sqrt(2272504), sqrt(4040009)),
               tolerance = tol)
  
  x <- c(1+1i, 2+3i, 3+5i, 4+7i, 5+9i)
  expect_equal(rssq(x), 14.8324, tolerance = tol)
  
})

