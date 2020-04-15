# gsignal filter analysis functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# freqs()

test_that("parameters to freqs() are correct", {
  expect_error(freqs())
  expect_error(freqs(1))
  expect_error(freqs(1, 2))
  expect_error(freqs(1, 2, 3, plot = 'invalid'))
})

test_that("freqs() tests are correct", {
  h <- freqs(1, 1, 1, plot = FALSE)
  expect_that(h, equals(1 + 0i))

  h <- freqs(c(1, 2), c(1, 1), 1:4, plot = FALSE)
  expect_that(h, equals(c(1.5-0.5i, 1.2-0.4i, 1.1-0.3i, 1.058824-0.235294i), tolerance = 1e-6))
  
})

# -----------------------------------------------------------------------
# fwhm()

test_that("parameters to fwhm() are correct", {
  expect_error(fwhm())
  expect_error(fwhm(1))
  expect_error(fwhm(1, 2, 3, 4, 5))
  expect_error(fwhm(array(1)))
  expect_error(fwhm(1, c(1,2)))
  expect_error(fwhm(1, array(1:12, dim = c(2, 3, 2))))
  expect_error(fwhm(1, 2, ref = 'invalid'))
  expect_error(fwhm(1, 2, level = 'invalid'))
})

test_that("fwhm() tests are correct", {

  x <- seq(-pi, pi, 0.001)
  y <- cos(x)
  expect_that(fwhm(x, y), equals(2 * pi / 3))
  
  expect_that(fwhm(y = -10:10), equals(0L))
  expect_that(fwhm(y = rep(1L, 50)), equals(0L))

  x <- seq(-20, 20, 1)
  y1 <- -4 + rep(0L, length(x)); y1[4:10] <- 8
  y2 <- -2 + rep(0L, length(x)); y2[4:11] <- 2
  y3 =  2 + rep(0L, length(x)); y3[5:13] <- 10
  expect_that(fwhm(x, cbind(y1, y2, y3)), equals(c(20/3, 7.5, 9.25)))

  x <- 1:3
  y <- c(-1, 3, -1)
  expect_that(fwhm(x, y), equals(0.75))
  expect_that(fwhm(x, y, 'max'), equals(0.75))
  expect_that(fwhm(x, y, 'zero'), equals(0.75))
  expect_that(fwhm(x, y, 'middle'), equals(1L))
  expect_that(fwhm(x, y, 'min'), equals(1L))

  x <- 1:3
  y <- c(-1, 3, -1)
  expect_that(fwhm(x, y, level = 0.1), equals(1.35))
  expect_that(fwhm(x, y, ref = 'max', level = 0.1), equals(1.35))
  expect_that(fwhm(x, y, ref = 'min', level = 0.1), equals(1.40))
  expect_that(fwhm(x, y, ref = 'abs', level = 2.5), equals(0.25))
  expect_that(fwhm(x, y, ref = 'abs', level = -0.5), equals(1.75))
  
  x <- -5:5
  y <- 18 - x * x
  expect_that(fwhm(y = y), equals(6))
  expect_that(fwhm(x, y), equals(6))
  expect_that(fwhm(x, y, 'min'), equals(7))
  
})
