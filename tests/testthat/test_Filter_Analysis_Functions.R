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

# -----------------------------------------------------------------------
# freqz()

test_that("parameters to freqz() are correct", {
  expect_error(freqz())
  expect_error(freqz('invalid'))
})

test_that("freqz() tests are correct", {
  
  # test correct values and fft-polyval consistency
  # butterworth filter, order 2, cutoff pi/2 radians
  b <- c(0.292893218813452, 0.585786437626905, 0.292893218813452)
  a <- c(1, 0, 0.171572875253810)
  hw <- freqz(b, a, 32)
  expect_that(Re(hw$h[1]), equals(1))
  expect_that(abs(hw$h[17])^2, equals(0.5))
  expect_that(hw$h, equals(freqz(b, a, hw$w)$h))  # fft should be consistent with polyval
  
  # test whole-half consistency
  b <- c(1, 1, 1)/3  # 3-sample average
  hw <- freqz(b, 1, 32, whole = TRUE)
  expect_that(hw$h[2:16], equals(Conj(hw$h[32:18])))
  hw2 <- freqz(b, 1, 16, whole = FALSE)
  expect_that(hw$h[1:16], equals(hw2$h))
  expect_that(hw$w[1:16], equals(hw2$w))
  
  # test sampling frequency properly interpreted
  b <- c(1, 1, 1) / 3; a <- c(1, 0.2)
  hw <- freqz(b, a, 16, fs = 320)
  expect_that(hw$w, equals((0:15) * 10))
  hw2 <- freqz(b, a, (0:15) * 10, fs = 320)
  expect_that(hw2$w, equals((0:15) * 10))
  expect_that(hw$h, equals(hw2$h))
  hw3 <- freqz(b, a, 32, whole = TRUE, fs = 320)
  expect_that(hw3$w, equals((0:31) * 10))
})


# -----------------------------------------------------------------------
# grpdelay()

test_that("parameters to grpdelay() are correct", {
  expect_error(grpdelay())
  expect_error(grpdelay('invalid'))
})

test_that("grpdelay() tests are correct", {
  
  gd1 <- grpdelay(c(0, 1))
  gd2 <- grpdelay(c(0, 1), 1)
  expect_equal(gd1$gd, gd2$gd)

  gd <- grpdelay(c(0, 1), 1, 4)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, pi/4 * 0:3)

  gd <- grpdelay(c(0, 1), 1, 4, whole = TRUE)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, pi/2 * 0:3)

  gd <- grpdelay(c(0, 1), 1, 4, fs = 0.5)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, 1/16 * 0:3)

  gd <- grpdelay(c(0, 1), 1, 4, TRUE, 1)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, 1/4 * 0:3)

  gd <- grpdelay(c(1, -0.9i), 1, 4, TRUE, 1)
  gd0 = 0.447513812154696; gdm1 =0.473684210526316;
  expect_equal(gd$gd, c(gd0, -9, gd0, gdm1))
  expect_equal(gd$w, 1/4 * 0:3)
  
  gd <- grpdelay(1, c(1, .9), n = 2 * pi * c(0, 0.125, 0.25, 0.375))
  expect_equal(gd$gd, c(-0.47368, -0.46918, -0.44751, -0.32316), tolerance = 1e-5)
  
  gd <- grpdelay(1, c(1, .9), c(0, 0.125, 0.25, 0.375), fs = 1)
  expect_equal(gd$gd, c(-0.47368, -0.46918, -0.44751, -0.32316), tolerance = 1e-5)
  
  gd <- grpdelay(c(1, 2), c(1, 0.5, .9), 4)
  expect_equal(gd$gd, c(-0.29167, -0.24218, 0.53077, 0.40658), tolerance = 1e-5)

  b1 <- c(1, 2); a1f <- c(0.25, 0.5, 1); a1 <- rev(a1f);
  gd1 <- grpdelay(b1, a1, 4)$gd
  gd <- grpdelay(conv(b1, a1f), 1, 4)$gd - 2
  expect_equal(gd, gd1)
  expect_equal(gd, c(0.095238, 0.239175, 0.953846, 1.759360), tolerance = 1e-5)
  
  a <- c(1, 0, 0.9)
  b <- c(0.9, 0, 1)
  dh <- grpdelay(b, a, 512, 'whole')$gd
  da <- grpdelay(1, a, 512, 'whole')$gd
  db <- grpdelay(b, 1, 512, 'whole')$gd
  expect_equal(dh, db + da, ttolerance = 1e-5)

})

# -----------------------------------------------------------------------
# impz()

test_that("parameters to impz() are correct", {
  expect_error(impz())
  expect_error(impz('invalid'))
})

test_that("impz() tests are correct", {
  
  xt <- impz(1, c(1, -1, 0.9), 100)
  expect_equal(length(xt$t), 100L)
  expect_equal(xt$t, 0:99)
  
  xt <- impz(1, c(1, -1, 0.9), 0:101)
  expect_equal(length(xt$t), 102L)
  expect_equal(xt$t, 0:101)
  
})

