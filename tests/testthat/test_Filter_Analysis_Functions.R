# gsignal filter analysis functions
library(gsignal)
library(testthat)

tol <- 1e-6
# -----------------------------------------------------------------------
# filternorm()

test_that("parameters to filternorm() are correct", {
  expect_error(filternorm())
  expect_error(filternorm(1, 1, 1))
  expect_error(filternorm(1, 1, 2, -1))
  expect_error(filternorm(matrix(1:10, 5, 2), 1))
  expect_error(filternorm(1, matrix(1:10, 5, 2)))
})

test_that("filternorm() tests are correct", {
  ba <- butter (5, 0.5)
  expect_equal(filternorm(ba), sqrt(2) / 2)
  expect_equal(filternorm (ba, Inf), 1)
  
  # identify filter
  expect_equal(filternorm(1, 1), 1)
  
  # known FIR filter
  ba <- Arma(c(1, 2, 1), 1)
  expect_equal(filternorm(ba), norm(ba$b, '2'))
  
  # gain scaling
  b <- c(1, -0.5, 0.25)
  a <- c(1, -0.2)
  n1 <- filternorm(b, a)
  n2 <- filternorm(3 * b, a)
  expect_equal(n2, 3 * n1)
  
  # zero numerator
  ba <- Arma(c(0, 0, 0), c(1, -0.5))
  expect_equal(filternorm(ba), 0)
  
  # denominator normalization
  b <- c(1, 2, 1)
  a <- c(1, -0.3)
  n1 <- filternorm(b, a)
  n2 <- filternorm(5 * b, 5 * a)
  expect_equal(n1, n2)
  
  # complex coefficients
  b <- c(1, 1i)
  a <- c(1, -0.2)
  L <- filternorm(b, a)
  expect_true(isPosscal(L))
  
  expect_equal(filternorm(c(1, 1i), 1), sqrt (2))

  ## Complex IIR impulse responses have real, nonnegative energy.
  b <- c(1, 1i)
  a <- c(1, -0.2)
  L <- filternorm(b, a)
  expect_true(is.numeric(L))
  expect_equal(filternorm(b, a, 2), L)
  
  # long filter
  b <- runif(100)
  a <- c(1, runif(20) * 0.01)
  L = filternorm(b, a)
  expect_true(is.finite(L))
  expect_true(isPosscal(L))

  # unstable filter
  zpg <- Zpg(1, c(1, -1.2), 1) # Pole outside the unit circle
  L <- filternorm(zpg)
  expect_true(is.na(L) || is.infinite(L) || isPosscal(L))

})

# -----------------------------------------------------------------------
# filtord()

test_that("parameters to filtord() are correct", {
  expect_error(filtord())
})

test_that("filtord() tests are correct", {
  b <- c(1, 0, 0)
  a <- c(1, 0, 0, 0)
  expect_equal(filtord(b, a), 3)

  ba <- butter(5, .5)
  expect_equal(filtord(ba), 5)
  
  ba <- butter(6, .5)
  expect_equal(filtord(ba), 6)

  sos <- tf2sos(c(1, 0, 0, 0, 0, 0, 0, 1), c(1, 0, 0, 0, 0, 0, 0, .5))
  expect_equal(filtord(sos), 7)  

  zpg <- tf2zp(c(1, 0, 0, 0, 0, 0, 1), c(1, 0, 0, 0, 0, 0, .5))
  expect_equal(filtord(zpg), 6)  
  
})

# -----------------------------------------------------------------------
# freqs()

test_that("parameters to freqs() are correct", {
  expect_error(freqs())
  expect_error(freqs(1))
  expect_error(freqs(1, 2))
})

test_that("freqs() tests are correct", {
  h <- freqs(1, 1, 1)
  expect_equal(h$h, 1 + 0i)

  h <- freqs(c(1, 2), c(1, 1), 1:4)
  expect_equal(h$h, c(1.5-0.5i, 1.2-0.4i, 1.1-0.3i, 1.058824-0.235294i), tolerance = tol)
  
})

# -----------------------------------------------------------------------
# freqspace()

test_that("parameters to freqspace() are correct", {
  expect_error(freqspace())
  expect_error(freqspace(-1))
  expect_error(freqspace(1, 2))
  expect_error(freqspace(c(1, 2, 3)))
  expect_error(freqspace(1, 'invalid'))
})

test_that("freqspace() tests are correct", {
  # 1D: freqspace(n) returns upper half-circle points
  expect_equal(freqspace(4), c(0, 0.5, 1))
  expect_equal(freqspace(5), c(0, 0.4, 0.8))
  expect_equal(freqspace(1), 0)

  # 1D: freqspace(n, 'whole') returns full unit circle
  expect_equal(freqspace(4, 'whole'), c(0, 0.5, 1, 1.5))
  expect_equal(freqspace(5, 'whole'), c(0, 0.4, 0.8, 1.2, 1.6))
  expect_equal(freqspace(2, 'whole'), c(0, 1))

  # 2D: f < freqspace(n, '2d')
  f <- freqspace(5, '2d')
  expect_equal(f$x, c(-0.8, -0.4, 0, 0.4, 0.80))
  expect_equal(f$y, c(-0.8, -0.4, 0, 0.4, 0.8))
        
  # 2D: f <- freqspace(c(m,n))
  f <- freqspace(c(4, 3))
  expect_equal(f$x, c(-2 / 3, 0, 2 / 3))
  expect_equal(f$y, c(-1, -0.5, 0, 0.5))
          
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
  expect_equal(fwhm(x, y), 2 * pi / 3, tolerance = tol)
  
  expect_equal(fwhm(y = -10:10), 0L)
  expect_equal(fwhm(y = rep(1L, 50)), 0L)

  x <- seq(-20, 20, 1)
  y1 <- -4 + rep(0L, length(x)); y1[4:10] <- 8
  y2 <- -2 + rep(0L, length(x)); y2[4:11] <- 2
  y3 <-  2 + rep(0L, length(x)); y3[5:13] <- 10
  expect_equal(fwhm(x, cbind(y1, y2, y3)), c(20 / 3, 7.5, 9.25),
               tolerance = tol)

  x <- 1:3
  y <- c(-1, 3, -1)
  expect_equal(fwhm(x, y), 0.75)
  expect_equal(fwhm(x, y, 'max'), 0.75)
  expect_equal(fwhm(x, y, 'zero'), 0.75)
  expect_equal(fwhm(x, y, 'middle'), 1L)
  expect_equal(fwhm(x, y, 'min'), 1L)

  x <- 1:3
  y <- c(-1, 3, -1)
  expect_equal(fwhm(x, y, level = 0.1), 1.35)
  expect_equal(fwhm(x, y, ref = 'max', level = 0.1), 1.35)
  expect_equal(fwhm(x, y, ref = 'min', level = 0.1), 1.40)
  expect_equal(fwhm(x, y, ref = 'abs', level = 2.5), 0.25)
  expect_equal(fwhm(x, y, ref = 'abs', level = -0.5), 1.75)
  
  x <- -5:5
  y <- 18 - x * x
  expect_equal(fwhm(y = y), 6)
  expect_equal(fwhm(x, y), 6)
  expect_equal(fwhm(x, y, 'min'), 7)
  
})

# -----------------------------------------------------------------------
# freqz()

test_that("parameters to freqz() are correct", {
  expect_error(freqz())
  expect_error(freqz('invalid'))
  expect_error(freqz(NA, 1))
  expect_error(freqz(1, NA))
  expect_error(freqz(1, 1, -1))
  expect_error(freqz(1, 1, 1, FALSE, 0))
})

test_that("freqz() tests are correct", {
  
  # test correct values and fft-polyval consistency
  # butterworth filter, order 2, cutoff pi/2 radians
  b <- c(0.292893218813452, 0.585786437626905, 0.292893218813452)
  a <- c(1, 0, 0.171572875253810)
  hw <- freqz(b, a, 32)
  expect_equal(Re(hw$h[1]), 1)
  expect_equal(abs(hw$h[17])^2, 0.5)
  expect_equal(hw$h, freqz(b, a, hw$w)$h,
               tolerance = tol)  # fft should be consistent with polyval
  
  # test whole-half consistency
  b <- c(1, 1, 1)/3  # 3-sample average
  hw <- freqz(b, 1, 32, whole = TRUE)
  expect_equal(hw$h[2:16], Conj(hw$h[32:18]), tolerance = tol)
  hw2 <- freqz(b, 1, 16, whole = FALSE)
  expect_equal(hw$h[1:16], hw2$h, tolerance = tol)
  expect_equal(hw$w[1:16], hw2$w, tolerance = tol)
  
  # test sampling frequency properly interpreted
  b <- c(1, 1, 1) / 3; a <- c(1, 0.2)
  hw <- freqz(b, a, 16, fs = 320)
  expect_equal(hw$w, (0:15) * 10)
  hw2 <- freqz(b, a, (0:15) * 10, fs = 320)
  expect_equal(hw2$w, (0:15) * 10)
  expect_equal(hw$h, hw2$h, tolerance = tol)
  hw3 <- freqz(b, a, 32, whole = TRUE, fs = 320)
  expect_equal(hw3$w, (0:31) * 10)
  
  # Github Issue #9
  sys <- Arma(b = c(0.000731569700125585, 0, -0.00292627880050234, 0, 
                    0.00438941820075351, 0, -0.00292627880050234,
                    0, 0.000731569700125585),
              a = c(1, -7.04203950456817, 21.7283807572297,
                    -38.3801847495478, 42.4586614063362, -30.1289923970756,
                    13.3939745471932, -3.41070179210893, 0.380901732541468)
              )
  hw <- freqz(sys, 5L)
  expect_equal(hw$h, c(3.255208e-06+0.000000e+00i, -9.759891e-04+1.062529e-01i,
                       3.328587e-03+2.658062e-03i,  3.101129e-04+1.143374e-04i,
                       1.305432e-05+2.075128e-06i), tol = 1e-3)
  
  hw <- freqz(sys, c(0, 1))
  expect_equal(hw$w, c(0, 1))
  expect_equal(hw$h, c(3e-06+0i, 8.254e-03+0.01047173i), tol = 1e-3)
  
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
  expect_equal(gd$w, pi/4 * 0:3, tolerance = tol)

  gd <- grpdelay(c(0, 1), 1, 4, whole = TRUE)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, pi/2 * 0:3, tolerance = tol)

  gd <- grpdelay(c(0, 1), 1, 4, fs = 0.5)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, 1/16 * 0:3, tolerance = tol)

  gd <- grpdelay(c(0, 1), 1, 4, TRUE, 1)
  expect_equal(gd$gd, rep(1L, 4))
  expect_equal(gd$w, 1/4 * 0:3)

  gd <- grpdelay(c(1, -0.9i), 1, 4, TRUE, 1)
  gd0 <- 0.447513812154696; gdm1 <- 0.473684210526316
  expect_equal(gd$gd, c(gd0, -9, gd0, gdm1), tolerance = tol)
  expect_equal(gd$w, 1/4 * 0:3)
  
  gd <- grpdelay(1, c(1, .9), n = 2 * pi * c(0, 0.125, 0.25, 0.375))
  expect_equal(gd$gd, c(-0.47368, -0.46918, -0.44751, -0.32316), tolerance = 1e-5)
  
  gd <- grpdelay(1, c(1, .9), c(0, 0.125, 0.25, 0.375), fs = 1)
  expect_equal(gd$gd, c(-0.47368, -0.46918, -0.44751, -0.32316), tolerance = 1e-5)
  
  gd <- grpdelay(c(1, 2), c(1, 0.5, .9), 4)
  expect_equal(gd$gd, c(-0.29167, -0.24218, 0.53077, 0.40658), tolerance = 1e-5)

  b1 <- c(1, 2); a1f <- c(0.25, 0.5, 1); a1 <- rev(a1f)
  gd1 <- grpdelay(b1, a1, 4)$gd
  gd <- grpdelay(conv(b1, a1f), 1, 4)$gd - 2
  expect_equal(gd, gd1, tolerance = 1e-5)
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
  
  ## FIR filter with scalar n
  xt <- impz(1, c(1, -0.5), 5)
  expect_equal(xt$x, c(1, 0.5, 0.25, 0.125, 0.0625))
  expect_equal(xt$t, 0:4)
  
  ## FIR filter with vector n
  xt <- impz(1, c(1, -0.5), c(0, 2, 4), 2)
  expect_equal(xt$x, c(1, 0.25, 0.0625))
  expect_equal(xt$t, c(0, 2, 4) / 2)
  
  ## With fs provided
  xt <- impz(1, c(1, -1, 0.9), 10, 1000)
  expect_equal(length(xt$x), 10)
  expect_equal(length(xt$t), 10)
  expect_equal(xt$t[2] - xt$t[1], 0.001)
  
})

# -----------------------------------------------------------------------
# impzlength()

test_that("parameters to impzlength() are correct", {
  expect_error(impzlength())
  expect_error(impzlength('invalid'))
  expect_error(impzlength(1, 'invalid'))
  expect_error(impzlength(1, 1, 'invalid'))
  expect_error(impzlength(1, 1, -1))
})

test_that("impzlength() tests are correct", {
  
  ## tests below shows different pole cases

  ## FIR filter
  expect_equal(impzlength(c(1, 2, 3, 4, 5)), 5)
  
  ## Stable IIR filter and with different tolerance
  b <- 1
  a <- c(1, -0.95)
  len1 <- impzlength(b, a)
  len2 <- impzlength(b, a, 1e-3)
  len3 <- impzlength(b, a, 1e-7)
  expect_equal(len1, 193)
  expect_equal(len2, 134)
  expect_equal(len3, 314)
    
  ## Unstable IIR filter
  b <- 1
  a <- c(1, -1.1)
  expect_equal(impzlength(b, a), 144)

  ## Oscillatory IIR filter
  b <- 1
  a <- c(1, 0.5, 1)
  expect_equal(impzlength(b, a), 17)

  ## IIR filter with repeated poles
  b <- 1
  a <- conv(c(1, -0.9), c(1, -0.9))
  expect_equal(impzlength(b, a), 187)
  
  ## Oscillatory IIR filter with damped component
  b <- 1
  a1 <- c(1, 0.5, 1)  ## oscillatory,len = 17
  a2 <- c(1, -0.1)    ## damped,len=4
  a3 <- c(1, -0.9)    ## damped,len=93
  len1 <- impzlength(b, conv(a1, a2))  ## damped component shorter than oscillatory
  len2 <- impzlength(b, a2)            ## damped
  len3 <- impzlength(b, conv(a1, a3))  ## damped component longer than oscillatory
  expect_equal(len1, 17)
  expect_equal(len2, 4)
  expect_equal(len3, 93)
  
  ## IIR filter with delay
  b <- c(0, 0, 1) 
  a <- c(1, -0.95)
  expect_equal(impzlength(b, a), 195)  ## 193 + 2 delay
              
  ## Special case: this is why r(i) = -r(i) to avoid case with arg(r) near 0
  b <- 1
  a <- c(1, -1)
  expect_equal(impzlength(b, a), 10)
  
  # zero numerator
  ba <- Arma(c(0, 0, 0), c(1, -0.5))
  expect_equal(impzlength(ba), 14)
  
  ## Second-order sections input
  sos <- ellip (4, 1, 60, 0.4, output = "Sos")
  expect_equal(impzlength(sos), 80)
})

# -----------------------------------------------------------------------
# isallpass()

test_that("parameters to isallpass() are correct", {
  expect_error(isallpass())
  expect_error(isallpass('invalid'))
})

test_that("isallpass() tests are correct", {
  
  b <- c((1 + 1i) / 2, -1)
  a <- c(1, -(1 - 1i) / 2)
  expect_true(isallpass(b, a))
  
  b <- c((1 + 1i) / 2, -1)
  a <- c(-1, (1 - 1i) / 2)
  expect_true(isallpass(b, a))
  
  ba <- butter(1, 0.5)
  expect_false(isallpass(ba))

  sos <- butter(1, 0.5, output = "Sos")
  expect_false(isallpass(sos))
  
  b1 <- 0.5 * (1 + 1i)
  b2 <- 0.7 * (cos(pi / 6) + 1i * sin(pi / 6))
  b <- conv(c(b1, -1), c(b2, -1))
  a <- conv(c(1, -Conj(b1)), c(1, -Conj(b2)))
  expect_true(isallpass(b, a))

})

# -----------------------------------------------------------------------
# ismaxphase()

test_that("parameters to ismaxphase() are correct", {
  expect_error(ismaxphase())
  expect_error(ismaxphase('invalid'))
  expect_error(ismaxphase(1, 1, tol = -1))
})

test_that("ismaxphase() tests are correct", {
  
  z1 <- c(0.9 * exp(1i * 0.6 * pi), 0.9 * exp(-1i * 0.6 * pi))
  z2 <- c(0.8 * exp(1i * 0.8 * pi), 0.8 * exp(-1i * 0.8 * pi))
  b <- gsignal::poly(c(z1, z2))
  a <- 1
  expect_false(ismaxphase(b, a))
  
  z1 <- c(0.9 * exp(1i * 0.6 * pi), 0.9 * exp(-1i * 0.6 * pi))
  z2 <- c(0.8 * exp(1i * 0.8 * pi), 0.8 * exp(-1i * 0.8 * pi))
  b <- gsignal::poly(c(1 / z1, 1 / z2))
  a <- 1
  expect_true(ismaxphase(b, a))
  
  z1 <- c(0.9 * exp(1i * 0.6 * pi), 0.9 * exp(-1i * 0.6 * pi))
  z2 <- c(0.8 * exp(1i * 0.8 * pi), 0.8 * exp(-1i * 0.8 * pi))
  b <- gsignal::poly(c(z1, 1 / z2))
  a <- 1
  expect_false(ismaxphase(b, a))
  
  z1 <- c(0.9 * exp(1i * 0.6 * pi), 0.9 * exp(-1i * 0.6 * pi))
  z2 <- c(0.8 * exp(1i * 0.8 * pi), 0.8 * exp(-1i * 0.8 * pi))
  b <- gsignal::poly(c(1 / z1, z2))
  a <- 1
  expect_false(ismaxphase(b, a))

  ba <- butter(1, 0.5)
  expect_false(ismaxphase(ba))

  sos <- butter(1, 0.5, output = "Sos")
  expect_false(ismaxphase(sos))
  
  zpg <- butter (8, .5);
  expect_false(ismaxphase(zpg))
})

# -----------------------------------------------------------------------
# isminphase()

test_that("parameters to isminphase() are correct", {
  expect_error(isminphase())
  expect_error(isminphase('invalid'))
  expect_error(isminphase(1, 1, tol = -1))
})

test_that("isminphase() tests are correct", {
  
  b <- c(3, 1)
  a <- c(1, 0.5)
  expect_true(isminphase(b, a))
  
  ba <- butter(1, 0.5)
  expect_false(isminphase(ba))

  zp <- butter(8, 0.5, output = "Zpg")
  expect_false(isminphase(zp))
  
  b <- 1.25^2 * conv(conv(conv(c(1, -0.9 * exp(-1i * 0.6 * pi)),
                               c(1, -0.9 * exp(1i * 0.6 * pi))),
                          c(1, -0.8 * exp(-1i * 0.8 * pi))),
                     c(1, -0.8 * exp(1i * 0.8 * pi)))
  a <- 1
  expect_true(isminphase(b, a))
})

# -----------------------------------------------------------------------
# isstable()

test_that("parameters to isstable() are correct", {
  expect_error(isstable())
  expect_error(isstable('invalid'))
  expect_error(isstable(1, NULL))
  
})

test_that("isstable() tests are correct", {
  
  b <- c(1, 2, 3, 4, 5, 5, 1, 2)
  a <- 1
  expect_true(isstable(b, a))
  
  b <- c(1, 2, 3, 4, 5, 5, 1, 2)
  a <- c(4, 5, 6, 7, 9, 10, 4, 6)
  expect_false(isstable(b, a))
  
  a <- gsignal::polystab(a)
  expect_true(isstable(b, a))
  
  zpg <- butter(6, 0.7, 'high', output = 'Zpg')
  expect_true(isstable(zpg))
  
  sos <- as.Sos(zpg)
  expect_true(isstable(sos))
  
  ba <- Arma(c(1, -0.5), c(1, -1))
  expect_false(isstable(ba))

})

