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

# -----------------------------------------------------------------------
# medfilt1()

test_that("parameters to medfilt1() are correct", {
  expect_error(medfilt1())
  expect_error(medfilt1(1, 2))
  expect_error(medfilt1(1, -1))
  expect_error(medfilt1(cbind(1:10, 1:10), 3, 3))
  expect_error(medfilt1('invalid'))
  expect_error(medfilt1(1:10, endrule = 'invalid'))
  expect_error(medfilt1(1:10, algorithm = 'invalid'))
  expect_error(medfilt1(1:10, printy.level = 'invalid'))
})

test_that("medfilt1() tests are correct", {
  expect_that(medfilt1(1:10), equals(1:10))
  expect_that(medfilt1(c(1, 1, 2, 3, 3, 4, 4, 4, 5)), equals(c(1, 1, 2, 3, 3, 4, 4, 4, 4)))
  expect_that(medfilt1(c(1, 1, 2, 3, NA, 4, 4, 4, 5)), equals(c(1, 1, 2, 3, 3.676871, 4, 4, 4, 4), tolerance = 1e-7))
  expect_that(medfilt1(c(1, 1, 2, 3, NA, 4, 4, 4, 5), na.omit = TRUE), equals(c(1, 1, 2, 3, 4, 4, 4, 4)))
  expect_that(medfilt1(cbind(1:5, 1:5)), equals(cbind(1:5, 1:5)))
  expect_that(medfilt1(cbind(1:5, 1:5), n = 1, dim = 1), equals(rbind(1:5, 1:5)))
})

# -----------------------------------------------------------------------
# movingrms()

test_that("parameters to movingrms() are correct", {
  expect_error(movingrms())
  expect_error(movingrms(1, -1))
  expect_error(movingrms(1, 1, -1))
  expect_error(movingrms(1, 1, 1, -1))
  expect_error(movingrms('invalid'))
  expect_error(movingrms(1, 2, 3, 4, 5))
})

test_that("movingrms() tests are correct", {
  r <- movingrms(1, 1)
  expect_that(r$rmsx, equals(Inf))
  expect_that(r$w, equals(1))

  r <- movingrms(matrix(1:100, 50), 1)
  expect_that(ncol(r$rmsx), equals(2))
  expect_that(nrow(r$rmsx), equals(50))
  expect_that(r$w, equals(c(rep(0, 23), 0.5, 1, 0.5, rep(0, 24))))
  
})

# tests for sgolayfilt are in test_FIR_Filter_design_functions.R
# together with the sgolay() function

# -----------------------------------------------------------------------
# sosfilt()

test_that("parameters to sosfilt() are correct", {
  expect_error(sosfilt())
  expect_error(sosfilt(1, -1))
  expect_error(sosfilt(rep(1, 6), 'invalid'))
  expect_error(sosfilt(1, 1, 1))
})

test_that("sosfilt() tests are correct", {
  expect_that(sosfilt(c(0,0,0,0,0,0), 1), equals(NA))
  expect_that(sosfilt(c(0,0,0,0,0,0), c(1, 1)), equals(c(NA, NA)))
  expect_that(sosfilt(c(0, 0, 0, 1, 0, 0), 1), equals(0))
  expect_that(sosfilt(c(0, 0, 0, 1, 0, 0), c(1, 1)), equals(c(0, 0)))

  sos <- rbind(c(0,1,0,1,-1,0),c(1,2,1,1,-2,1))
  x=1:10
  y=sosfilt(sos,x)
  expect_that(y, equals(c(0, 1, 7, 26, 70, 155, 301, 532, 876, 1365)))
})

# -----------------------------------------------------------------------
# fftfilt()

test_that("parameters to fftfilt() are correct", {
  expect_error(fftfilt())
  expect_error(fftfilt(1))
  expect_error(fftfilt(1, 2, 3, 4))
  expect_error(fftfilt(matrix(rep(1L, 4), 2), 1))
  expect_error(fftfilt(2, array(rep(1L, 12), dim = c(2, 3, 2))))
  expect_error(fftfilt(2, 1, matrix(rep(1L, 4), 2)))
})

test_that("fftfilt() tests are correct", {

  b <- c(1, 1)
  x <- c(1L, rep(0L, 9))
  res <- c(rep(1L, 2), rep(0L, 8))
  expect_that(fftfilt(b, x), equals(res))
  expect_that(fftfilt(b, replicate(2, x)), equals(replicate(2,res)))
  expect_that(fftfilt(b, replicate(2, x + 2 *.Machine$double.eps)), equals(replicate(2,res)))
  
  r <- sqrt (1/2) * (1+1i)
  b <-  c(1, 1) * r
  x <- c(1L, rep(0L, 9))
  res <- c(rep(1L, 2), rep(0L, 8))
  expect_that(fftfilt(b, x), equals(r * res))
  expect_that(fftfilt(b, r * x), equals(r * r * res))

  b  <- c(1, 1)
  x  <- matrix(rep(0L, 30), 10, 3); x[1, 1] <--1; x[1, 2] <- 1
  y0 <- matrix(rep(0L, 30), 10, 3); y0[1:2, 1] = -1; y0[1:2, 2] <- 1
  y  <- fftfilt(b, x)
  expect_that(y0, equals(y))
  y  <- fftfilt(b * 1i, x)
  expect_that(y0 * 1i, equals(y))
  y  <- fftfilt(b, x * 1i)
  expect_that(y0 * 1i, equals(y))
  y  <- fftfilt(b * 1i, x * 1i)
  expect_that(-y0, equals(y))
  x  <- runif(10)
  y  <- fftfilt(b, cbind(x, x * 1i))
  expect_that(all(abs(Im(y[, 1])) < .Machine$double.eps), equals(TRUE))
  expect_that(all(abs(Re(y[, 2])) < .Machine$double.eps), equals(TRUE))
  
  b  <- runif(10)
  x  <- runif(10)
  y0 <- filter(b, 1, x)
  y  <- fftfilt(b, x)
  expect_that(y0, equals(y))
  
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
