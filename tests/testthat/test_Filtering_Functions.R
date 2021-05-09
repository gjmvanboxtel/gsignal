# gsignal filtering functions
library(gsignal)
library(testthat)

tol <- 1e-6

# -----------------------------------------------------------------------
# filtfilt()

test_that("parameters to filtfilt() are correct", {
  expect_error(filtfilt())
  expect_error(filtfilt(0, 0, 1:10))
  expect_error(filtfilt(1, 2, array(1:8, c(2, 2, 2))))
  expect_error(filtfilt(1, 1, c('invalid', 'invalid')))
})

test_that("filtfilt() tests are correct", {
  expect_equal(filtfilt(1, 1, 1:2), c(1,2))
  expect_equal(filtfilt(1, 2, 1:2), c(0.25, 0.50))
  expect_equal(filtfilt(2, 1, 1:2), c(4, 8))
  x <- runif(100)
  y <- filtfilt(1, 1, x)
  expect_equal(length(y), length(x))
  x <- matrix(runif(200), 100, 2)
  y <- filtfilt(1, 1, x)
  expect_equal(ncol(y), ncol(x))
  expect_equal(nrow(y), nrow(x))
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
  expect_equal(filtic(b, a, 1, 1), 0.75)
  
  # Simple high pass filter
  b <- c(0.25, -0.25)
  a <- c(1.0, 0.5)
  expect_equal(filtic(b, a, 0, 1), -0.25)

  # Second order cases
  # bs <- butter(2, 0.4)
  b <- c(0.2065721, 0.4131442, 0.2065721)
  a <- c(1.0000000, -0.3695274,  0.1958157)
  x <- y <- c(1, 1)
  expect_equal(filtic(b, a, y, x), c(0.7934280, 0.0107564), tolerance = tol)
  N <- 1000
  xx <- cos(2 * pi * seq(0, N-1, length.out = N)/8)
  yy <- filter(b, a, xx)
  x <- xx[seq(N, N - 1, -1)]
  y <- yy[seq(N, N - 1, -1)]
  zf <- filtic(b, a, y, x)
  expect_equal(filtic(b, a, y, x), c( 0.4039015, 0.1625113), tolerance = tol)
  
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
  expect_equal(medfilt1(1:10), 1:10)
  expect_equal(medfilt1(c(1, 1, 2, 3, 3, 4, 4, 4, 5)), 
               c(1, 1, 2, 3, 3, 4, 4, 4, 4))
  expect_equal(medfilt1(c(1, 1, 2, 3, NA, 4, 4, 4, 5)),
               c(1, 1, 2, 3, 3.676871, 4, 4, 4, 4), tolerance = tol)
  expect_equal(medfilt1(c(1, 1, 2, 3, NA, 4, 4, 4, 5), na.omit = TRUE),
               c(1, 1, 2, 3, 4, 4, 4, 4))
  expect_equal(medfilt1(cbind(1:5, 1:5)), cbind(1:5, 1:5))
  expect_equal(medfilt1(cbind(1:5, 1:5), n = 1, MARGIN = 1),
              rbind(1:5, 1:5))
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
  expect_equal(r$rmsx, Inf)
  expect_equal(r$w, 1)

  r <- movingrms(matrix(1:100, 50), 1)
  expect_equal(ncol(r$rmsx), 2)
  expect_equal(nrow(r$rmsx), 50)
  expect_equal(r$w, c(rep(0, 23), 0.5, 1, 0.5, rep(0, 24)))
  
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
  expect_error(sosfilt(c(0,0,0,0,0,0), 1))
  expect_error(sosfilt(rep(1, 6), 1, 'invalid'))
})

test_that("sosfilt() tests are correct", {
  expect_equal(sosfilt(c(0, 0, 0, 1, 0, 0), 1), 0)
  expect_equal(sosfilt(c(0, 0, 0, 1, 0, 0), c(1, 1)), c(0, 0))

  sos <- rbind(c(0,1,0,1,-1,0),c(1,2,1,1,-2,1))
  x <- 1:10
  y <- sosfilt(sos,x)
  expect_equal(y, c(0, 1, 7, 26, 70, 155, 301, 532, 876, 1365))

  # initial conditions  
  sos <- rbind(c(0,1,0,1,-1,0), c(1,2,1,1,-2,1))
  x1 <- 1:10
  y1 <- sosfilt(sos, x1, "zf")
  expect_equal(y1$y, c(0, 1, 7, 26, 70, 155, 301, 532, 876, 1365))
  expect_equal(y1$zf, matrix(c(55, 1980, 0, -1320), ncol = 2))
  x2 <- 11:20
  y2 <- sosfilt(sos, x2, y1$zf)
  expect_equal(y2$y, c(2035,2926,4082,5551,7385,9640,12376,15657,19551,24130))
  expect_equal(y2$zf, matrix(c(210, 29260, 0, -23940), ncol = 2))
  x <- 1:20
  y <- sosfilt(sos, x)
  expect_equal(y, c(y1$y, y2$y))

  # multidimensional
  sos <- rbind(c(0,1,0,1,-1,0), c(1,2,1,1,-2,1))
  x <- cbind(1:10, 11:20)
  y <- sosfilt(sos, x, "zf")
  expect_equal(y$y, cbind(c(0,1,7,26,70,155,301,532,876,1365),
                          c(0,11,67,216,510,1005,1761,2842,4316,6255)))
  expect_equal(y$zf, array(c(55,1980,0,-1320,155,8580,0,-6120), c(2,2,2)))
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
  expect_equal(fftfilt(b, x), res)
  expect_equal(fftfilt(b, replicate(2, x)), replicate(2,res))
  expect_equal(fftfilt(b, replicate(2, x + 2 *.Machine$double.eps)),
               replicate(2,res), tolerance = tol)
  
  r <- sqrt (1/2) * (1+1i)
  b <-  c(1, 1) * r
  x <- c(1L, rep(0L, 9))
  res <- c(rep(1L, 2), rep(0L, 8))
  expect_equal(fftfilt(b, x), r * res, tolerance = tol)
  expect_equal(fftfilt(b, r * x), r * r * res, tolerance = tol)

  b  <- c(1, 1)
  x  <- matrix(rep(0L, 30), 10, 3); x[1, 1] <--1; x[1, 2] <- 1
  y0 <- matrix(rep(0L, 30), 10, 3); y0[1:2, 1] <- -1; y0[1:2, 2] <- 1
  y  <- fftfilt(b, x)
  expect_equal(y0, y)
  y  <- fftfilt(b * 1i, x)
  expect_equal(y0 * 1i, y)
  y  <- fftfilt(b, x * 1i)
  expect_equal(y0 * 1i, y)
  y  <- fftfilt(b * 1i, x * 1i)
  expect_equal(-y0, y)
  x  <- runif(10)
  y  <- fftfilt(b, cbind(x, x * 1i))
  expect_equal(all(abs(Im(y[, 1])) < tol), TRUE)
  expect_equal(all(abs(Re(y[, 2])) < tol), TRUE)
  
  b  <- runif(10)
  x  <- runif(10)
  y0 <- filter(b, 1, x)
  y  <- fftfilt(b, x)
  expect_equal(y0, y, tolerance =  tol)
  
})

# -----------------------------------------------------------------------
# filter_zi()

test_that("parameters to filter_zi() are correct", {
  expect_error(filter_zi())
  expect_error(filter_zi(1))
  expect_error(filter_zi(1, 2))
  expect_error(filter_zi(1, 2, 3, 4, 5))
  expect_error(filter_zi(0, 0, 'invalid'))
})

test_that("filter_zi() tests are correct", {
  
  h <- butter(2, 0.4)
  l <- max(length(h$b), length(h$a)) - 1
  x <- y <- rep(1, l)
  expect_equal(filtic(h, y, x), filter_zi(h), tolerance = tol)

})
