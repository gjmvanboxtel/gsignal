# gsignal Filter Conversion functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# sos2tf()

test_that("parameters to sos2tf() are correct", {
  expect_error(sos2tf())
  expect_error(sos2tf(1, 2, 3))
  expect_error(sos2tf(1, 1))
})

test_that("sos2tf() tests are correct", {
  sos <- rbind(c(1, 1, 1, 1, 0, -1), c(-2, 3, 1, 1, 10, 1))
  ba <- sos2tf(sos)
  expect_equal(ba$b, c(-2, 1, 2, 4, 1))
  expect_equal(ba$a, c(1, 10, 0, -10, -1))
  
  ba <- sos2tf(sos, 2)
  expect_equal(ba$b, c(-4, 2, 4, 8, 2))
  expect_equal(ba$a, c(1, 10, 0, -10, -1))
  
  ba <- sos2tf(sos, c(2, 2, 2))
  expect_equal(ba$b, c(-16, 8, 16, 32, 8))
  expect_equal(ba$a, c(1, 10, 0, -10, -1))
  
})

# -----------------------------------------------------------------------
# sos2zp()

test_that("parameters to sos2zp() are correct", {
  expect_error(sos2zp())
  expect_error(sos2zp(1, 2, 3))
  expect_error(sos2zp(1, 1))
})

test_that("sos2zp() tests are correct", {
  sos <- rbind(c(1, 2, 3, 1, 0.2, 0.3), c(4, 5, 6, 1, 0.4, 0.5))
  zref <- c(-1-1.41421356237310i, -1+1.41421356237310i, -0.625-1.05326872164704i, -0.625+1.05326872164704i)
  pref <- c(-0.2-0.678232998312527i, -0.2+0.678232998312527i, -0.1-0.538516480713450i, -0.1+0.538516480713450i)
  kref <- 4
  zpg <- sos2zp(sos, 1)
  expect_equal(cplxpair(zpg$z, 1e-7), as.vector(zref))
  expect_equal(cplxpair(zpg$p, 1e-7), as.vector(pref))
  expect_equal(zpg$k, 4)
})

# -----------------------------------------------------------------------
# tf2zp()

test_that("parameters to tf2zp() are correct", {
  expect_error(tf2zp())
  expect_error(tf2zp(1, 2, 3))
  expect_error(tf2zp('invalid', 'invalid'))
})

test_that("tf2zp() tests are correct", {
  b <- c(2, 3)
  a <- c(1, 1/sqrt(2), 1/4)
  zpk <- tf2zp(b, a)
  expect_equal(zpk$z, pracma::roots(b))
  expect_equal(zpk$p, pracma::roots(a))
  expect_equal(zpk$k, 2)
  
})

# -----------------------------------------------------------------------
# zp2sos()

test_that("parameters to zp2sos() are correct", {
  expect_error(zp2sos())
  expect_error(zp2sos(1, 2, 3, 4))
  expect_error(zp2sos('invalid', 'invalid'))
})

test_that("zp2sos() tests are correct", {
  sosg <- zp2sos(c(0+1i, 0-1i), c(0+1i, 0-1i))
  expect_equal(sosg$sos, c(1, 0, 1, 1, 0, 1))
  expect_equal(sosg$g, 1)

  sosg <- zp2sos(c(1+1i, 1-1i), c(1+1i, 1-1i))
  expect_equal(sosg$sos, c(1, -2, 2, 1, -2, 2))
  expect_equal(sosg$g, 1)

  sosg <- zp2sos(c(1+1i, 1-1i), c(1+1i, 1-1i), 3)
  expect_equal(sosg$sos, c(1, -2, 2, 1, -2, 2))
  expect_equal(sosg$g, 3)
  
  # these are slightly different in Matlab (b[0] and b[1] swapped),
  # and produce errors in Octave
  expect_equal(as.vector(zp2sos(NULL, 0, 0)$sos), c(1, 0, 0, 1, 0, 0))
  expect_equal(as.vector(zp2sos(NULL, 1, 0)$sos), c(1, 0, 0, 1, -1, 0))
  expect_equal(as.vector(zp2sos(NULL, -1, 1)$sos), c(1, 0, 0, 1, 1, 0))
  
})

# -----------------------------------------------------------------------
# tf2sos()

test_that("parameters to tf2sos() are correct", {
  expect_error(tf2sos())
  expect_error(tf2sos(1, 2, 3))
  expect_error(tf2sos('invalid', 'invalid'))
})

test_that("tf2sos() tests are correct", {
  
  b <- c(1, 0, 0, 0, 0, 1)
  a <- c(1, 0, 0, 0, 0, .9)
  sosg <- tf2sos (b, a)
  sec1 <- c(1, 1, 0, 1,  0.9791484, 0)
  sec2 <- c(1, -1.618034, 1, 1, -1.5842953, 0.9587315)
  sec3 <- c(1, 0.618034, 1, 1, 0.6051470, 0.9587315)
  expect_equal(sosg$sos, rbind(sec1, sec2, sec3, deparse.level = 0))
  
  # these are slightly different in Matlab (b[0] and b[1] swapped),
  # and produce errors in Octave
  sosg <- tf2sos(c(0, 0), c(1,1))
  expect_equal(sosg$sos, c(1, 0, 0, 1, 1, 0))
  expect_equal(sosg$g, 1)

})

# -----------------------------------------------------------------------
# zp2tf()

test_that("parameters to zp2tf() are correct", {
  expect_error(zp2tf())
  expect_error(zp2tf(1, 2, 3, 4))
  expect_error(zp2tf('invalid', 'invalid'))
})

test_that("zp2tf() tests are correct", {
  # Matlab returns a = [1 0.01 1] - bug?
  # Octave gives an error
  ba <- zp2tf(c(0, 0), pracma::roots(c(1, 0.01, 1)), 1)
  expect_equal(ba$b, c(1, 0, 0))
  expect_equal(ba$a, c(1, 0.01, 1))
  
  # design 2-pole notch filter at pi/4 radians = 0.5/4 = 0.125 * fs
  w = pi/4
  # zeroes at r = 1
  r <- 1
  z1 <- r * exp(1i * w)
  z2 <- r * exp(1i * -w)
  # poles at r = 0.9
  r = 0.9
  p1 <- r * exp(1i * w)
  p2 <- r * exp(1i * -w)
  
  zeros <- c(z1, z2)
  poles <- c(p1, p2)
  ba <- zp2tf(zeros, poles, 1)
  inv <- tf2zp(ba$b, ba$a)
  expect_equal(inv$z, zeros)
  expect_equal(inv$p, poles)
  expect_equal(inv$k, 1)
})

