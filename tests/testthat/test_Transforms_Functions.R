# gsignal Transforms Functions
library(gsignal)
library(testthat)

tol <- 1e-6

# -----------------------------------------------------------------------
# cplxreal()

test_that("parameters to cplxreal() are correct", {
  expect_error(cplxreal())
  expect_error(cplxreal(1, 2, 3, 4))
  expect_error(cplxreal(1, matrix(1L, 2, 3)))
  expect_error(cplxreal(1, -1))
  expect_error(cplxreal(1, dim = 3))
})

test_that("cplxreal() tests are correct", {
  ret <- cplxreal(1)
  expect_equal(length(ret$zc), 0)
  expect_equal(ret$zr, 1)
  
  ret <- cplxreal(c(1 + 2i, 1 - 2i))
  expect_equal(ret$zc, 1 + 2i)
  expect_equal(length(ret$zr), 0)
  
  ret <- cplxreal(polyroot(c(1, 0, 0, 1)))
  expect_equal(ret$zc, complex(real = 0.5, imag = sinpi(1 / 3)), tolerance = tol)
  expect_equal(ret$zr, -1)
})

# -----------------------------------------------------------------------
# digitrevorder()

test_that("parameters to digitrevorder() are correct", {
  expect_error(digitrevorder())
  expect_error(digitrevorder(1))
  expect_error(digitrevorder(1, 2, 3))
  expect_error(digitrevorder(1, 1))
  expect_error(digitrevorder(1, 37))
  expect_error(digitrevorder(0:3, 8))
})

test_that("digitrevorder() tests are correct", {
  expect_equal(digitrevorder(0, 2), 0)
  expect_equal(digitrevorder(0, 36), 0)
  expect_equal(digitrevorder(0:3, 4), 0:3)
  expect_equal(digitrevorder(0:7, 2), c(0, 4, 2, 6, 1, 5, 3, 7))
  expect_equal(digitrevorder(0:7 * 1i, 2), c(0, 4, 2, 6, 1, 5, 3, 7) * 1i)
  expect_equal(digitrevorder(0:15, 2), c(0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15))
  expect_equal(digitrevorder(0:15, 4), c(0, 4, 8, 12, 1, 5, 9, 13, 2, 6, 10, 14, 3, 7, 11, 15))
})

# -----------------------------------------------------------------------
# bitrevorder()

test_that("parameters to bitrevorder() are correct", {
  expect_error(bitrevorder())
  expect_error(bitrevorder(1, 2))
  expect_error(bitrevorder(1, 2, 3))
  expect_error(bitrevorder(NULL))
  expect_error(bitrevorder(0:2))
})

test_that("bitrevorder() tests are correct", {
  expect_equal(bitrevorder(0), 0)
  expect_equal(bitrevorder(0:1), 0:1)
  expect_equal(bitrevorder(0:7), c(0, 4, 2, 6, 1, 5, 3, 7))
  expect_equal(bitrevorder(0:7 * 1i), c(0, 4, 2, 6, 1, 5, 3, 7) * 1i)
  expect_equal(bitrevorder(0:15), c(0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15))
})

# -----------------------------------------------------------------------
# fftshift()

test_that("parameters to fftshift() are correct", {
  expect_error(fftshift())
  expect_error(fftshift(1, 2, 3))
  #expect_error(fftshift(matrix(1:4, 2, 2), -1))
  expect_error(fftshift(NULL))
  expect_error(fftshift(array(1:8, c(2, 2, 2))))
})

test_that("fftshift() tests are correct", {
  
  expect_equal(fftshift(1), 1)
  
  x <- 0:7
  y <- fftshift(x)
  expect_equal(y, c(4, 5, 6, 7, 0, 1, 2, 3))
  expect_equal(fftshift(y), x)

  x <- 0:6
  y <- fftshift(x)
  expect_equal(y, c(4, 5, 6, 0, 1, 2, 3))
  expect_equal(fftshift(y), c(1, 2, 3, 4, 5, 6, 0))
  
  x <- 0:3
  x <- matrix(c(x, 2 * x, 3 * x + 1, 4 * x + 1), 4, byrow = TRUE)
  y <- fftshift(x, 1)
  expect_equal(y, matrix(c(1, 4, 7, 10, 1, 5, 9, 13, 0, 1, 2, 3, 0, 2, 4, 6), 4, byrow = TRUE))
  y <- fftshift(x, 2)
  expect_equal(y, matrix(c(2, 3, 0, 1, 4, 6, 0, 2, 7, 10, 1, 4, 9, 13, 1, 5), 4, byrow = TRUE))
  y <- fftshift(x, c(1, 2))
  expect_equal(y, matrix(c(7, 10, 1, 4, 9, 13, 1, 5, 2, 3, 0, 1, 4, 6, 0, 2), 4, byrow = TRUE))
  
})

# -----------------------------------------------------------------------
# ifftshift()

test_that("parameters to ifftshift() are correct", {
  expect_error(ifftshift())
  expect_error(ifftshift(1, 2, 3))
  #expect_error(ifftshift(matrix(1:4, 2, 2), -1))
  expect_error(ifftshift(NULL))
  expect_error(ifftshift(array(1:8, c(2, 2, 2))))
})

test_that("ifftshift() tests are correct", {
  
  expect_equal(ifftshift(1), 1)
  
  x <- 0:7
  y <- ifftshift(x)
  expect_equal(y, c(4, 5, 6, 7, 0, 1, 2, 3))
  expect_equal(ifftshift(y), x)
  
  x <- 0:6
  y <- ifftshift(x)
  expect_equal(y, c(3, 4, 5, 6, 0, 1, 2))
  expect_equal(ifftshift(y), c(6, 0, 1, 2, 3, 4, 5))
  
  x <- 0:3
  x <- matrix(c(x, 2 * x, 3 * x + 1, 4 * x + 1), 4, byrow = TRUE)
  y <- ifftshift(x, 1)
  expect_equal(y, matrix(c(1, 4, 7, 10, 1, 5, 9, 13, 0, 1, 2, 3, 0, 2, 4, 6), 4, byrow = TRUE))
  expect_equal(ifftshift(y, 1), x)
  y <- ifftshift(x, 2)
  expect_equal(y, matrix(c(2, 3, 0, 1, 4, 6, 0, 2, 7, 10, 1, 4, 9, 13, 1, 5), 4, byrow = TRUE))
  expect_equal(ifftshift(y, 2), x)
  y <- ifftshift(x, c(1, 2))
  expect_equal(y, matrix(c(7, 10, 1, 4, 9, 13, 1, 5, 2, 3, 0, 1, 4, 6, 0, 2), 4, byrow = TRUE))
  expect_equal(ifftshift(y, c(1, 2)), x)
  
})

# -----------------------------------------------------------------------
# cceps()

test_that("parameters to cceps() are correct", {
  expect_error(cceps())
  expect_error(cceps(1, 2))
  expect_error(cceps(matrix(1:4, 2, 2)))
  expect_error(cceps(TRUE))
  expect_error(cceps(1:10 * 1i))
})

test_that("cceps() tests are correct", {

  expect_error(cceps(rep(1L, 4)))
  expect_error(cceps(0))
  
  x <- runif (256)
  cps <- cceps(x)
  expect_equal(length(x), length(cps))
  
})

# -----------------------------------------------------------------------
# rceps()

test_that("parameters to rceps() are correct", {
  expect_error(rceps())
  expect_error(rceps(1, 2))
  expect_error(rceps(1, TRUE, 3))
  expect_error(rceps(matrix(1:4, 2, 2)))
  expect_error(rceps(TRUE))
  expect_error(rceps(1:10 * 1i))
})

test_that("rceps() tests are correct", {
  
  # Test that an odd-length input produces an odd-length output
  x <- runif(33)
  rc <- rceps(x, TRUE)
  expect_equal(length(rc$y), length(x))
  expect_equal(length(rc$ym), length(x))
})

# -----------------------------------------------------------------------
# czt()

test_that("parameters to czt() are correct", {
  expect_error(czt())
  expect_error(czt('a'))
  expect_error(czt(1, -1))
  expect_error(czt(1, 1, 'a'))
  expect_error(czt(1, 1, 1, 'a'))
  expect_error(czt(1, 2, 3, 4, 5))
})

test_that("czt() tests are correct", {
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5, 1)
  expect_equal(stats::fft(x), czt(x), tolerance = 1e-6)
  expect_equal(stats::mvfft(cbind(x, x)), czt(cbind(x, x)), tolerance = 1e-6)
})

# -----------------------------------------------------------------------
# dct() and idct

test_that("parameters to dct() and idct() are correct", {
  expect_error(dct())
  expect_error(dct('a'))
  expect_error(dct(1, -1))
  expect_error(dct(1, 1, 1))
  expect_error(idct())
  expect_error(idct('a'))
  expect_error(idct(1, -1))
  expect_error(idct(1, 1, 1))
})

test_that("dct() and idct() tests are correct", {
  # even-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5, 1)
  expect_equal(x, idct(dct(x)))
  expect_equal(cbind(x, x), idct(dct(cbind(x, x))))
  #uneven-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5)
  expect_equal(x, idct(dct(x)))
  expect_equal(cbind(x, x), idct(dct(cbind(x, x))))
  
})

# -----------------------------------------------------------------------
# dct2() and idct2

test_that("parameters to dct2() and idct2() are correct", {
  expect_error(dct2())
  expect_error(dct2('a'))
  expect_error(dct2(matrix(1:4, 2, 2), -1))
  expect_error(dct2(matrix(1:4, 2, 2), 1, -1))
  expect_error(dct2(1, 1, 1, 1))
  expect_error(idct2())
  expect_error(idct2('a'))
  expect_error(idct2(matrix(1:4, 2, 2), -1))
  expect_error(idct2(matrix(1:4, 2, 2), 1, -1))
  expect_error(idct2(1, 1, 1, 1))
})

test_that("dct2() and idct2() tests are correct", {
  # even-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5, 1)
  expect_equal(cbind(x, x), idct2(dct2(cbind(x, x))))
  #uneven-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5)
  expect_equal(cbind(x, x), idct2(dct2(cbind(x, x))))
  
})

# -----------------------------------------------------------------------
# dst() and idst

test_that("parameters to dst() and idst() are correct", {
  expect_error(dst())
  expect_error(dst('a'))
  expect_error(dst(1, -1))
  expect_error(dst(1, 1, 1))
  expect_error(idst())
  expect_error(idst('a'))
  expect_error(idst(1, -1))
  expect_error(idst(1, 1, 1))
})

test_that("dst() and idst() tests are correct", {
  # even-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5, 1)
  expect_equal(x, idst(dst(x)))
  expect_equal(cbind(x, x), idst(dst(cbind(x, x))))
  #uneven-length series
  x <- c(1, 2, 4, 1, 2, 3, 5, 2, 3, 5, 6, 7, 8, 4, 3, 6, 3, 2, 5)
  expect_equal(x, idst(dst(x)))
  expect_equal(cbind(x, x), idst(dst(cbind(x, x))))
  
})

# -----------------------------------------------------------------------
# fwht() and ifwht

test_that("parameters to fwht() and ifwht() are correct", {
  expect_error(fwht())
  expect_error(fwht('a'))
  expect_error(fwht(1, -1))
  expect_error(fwht(1, 1, 1))
  expect_error(fwht(1, 1, 'invalid'))
  expect_error(ifwht())
  expect_error(ifwht('a'))
  expect_error(ifwht(1, -1))
  expect_error(ifwht(1, 1, 1))
  expect_error(ifwht(1, 1, 'invalid'))
})

test_that("fwht() and ifwht() tests are correct", {
  expect_equal(fwht(rep(0L, 16)), rep(0L, 16))
  expect_equal(fwht(rep(1L, 16)), c(1L, rep(0L, 15)))
  expect_equal(fwht(rep(0L, 17)), rep(0L, 32))
  expect_equal(fwht(c(1, -1, 1, -1, 1, -1, 1, -1)), c(0, 0, 0, 0, 0, 0, 0, 1))
  
  expect_equal(ifwht(rep(0L, 16)), rep(0L, 16))
  expect_equal(ifwht(c(1L, rep(0L, 15))), rep(1L, 16))
  expect_equal(ifwht(rep(0L, 17)), rep(0L, 32))
  expect_equal(ifwht(c(0, 0, 0, 0, 0, 0, 0, 1)), c(1, -1, 1, -1, 1, -1, 1, -1))
  
  x <- matrix(round(runif(256) * 100), ncol = 16)
  expect_equal(ifwht(fwht(x)), x)
  expect_equal(ifwht(fwht(x, ordering = "sequency"), ordering = "sequency"), x)
  expect_equal(ifwht(fwht(x, ordering = "hadamard"), ordering = "hadamard"), x)
  expect_equal(ifwht(fwht(x, ordering = "dyadic"), ordering = "dyadic"), x)
})

# -----------------------------------------------------------------------
# hilbert()

test_that("parameters to hilbert() are correct", {
  expect_error(hilbert())
  expect_error(hilbert('a'))
  expect_warning(hilbert(1 + 1i))
  expect_error(hilbert(1, 'a'))
  expect_error(hilbert(1, -1))
})

test_that("hilbert() tests are correct", {
  x <- 1:4
  i <- c(1, -1, -1, 1)
  expect_equal(Re(hilbert(x)), x)
  expect_equal(Im(hilbert(x)), i)
  
  # test with small numbers (Github bug #4)
  expect_equal(Re(hilbert(1e-12 * x)), 1e-12 * x)
  expect_equal(Im(hilbert(1e-12 * x)), 1e-12 * i)
  
})
