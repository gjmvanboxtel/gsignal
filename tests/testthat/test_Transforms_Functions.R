# gsignal Transforms Functions
library(gsignal)
library(testthat)

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
  expect_equal(ret$zc, complex(real = 0.5, imag = sinpi(1 / 3)))
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

