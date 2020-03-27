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
  expect_equal(ret$zc, c(1 - 2i, 1 + 2i))
  expect_equal(length(ret$zr), 0)
  
  ret <- cplxreal(polyroot(c(1, 0, 0, 1)))
  expect_equal(ret$zc, c(complex(real = 0.5, imag = -sin(pi / 3)),
                         complex(real = 0.5, imag = +sin(pi / 3))))
  expect_equal(ret$zr, -1)
})
