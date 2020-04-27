# gsignal IIR filter design functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# cheb()

test_that("parameters to cheb() are correct", {
  expect_error(cheb())
  expect_error(cheb(0.5))
  expect_error(cheb(-1L))
  expect_error(cheb(array(1L, c(1, 4))))
})

test_that("cheb() tests are correct", {
  expect_that(cheb(1, 1), equals(1))
  expect_that(cheb(2, 1), equals(1))
  expect_that(cheb(5, 2), equals(362))
  expect_that(cheb(5, c(2,3)), equals(c(362, 3363)))
})
  
# -----------------------------------------------------------------------
# besselap()

test_that("parameters to besselap() are correct", {
  expect_error(besselap())
  expect_error(besselap(0.5))
  expect_error(besselap(-1L))
  expect_error(besselap(array(1L, c(1, 4))))
  expect_error(besselap(1, 2))
})

test_that("besselap() tests are correct", {
  expect_equal(besselap(1)$p, -1)
  expect_equal(besselap(2)$p, c(-0.8660254+0.5i, -0.8660254-0.5i))
  expect_equal(besselap(3)$p, c(-0.7456404+0.7113666i, -0.7456404-0.7113666i, -0.9416000+0.0000000i),
               tolerance = 1e-7)
})
