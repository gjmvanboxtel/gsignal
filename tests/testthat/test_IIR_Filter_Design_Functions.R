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
  