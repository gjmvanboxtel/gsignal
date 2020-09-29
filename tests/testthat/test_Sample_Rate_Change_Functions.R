# gsignal Sample Rate Change Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# upfirdn()

test_that("parameters to upfirdn() are correct", {
  expect_error(upfirdn())
  expect_error(upfirdn(1))
  expect_error(upfirdn(1, 2, 3, 4, 5))
  expect_error(upfirdn(0:10, matrix(1:4, 2, 2), 1, 1))
})

test_that("upfirdn() tests are correct", {
  expect_equal(upfirdn(1:100, 1, 1, 1), seq(1, 100, 1))
  expect_equal(upfirdn(1:100, 1, 1, 2), seq(1, 100, 2))
})

# -----------------------------------------------------------------------
# resample()

test_that("parameters to resample() are correct", {
  expect_error(resample())
  expect_error(resample(1))
  expect_error(resample(1, 2))
  expect_error(resample(1, 2, 3, 4, 5))
  expect_error(resample(1, 1, 0.1))
  expect_error(resample(1, 0.1, 1))
})

test_that("resample() tests are correct", {
  expect_equal(resample(1:100, 1, 1), seq(1, 100, 1))
  expect_equal(length(resample(1:100, 1, 2)), 50)
  expect_equal(length(resample(1:100, 2, 1)), 200)
})
