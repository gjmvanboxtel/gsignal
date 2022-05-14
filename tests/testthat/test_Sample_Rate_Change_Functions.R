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
  expect_equal(upfirdn(cbind(1:100, 1:100), 1, 1, 1),
               cbind(seq(1, 100, 1), seq(1, 100, 1)))
  expect_equal(upfirdn(cbind(1:100, 1:100), 1, 1, 2),
               cbind(seq(1, 100, 2), seq(1, 100, 2)))
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
  expect_equal(resample(cbind(1:100, 1:100), 1, 1),
               cbind(seq(1, 100, 1), seq(1, 100, 1)))
  expect_equal(length(resample(1:100, 1, 2)), 50)
  expect_equal(nrow(resample(cbind(1:100, 1:100), 1, 2)), 50)
  expect_equal(length(resample(1:100, 2, 1)), 200)
  expect_equal(nrow(resample(cbind(1:100, 1:100), 2, 1)), 200)
})

# -----------------------------------------------------------------------
# downsample()

test_that("parameters to downsample() are correct", {
  expect_error(downsample())
  expect_error(downsample(1))
  expect_error(downsample(1, -1))
  expect_error(downsample(1, 1, -1))
  expect_error(downsample(1, 1, 1))
  expect_error(downsample(1, 2, 3, 4))
})

test_that("downsample() tests are correct", {
  expect_equal(downsample(1:5, 2), c(1, 3, 5))
  expect_equal(downsample(matrix(1:10, 5, byrow = TRUE), 2),
               matrix(c(1, 2, 5, 6, 9, 10), 3, byrow = TRUE))
  expect_equal(downsample(1:5, 2, 1), c(2, 4))
  expect_equal(downsample(matrix(1:10, 5, byrow = TRUE), 2, 1),
               matrix(c(3, 4, 7, 8), 2, byrow = TRUE))
})

# -----------------------------------------------------------------------
# upsample()

test_that("parameters to upsample() are correct", {
  expect_error(upsample())
  expect_error(upsample(1))
  expect_error(upsample(1, -1))
  expect_error(upsample(1, 1, -1))
  expect_error(upsample(1, 1, 1))
  expect_error(upsample(1, 2, 3, 4))
})

test_that("upsample() tests are correct", {
  expect_equal(upsample(c(1, 3, 5), 2), c(1, 0, 3, 0, 5, 0))
  expect_equal(upsample(matrix(c(1, 2, 5, 6, 9, 10), 3, byrow = TRUE), 2),
               matrix(c(1, 2, 0, 0, 5, 6, 0, 0, 9, 10, 0, 0), 6, byrow = TRUE))
  expect_equal(upsample(c(2, 4), 2, 1), c(0, 2, 0, 4))
  expect_equal(upsample(matrix(c(3, 4, 7, 8), 2, byrow = TRUE), 2, 1),
               matrix(c(0, 0, 3, 4, 0, 0, 7, 8), 4, byrow = TRUE))
})

# -----------------------------------------------------------------------
# decimate()

test_that("parameters to upsample() are correct", {
  expect_error(decimate())
  expect_error(decimate(1))
  expect_error(decimate(1, -1))
  expect_error(decimate(1, 1, -1))
  expect_error(decimate(1, 2, 3, 4))
  expect_error(decimate(1, 2, 3, "error"))
})

test_that("upsample() tests are correct", {
  expect_equal(round(decimate(1:10, 2), 3), c(1.306, 2.833, 5.037, 6.789, 9.142))
  expect_equal(round(decimate(1:10, 2, ftype = "fir"), 3),
               c(-0.002, -0.002, -0.007, -0.003, -0.014))
  expect_equal(round(decimate(matrix(c(1:10, 1:10), ncol = 2), 2), 3),
               matrix(c(1.306, 2.833, 5.037, 6.789, 9.142,
                        1.306, 2.833, 5.037, 6.789, 9.142), ncol = 2))
})
