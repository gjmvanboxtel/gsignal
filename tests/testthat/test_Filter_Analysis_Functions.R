# gsignal filter analysis functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# freqs()

test_that("parameters to freqs() are correct", {
  expect_error(freqs())
  expect_error(freqs(1))
  expect_error(freqs(1, 2))
  expect_error(freqs(1, 2, 3, plot = 'invalid'))
})

test_that("freqs() tests are correct", {
  h <- freqs(1, 1, 1, plot = FALSE)
  expect_that(h, equals(1 + 0i))

  h <- freqs(c(1, 2), c(1, 1), 1:4, plot = FALSE)
  expect_that(h, equals(c(1.5-0.5i, 1.2-0.4i, 1.1-0.3i, 1.058824-0.235294i), tolerance = 1e-6))
  
})

