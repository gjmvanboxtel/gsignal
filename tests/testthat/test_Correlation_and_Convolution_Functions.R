# gsignal Standard Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# cconv()

test_that("parameters to cconv() are correct", {
  expect_error(cconv())
  expect_error(cconv('invalid'))
  expect_error(cconv(1, 1, c(1,1)))
  expect_error(cconv(1, 1, -1))
  expect_error(cconv(1, 1, 'invalid'))
})

test_that("cconv() tests are correct", {
  x <- 1:5
  expect_equal(cconv(x, 1), 1:5)
  expect_equal(cconv(x, c(1, 1)), c(1, 3, 5, 7, 9, 5))
  expect_equal(cconv(x, c(1, 1), 3), c(8, 12, 10))
  
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4)), c(2, 5, 10, 16, 12, 11, 4))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 4), c(14, 16, 14, 16))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 3), c(22, 17, 21))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 2), c(28, 32))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 1), 60)
  
  expect_equal(cconv(x*1i, 1), c(0+1i, 0+2i, 0+3i, 0+4i, 0+5i))
})
