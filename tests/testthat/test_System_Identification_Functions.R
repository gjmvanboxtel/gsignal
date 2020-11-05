# gsignal System Identification Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# arburg()

test_that("parameters to arburg() are correct", {
  expect_error(arburg())
  expect_error(arburg('a'))
  expect_error(arburg(c(0,0)))
  expect_error(arburg(1:10, -1))
  expect_error(arburg(1:10, 9))
  expect_error(arburg(1:10, 7, 'invalid'))
})

test_that("arburg() tests are correct", {
  
  cf <- arburg(rep(0L, 5), 1)
  expect_equal(cf$a, c(1, NaN))
  expect_equal(cf$e, NaN)
  expect_equal(cf$k, NaN)
  
  cf <- arburg(c(1L, rep(0L, 4)), 1)
  expect_equal(cf$a, c(1, 0))
  expect_equal(cf$e, 0.2)
  expect_equal(cf$k, 0)

  cf <- arburg(c(1L, 0L, 1L, 0L, 0L), 1)
  expect_equal(cf$a, c(1, 0))
  expect_equal(cf$e, 0.4)
  expect_equal(cf$k, 0)

  cf <- arburg(c(1L, 0L, 1L, 0L, 0L), 2)
  expect_equal(cf$a, c(1, 0, -2 / 3))
  expect_equal(cf$e, 2 / 9)
  expect_equal(cf$k, c(0, -2 / 3))

  x <- filter(1, c(1, -0.75, 0.5), 0.2 * rnorm(1024))
  y <- cbind(x, x)
  cf <- arburg(y, 2)
  expect_equal(ncol(cf$a), 3)
  expect_equal(nrow(cf$a), 2)
  expect_equal(length(cf$e), 2)
  expect_equal(ncol(cf$k), 2)
  expect_equal(nrow(cf$k), 2)
  expect_equal(cf$a[1, ], cf$a[2, ])
  expect_equal(cf$e[1], cf$e[2])
  expect_equal(cf$k[, 1], cf$k[, 1])
  
})

