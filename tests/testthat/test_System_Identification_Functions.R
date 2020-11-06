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

# -----------------------------------------------------------------------
# levinson()

test_that("parameters to levinson() are correct", {
  expect_error(levinson())
  expect_error(levinson('invalid'))
  expect_error(levinson(1:10, -1))
  expect_error(levinson(1:10, 'invalid'))
  expect_error(levinson(1:10, 7, 'invalid'))
})

test_that("levinson() tests are correct", {
  
  cf <- levinson(rep(0L, 5), 1)
  expect_equal(cf$a, c(1, NaN))
  expect_equal(cf$e, NaN)
  expect_equal(cf$k, NaN)
  
  cf <- levinson(c(1L, rep(0L, 4)), 1)
  expect_equal(cf$a, c(1, 0))
  expect_equal(cf$e, 1)
  expect_equal(cf$k, 0)
  

  cf <- levinson(c(1L, 0L, 1L, 0L, 0L), 2)
  expect_equal(cf$a, c(1, 0, -1))
  expect_equal(cf$e, 0)
  expect_equal(cf$k, c(0, -1))
  
  x <- filter(1, c(1, -0.75, 0.5), 0.2 * rnorm(1024))
  y <- cbind(x, x)
  cf <- levinson(y, 2)
  expect_equal(ncol(cf$a), 3)
  expect_equal(nrow(cf$a), 2)
  expect_equal(length(cf$e), 2)
  expect_equal(ncol(cf$k), 2)
  expect_equal(nrow(cf$k), 2)
  expect_equal(cf$a[1, ], cf$a[2, ])
  expect_equal(cf$e[1], cf$e[2])
  expect_equal(cf$k[, 1], cf$k[, 1])
  
})

# -----------------------------------------------------------------------
# aryule()

test_that("parameters to aryule() are correct", {
  expect_error(aryule())
  expect_error(aryule('invalid'))
  expect_error(aryule(c(0,0)))
  expect_error(aryule(1:10, -1))
  expect_error(aryule(1:10, 'invalid'))
  expect_error(aryule(1:10, 7, 'invalid'))
})

test_that("aryule() tests are correct", {
  
  cf <- aryule(rep(0L, 5), 1)
  expect_equal(cf$a, c(1, NaN))
  expect_equal(cf$e, NaN)
  expect_equal(cf$k, NaN)
  
  cf <- aryule(c(1L, rep(0L, 4)), 1)
  expect_equal(cf$a, c(1, 0))
  expect_equal(cf$e, 0.2)
  expect_equal(cf$k, 0)
  
  cf <- aryule(c(1L, 0L, 1L, 0L, 0L), 2)
  expect_equal(cf$a, c(1, 0, -0.5))
  expect_equal(cf$e, 0.3)
  expect_equal(cf$k, c(0, -0.5))
  
  x <- filter(1, c(1, -0.75, 0.5), 0.2 * rnorm(1024))
  y <- cbind(x, x)
  cf <- aryule(y, 2)
  expect_equal(ncol(cf$a), 3)
  expect_equal(nrow(cf$a), 2)
  expect_equal(length(cf$e), 2)
  expect_equal(ncol(cf$k), 2)
  expect_equal(nrow(cf$k), 2)
  expect_equal(cf$a[1, ], cf$a[2, ])
  expect_equal(cf$e[1], cf$e[2])
  expect_equal(cf$k[, 1], cf$k[, 1])
  
})

