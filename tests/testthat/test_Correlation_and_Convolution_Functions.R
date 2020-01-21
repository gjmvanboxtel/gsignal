# gsignal Standard Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# cconv()

# test_that("parameters to cconv() are correct", {
#   expect_error(cconv())
#   expect_error(cconv('invalid'))
#   expect_error(cconv(1:10, -1))
#   expect_error(cconv(1:10, 'invalid'))
#   expect_error(cconv(1:10, 0, 1))
# })
# 
# test_that("cconv() tests are correct", {
#   N <- 32
#   x <- seq(0, N - 1, 1) / N + 2
#   y <- cconv (x)
#   expect_true(all(abs(y) < 1e-10))
#   
#   N <- 32
#   t <- seq(0, N - 1, 1) / N
#   x <- t * t + 2
#   y <- cconv (x, 2)
#   expect_true(all(abs(y) < 1e-10))
#   
#   N <- 32
#   t <- seq(0, N - 1, 1) / N
#   x <- cbind(t, 4 * t - 3)
#   y <- cconv (x)
#   expect_true(all(abs(y) < 1e-10))
# })
