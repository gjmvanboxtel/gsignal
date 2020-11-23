# gsignal Utility Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# fracshift()

test_that("parameters to fracshift() are correct", {
  expect_error(fracshift())
  expect_error(fracshift('invalid'))
  expect_error(fracshift(array()))
  expect_error(fracshift(1:10, c(1, 1)))
  expect_error(fracshift(1:10, 'invalid'))
  expect_error(fracshift(1:10, 7, 'invalid'))
})

test_that("fracshift() tests are correct", {
  
  N  <- 1024
  p  <- 6
  q  <- 7
  d1 <- 64
  d2 <- d1 * p/q
  t  <- 128
  ba <- butter (10, .25)
  n <- rep(0, N)
  n[N / 2 + (-t:t)] <- rnorm(2 * t + 1)
  n  <-  filter(ba, n)
  n1 <- fracshift(n, d1)
  n1 <- resample(n1, p, q)
  n2 <- resample(n, p, q)
  n2 <- fracshift(n2, d2)
  err <- abs(n2 - n1)
  expect_equal(max(err), 0, tolerance = 1e-3)

  # test #integer shift similar similar to non-integer
  N <- 1024
  t <- seq(0, 1, length.out = N)
  x <- exp(-t^2 / 2 / 0.25^2) * sin(2 * pi * 10 * t)
  d  <- 10
  y <- fracshift(x, as.integer(d))
  yh <- fracshift(x, as.double(d) + 1e-8)
  expect_equal(y, yh, tolerance = 1e-8)
  
  # test Octave bug #52758
  x <- c(0, 1, 0, 0, 0, 0, 0, 0)
  y <- fracshift(x, 1)
  expect_equal(length(x), length(y))
  
})

