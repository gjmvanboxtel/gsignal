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

# -----------------------------------------------------------------------
# clustersegment()

test_that("parameters to clustersegment() are correct", {
  expect_error(clustersegment())
  expect_error(clustersegment('invalid'))
  expect_error(clustersegment(array()))
  expect_error(clustersegment(1, 2))
})

test_that("clustersegment() tests are correct", {
  
  x <- rep(0L, 5)
  rng <- clustersegment(x)
  expect_equal(rng, matrix(NA, 1))
  
  x <- rep(10L, 5)
  rng <- clustersegment(x)
  expect_equal(rng, matrix(NA, 1))
  
  x <- c(rep(0L, 5), 2)
  rng <- clustersegment(x)
  expect_equal(rng, matrix(c(6, 6), nrow = 2))
  
  x <- c(2, rep(0L, 5))
  rng <- clustersegment(x)
  expect_equal(rng, matrix(c(1, 1), nrow = 2))
  
  x <- matrix(c(1, 1, 1, 0, 0, 1, 1, 1, 1, 1,
                1, 1, 1, 0, 1, 0, 0, 1, 1, 1,
                1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
              nrow = 3, byrow = TRUE)
  rng <- clustersegment(x)
  expect_equal(rng[[1]], matrix(c(1, 3, 6, 10), nrow = 2))
  expect_equal(rng[[2]], matrix(c(1, 3, 5, 5, 8, 10), nrow = 2))
  expect_equal(rng[[3]], matrix(c(1, 1, 9, 9), nrow = 2))
  
})

# -----------------------------------------------------------------------
# schtrig()

test_that("parameters to schtrig() are correct", {
  expect_error(schtrig())
  expect_error(schtrig('invalid'))
  expect_error(schtrig(array()))
  expect_error(schtrig(1, 2, 3, 4))
})

test_that("schtrig() tests are correct", {
  
  x <- c(0, 0.5, 1, 1.5, 2, 1.5, 1.5, 1.2, 1, 0, 0)
  y <- schtrig(x, c(1.3, 1.6))
  expect_equal(y$v, c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0))
  expect_equal(y$rng, matrix(c(4, 4, 6, 7), nrow = 2))
  expect_equal(y$st, 0)

  x <- c(0, 0.5, 1, 1.5, 2, 1.5, 1.5, 1.2, 1, 0, 0)
  y <- schtrig(x, c(1.3, 1.6, 1.8))
  expect_equal(y$v, c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0))
  expect_equal(y$rng, matrix(NA, nrow = 1))
  expect_equal(y$st, 0)
  
  expect_equal(schtrig(x, c(1.3, 1.6)), schtrig(x, c(1.3, 1.6, 1.3)))
})

# -----------------------------------------------------------------------
# upsamplefill()

test_that("parameters to upsamplefill() are correct", {
  expect_error(upsamplefill())
  expect_error(upsamplefill('invalid'))
  expect_error(upsamplefill(array()))
  expect_error(upsamplefill(1, 'invalid'))
  expect_error(upsamplefill(1, -1, 'invalid'))
  expect_error(upsamplefill(1, -1, TRUE))
  
  expect_error(upsamplefill(1, 2, 3, 4))
})

test_that("upsamplefill() tests are correct", {
  expect_equal(upsamplefill(c(1, 3, 5), 2), c(1, 2, 3, 2, 5, 2))
  expect_equal(upsamplefill(c(1, 2, 5), c(2, -2)), c(1, 2, -2, 2, 2, -2, 5, 2, -2))
  expect_equal(upsamplefill(diag(2), 2, TRUE), matrix(c(rep(1, 3), rep(0, 6), rep(1, 3)), ncol = 2))
  expect_equal(upsamplefill(c(1, 3, 5), 2, TRUE), c(1, 1, 1, 3, 3, 3, 5, 5, 5))
})

# -----------------------------------------------------------------------
# wkeep()

test_that("parameters to wkeep() are correct", {
  expect_error(wkeep())
  expect_error(wkeep('invalid'))
  expect_error(wkeep(1, 'invalid'))
  expect_error(wkeep(1, 1, 'invalid'))
  expect_error(wkeep(1, -1, 'invalid'))
  expect_error(wkeep(1:10, 11))
  expect_error(wkeep(1:10, 6, 7))
  expect_error(wkeep(1:10, 6, -1))
  expect_error(wkeep(matrix(1:10), 1))
  expect_error(wkeep(matrix(1:10), c(1,2), 1))
  expect_error(wkeep(1, 2, 3, 4))
})

test_that("wkeep() tests are correct", {
  expect_equal(wkeep(1:10, 2), c(5, 6))
  expect_equal(wkeep(1:10, 2, "l"), c(1, 2))
  expect_equal(wkeep(1:10, 2, "r"), c(9, 10))
  m <- matrix(c(17, 23, 4, 10, 11, 24, 5, 6, 12, 18, 1, 7, 13,
                19, 25, 8, 14, 20, 21, 2, 15, 16, 22, 3, 9 ), 5, 5)
  expect_equal(wkeep(m, c(3, 2)), matrix(c(5, 6, 12, 7, 13, 19), 3))
  expect_equal(wkeep(m, c(2, 4), c(3, 1)), matrix(c(4, 10, 6, 12, 13, 19, 20, 21), 2))
})

# -----------------------------------------------------------------------
# zerocrossing()

test_that("parameters to zerocrossing() are correct", {
  expect_error(zerocrossing())
  expect_error(zerocrossing(1))
  expect_error(zerocrossing('invalid', 1))
  expect_error(zerocrossing(1, 'invalid'))
  expect_error(zerocrossing(1:2, 1:3))
  expect_error(zerocrossing(1, 1, 3))
})

test_that("zerocrossing() tests are correct", {
  x <- 1:10
  y <- sawtooth(x)
  expect_equal(floor(zerocrossing(x, y)), c(3, 6, 9))
  y <- sin(x)
  expect_equal(floor(zerocrossing(x, y)), c(3, 6, 9))
  y <- cos(x)
  expect_equal(floor(zerocrossing(x, y)), c(1, 4, 7))
})

