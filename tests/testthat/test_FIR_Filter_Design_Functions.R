# gsignal filtering functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# sgolay()

test_that("parameters to sgolay() are correct", {
  expect_error(sgolay())
  expect_error(sgolay(-1))
  expect_error(sgolay(0.5))
  expect_error(sgolay(1))
  expect_error(sgolay(1, -1))
  expect_error(sgolay(1, 0.5))
  expect_error(sgolay(1, 2))
  expect_error(sgolay(1, 2, 3, 4, 5))
})

test_that("sgolay() tests are correct", {
  N <- 2^12
  t <- seq(0, N-1) / N
  dt <- t[2] - t[1]
  w <- 2 * pi * 50
  offset <- 0.5 # 50 Hz carrier
  # exponential modulation and its derivatives
  d <- 1 + exp(-3 * (t - offset))
  dd <- -3 * exp(-3 * (t - offset))
  d2d <- 9 * exp(-3 * (t - offset))
  d3d <- -27 * exp(-3 * (t -offset))
  # modulated carrier and its derivatives
  x <- d * sin(w * t)
  dx <- dd * sin(w * t) + w * d * cos(w * t)
  d2x <- (d2d - w^2 * d) * sin(w * t) + 2 * w * dd * cos(w * t)
  d3x <- (d3d - 3 * w^2 * dd) * sin(w * t) + (3 * w * d2d - w^3 * d) * cos(w * t)

  y <- sgolayfilt(x, sgolay(8, 41, 0, dt))
  expect_that(norm(y - x, '2') / norm(x, '2'), equals(0, tolerance = 5e-6))

  y <- sgolayfilt(x, sgolay(8, 41, 1, dt))
  expect_that(norm(y - dx, '2') / norm(dx, '2'), equals(0, tolerance = 5e-6))

  y <- sgolayfilt(x,sgolay(8, 41, 2, dt))
  expect_that(norm(y - d2x, '2') / norm(d2x, '2'), equals(0, tolerance = 1e-5))

  y <- sgolayfilt(x, sgolay(8, 41, 3, dt))
  expect_that(norm(y - d3x, '2') / norm(d3x, '2'), equals(0, tolerance = 1e-4))
})
