# gsignal filtering functions
library(gsignal)
library(testthat)

tol <- 1e-6

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
  expect_equal(norm(y - x, '2') / norm(x, '2'), 0, tolerance = 5e-6)

  y <- sgolayfilt(x, sgolay(8, 41, 1, dt))
  expect_equal(norm(y - dx, '2') / norm(dx, '2'), 0, tolerance = 5e-6)

  y <- sgolayfilt(x,sgolay(8, 41, 2, dt))
  expect_equal(norm(y - d2x, '2') / norm(d2x, '2'), 0, tolerance = 1e-5)

  y <- sgolayfilt(x, sgolay(8, 41, 3, dt))
  expect_equal(norm(y - d3x, '2') / norm(d3x, '2'), 0, tolerance = 1e-4)
})

# -----------------------------------------------------------------------
# fir2()

test_that("parameters to fir2() are correct", {
  expect_error(fir2())
  expect_error(fir2(-1))
  expect_error(fir2(0.5))
  expect_error(fir2(1))
  expect_error(fir2(1, -1))
  expect_error(fir2(1, 0.5))
  expect_error(fir2(1, 2))
  expect_error(fir2(1, c(0,1), c(1,0), 4, 5, hamming(-1)))
  expect_error(fir2(1, c(0,1), c(1,0), 4, 5, hamming(1), 7))
})

test_that("fir2() tests are correct", {

  # Test that the grid size is rounded up to the next power of 2
  f <- c(0, 0.6, 0.6, 1); m <- c(1, 1, 0, 0)
  b9  <- fir2 (30, f, m, 9)
  b16 <- fir2 (30, f, m, 16)
  b17 <- fir2 (30, f, m, 17)
  b32 <- fir2 (30, f, m, 32)
  expect_equal(b9,  b16, tolerance = tol)
  expect_equal(b17, b32, tolerance = tol)
  expect_false(isTRUE(all.equal(b16, b17, tolerance = tol)))
  
  # Test expected magnitudes of passbands, stopbands, and cutoff frequencies
  f <- c(0, 0.7, 0.7, 1); m <- c(0, 0, 1, 1)
  b <- fir2 (50, f, m)
  h <- abs(freqz (b, c(0, 0.7, 1), fs = 2)$h)
  expect_lte(h[1], 3e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1, tolerance = 2e-3)
  
  f <- c(0, 0.25, 0.25, 0.75, 0.75, 1); m <- c(0, 0, 1, 1, 0, 0)
  b <- fir2 (50, f, m)
  h <- abs (freqz (b, c(0, 0.25, 0.5, 0.75, 1), fs = 2)$h)
  expect_lte(h[1], 3e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1, tolerance = 2e-3)
  expect_lte(h[4], 1 / sqrt (2))
  expect_lte(h[5], 3e-3)
  
  f <- c(0, 0.45, 0.45, 0.55, 0.55, 1); m <- c(1, 1, 0, 0, 1, 1)
  b <- fir2 (50, f, m)
  h <- abs (freqz (b, c(0, 0.45, 0.5, 0.55, 1), fs = 2)$h)
  expect_equal(h[1], 1, tolerance = 2e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1e-1, tolerance = 2e-2)
  expect_lte(h[4], 1 / sqrt (2))
  expect_equal(h[5], 1, tolerance = 2e-3)
  
})

# -----------------------------------------------------------------------
# fir1()

test_that("parameters to fir1() are correct", {
  expect_error(fir1())
  expect_error(fir1(-1))
  expect_error(fir1(0.5))
  expect_error(fir1(1))
  expect_error(fir1(1, -1))
  expect_error(fir1(1, 2))
  expect_error(fir1(1, 0.5, 'invalid'))
  expect_error(fir1(1, 0.5, "low", 'invalid(1)'))
  expect_error(fir1(1, 0.5, "low", hamming(2), 'invalid'))
  expect_error(fir1(1, 0.5, "low", hamming(2), 'scale', 6))
})

test_that("fir1() tests are correct", {
  
  b <- fir1(30, 0.3)
  h <- abs (freqz (b, c(0, 0.3, 1), fs = 2)$h)
  expect_equal(h[1], 1, tolerance = 1e-2)
  expect_true(all(h[2:3] <= 1 / sqrt(2)))
  
  b <- fir1(30, 0.7, "high")
  h <- abs (freqz (b, c(0, 0.7, 1), fs = 2)$h)
  expect_equal(h[3], 1, tolerance = 1e-2)
  expect_true(all(h[1:2] <= c(3e-3, 1 / sqrt(2))))
  
  b <- fir1 (30, c(0.3, 0.7), 'pass')
  h <- abs (freqz (b, c(0, 0.3, 0.5, 0.7, 1), fs = 2)$h)
  expect_equal(h[3], 1, tolerance = 1e-3)
  expect_true(all(h[-3] <= c(3e-3, 1 / sqrt(2), 1 / sqrt(2), 3e-3)))
  
  b <- fir1(50, c(0.3, 0.7), "stop")
  h <- abs (freqz (b, c(0, 0.3, 0.5, 0.7, 1), fs = 2)$h)
  expect_equal(h[c(1, 5)], c(1, 1), tolerance = 1e-3)
  expect_true(all(h[2:4] <= c(1 / sqrt(2), 3e-3, 1 / sqrt(2))))
})

# -----------------------------------------------------------------------
# firls()

test_that("parameters to firls() are correct", {
  expect_error(firls())
  expect_error(firls(-1))
  expect_error(firls(0.5))
  expect_error(firls(1))
  expect_error(firls(1, 2))
  expect_error(firls(1, 2, 3))
  expect_error(firls(1, 0.5, c(1, 1)))
  expect_error(firls(1, c(0.1, 0.5), c(1, 0), c(3, 3)))
  expect_error(firls(1, 2, 3, 4, 5))
})

