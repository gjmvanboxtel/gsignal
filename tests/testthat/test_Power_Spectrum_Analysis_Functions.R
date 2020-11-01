# gsignal Power Spectrum Analysis Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# pwelch()

test_that("parameters to pwelch() are correct", {
  expect_error(pwelch())
  expect_error(pwelch('a'))
  expect_error(pwelch(1:10, -1))
  expect_error(pwelch(1:10, 4, -1))
  expect_error(pwelch(1:10, 4, 0.5, -1))
  expect_error(pwelch(1:10, 4, 0.5, 12, -1))
  expect_error(pwelch(1:10, 4, 0.5, 12, pi, 'invalid'))
  expect_error(pwelch(1:10, 4, 0.5, 12, pi, 'none', 'unused'))
})

test_that("pwelch() tests are correct", {
  fs <- 1000
  secs <- 10
  freq <- 30
  A <- 1
  t <- seq(0, secs, length.out = fs * secs)
  x <- A * cos(freq * 2 * pi * t)
  Pxx <- pwelch(x, fs = fs)
  expect_equal(length(Pxx$freq), 65L)
  expect_equal(length(Pxx$spec), 65L)
  expect_equal(round(Pxx$freq[which.max(Pxx$spec)], -1), 30L)
  expect_null(Pxx$cross)
  expect_null(Pxx$phase)
  expect_null(Pxx$coh)
  expect_null(Pxx$trans)
  expect_equal(Pxx$x_len, 9984L)
  expect_equal(Pxx$seg_len, 128L)
  expect_equal(Pxx$psd_len, 65L)
  expect_equal(Pxx$nseries, 1L)
  
  y <- A * sin(freq * 2 * pi * t)
  Pxy <- pwelch(cbind(x, y), fs = fs)
  
  expect_equal(length(Pxy$freq), 65L)
  expect_equal(dim(Pxy$spec), c(65L, 2L))

  expect_equal(dim(Pxy$cross), c(65L, 1L))
  expect_equal(round(Pxy$freq[which.max(Pxy$cross[, 1])], -1), 30L)
  expect_equal(colnames(Pxy$cross), "x-y")
  
  expect_equal(dim(Pxy$phase), c(65L, 1L))
  expect_equal(unname(Pxy$phase[which(Pxy$freq>30)[1], 1]), pi/2, tolerance = 1e-4)
  expect_equal(colnames(Pxy$phase), "x-y")
  
  expect_equal(dim(Pxy$coh), c(65L, 1L))
  expect_equal(unname(Pxy$coh[which(Pxy$freq>30)[1], 1]), 1, tolerance = 1e-4)
  expect_equal(colnames(Pxy$coh), "x-y")
  
  expect_equal(dim(Pxy$trans), c(65L, 1L))
  expect_equal(unname(abs(Pxy$trans[which(Pxy$freq>30)[1], 1])), 1, tolerance = 1e-4)
  expect_equal(colnames(Pxy$trans), "x-y")
  
  expect_equal(Pxx$x_len, 9984L)
  expect_equal(Pxy$seg_len, 128L)
  expect_equal(Pxy$psd_len, 65L)
  expect_equal(Pxy$nseries, 2L)

})

# -----------------------------------------------------------------------
# ar_psd()

test_that("parameters to ar_psd() are correct", {
  expect_error(ar_psd())
  expect_error(ar_psd('a'))
  expect_error(ar_psd(c(0,0)))
  expect_error(ar_psd(1:10, -1))
  expect_error(ar_psd(1:10, 4, -1))
  expect_error(ar_psd(1:10, 4, 2, -1))
  expect_error(ar_psd(1:10, 4, 2, 1, 'invalid'))
  expect_error(ar_psd(1:10, 4, 2, 1, 'whole', 'invalid'))
  expect_error(ar_psd(1:10, 4, 2, 1, 'whole', 'fft', 7))
})

test_that("ar_psd() tests are correct", {
  
  psd <- ar_psd(c(1,0), 1)
  expect_equal(psd$freq, (1 / 2 / 256) * seq(0, 255))
  expect_equal(psd$psd, rep(2L, 256))

  n <- 64
  psd <- ar_psd(c(1,0, 0), 1, n)
  expect_equal(psd$freq, (1 / 2 / n) * seq(0, n - 1))
  expect_equal(psd$psd, rep(2L, n))

  psd <- ar_psd(c(1,0, 1), 1, n)
  expect_equal(which.max(psd$psd), (n / 2) + 1)

  psd <- ar_psd(c(1,0, 1), 1, n, range = "whole")
  expect_equal(which.max(psd$psd), (n / 4) + 1)

  psd <- ar_psd(c(1,0, 1), 1, n, range = "centerdc")
  expect_equal(which(is.infinite(psd$psd)), c(17, 49))
  
})


