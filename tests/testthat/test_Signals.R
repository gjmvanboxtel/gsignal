# gsignal Signals functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# buffer()

test_that("parameters to buffer() are correct", {
  expect_error(buffer())
  expect_error(buffer(x = 1:10, n = 4.1))
  expect_error(buffer(x = 1:10, n = 4, p = 3.1))
  expect_error(buffer(x = 1:10, n = 4, p = 4))
  expect_error(buffer(x = 1:10, n = 4, p = 1, opt = 10:11))
  expect_error(buffer(x = 1:10, n = 4, p = 1, opt = 'badstring'))
  expect_error(buffer(x = 1:10, n = 3, p = -2, opt = 4))
  expect_error(buffer(x = 1:10, n = 4, zopt = 5))
})

test_that("buffer() tests returning only y are correct", {
  expect_that(buffer(1:10, 4), equals(matrix(c(1:10, 0, 0), 4, 3)))
  expect_that(buffer(1:10, 4, 1), equals(matrix(c(0:3, 3:6, 6:9, 9, 10, 0, 0), 4, 4)))
  expect_that(buffer(1:10, 4, 2), equals(matrix(c(0, 0:2, 1:4, 3:6, 5:8, 7:10), 4, 5)))
  expect_that(buffer(1:10, 4, 3), equals(rbind(c(0, 0, 0:7), c(0, 0:8), 0:9, 1:10)))
  expect_that(buffer(1:10, 4, -1), equals(matrix(c(1:4, 6:9), 4, 2)))
  expect_that(buffer(1:10, 4, -2), equals(matrix(c(1:4, 7:10), 4, 2)))
  expect_that(buffer(1:10, 4, -3), equals(matrix(c(1:4, 8:10, 0), 4, 2)))
  expect_that(buffer(1:10, 4, 1, 11), equals(matrix(c(11,1:3,3:6,6:9,9,10,0,0), 4, 4)))
  expect_that(buffer(1:10, 4, 1, 'nodelay'), equals(matrix(c(1:4,4:7,7:10), 4, 3)))
  expect_that(buffer(1:10, 4, 2, 'nodelay'), equals(matrix(c(1:4,3:6,5:8,7:10), 4, 4)))
  expect_that(buffer(1:10, 4, 3, c(11, 12, 13)), equals(rbind(c(11:13, 1:7), c(12:13, 1:8), c(13, 1:9), 1:10)))
  expect_that(buffer(1:10, 4, 3, 'nodelay'), equals(rbind(1:8, 2:9, 3:10, c(4:10, 0))))
  expect_that(buffer(1:11, 4, -2, 1), equals(matrix(c(2:5, 8:11), 4, 2)))
})

test_that("buffer() tests returning y, and z are correct", {
  buf <- buffer(1:12, 4, zopt = TRUE)
  expect_that(buf$y, equals(matrix(1:12, 4, 3)))
  expect_that(buf$z, equals(NULL))
  
  buf <- buffer(1:11, 4, zopt = TRUE)
  expect_that(buf$y, equals(matrix(1:8, 4, 2)))
  expect_that(buf$z, equals(9:11))
  
  buf <- buffer(t(1:12), 4, zopt = TRUE)
  expect_that(buf$y, equals(matrix(1:12, 4, 3)))
  expect_that(buf$z, equals(NULL))
  
  # slightly different from Matlab implementation (column vector)
  # not sure if this matters - find field tests for this situation
  buf <- buffer(t(1:11), 4, zopt = TRUE)
  expect_that(buf$y, equals(matrix(1:8, 4, 2)))
  expect_that(buf$z, equals(9:11))
})

test_that("buffer() tests returning y, z, and opt are correct", {
  buf <- buffer(1:15, 4, -2, 1, zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(2:5,8:11), 4, 2)))
  expect_that(buf$z, equals(c(14,15)))
  expect_that(buf$opt, equals(0L))
  
  buf <- buffer(1:11, 4, -2, 1, zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(2:5,8:11), 4, 2)))
  expect_that(buf$z, equals(NULL))
  expect_that(buf$opt, equals(2))
  
  # slightly different from Matlab implementation (column vector)
  # not sure if this matters - find field tests for this situation
  buf <- buffer(t(1:15), 4, -2, 1, zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(2:5,8:11), 4, 2)))
  expect_that(buf$z, equals(c(14,15)))
  expect_that(buf$opt, equals(0L))
  
  buf <- buffer(t(1:11), 4, -2, 1, zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(2:5,8:11), 4, 2)))
  expect_that(buf$z, equals(NULL))
  expect_that(buf$opt, equals(2))
  
  buf <- buffer(1:11, 5, 2, c(-1,0), zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(-1:3,2:6,5:9), 5, 3)))
  expect_that(buf$z, equals(c(10, 11)))
  expect_that(buf$opt, equals(c(8, 9)))
  
  buf <- buffer(t(1:11), 5, 2, c(-1,0), zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(-1:3,2:6,5:9), 5, 3)))
  expect_that(buf$z, equals(c(10, 11)))
  expect_that(buf$opt, equals(c(8, 9)))
  
  buf <- buffer(t(1:10), 6, 4, zopt = TRUE)
  expect_that(buf$y, equals(matrix(c(rep(0, 4), 1:2, rep(0, 2), 1:4, 1:6, 3:8, 5:10), 6, 5)))
  expect_that(buf$z, equals(NULL))
  expect_that(buf$opt, equals(7:10))
  
})

test_that("buffer() works correctly with continuous buffering", {
  
  # overlap
  data <- buffer(1:1100, 11)
  n <- 4
  p <- 1
  buf <- list(y = NULL, z = NULL, opt = -5)
  for (i in 1:ncol(data)) {
    x <- data[,i]
    buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
  }
  expect_that(buf$y, equals(matrix(c(1089:1092, 1092:1095, 1095:1098), 4, 3)))
  expect_that(buf$z, equals(c(1099, 1100)))
  expect_that(buf$opt, equals(1098))
  
  # underlap
  data <- buffer(1:1100, 11)
  n <- 4
  p <- -2
  buf <- list(y = NULL, z = NULL, opt = 1)
  for (i in 1:ncol(data)) {
    x <- data[,i]
    buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
  }
  expect_that(buf$y, equals(matrix(c(1088:1091, 1094:1097), 4, 2)))
  expect_that(buf$z, equals(1100))
  expect_that(buf$opt, equals(0))
  
})

# -----------------------------------------------------------------------
# chirp()

test_that("parameters to chirp() are correct", {
  expect_error(chirp())
  expect_error(chirp(1, 2, 3, 4, 5, 6, 7))
  expect_error(chirp(0, shape = "foo"))
})

test_that("chirp() works for linear, quadratic and logarithmic shapes", {
  
  t <- seq(0, 5, 0.001)
  y <- chirp (t)
  expect_that(sum(head(y)), equals(5.999952, tolerance = 1e-6))
  expect_that(sum(tail(y)), equals(2.146626e-05, tolerance = 1e-6))
  
  t <- seq(-2, 15, 0.001)
  y <- chirp (t, 400, 10, 100, "quadratic")
  expect_that(sum(head(y)), equals(0.8976858, tolerance = 1e-6))
  expect_that(sum(tail(y)), equals(0.4537373, tolerance = 1e-6))
  
  t <- seq(0, 5, 1/8000)
  y <- chirp (t, 200, 2, 500, "logarithmic")
  expect_that(sum(head(y)), equals(-4.56818, tolerance = 1e-6))
  expect_that(sum(tail(y)), equals(0.8268064, tolerance = 1e-6))
  
})

# -----------------------------------------------------------------------
# cmorwavf()

test_that("parameters to cmorwavf() are correct", {
  expect_error(cmorwavf(n = -1))
  expect_error(cmorwavf(n = 2.5))
  expect_error(cmorwavf(fb = -1))
  expect_error(cmorwavf(fb = 0))
  expect_error(cmorwavf(fc = -1))
  expect_error(cmorwavf(fc = 0))
})

test_that("cmorwavf() works correctly", {
  expect_that(round(mean(Re(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 4), equals(0))
  expect_that(round(mean(Im(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 4), equals(0))
  expect_lt(max(Re(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 1L)
  expect_lt(max(Im(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 1L)
  expect_gt(min(Re(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), -1L)
  expect_gt(min(Im(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), -1L)
})

# -----------------------------------------------------------------------
# diric()

test_that("parameters to diric() are correct", {
  expect_error(diric())
  expect_error(diric(seq(-2*pi, 2*pi, len = 301)))
  expect_error(diric(seq(-2*pi, 2*pi, len = 301), 0))
  expect_error(diric(seq(-2*pi, 2*pi, len = 301), -1))
  expect_error(diric(seq(-2*pi, 2*pi, len = 301), 2.5))
})

# -----------------------------------------------------------------------
# gauspuls()

test_that("parameters to gauspuls() are correct", {
  expect_error(gauspuls())
  expect_error(gauspuls(seq(-2*pi, 2*pi, len = 301), -1))
  expect_error(gauspuls(seq(-2*pi, 2*pi, len = 301), 2, 0))
  expect_error(gauspuls(seq(-2*pi, 2*pi, len = 301), 2, -1))
})

# -----------------------------------------------------------------------
# gmonopuls()

test_that("parameters to gmonopuls() are correct", {
  expect_error(gmonopuls())
  expect_error(gmonopuls(seq(-2*pi, 2*pi, len = 301), -1))
})

# -----------------------------------------------------------------------
# mexihat()

test_that("parameters to mexihat() are correct", {
  expect_error(mexihat(n = -1))
  expect_error(mexihat(n = 2.5))
})

# -----------------------------------------------------------------------
# meyeraux()

test_that("parameters to meyeraux() are correct", {
  expect_error(meyeraux())
})

# -----------------------------------------------------------------------
# morlet()

test_that("parameters to morlet() are correct", {
  expect_error(morlet(n = -1))
  expect_error(morlet(n = 2.5))
})

# -----------------------------------------------------------------------
# pulstran()

test_that("parameters to pulstran() are correct", {
  expect_error(pulstran())
  expect_error(pulstran(NULL))
  expect_error(pulstran(1, 2, 3, 4, 5, 6))
  expect_error(pulstran(d = seq(0, 0.1, 0.01)))
})

test_that("rectpuls() works correctly", {
  t <- seq(0, 1, 0.01)
  d <- seq(0, 1, 0.1)
  expect_that(pulstran(NA, d, 'sin'), equals(NA_integer_))
  expect_that(pulstran(t, NULL, 'sin'), equals(rep(0L, length(t))))
  expect_that(pulstran(seq(0, 0.1, 0.001)), equals(rep(0L, length(seq(0, 0.1, 0.001)))))
  expect_that(length(pulstran(t, d, 'sin')), equals(length(t)))
})

# -----------------------------------------------------------------------
# rectpuls()

test_that("parameters to rectpuls() are correct", {
  expect_error(rectpuls())
  expect_error(rectpuls(NULL, 0.1))
  expect_error(rectpuls(seq(-2*pi, 2*pi, len = 301), -1))
  expect_error(rectpuls(seq(-2*pi, 2*pi, len = 301), 1, 3))
  expect_error(rectpuls(seq(-2*pi, 2*pi, len = 301), 1i))
})

test_that("rectpuls() works correctly", {
  expect_that(rectpuls(0, 0), equals(0))
  expect_that(rectpuls(0, 0.1), equals(1))
  expect_that(rectpuls(rep(0L, 10)), equals(rep(1L, 10)))
  expect_that(rectpuls(-1:1), equals(c(0, 1, 0)))
  expect_that(rectpuls(-5:5, 9), equals(c(0, rep(1L, 9), 0)))
})

# -----------------------------------------------------------------------
# sawtooth()

test_that("parameters to sawtooth() are correct", {
  expect_error(sawtooth())
  expect_error(sawtooth(NULL, 0.1))
  expect_error(sawtooth(0:10, -1))
  expect_error(sawtooth(0:10, 2))
  expect_error(sawtooth(0:10, 1, 3))
  expect_error(sawtooth(0:10, 1i))
})

test_that("sawtooth() works correctly", {
  expect_that(sawtooth(0, 0), equals(1))
  expect_that(sawtooth(0, 1), equals(-1))
  expect_that(sawtooth(rep(0L, 10)), equals(rep(-1L, 10)))
})

# -----------------------------------------------------------------------
# square()

test_that("parameters to square() are correct", {
  expect_error(square())
  expect_error(square(NULL, 1))
  expect_error(square(0:10, -1))
  expect_error(square(0:10, 150))
  expect_error(square(0:10, 1, 3))
  expect_error(square(0:10, 1i))
})

test_that("square() works correctly", {
  expect_that(square(0, 0), equals(-1))
  expect_that(square(0, 1), equals(1))
  expect_that(square(rep(0L, 10)), equals(rep(1L, 10)))
  expect_that(square(1:12, 50), equals(rep(c(rep(1,3), rep(-1, 3)), 2)))
})

# -----------------------------------------------------------------------
# tripuls()

test_that("parameters to tripuls() are correct", {
  expect_error(tripuls())
  expect_error(tripuls(NULL, 1))
  expect_error(tripuls(0:10, c(0,1)))
  expect_error(tripuls(0:10, 1, -2))
  expect_error(tripuls(0:10, 1, 2))
  expect_error(tripuls(0:10, 1i))
})

test_that("tripuls() works correctly", {
  expect_that(tripuls(0, 1), equals(1))
  expect_that(tripuls(rep(0L, 10)), equals(rep(1L, 10)))
})

# -----------------------------------------------------------------------
# shanwavf()

test_that("parameters to shanwavf() are correct", {
  expect_error(shanwavf(n = -1))
  expect_error(shanwavf(n = 2.5))
  expect_error(shanwavf(fb = -1))
  expect_error(shanwavf(fb = 0))
  expect_error(shanwavf(fc = -1))
  expect_error(shanwavf(fc = 0))
})

test_that("shanwavf() works correctly", {
  expect_that(round(mean(Re(shanwavf(-20, 20, 1000, 1.5, 1)$psi)), 4), equals(0))
  expect_that(round(mean(Im(shanwavf(-20, 20, 1000, 1.5, 1)$psi)), 4), equals(0))
})

# -----------------------------------------------------------------------
# shiftdata()

test_that("parameters to shiftdata() are correct", {
  expect_error(shiftdata())
  expect_error(shiftdata(1, 2, 3))
  expect_error(shiftdata(1, 2.5))
  expect_error(shiftdata(1, 2i))
  expect_error(shiftdata(1:5, 2))
  expect_error(shiftdata(array(1:24, c(2,3)), 3))
})

test_that("shiftdata() works correctly", {
  sd <- shiftdata(matrix(1:9, 3, 3, byrow = TRUE), 2)
  expect_that(sd$x, equals(matrix(c(1, 4, 7, 2, 5, 8, 3, 6, 9), 3, 3, byrow = TRUE)))
  expect_that(sd$perm, equals(c(2,1)))
  expect_that(sd$nshifts, equals(NA))
  
  sd <- shiftdata(array(c(27, 63, 67, 42, 48, 74, 11, 5, 93, 15, 34, 70, 23, 60, 54, 81, 28, 38), c(3, 3, 2)), 2)
  expect_that(sd$x, equals(array(c(27, 42, 11, 63, 48, 5, 67, 74, 93, 15, 23, 81, 34, 60, 28, 70, 54, 38), c(3, 3, 2))))
  expect_that(sd$perm, equals(c(2, 1, 3)))
  expect_that(sd$nshifts, equals(NA))

  X <- array(round(runif(4 * 4 * 4 * 4) * 100), c(4, 4, 4, 4))
  Y <- shiftdata(X, 3)
  T <- NULL
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:2) {
        for (l in 1:2) {
          T <- c(T, Y$x[k, i, j, l] - X[i, j, k ,l])
        }
      }
    }
  }
  expect_that(T, equals(rep(0L, length(T))))
})

# -----------------------------------------------------------------------
# unshiftdata()

test_that("parameters to unshiftdata() are correct", {
  expect_error(unshiftdata())
  expect_error(unshiftdata(1, 2, 3))
  expect_error(unshiftdata(1))
  expect_error(unshiftdata(2i))
  expect_error(unshiftdata(list(x = 1:5, perm = 1, nshifts = 0)))
  expect_error(unshiftdata(list(x=array(1:5), perm = 2i, nshifts = 0)))
  expect_error(unshiftdata(list(x=array(1:5), perm = NULL, nshifts = NULL)))
})

test_that("unshiftdata() works correctly", {
  x <- 1:5
  sd <- shiftdata(x)
  x2 <- unshiftdata(sd)
  expect_that(array(x), equals(x2))
  
  x <- array(round(runif(3 * 3) * 100), c(3, 3))
  sd <- shiftdata(x, 2)
  x2 <- unshiftdata(sd)
  expect_that(x, equals(x2))
  
  x <- array(round(runif(4 * 4 * 4 * 4) * 100), c(4, 4, 4, 4))
  sd <- shiftdata(x, 3)
  x2 <- unshiftdata(sd)
  expect_that(x, equals(x2))

  x <- array(round(runif(1 * 1 * 3 * 4) * 100), c(1, 1, 3, 4))
  sd <- shiftdata(x)
  x2 <- unshiftdata(sd)
  expect_that(x, equals(x2))
  
})

# -----------------------------------------------------------------------
# sigmoid_train()

test_that("parameters to sigmoid_train() are correct", {
  expect_error(sigmoid_train())
  expect_error(sigmoid_train(1:10, NULL, NULL))
  expect_error(sigmoid_train(1:10, c(1,2), NULL))
  expect_error(sigmoid_train(1:10, rbind(c(1,2),1), NULL))
  expect_error(sigmoid_train(1:10, rbind(c(1,2),1), 2i))
})

test_that("sigmoid_train() works correctly", {
  st <- sigmoid_train(1:10, rbind(c(2,3)), 1)
  expect_that(st$y, equals(st$s))
})

# -----------------------------------------------------------------------
# specgram()

test_that("parameters to specgram() are correct", {
  expect_error(specgram())
  expect_error(specgram(matrix(1:10, 2, 5)))
  expect_error(specgram(x = 1:10, n = 4.1))
  expect_warning(specgram(x = 1:10, n = 11))
  expect_warning(specgram(x = 1:10, n = 2, window = 1:11))
  expect_error(specgram(x = 1:10, n = 2, overlap = 3))
})

test_that("buffer() tests returning only y are correct", {
  sp <- specgram(chirp(seq(-2, 15, by = 0.001), 400, 10, 100, 'quadratic'), plot = FALSE)
  expect_that(length(sp$f), equals(128))
  expect_that(length(sp$t), equals(131))
  expect_that(nrow(sp$S), equals(length(sp$f)))
  expect_that(ncol(sp$S), equals(length(sp$t)))
})
