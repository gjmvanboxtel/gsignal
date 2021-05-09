# gsignal Signals functions
library(gsignal)
library(testthat)

tol <- 1e-6

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
  expect_equal(buffer(1:10, 4), matrix(c(1:10, 0, 0), 4, 3))
  expect_equal(buffer(1:10, 4, 1), matrix(c(0:3, 3:6, 6:9, 9, 10, 0, 0), 4, 4))
  expect_equal(buffer(1:10, 4, 2), matrix(c(0, 0:2, 1:4, 3:6, 5:8, 7:10), 4, 5))
  expect_equal(buffer(1:10, 4, 3), rbind(c(0, 0, 0:7), c(0, 0:8), 0:9, 1:10))
  expect_equal(buffer(1:10, 4, -1), matrix(c(1:4, 6:9), 4, 2))
  expect_equal(buffer(1:10, 4, -2), matrix(c(1:4, 7:10), 4, 2))
  expect_equal(buffer(1:10, 4, -3), matrix(c(1:4, 8:10, 0), 4, 2))
  expect_equal(buffer(1:10, 4, 1, 11), matrix(c(11,1:3,3:6,6:9,9,10,0,0), 4, 4))
  expect_equal(buffer(1:10, 4, 1, 'nodelay'), matrix(c(1:4,4:7,7:10), 4, 3))
  expect_equal(buffer(1:10, 4, 2, 'nodelay'), matrix(c(1:4,3:6,5:8,7:10), 4, 4))
  expect_equal(buffer(1:10, 4, 3, c(11, 12, 13)),
               rbind(c(11:13, 1:7), c(12:13, 1:8), c(13, 1:9), 1:10))
  expect_equal(buffer(1:10, 4, 3, 'nodelay'), rbind(1:8, 2:9, 3:10, c(4:10, 0)))
  expect_equal(buffer(1:11, 4, -2, 1), matrix(c(2:5, 8:11), 4, 2))
})

test_that("buffer() tests returning y, and z are correct", {
  buf <- buffer(1:12, 4, zopt = TRUE)
  expect_equal(buf$y, matrix(1:12, 4, 3))
  expect_equal(buf$z, NULL)
  
  buf <- buffer(1:11, 4, zopt = TRUE)
  expect_equal(buf$y, matrix(1:8, 4, 2))
  expect_equal(buf$z, 9:11)
  
  buf <- buffer(t(1:12), 4, zopt = TRUE)
  expect_equal(buf$y, matrix(1:12, 4, 3))
  expect_equal(buf$z, NULL)
  
  # slightly different from Matlab implementation (column vector)
  # not sure if this matters - find field tests for this situation
  buf <- buffer(t(1:11), 4, zopt = TRUE)
  expect_equal(buf$y, matrix(1:8, 4, 2))
  expect_equal(buf$z, 9:11)
})

test_that("buffer() tests returning y, z, and opt are correct", {
  buf <- buffer(1:15, 4, -2, 1, zopt = TRUE)
  expect_equal(buf$y, matrix(c(2:5,8:11), 4, 2))
  expect_equal(buf$z, c(14,15))
  expect_equal(buf$opt, 0L)
  
  buf <- buffer(1:11, 4, -2, 1, zopt = TRUE)
  expect_equal(buf$y, matrix(c(2:5,8:11), 4, 2))
  expect_equal(buf$z, NULL)
  expect_equal(buf$opt, 2)
  
  # slightly different from Matlab implementation (column vector)
  # not sure if this matters - find field tests for this situation
  buf <- buffer(t(1:15), 4, -2, 1, zopt = TRUE)
  expect_equal(buf$y, matrix(c(2:5,8:11), 4, 2))
  expect_equal(buf$z, c(14,15))
  expect_equal(buf$opt, 0L)
  
  buf <- buffer(t(1:11), 4, -2, 1, zopt = TRUE)
  expect_equal(buf$y, matrix(c(2:5,8:11), 4, 2))
  expect_equal(buf$z, NULL)
  expect_equal(buf$opt, 2)
  
  buf <- buffer(1:11, 5, 2, c(-1,0), zopt = TRUE)
  expect_equal(buf$y, matrix(c(-1:3,2:6,5:9), 5, 3))
  expect_equal(buf$z, c(10, 11))
  expect_equal(buf$opt, c(8, 9))
  
  buf <- buffer(t(1:11), 5, 2, c(-1,0), zopt = TRUE)
  expect_equal(buf$y, matrix(c(-1:3,2:6,5:9), 5, 3))
  expect_equal(buf$z, c(10, 11))
  expect_equal(buf$opt, c(8, 9))
  
  buf <- buffer(t(1:10), 6, 4, zopt = TRUE)
  expect_equal(buf$y, matrix(c(rep(0, 4), 1:2, rep(0, 2), 1:4, 1:6, 3:8, 5:10), 6, 5))
  expect_equal(buf$z, NULL)
  expect_equal(buf$opt, 7:10)
  
})

test_that("buffer() works correctly with continuous buffering", {
  
  # overlap
  data <- buffer(1:1100, 11)
  n <- 4
  p <- 1
  buf <- list(y = NULL, z = NULL, opt = -5)
  for (i in seq_len(ncol(data))) {
    x <- data[,i]
    buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
  }
  expect_equal(buf$y, matrix(c(1089:1092, 1092:1095, 1095:1098), 4, 3))
  expect_equal(buf$z, c(1099, 1100))
  expect_equal(buf$opt, 1098)
  
  # underlap
  data <- buffer(1:1100, 11)
  n <- 4
  p <- -2
  buf <- list(y = NULL, z = NULL, opt = 1)
  for (i in seq_len(ncol(data))) {
    x <- data[,i]
    buf <- buffer(x = c(buf$z,x), n, p, opt=buf$opt, zopt = TRUE)
  }
  expect_equal(buf$y, matrix(c(1088:1091, 1094:1097), 4, 2))
  expect_equal(buf$z, 1100)
  expect_equal(buf$opt, 0)
  
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
  expect_equal(sum(head(y)), 5.999952, tolerance = tol)
  expect_equal(sum(tail(y)), 2.146626e-05, tolerance = tol)
  
  t <- seq(-2, 15, 0.001)
  y <- chirp (t, 400, 10, 100, "quadratic")
  expect_equal(sum(head(y)), 0.8976858, tolerance = tol)
  expect_equal(sum(tail(y)), 0.4537373, tolerance = tol)
  
  t <- seq(0, 5, 1/8000)
  y <- chirp (t, 200, 2, 500, "logarithmic")
  expect_equal(sum(head(y)), -4.56818, tolerance = tol)
  expect_equal(sum(tail(y)), 0.8268064, tolerance = tol)
  
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
  expect_equal(round(mean(Re(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 4), 0)
  expect_equal(round(mean(Im(cmorwavf(-8, 8, 1000, 1.5, 1)$psi)), 4), 0)
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
  expect_equal(pulstran(NA, d, 'sin'), NA_integer_)
  expect_equal(pulstran(t, NULL, 'sin'), rep(0L, length(t)))
  expect_equal(pulstran(seq(0, 0.1, 0.001)), rep(0L, length(seq(0, 0.1, 0.001))))
  expect_equal(length(pulstran(t, d, 'sin')), length(t))
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
  expect_equal(rectpuls(0, 0), 0L)
  expect_equal(rectpuls(0, 0.1), 1L)
  expect_equal(rectpuls(rep(0L, 10)), rep(1L, 10))
  expect_equal(rectpuls(-1:1), c(0, 1, 0))
  expect_equal(rectpuls(-5:5, 9), c(0, rep(1L, 9), 0))
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
  expect_equal(sawtooth(0, 0), 1L)
  expect_equal(sawtooth(0, 1), -1L)
  expect_equal(sawtooth(rep(0L, 10)), rep(-1L, 10))
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
  expect_equal(square(0, 0), -1L)
  expect_equal(square(0, 1), 1L)
  expect_equal(square(rep(0L, 10)), rep(1L, 10))
  expect_equal(square(1:12, 50), rep(c(rep(1,3), rep(-1, 3)), 2))
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
  expect_equal(tripuls(0, 1), 1L)
  expect_equal(tripuls(rep(0L, 10)), rep(1L, 10))
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
  expect_equal(mean(Re(shanwavf(-20, 20, 1000, 1.5, 1)$psi)), 0, tolerance = 1e-3)
  expect_equal(mean(Im(shanwavf(-20, 20, 1000, 1.5, 1)$psi)), 0, tolerance = 1e-3)
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
  expect_equal(sd$x, matrix(c(1, 4, 7, 2, 5, 8, 3, 6, 9), 3, 3, byrow = TRUE))
  expect_equal(sd$perm, c(2,1))
  expect_equal(sd$nshifts, NA)
  
  sd <- shiftdata(array(c(27, 63, 67, 42, 48, 74, 11, 5, 93, 15, 34, 70, 23, 60, 54, 81, 28, 38), c(3, 3, 2)), 2)
  expect_equal(sd$x, array(c(27, 42, 11, 63, 48, 5, 67, 74, 93, 15, 23, 81, 34, 60, 28, 70, 54, 38), c(3, 3, 2)))
  expect_equal(sd$perm, c(2, 1, 3))
  expect_equal(sd$nshifts, NA)

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
  expect_equal(T, rep(0L, length(T)))
})

# -----------------------------------------------------------------------
# unshiftdata()

test_that("parameters to unshiftdata() are correct", {
  expect_error(unshiftdata())
  expect_error(unshiftdata(1, 2, 3))
  expect_error(unshiftdata(1))
  expect_error(unshiftdata(2i))
  #expect_error(unshiftdata(list(x = 1:5, perm = 1, nshifts = 0)))
  expect_error(unshiftdata(list(x=array(1:5), perm = 2i, nshifts = 0)))
  expect_error(unshiftdata(list(x=array(1:5), perm = NULL, nshifts = NULL)))
})

test_that("unshiftdata() works correctly", {
  x <- 1:5
  sd <- shiftdata(x)
  x2 <- unshiftdata(sd)
  expect_equal(array(x), x2)
  
  x <- array(round(runif(3 * 3) * 100), c(3, 3))
  sd <- shiftdata(x, 2)
  x2 <- unshiftdata(sd)
  expect_equal(x, x2)
  
  x <- array(round(runif(4 * 4 * 4 * 4) * 100), c(4, 4, 4, 4))
  sd <- shiftdata(x, 3)
  x2 <- unshiftdata(sd)
  expect_equal(x, x2)

  x <- array(round(runif(1 * 1 * 3 * 4) * 100), c(1, 1, 3, 4))
  sd <- shiftdata(x)
  x2 <- unshiftdata(sd)
  expect_equal(x, x2)
  
})

# -----------------------------------------------------------------------
# sigmoid_train()

test_that("parameters to sigmoid_train() are correct", {
  expect_error(sigmoid_train())
  expect_error(sigmoid_train(1:10, NULL, NULL))
  expect_error(sigmoid_train(1:10, rbind(c(1,2),1), NULL))
  expect_error(sigmoid_train(1:10, rbind(c(1,2),1), 2i))
})

test_that("sigmoid_train() works correctly", {
  st <- sigmoid_train(1:10, rbind(c(2,3)), 1)
  expect_equal(st$y, st$s, tolerance = tol)
  st <- sigmoid_train(1:10, c(2,3), 1)
  expect_equal(st$y, st$s, tolerance = tol)
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

test_that("specgram() works correctly", {
  sp <- specgram(chirp(seq(-2, 15, by = 0.001), 400, 10, 100, 'quadratic'))
  expect_equal(length(sp$f), 128L)
  expect_equal(length(sp$t), 131L)
  expect_equal(nrow(sp$S), length(sp$f))
  expect_equal(ncol(sp$S), length(sp$t))
})

# -----------------------------------------------------------------------
# uencode()

test_that("parameters to uencode() are correct", {
  expect_error(uencode())
  expect_error(uencode(1))
  expect_error(uencode(1, 2, 3, 4, 5))
  expect_error(uencode(1, 100))
  expect_error(uencode(1, 4, 0))
  expect_error(uencode(1, 4, -1))
  expect_error(uencode(1, 4, 2, 'invalid'))
})

test_that("uencode() works correctly", {
  expect_equal(uencode(seq(-3, 3, 0.5), 2), 
               c(0, 0, 0, 0, 0, 1, 2, 3, 3, 3, 3, 3, 3))
  expect_equal(uencode(seq(-4, 4, 0.5), 3, 4), 
               c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7))
  expect_equal(uencode(seq(-8, 8, 0.5), 4, 8, FALSE),
              c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10,
                10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, 15))
  expect_equal(uencode(seq(-8, 8, 0.5), 4, 8, TRUE),
              c(-8, -8, -7, -7, -6, -6, -5, -5, -4, -4, -3, -3, -2, -2, -1, -1, 0,
                0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7))
  expect_equal(uencode(matrix(c(-2, 1, -1, 2), 2, 2), 2),
               matrix(c(0, 3, 0, 3), 2, 2))
  expect_equal(uencode(matrix(c(1+1i, 2+1i, 3+1i, 4+2i, 5+2i, 6+2i, 7+3i, 8+3i, 9+3i), 3, 3, byrow = TRUE), 2),
              matrix(rep(3, 9), 3, 3))
})

# -----------------------------------------------------------------------
# udecode()

test_that("parameters to udecode() are correct", {
  expect_error(udecode())
  expect_error(udecode(1))
  expect_error(udecode(1, 2, 3, 4, 5))
  expect_error(udecode(1, 100))
  expect_error(udecode(1, 4, 0))
  expect_error(udecode(1, 4, -1))
  expect_error(udecode(1, 4, 2, 'invalid'))
})

test_that("udecode() works correctly", {
  expect_equal(udecode(c(rep(0, 5), 1, 2, rep(3, 6)), 2),
               c(-1, -1, -1, -1, -1, -0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
  expect_equal(udecode(0:10, 2, 1, TRUE),
               c(-1, -0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
  expect_equal(udecode(0:10, 2, 1, FALSE),
               c(-1, -0.5, 0, 0.5, -1, -0.5, 0, 0.5, -1, -0.5, 0))
  expect_equal(udecode(-4:3, 3, 2),  c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5))
  expect_equal(udecode(-7:7, 3, 2, TRUE),
               c(-2, -2, -2, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 1.5, 1.5, 1.5, 1.5))
  expect_equal(udecode(-7:7, 3, 2, FALSE),
               c(0.5, 1, 1.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, -2, -1.5, -1, -0.5))
  expect_equal(udecode(matrix(c(-2, 1, -1, 2), 2, 2), 2),
               matrix(c(-1, 0.5, -0.5, 0.5), 2, 2))
  expect_equal(udecode(matrix(c(1+1i, 2+1i, 3+1i, 4+2i, 5+2i, 6+2i, 7+3i, 8+3i, 9+3i), 3, 3, byrow = TRUE), 2),
               matrix(complex(real = c(-0.5, 0.0, rep(0.5, 7)), imaginary = c(rep(-0.5, 3), rep(0, 3), rep(0.5,3))), 3, 3))
})

# -----------------------------------------------------------------------
# sinetone()

test_that("parameters to sinetone() are correct", {
  expect_error(sinetone())
  expect_error(sinetone('invalid'))
  expect_error(sinetone(-1))
  expect_error(sinetone(1, 'invalid'))
  expect_error(sinetone(1, 0))
  expect_error(sinetone(1, 1, 'invalid'))
  expect_error(sinetone(1, 1, 0))
  expect_error(sinetone(1, 1, 1, 'invalid'))
  expect_error(sinetone(1, 1, 1, 1, 1))
  
})

test_that("sinetone() works correctly", {
  y <- sinetone(0)
  expect_equal(length(y), 8000)
  expect_equal(y, rep(0, 8000))
  y <-sinetone (18e6, 150e6, 19550/150e6, 1)
  expect_equal(length(y), 19550)
})

# -----------------------------------------------------------------------
# sinewave()

test_that("parameters to sinewave() are correct", {
  expect_error(sinewave())
  expect_error(sinewave(1, 'invalid'))
  expect_error(sinewave(1, 1, 'invalid'))
  expect_error(sinewave(1, 2, 3, 4))
})

test_that("sinetone() works correctly", {
  expect_equal(sinewave(1), 0)
  expect_equal(sinewave(1, 4, 1), 1)
  expect_equal(sinewave(1, 12, 1), 1 / 2, tolerance = tol)
  expect_equal(sinewave(1, 12, 2), sqrt(3) / 2, tolerance = tol)
  expect_equal(sinewave(1, 20, 1), (sqrt(5) - 1) / 4, tolerance = tol)
  expect_equal(sinewave(1), sinewave(1, 1, 0), tolerance = tol)
  expect_equal(sinewave(3, 4), sinewave(3, 4, 0), tolerance = tol)
})

