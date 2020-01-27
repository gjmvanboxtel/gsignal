# gsignal Miscellaneous Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# ifft() and imvfft()

test_that("parameters to ifft() are correct", {
  expect_error(ifft())
  expect_error(ifft('invalid'))
  expect_error(ifft(1, -2))
  expect_error(ifft(1, 2, 3, 4, 5))
})

test_that("ifft() tests are correct", {
  expect_equal(ifft(stats::fft(1:10)), 1:10)
  expect_equal(ifft(stats::fft(c(1+5i, 2+3i, 3+2i, 4+6i, 5+2i))), c(1+5i, 2+3i, 3+2i, 4+6i, 5+2i))
  expect_equal(imvfft(stats::mvfft(matrix(1:20, 4, 5))), matrix(1:20, 4, 5))
})

# -----------------------------------------------------------------------
# pad(), prepad(), postpad()

test_that("parameters to pad() are correct", {
  expect_error(pad())
  expect_error(pad('invalid'))
  expect_error(pad(1, -2))
  expect_error(pad(1, 2, 3, 4, 5, 6))
})

test_that("pad() tests are correct", {
  v <- 1:24
  expect_equal(postpad(v, 30), c(1:24, rep(0, 6)))
  expect_equal(postpad(v, 20), 1:20)
  expect_equal(prepad(v, 30), c(rep(0, 6), 1:24))
  expect_equal(prepad(v, 20), 5:24)
  
  m <- matrix(1:24, 4, 6)
  expect_equal(postpad(m, 8, 100), matrix(c(1:4, rep(100, 4), 5:8, rep(100, 4), 9:12, rep(100, 4),
                                            13:16, rep(100, 4), 17:20, rep(100, 4), 21:24, rep(100, 4)),
                                          8, 6, byrow = FALSE))
  expect_equal(postpad(m, 8, 100, MARGIN = 1), matrix(c(1:24, rep(100, 8)), 4, 8))
  expect_equal(prepad(m, 8, 100), matrix(c(rep(100, 4), 1:4, rep(100, 4), 5:8, rep(100, 4), 9:12, rep(100, 4),
                                            13:16, rep(100, 4), 17:20, rep(100, 4), 21:24),
                                          8, 6, byrow = FALSE))
  expect_equal(prepad(m, 8, 100, MARGIN = 1), matrix(c(rep(100, 8), 1:24), 4, 8))
  
  expect_equal(postpad(m, 2), matrix(c(1, 2, 5, 6, 9, 10, 13, 14, 17, 18, 21, 22), 2, 6))
  expect_equal(postpad(m, 2, MARGIN = 1), matrix(1:8, 4, 2))
  expect_equal(prepad(m, 2), matrix(c(3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24), 2, 6))
  expect_equal(prepad(m, 2, MARGIN = 1), matrix(17:24, 4, 2))
})

# -----------------------------------------------------------------------
# poly() 

test_that("parameters to poly() are correct", {
  expect_error(poly())
  expect_error(poly('invalid'))
  expect_error(poly(1, 2))
  expect_error(poly(matrix(1:6, 2, 3)))
})

test_that("poly() tests are correct", {
  expect_equal(poly(0), c(1, 0))
  expect_equal(poly(1), c(1, -1))
  expect_equal(poly(-1), c(1, 1))
  expect_equal(poly(c(1, 2, 3)), c(1, -6, 11, -6))
  expect_equal(poly(matrix(1:4, 2, 2, byrow = TRUE)), c(1, -5, -2))
  expect_equal(poly(c(-1 + 1i)), c(1 + 0i, 1 - 1i))
})

# -----------------------------------------------------------------------
# roots() 

test_that("parameters to roots() are correct", {
  expect_error(roots())
  expect_error(roots('invalid'))
  expect_error(roots(1, 2))
})

test_that("roots() tests are correct", {
  expect_equal(roots(0), NULL)
  expect_equal(roots(0, "eigen"), NULL)
  expect_equal(roots(1), numeric(0))
  expect_equal(roots(1, "eigen"), numeric(0))
  
  p <- c(poly(rep(3, 4)), rep(0, 4))
  r <- sort(roots (p))
  expect_equal(r, c(rep(0, 4), rep(3, 4)))
  
  expect_equal(roots(c(1e-200, -1e200, 1)), 1e-200)
  expect_equal(roots(c(1e-200, -1e200 * 1i, 1)), 1e-200 * 1i)
})

