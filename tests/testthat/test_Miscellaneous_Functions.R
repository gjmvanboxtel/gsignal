# gsignal Miscellaneous Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# ifft() and imvfft()

test_that("parameters to ifft() are correct", {
  expect_error(pad())
  expect_error(pad('invalid'))
  expect_error(pad(1, -2))
  expect_error(pad(1, 2, 3, 4, 5))
  expect_error(pad(1, 2, direction = 'invalid'))
  expect_error(pad(matrix(1:6, 2, 3), 2, MARGIN = 3))
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

