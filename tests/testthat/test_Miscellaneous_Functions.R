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

# -----------------------------------------------------------------------
# filter() 

test_that("parameters to filter() are correct", {
  expect_error(filter())
  expect_error(filter(1, 2))
  expect_warning(filter(1, 2, 'invalid'))
  expect_error(filter(1, 1, 1:10, init.x = 1))
})

test_that("filter() tests are correct", {
  a <- c(1, 1)
  b <- c(1, 1)
  x <- c(1, rep(0L, 9))
  expect_equal(filter(b, 1, x), c(rep(1L, 2), rep(0L, 8)))
  filt <- Ma(b)
  expect_equal(filter(filt, x), c(rep(1L, 2), rep(0L, 8)))
  expect_equal(filter(1, a, x), rep(c(1L, -1L), 5))
  filt <- Arma(b, a)
  expect_equal(filter(filt, x), c(1L, rep(0L, 9)))

  # complex input  
  r <- sqrt (1/2) * (1 + 1i)
  a <- a * r
  b <- b * r
  expect_equal(suppressWarnings(filter (b, 1, x)), Re(r * c(rep(1L, 2), rep(0L, 8))))
  expect_equal(suppressWarnings(filter (b, a, x)), c(1L, rep(0L, 9)))
  
  # initial conditions
  expect_equal(filter (c(1,1,1), c(1,1), c(1,2), init.x=c(1,1), init.y=1), c(2, 2))
})

# -----------------------------------------------------------------------
# conv()

test_that("parameters to conv() are correct", {
  expect_error(conv())
  expect_error(conv(1))
  expect_error(conv(1, 2, 3, 4))
  expect_error(conv(1, 2, 'invalid'))
})

test_that("conv() tests are correct", {
  x <- rep(1L, 3); b <- 2; c <- 3
  expect_equal(conv(x, x), c(1, 2, 3, 2, 1))
  expect_equal(conv(x, b), rep(2L, 3))
  expect_equal(conv(b, x), rep(2L, 3))
  expect_equal(conv(x, c), rep(3L, 3))
  expect_equal(conv(c, x), rep(3L, 3))
  expect_equal(conv(b, c), 6)

  a <- 1:10; b <- 1:3
  expect_equal(length(conv (a,b)), length(a) + length(b) - 1)
  expect_equal(length(conv (b,a)), length(a) + length(b) - 1)
  expect_equal(conv(a, b, "full"), conv (a,b))
  expect_equal(conv(b, a, "full"), conv (b,a))
  expect_equal(conv(a, b, "same"), c(4, 10, 16, 22, 28, 34, 40, 46, 52, 47))
  expect_equal(conv(b, a, "same"), c(28, 34, 40))
  expect_equal(conv(a, b, "valid"), c(10, 16, 22, 28, 34, 40, 46, 52))
  expect_equal(conv(b, a, "valid"), NULL)
  expect_equal(conv(a, a, "valid"), 220L)
  expect_equal(conv(b, b, "valid"), 10L)
})

# -----------------------------------------------------------------------
# conv2()

test_that("parameters to conv2() are correct", {
  expect_error(conv2())
  expect_error(conv2(1))
  expect_error(conv2(1, 2, 3, 4))
  expect_error(conv2(matrix(1,1), matrix(2,1), 'invalid'))
})

test_that("conv2() tests are correct", {
  a <- matrix(1:16, 4, 4)
  b <- matrix(1:9, 3,3)
  ans <- matrix(c(1, 9, 36, 84, 115, 91,
                  4, 29, 99, 207, 263, 202,
                  10, 62, 192, 372, 446, 334,
                  16, 83, 237, 417, 485, 358,
                  17, 75, 198, 330, 365, 263,
                  12, 48, 120, 192, 204, 144),
                6, 6, byrow = TRUE)
  expect_equal(conv2(a, b), ans)
  expect_equal(conv2(a, b, 'same'), ans[2:5, 2:5])
  expect_equal(conv2(a, b, 'valid'), ans[3:4, 3:4])

  a <- matrix(c(1:5, 1:5), 2, 5, byrow = TRUE)
  b <- matrix(1:2, 1, 2)
  ans <- matrix(rep(c(1,4,7,10,13,10),2),2,6, byrow=T)
  expect_equal(conv2(a, b), ans)
  expect_equal(conv2(a, b, 'same'), ans[1:2, 2:6])
  expect_equal(conv2(a, b, 'valid'), ans[1:2, 2:5])
})
