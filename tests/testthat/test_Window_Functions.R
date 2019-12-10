# gsignal Window functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# bartlett()

test_that("parameters to bartlett() are correct", {
  expect_error(bartlett())
  expect_error(bartlett(0.5))
  expect_error(bartlett(-1L))
  expect_error(bartlett(array(1L, c(1, 4))))
})

test_that("bartlett() tests are correct", {
  expect_that(bartlett(1), equals(1))
  expect_that(bartlett(2), equals(c(0, 0)))
  expect_that(rev(bartlett(15)), equals(bartlett(15)))
  expect_that(rev(bartlett(16)), equals(bartlett(16)))
  N <- 9
  A <- bartlett(N)
  expect_that(A[ceiling(N / 2)], equals(1L))
})

# -----------------------------------------------------------------------
# hamming()

test_that("parameters to hamming() are correct", {
  expect_error(hamming())
  expect_error(hamming(0.5))
  expect_error(hamming(-1L))
  expect_error(hamming(array(1L, c(1, 4))))
  expect_error(hamming(1, 'invalid'))
})

test_that("hamming() tests are correct", {
  expect_that(hamming(1), equals(1))
  expect_that(hamming(2), equals(25/46 - 21/46 * rep(1L, 2)))
  expect_that(rev(hamming(15)), equals(hamming(15)))
  expect_that(rev(hamming(16)), equals(hamming(16)))
  
  N <- 15
  A <- hamming(N)
  expect_that(A[ceiling(N / 2)], equals(1L))
  
  expect_that(hamming(15), equals(hamming(15, 'symmetric')))
  expect_that(hamming(16)[1:15], equals(hamming(15, 'periodic')))
  
  N <- 16
  A <- hamming(N, 'periodic')
  expect_that(A[N / 2 + 1], equals(1L))
  
})

# -----------------------------------------------------------------------
# hann()

test_that("parameters to hann() are correct", {
  expect_error(hann())
  expect_error(hann(0.5))
  expect_error(hann(-1L))
  expect_error(hann(array(1L, c(1, 4))))
  expect_error(hann(1, 'invalid'))
})

test_that("hann() tests are correct", {
  expect_that(hann(1), equals(1))
  expect_that(hann(2), equals(0.5 - 0.5 * rep(1L, 2)))
  expect_that(rev(hann(15)), equals(hann(15)))
  expect_that(rev(hann(16)), equals(hann(16)))
  
  N <- 15
  A <- hann(N)
  expect_that(A[ceiling(N / 2)], equals(1L))
  
  expect_that(hann(15), equals(hann(15, 'symmetric')))
  expect_that(hann(16)[1:15], equals(hann(15, 'periodic')))
  
  N <- 16
  A <- hann(N, 'periodic')
  expect_that(A[N / 2 + 1], equals(1L))
  
})

# -----------------------------------------------------------------------
# triang()

test_that("parameters to triang() are correct", {
  expect_error(triang())
  expect_error(triang(0.5))
  expect_error(triang(-1L))
  expect_error(triang(array(1L, c(1, 4))))
})

test_that("triang() tests are correct", {
  expect_that(triang(1), equals(1))
  expect_that(triang(2), equals(c(1, 1) / 2))
  expect_that(triang(3), equals(c(1, 2, 1) / 2))
  expect_that(triang(4), equals(c(1, 3, 3, 1) / 4))
  x <- bartlett(5)
  expect_that(triang(3), equals(x[2:4]))
})

# -----------------------------------------------------------------------
# blackman()

test_that("parameters to blackman() are correct", {
  expect_error(blackman())
  expect_error(blackman(0.5))
  expect_error(blackman(-1L))
  expect_error(blackman(array(1L, c(1, 4))))
  expect_error(blackman(1, 'invalid'))
})

test_that("blackman() tests are correct", {
  expect_that(blackman(1), equals(1))
  expect_that(blackman(2), equals(c(0, 0)))
  expect_that(rev(blackman(15)), equals(blackman(15)))
  expect_that(rev(blackman(16)), equals(blackman(16)))
  
  N <- 9
  A <- blackman(N)
  expect_that(A[ceiling(N / 2)], equals(1L))
  
  expect_that(blackman(15), equals(blackman(15, 'symmetric')))
  expect_that(blackman(16)[1:15], equals(blackman(15, 'periodic')))
  
  N <- 16
  A <- blackman(N, 'periodic')
  expect_that(A[N / 2 + 1], equals(1L))
  
})
