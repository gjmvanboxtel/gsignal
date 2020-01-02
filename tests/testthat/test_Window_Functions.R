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

# -----------------------------------------------------------------------
# barthannwin()

test_that("parameters to barthannwin() are correct", {
  expect_error(barthannwin())
  expect_error(barthannwin(0.5))
  expect_error(barthannwin(-1L))
  expect_error(barthannwin(array(1L, c(1, 4))))
  expect_error(barthannwin(1, 2))
})

test_that("barthannwin() tests are correct", {
  expect_that(barthannwin(1), equals(1))
  expect_that(barthannwin(2), equals(c(0, 0)))
  expect_that(rev(barthannwin(15)), equals(barthannwin(15)))
  expect_that(rev(barthannwin(16)), equals(barthannwin(16)))
})

# -----------------------------------------------------------------------
# blackmanharris()

test_that("parameters to blackmanharris() are correct", {
  expect_error(blackmanharris())
  expect_error(blackmanharris(0.5))
  expect_error(blackmanharris(-1L))
  expect_error(blackmanharris(array(1L, c(1, 4))))
  expect_error(blackmanharris(1, 'invalid'))
})

test_that("blackmanharris() tests are correct", {
  expect_that(blackmanharris(1), equals(1))
  expect_that(blackmanharris(2), equals(c(6e-5, 6e-5)))
  expect_that(rev(blackmanharris(15)), equals(blackmanharris(15)))
  expect_that(rev(blackmanharris(16)), equals(blackmanharris(16)))
  expect_that(blackmanharris(15), equals(blackmanharris(15, 'symmetric')))
  expect_that(blackmanharris(16)[1:15], equals(blackmanharris(15, 'periodic')))
})

# -----------------------------------------------------------------------
# blackmannuttall()

test_that("parameters to blackmannuttall() are correct", {
  expect_error(blackmannuttall())
  expect_error(blackmannuttall(0.5))
  expect_error(blackmannuttall(-1L))
  expect_error(blackmannuttall(array(1L, c(1, 4))))
  expect_error(blackmannuttall(1, 'invalid'))
})

test_that("blackmannuttall() tests are correct", {
  expect_that(blackmannuttall(1), equals(1))
  expect_that(blackmannuttall(2), equals(c(0.0003628, 0.0003628)))
  expect_that(rev(blackmannuttall(15)), equals(blackmannuttall(15)))
  expect_that(rev(blackmannuttall(16)), equals(blackmannuttall(16)))
  expect_that(blackmannuttall(15), equals(blackmannuttall(15, 'symmetric')))
  expect_that(blackmannuttall(16)[1:15], equals(blackmannuttall(15, 'periodic')))
})

# -----------------------------------------------------------------------
# barthannwin()

test_that("parameters to barthannwin() are correct", {
  expect_error(barthannwin())
  expect_error(barthannwin(0.5))
  expect_error(barthannwin(-1L))
  expect_error(barthannwin(array(1L, c(1, 4))))
  expect_error(barthannwin(1, 2))
})

test_that("barthannwin() tests are correct", {
  expect_that(barthannwin(1), equals(1))
  expect_that(barthannwin(2), equals(c(0, 0)))
  expect_that(rev(barthannwin(15)), equals(barthannwin(15)))
  expect_that(rev(barthannwin(16)), equals(barthannwin(16)))
})

# -----------------------------------------------------------------------
# bohmanwin()

test_that("parameters to bohmanwin() are correct", {
  expect_error(bohmanwin())
  expect_error(bohmanwin(0.5))
  expect_error(bohmanwin(-1L))
  expect_error(bohmanwin(array(1L, c(1, 4))))
  expect_error(bohmanwin(1, 2))
})

test_that("bohmanwin() tests are correct", {
  expect_that(bohmanwin(1), equals(1))
  expect_that(bohmanwin(2), equals(rep(0, 2)))
  expect_that(rev(bohmanwin(15)), equals(bohmanwin(15)))
  expect_that(rev(bohmanwin(16)), equals(bohmanwin(16)))
  expect_that(bohmanwin(15)[1], equals(0L))
  expect_that(bohmanwin(15)[15], equals(0L))
})

# -----------------------------------------------------------------------
# boxcar()

test_that("parameters to boxcar() are correct", {
  expect_error(boxcar())
  expect_error(boxcar(0.5))
  expect_error(boxcar(-1L))
  expect_error(boxcar(array(1L, c(1, 4))))
  expect_error(boxcar(1, 2))
})

test_that("boxcar() tests are correct", {
  expect_that(boxcar(1), equals(1L))
  expect_that(boxcar(2), equals(rep(1L, 2)))
  expect_that(rev(boxcar(100)), equals(rep(1L, 100)))
})

# -----------------------------------------------------------------------
# chebwin()

test_that("parameters to chebwin() are correct", {
  expect_error(chebwin())
  expect_error(chabwin(0.5))
  expect_error(chebwin(-1L))
  expect_error(chebwin(array(1L, c(1, 4))))
})

test_that("boxcar() tests are correct", {
  expect_that(chebwin(1), equals(1L))
  expect_that(chebwin(2), equals(rep(1L, 2)))
  expect_that(rev(chebwin(15)), equals(chebwin(15)))
  expect_that(rev(chebwin(16)), equals(chebwin(16)))
})

# -----------------------------------------------------------------------
# flattopwin()

test_that("parameters to flattopwin() are correct", {
  expect_error(flattopwin())
  expect_error(flattopwin(0.5))
  expect_error(flattopwin(-1L))
  expect_error(flattopwin(array(1L, c(1, 4))))
  expect_error(flattopwin(1, 'invalid'))
})

test_that("flattopwin() tests are correct", {
  expect_that(flattopwin(1), equals(1))
  expect_that(flattopwin(2), equals(0.0042 / 4.6402 * rep(1L, 2)))
  expect_that(rev(flattopwin(15)), equals(flattopwin(15)))
  expect_that(rev(flattopwin(16)), equals(flattopwin(16)))
  expect_that(flattopwin(15), equals(flattopwin(15, 'symmetric')))
  expect_that(flattopwin(16)[1:15], equals(flattopwin(15, 'periodic')))
})

# -----------------------------------------------------------------------
# gaussian()

test_that("parameters to gaussian() are correct", {
  expect_error(gaussian())
  expect_error(gaussian(0.5))
  expect_error(gaussian(-1L))
  expect_error(gaussian(array(1L, c(1, 4))))
  expect_error(gaussian(1, 2, 3))
})

test_that("gaussian() tests are correct", {
  expect_that(gaussian(1), equals(1))
  expect_that(rev(gaussian(15)), equals(gaussian(15)))
  expect_that(rev(gaussian(16)), equals(gaussian(16)))
})

# -----------------------------------------------------------------------
# gausswin()

test_that("parameters to gausswin() are correct", {
  expect_error(gausswin())
  expect_error(gausswin(0.5))
  expect_error(gausswin(-1L))
  expect_error(gausswin(array(1L, c(1, 4))))
  expect_error(gausswin(1, 2, 3))
})

test_that("gausswin() tests are correct", {
  expect_that(gausswin(1), equals(1))
  expect_that(gausswin(2), equals(c(exp(-3.125), exp(-3.125))))
  expect_that(gausswin(3), equals(c(exp(-3.125), 1, exp(-3.125))))
  expect_that(rev(gausswin(15)), equals(gausswin(15)))
  expect_that(rev(gausswin(16)), equals(gausswin(16)))
})

# -----------------------------------------------------------------------
# kaiser()

test_that("parameters to kaiser() are correct", {
  expect_error(kaiser())
  expect_error(kaiser(0.5))
  expect_error(kaiser(-1L))
  expect_error(kaiser(array(1L, c(1, 4))))
  expect_error(kaiser(1, 2, 3))
})

test_that("kaiser() tests are correct", {
  expect_that(kaiser(1), equals(1))
  expect_that(round(kaiser(2), 4), equals(rep(0.9403, 2)))
  expect_that(rev(kaiser(15)), equals(kaiser(15)))
  expect_that(rev(kaiser(16)), equals(kaiser(16)))
})

# -----------------------------------------------------------------------
# nuttallwin()

test_that("parameters to nuttallwin() are correct", {
  expect_error(nuttallwin())
  expect_error(nuttallwin(0.5))
  expect_error(nuttallwin(-1L))
  expect_error(nuttallwin(array(1L, c(1, 4))))
  expect_error(nuttallwin(1, 2))
  expect_error(nuttallwin(1, 'invalid'))
})

test_that("nuttallwin() tests are correct", {
  expect_that(nuttallwin(1), equals(1))
  expect_that(nuttallwin(2), equals(c(0, 0)))
  expect_that(rev(nuttallwin(15)), equals(nuttallwin(15)))
  expect_that(rev(nuttallwin(16)), equals(nuttallwin(16)))
  expect_that(nuttallwin(15), equals(nuttallwin(15, 'symmetric')))
  expect_that(nuttallwin(16)[1:15], equals(nuttallwin(15, 'periodic')))
})

# -----------------------------------------------------------------------
# parzenwin()

test_that("parameters to parzenwin() are correct", {
  expect_error(parzenwin())
  expect_error(parzenwin(0.5))
  expect_error(parzenwin(-1L))
  expect_error(parzenwin(array(1L, c(1, 4))))
  expect_error(parzenwin(1, 2))
})

test_that("parzenwin() tests are correct", {
  expect_that(parzenwin(1), equals(1))
  expect_that(parzenwin(2), equals(0.25 * rep(1, 2)))
  expect_that(rev(parzenwin(15)), equals(parzenwin(15)))
  expect_that(rev(parzenwin(16)), equals(parzenwin(16)))
})

# -----------------------------------------------------------------------
# rectwin()

test_that("parameters to rectwin() are correct", {
  expect_error(rectwin())
  expect_error(rectwin(0.5))
  expect_error(rectwin(-1L))
  expect_error(rectwin(array(1L, c(1, 4))))
  expect_error(rectwin(1, 2))
})

test_that("rectwin() tests are correct", {
  expect_that(rectwin(1), equals(1L))
  expect_that(rectwin(2), equals(rep(1L, 2)))
  expect_that(rev(rectwin(100)), equals(rep(1L, 100)))
})

# -----------------------------------------------------------------------
# tukeywin()

test_that("parameters to tukeywin() are correct", {
  expect_error(tukeywin())
  expect_error(tukeywin(0.5))
  expect_error(tukeywin(-1L))
  expect_error(tukeywin(array(1L, c(1, 4))))
  expect_error(tukeywin(1, 2, 3))
})

test_that("tukeywin() tests are correct", {
  expect_that(tukeywin(1, 0), equals(1L))
  expect_that(tukeywin(1, 1), equals(1L))
  expect_that(tukeywin(2, 0), equals(rep(1L, 2)))
  expect_that(tukeywin(2, 1), equals(rep(0L, 2)))
  expect_that(tukeywin(3, 0), equals(rep(1L, 3)))
  expect_that(tukeywin(3, 1), equals(c(0, 1, 0)))
  expect_that(tukeywin(4, 0), equals(rep(1L, 4)))
  expect_that(tukeywin(4, 1), equals(c(0, 0.75, 0.75, 0)))
  expect_that(tukeywin(5, 0), equals(rep(1L, 5)))
  expect_that(tukeywin(5, 1), equals(c(0, 0.5, 1, 0.5, 0)))
  expect_that(tukeywin(16, 0), equals(rectwin(16)))
  expect_that(tukeywin(16, 1), equals(hann(16)))
})

# -----------------------------------------------------------------------
# welchwin()

test_that("parameters to welchwin() are correct", {
  expect_error(welchwin())
  expect_error(welchwin(0.5))
  expect_error(welchwin(1))
  expect_error(welchwin(2, "symmatric"))
  expect_error(welchwin(-1L))
  expect_error(welchwin(array(1L, c(1, 4))))
  expect_error(welchwin(1, 'invalid'))
})

test_that("welchwin() tests are correct", {
  expect_that(welchwin(2, 'periodic'), equals(c(0,1)))
  expect_that(welchwin(3, 'symmetric'), equals(c(0, 1, 0)))
  expect_that(rev(welchwin(15)), equals(welchwin(15)))
  expect_that(rev(welchwin(16)), equals(welchwin(16)))
  expect_that(welchwin(15), equals(welchwin(15, 'symmetric')))
  expect_that(welchwin(16)[1:15], equals(welchwin(15, 'periodic')))
})
