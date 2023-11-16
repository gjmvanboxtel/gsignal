# gsignal Window functions
library(gsignal)
library(testthat)

tol <- 1e-6

# -----------------------------------------------------------------------
# bartlett()

test_that("parameters to bartlett() are correct", {
  expect_error(bartlett())
  expect_error(bartlett(0.5))
  expect_error(bartlett(-1L))
  expect_error(bartlett(array(1L, c(1, 4))))
})

test_that("bartlett() tests are correct", {
  expect_equal(bartlett(1), 1L)
  expect_equal(bartlett(2), c(0, 0))
  expect_equal(rev(bartlett(15)), bartlett(15))
  expect_equal(rev(bartlett(16)), bartlett(16))
  N <- 9
  A <- bartlett(N)
  expect_equal(A[ceiling(N / 2)], 1L)
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
  expect_equal(hamming(1), 1L)
  expect_equal(hamming(2), 0.54 - 0.46 * rep(1L, 2), tolerance = tol)
  expect_equal(rev(hamming(15)), hamming(15))
  expect_equal(rev(hamming(16)), hamming(16))
  
  N <- 15
  A <- hamming(N)
  expect_equal(A[ceiling(N / 2)], 1L)
  
  expect_equal(hamming(15), hamming(15, 'symmetric'), tolerance = tol)
  expect_equal(hamming(16)[1:15], hamming(15, 'periodic'), tolerance = tol)
  
  N <- 16
  A <- hamming(N, 'periodic')
  expect_equal(A[N / 2 + 1], 1L)
  
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
  expect_equal(hann(1), 1L)
  expect_equal(hann(2), 0.5 - 0.5 * rep(1L, 2))
  expect_equal(rev(hann(15)), hann(15))
  expect_equal(rev(hann(16)), hann(16))
  
  N <- 15
  A <- hann(N)
  expect_equal(A[ceiling(N / 2)], 1L)
  
  expect_equal(hann(15), hann(15, 'symmetric'), tolerance = tol)
  expect_equal(hann(16)[1:15], hann(15, 'periodic'), tolerance = tol)
  
  N <- 16
  A <- hann(N, 'periodic')
  expect_equal(A[N / 2 + 1], 1L)
  
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
  expect_equal(triang(1), 1L)
  expect_equal(triang(2), c(1, 1) / 2)
  expect_equal(triang(3), c(1, 2, 1) / 2)
  expect_equal(triang(4), c(1, 3, 3, 1) / 4)
  x <- bartlett(5)
  expect_equal(triang(3), x[2:4])
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
  expect_equal(blackman(1), 1L)
  expect_equal(blackman(2), c(0, 0))
  expect_equal(rev(blackman(15)), blackman(15))
  expect_equal(rev(blackman(16)), blackman(16))
  
  N <- 9
  A <- blackman(N)
  expect_equal(A[ceiling(N / 2)], 1L)
  
  expect_equal(blackman(15), blackman(15, 'symmetric'), tolerance = tol)
  expect_equal(blackman(16)[1:15], blackman(15, 'periodic'), tolerance = tol)
  
  N <- 16
  A <- blackman(N, 'periodic')
  expect_equal(A[N / 2 + 1], 1L)
  
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
  expect_equal(barthannwin(1), 1L)
  expect_equal(barthannwin(2), c(0, 0))
  expect_equal(rev(barthannwin(15)), barthannwin(15))
  expect_equal(rev(barthannwin(16)), barthannwin(16))
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
  expect_equal(blackmanharris(1), 1L)
  expect_equal(blackmanharris(2), c(6e-5, 6e-5))
  expect_equal(rev(blackmanharris(15)), blackmanharris(15))
  expect_equal(rev(blackmanharris(16)), blackmanharris(16))
  expect_equal(blackmanharris(15), blackmanharris(15, 'symmetric'), tolerance = tol)
  expect_equal(blackmanharris(16)[1:15], blackmanharris(15, 'periodic'), tolerance = tol)
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
  expect_equal(blackmannuttall(1), 1L)
  expect_equal(blackmannuttall(2), c(0.0003628, 0.0003628), , tolerance = tol)
  expect_equal(rev(blackmannuttall(15)), blackmannuttall(15))
  expect_equal(rev(blackmannuttall(16)), blackmannuttall(16))
  expect_equal(blackmannuttall(15), blackmannuttall(15, 'symmetric'), tolerance = tol)
  expect_equal(blackmannuttall(16)[1:15], blackmannuttall(15, 'periodic'), tolerance = tol)
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
  expect_equal(barthannwin(1), 1L)
  expect_equal(barthannwin(2), c(0, 0))
  expect_equal(rev(barthannwin(15)), barthannwin(15))
  expect_equal(rev(barthannwin(16)), barthannwin(16))
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
  expect_equal(bohmanwin(1), 1L)
  expect_equal(bohmanwin(2), rep(0, 2))
  expect_equal(rev(bohmanwin(15)), bohmanwin(15))
  expect_equal(rev(bohmanwin(16)), bohmanwin(16))
  expect_equal(bohmanwin(15)[1], 0L)
  expect_equal(bohmanwin(15)[15], 0L)
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
  expect_equal(boxcar(1), 1L)
  expect_equal(boxcar(2), rep(1L, 2))
  expect_equal(rev(boxcar(100)), rep(1L, 100))
})

# -----------------------------------------------------------------------
# chebwin()

test_that("parameters to chebwin() are correct", {
  expect_error(chebwin())
  expect_error(chabwin(0.5))
  expect_error(chebwin(-1L))
  expect_error(chebwin(array(1L, c(1, 4))))
})

test_that("chebwin() tests are correct", {
  expect_equal(chebwin(1), 1L)
  expect_equal(chebwin(2), rep(1L, 2))
  expect_equal(rev(chebwin(15)), chebwin(15))
  expect_equal(rev(chebwin(16)), chebwin(16))
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
  expect_equal(flattopwin(1), 1L)
  expect_equal(flattopwin(2), 0.0042 / 4.6402 * rep(1L, 2), tolerance = tol)
  expect_equal(rev(flattopwin(15)), flattopwin(15))
  expect_equal(rev(flattopwin(16)), flattopwin(16))
  expect_equal(flattopwin(15), flattopwin(15, 'symmetric'), tolerance = tol)
  expect_equal(flattopwin(16)[1:15], flattopwin(15, 'periodic'), tolerance = tol)
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
  expect_equal(gaussian(1), 1L)
  expect_equal(rev(gaussian(15)), gaussian(15))
  expect_equal(rev(gaussian(16)), gaussian(16))
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
  expect_equal(gausswin(1), 1)
  expect_equal(gausswin(2), c(exp(-3.125), exp(-3.125)), tolerance = tol)
  expect_equal(gausswin(3), c(exp(-3.125), 1, exp(-3.125)), tolerance = tol)
  expect_equal(rev(gausswin(15)), gausswin(15))
  expect_equal(rev(gausswin(16)), gausswin(16))
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
  expect_equal(kaiser(1), 1)
  expect_equal(round(kaiser(2), 4), rep(0.9403, 2))
  expect_equal(rev(kaiser(15)), kaiser(15))
  expect_equal(rev(kaiser(16)), kaiser(16))
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
  expect_equal(nuttallwin(1), 1L)
  expect_equal(nuttallwin(2), c(0, 0))
  expect_equal(rev(nuttallwin(15)), nuttallwin(15))
  expect_equal(rev(nuttallwin(16)), nuttallwin(16))
  expect_equal(nuttallwin(15), nuttallwin(15, 'symmetric'), tolerance = tol)
  expect_equal(nuttallwin(16)[1:15], nuttallwin(15, 'periodic'), tolerance = tol)
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
  expect_equal(parzenwin(1), 1L)
  expect_equal(parzenwin(2), 0.25 * rep(1, 2))
  expect_equal(rev(parzenwin(15)), parzenwin(15))
  expect_equal(rev(parzenwin(16)), parzenwin(16))
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
  expect_equal(rectwin(1), 1L)
  expect_equal(rectwin(2), rep(1L, 2))
  expect_equal(rev(rectwin(100)), rep(1L, 100))
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
  expect_equal(tukeywin(1, 0), 1L)
  expect_equal(tukeywin(1, 1), 1L)
  expect_equal(tukeywin(2, 0), rep(1L, 2))
  expect_equal(tukeywin(2, 1), rep(0L, 2))
  expect_equal(tukeywin(3, 0), rep(1L, 3))
  expect_equal(tukeywin(3, 1), c(0, 1, 0))
  expect_equal(tukeywin(4, 0), rep(1L, 4))
  expect_equal(tukeywin(4, 1), c(0, 0.75, 0.75, 0))
  expect_equal(tukeywin(5, 0), rep(1L, 5))
  expect_equal(tukeywin(5, 1), c(0, 0.5, 1, 0.5, 0))
  expect_equal(tukeywin(16, 0), rectwin(16))
  expect_equal(tukeywin(16, 1), hann(16), tolerance = tol)
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
  expect_equal(welchwin(2, 'periodic'), c(0,1))
  expect_equal(welchwin(3, 'symmetric'), c(0, 1, 0))
  expect_equal(rev(welchwin(15)), welchwin(15))
  expect_equal(rev(welchwin(16)), welchwin(16))
  expect_equal(welchwin(15), welchwin(15, 'symmetric'), tolerance = tol)
  expect_equal(welchwin(16)[1:15], welchwin(15, 'periodic'), tolerance = tol)
})
