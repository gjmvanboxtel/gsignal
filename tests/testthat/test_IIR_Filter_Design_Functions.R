# gsignal IIR filter design functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# cheb()

test_that("parameters to cheb() are correct", {
  expect_error(cheb())
  expect_error(cheb(0.5))
  expect_error(cheb(-1L))
  expect_error(cheb(array(1L, c(1, 4))))
})

test_that("cheb() tests are correct", {
  expect_that(cheb(1, 1), equals(1))
  expect_that(cheb(2, 1), equals(1))
  expect_that(cheb(5, 2), equals(362))
  expect_that(cheb(5, c(2,3)), equals(c(362, 3363)))
})
  
# -----------------------------------------------------------------------
# besselap()

test_that("parameters to besselap() are correct", {
  expect_error(besselap())
  expect_error(besselap(0.5))
  expect_error(besselap(-1L))
  expect_error(besselap(array(1L, c(1, 4))))
  expect_error(besselap(1, 2))
})

test_that("besselap() tests are correct", {
  expect_equal(besselap(1)$z, complex(0))
  expect_equal(besselap(1)$p, -1)
  expect_equal(besselap(2)$p, c(-0.8660254+0.5i, -0.8660254-0.5i))
  expect_equal(besselap(3)$p, c(-0.7456404+0.7113666i, -0.7456404-0.7113666i, -0.9416000+0.0000000i),
               tolerance = 1e-7)
})

# -----------------------------------------------------------------------
# besself()

test_that("parameters to besself() are correct", {
  expect_error(besself())
  expect_error(besself(1))
  expect_error(besself(1, 2, 3, 4))
  expect_error(besself(0.5, 0.2))
  expect_error(besself(3, -1))
  expect_error(besself(3, 2, "invalid"))
})

test_that("besself() tests are correct", {
  zpg <- besself(1, 1, 'low')
  expect_equal(zpg$z, complex(0))
  expect_equal(zpg$p, -1L)
  expect_equal(zpg$g, 1L)

  zpg <- besself(1, 1, 'high')
  expect_equal(zpg$z, 0L)
  expect_equal(zpg$p, -1L)
  expect_equal(zpg$g, 1L)

  zpg <- besself(1, c(1, 2), 'stop')
  expect_equal(zpg$z, c(0-1.414214i, 0+1.414214i), tolerance = 1e-6)
  expect_equal(zpg$p, c(-0.5+1.322876i, -0.5-1.322876i), tolerance = 1e-6)
  expect_equal(zpg$g, 1L)

  zpg <- besself(1, c(1, 2), 'pass')
  expect_equal(zpg$z, 0L)
  expect_equal(zpg$p, c(-0.5+1.322876i, -0.5-1.322876i), tolerance = 1e-6)
  expect_equal(zpg$g, 1L)

  zpg <- besself(2, 1, 'low')
  expect_equal(zpg$z, complex(0))
  expect_equal(zpg$p, c(-0.8660254+0.5i, -0.8660254-0.5i))
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, 1, 'high')
  expect_equal(zpg$z, c(0L, 0L))
  expect_equal(zpg$p, c(-0.8660254-0.5i, -0.8660254+0.5i))
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, c(1, 2), 'stop')
  expect_equal(zpg$z, c(0-1.414214i, 0+1.414214i, 0-1.414214i, 0+1.414214i), tolerance = 1e-6)
  expect_equal(zpg$p, c(-0.354087+1.121579i, -0.354087-1.121579i, -0.511939-1.621579i, -0.511939+1.621579i), tolerance = 1e-6)
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, c(1, 2), 'pass')
  expect_equal(zpg$z, c(0L, 0L))
  expect_equal(zpg$p, c( -0.354087-1.121579i, -0.354087+1.121579i, -0.511939+1.621579i, -0.511939-1.621579i), tolerance = 1e-6)
  expect_equal(zpg$g, 1L)
  
})

# -----------------------------------------------------------------------
# bilinear()

test_that("parameters to bilinear() are correct", {
  expect_error(bilinear())
  expect_error(bilinear(1))
  expect_error(bilinear(1, 2))
  expect_error(bilinear(Zpg(c(1,1,1), 1, 1)))
})

test_that("bilinear() tests are correct", {
  
  res <- bilinear(1, 1, 1, 1)
  expect_equal(res$z, 3)
  expect_equal(res$p, 3)
  expect_equal(res$g, 1)

  res <- bilinear(1, 2, 1, 1)
  expect_equal(res$z, 3)
  expect_equal(res$p, Inf)
  expect_equal(res$g, Inf)
  
  res <- bilinear(1, 3, 1, 1)
  expect_equal(res$z, 3)
  expect_equal(res$p, -5)
  expect_equal(res$g, -1)
  
})

# -----------------------------------------------------------------------
# sftrans()

test_that("parameters to sftrans() are correct", {
  expect_error(sftrans())
  expect_error(sftrans(1))
  expect_error(sftrans(1, 2))
  expect_error(sftrans(Zpg(c(1,1,1), 1, 1)))
})

test_that("sftrans() tests are correct", {
  
  res <- sftrans(1, 1, 1, 1, TRUE)
  expect_equal(res$z, 1)
  expect_equal(res$p, 1)
  expect_equal(res$g, 1)
  
  res <- sftrans(1, 2, 1, 1, TRUE)
  expect_equal(res$z, 1)
  expect_equal(res$p, 0.5)
  expect_equal(res$g, 0.5)
  
  res <- sftrans(1, 3, 1, 1, TRUE)
  expect_equal(res$z, 1)
  expect_equal(res$p, 1/3)
  expect_equal(res$g, 1/3)

  res <- sftrans(1, 3, 1, 1, FALSE)
  expect_equal(res$z, 1)
  expect_equal(res$p, 3)
  expect_equal(res$g, 1)
  
})
