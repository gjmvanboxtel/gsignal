# gsignal Standard Functions
library(gsignal)
library(testthat)

# -----------------------------------------------------------------------
# cconv()

test_that("parameters to cconv() are correct", {
  expect_error(cconv())
  expect_error(cconv('invalid'))
  expect_error(cconv(1, 1, c(1,1)))
  expect_error(cconv(1, 1, -1))
  expect_error(cconv(1, 1, 'invalid'))
})

test_that("cconv() tests are correct", {
  x <- 1:5
  expect_equal(cconv(x, 1), 1:5)
  expect_equal(cconv(x, c(1, 1)), c(1, 3, 5, 7, 9, 5))
  expect_equal(cconv(x, c(1, 1), 3), c(8, 12, 10))
  
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4)), c(2, 5, 10, 16, 12, 11, 4))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 4), c(14, 16, 14, 16))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 3), c(22, 17, 21))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 2), c(28, 32))
  expect_equal(cconv(c(2, 1, 2, 1), c(1, 2, 3, 4), 1), 60)
  
  expect_equal(cconv(x*1i, 1), c(0+1i, 0+2i, 0+3i, 0+4i, 0+5i))
})

# -----------------------------------------------------------------------
# convmtx()

test_that("parameters to convmtx() are correct", {
  expect_error(convmtx())
  expect_error(convmtx(1, 'invalid'))
  expect_error(convmtx(1, -1))
})

test_that("convmtx() tests are correct", {
  expect_equal(convmtx(c(3, 4, 5), 3), matrix(c(3,4,5,0,0,0,3,4,5,0,0,0,3,4,5), 5, 3))
})

# -----------------------------------------------------------------------
# wconv()

test_that("parameters to wconv() are correct", {
  expect_error(wconv())
  expect_error(wconv('invalid'))
  expect_error(wconv(1, -1))
  expect_error(wconv('1', -1))
  expect_error(wconv('1', 1:5, 1:2, 'invalid'))
})

test_that("wconv() tests are correct", {
  a <- matrix(1:16, 4, 4)
  b <- matrix(1:9, 3,3)
  expect_equal(wconv('2', a, b), 
               matrix(c(1,4,10,16,17,12,9,29,62,83,75,48,36,99,192,237,198,120,84,207,
                        372,417,330,192,115,263,446,485,365,204,91,202,334,358,263,144), 6,6))
  expect_equal(wconv('1', a, b, 'same'), c(35,56,84,120,165,210,255,300,345,390,435,480,508,518,509,480))
  expect_equal(wconv('r', a, b), matrix(c(1,2,3,4,7,10,13,16,22,28,34,40,50,60,70,80,78,92,106,120,
                                          106,124,142,160,134,156,178,200,162,188,214,240,190,220,250,
                                          280,208,232,256,280,185,202,219,236,117,126,135,144), 4, 12))
  expect_equal(wconv('r', a, c(0,1), 'same'), matrix(1:16, 4, 4))
  expect_equal(wconv('c', a, c(0,1), 'valid'), matrix(c(1:3,5:7,9:11,13:15), 3, 4))
})
    
# -----------------------------------------------------------------------
# xcorr()

test_that("parameters to xcorr() are correct", {
  expect_error(xcorr())
  expect_error(xcorr('invalid'))
  expect_error(xcorr(array(1:12, dim = c(2, 2, 3))))
  expect_error(xcorr(1, 'invalid'))
  expect_error(xcorr(1, array(1:12, dim = c(2, 2, 3))))
  expect_error(xcorr(1, -1, maglag = -1))
  expect_error(xcorr(matrix(1:9, 3, 3), 1))
  expect_error(xcorr(1, -1, scale = 'invalid'))
  expect_error(xcorr(1:10, 1:10, 2, 'none', 'extra'))
})

test_that("xcorr() tests are correct", {
  rl <- xcorr(1, -1)
  expect_equal(rl$R, -1)
  expect_equal(rl$lags, 0)
  
  rl <- xcorr(c(1, 2))
  expect_equal(rl$R, c(2, 5, 2))
  expect_equal(rl$lags, c(-1, 0, 1))
  
  rl <- xcorr(1:10, 1:10, 2, 'none')
  expect_equal(rl$R, c(276, 330, 385, 330, 276))
  expect_equal(rl$lags, -2:2)

  rl <- xcorr(1:10, 1:10, 2, 'biased')
  expect_equal(rl$R, c(27.6, 33.0, 38.5, 33.0, 27.6))
  expect_equal(rl$lags, -2:2)

  rl <- xcorr(1:10, 1:10, 2, 'unbiased')
  expect_equal(rl$R, c(34.5, 36.666667, 38.5, 36.666667, 34.5))
  expect_equal(rl$lags, -2:2)

  rl <- xcorr(1:10, 1:10, 2, 'coeff')
  expect_equal(rl$R, c(0.7168831, 0.8571429, 1.0000000, 0.8571429, 0.7168831), tolerance = 1e-7)
  expect_equal(rl$lags, -2:2)
  
})

# -----------------------------------------------------------------------
# xcorr2()

test_that("parameters to xcorr2() are correct", {
  expect_error(xcorr2())
  expect_error(xcorr2('invalid'))
  expect_error(xcorr2(array(1:12, dim = c(2, 2, 3))))
  expect_error(xcorr2(1, 'invalid'))
  expect_error(xcorr2(1, array(1:12, dim = c(2, 2, 3))))
  expect_error(xcorr2(1, -1))
  expect_error(xcorr2(matrix(1:9, 3, 3), 1))
  expect_error(xcorr2(matrix(1:9, 3, 3), scale = 'invalid'))
  expect_error(xcorr2(matrix(1:9, 3, 3), matrix(1:9, 3, 3), 'none', 'extra'))
})

test_that("xcorr2() tests are correct", {

  a <- pracma::magic(3)
  b <- matrix(c(6, 13, 10, 18), 2, 2)
  R <- matrix(c(144,122,121, 78,
                134,187,257,127,
                102,282,253, 68,
                 40,114, 74, 12), 4, 4, byrow = TRUE)
  expect_equal(xcorr2(a, b, 'none'), R)
  expect_equal(xcorr2(a, b, 'biased'), R / 4)
  expect_equal(xcorr2(a, b, 'unbiased'), R / matrix(c(1,2,2,1,2,4,4,2,2,4,4,2,1,2,2,1), 4, 4))
  Rc <- matrix(c(0.71771, 0.60336, 0.79316, 0.51834,
                 0.62534, 0.74937, 0.97263, 0.54925,
                 0.81340, 0.98240, 0.80001, 0.37243,
                 0.39873, 0.46152, 0.32003, 0.23924), 4, 4, byrow = TRUE)
  expect_equal(xcorr2(a, b, 'coeff'), Rc, tolerance = 1e-5)

  row_shift <- 18
  col_shift <- 20
  a <- matrix(runif(900, 1, 255), 30, 30)
  b <- a[(row_shift - 10):row_shift, (col_shift - 7):col_shift]
  R <- xcorr2(a, b, "coeff")
  expect_equal(as.vector(which(R == max(R), arr.ind = TRUE)), c(row_shift, col_shift))
})

# -----------------------------------------------------------------------
# xcov()

test_that("parameters to xcov() are correct", {
  expect_error(xcov())
  expect_error(xcov('invalid'))
  expect_error(xcov(array(1:12, dim = c(2, 2, 3))))
  expect_error(xcov(1, 'invalid'))
  expect_error(xcov(1, array(1:12, dim = c(2, 2, 3))))
  expect_error(xcov(1, -1, maglag = -1))
  expect_error(xcov(matrix(1:9, 3, 3), 1))
  expect_error(xcov(1, -1, scale = 'invalid'))
  expect_error(xcov(1:10, 1:10, 2, 'none', 'extra'))
})

test_that("xcov() tests are correct", {
  cl <- xcov(1, -1)
  expect_equal(cl$C, 0)
  expect_equal(cl$lags, 0)
  
  cl <- xcov(c(1, 2))
  expect_equal(cl$C, c(-0.25, 0.50, -0.25))
  expect_equal(cl$lags, c(-1, 0, 1))
  
  cl <- xcov(1:10, 1:10, 2, 'none')
  expect_equal(cl$C, c(34, 57.75, 82.50, 57.75, 34))
  expect_equal(cl$lags, -2:2)
  
  cl <- xcov(1:10, 1:10, 2, 'biased')
  expect_equal(cl$C, c(3.4, 5.775, 8.25, 5.775, 3.4))
  expect_equal(cl$lags, -2:2)
  
  cl <- xcov(1:10, 1:10, 2, 'unbiased')
  expect_equal(cl$C, c(4.25, 6.4167, 8.25, 6.4167, 4.25), tolerance = 1e-5)
  expect_equal(cl$lags, -2:2)
  
  cl <- xcov(1:10, 1:10, 2, 'coeff')
  expect_equal(cl$C, c(0.4121212, 0.7, 1, 0.7, 0.4121212), tolerance = 1e-7)
  expect_equal(cl$lags, -2:2)
  
})
