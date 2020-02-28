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

