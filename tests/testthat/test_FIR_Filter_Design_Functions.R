# gsignal filtering functions
library(gsignal)
library(testthat)

tol <- 1e-6

# -----------------------------------------------------------------------
# sgolay()

test_that("parameters to sgolay() are correct", {
  expect_error(sgolay())
  expect_error(sgolay(-1))
  expect_error(sgolay(0.5))
  expect_error(sgolay(1))
  expect_error(sgolay(1, -1))
  expect_error(sgolay(1, 0.5))
  expect_error(sgolay(1, 2))
  expect_error(sgolay(1, 2, 3, 4, 5))
})

test_that("sgolay() tests are correct", {
  N <- 2^12
  t <- seq(0, N-1) / N
  dt <- t[2] - t[1]
  w <- 2 * pi * 50
  offset <- 0.5 # 50 Hz carrier
  # exponential modulation and its derivatives
  d <- 1 + exp(-3 * (t - offset))
  dd <- -3 * exp(-3 * (t - offset))
  d2d <- 9 * exp(-3 * (t - offset))
  d3d <- -27 * exp(-3 * (t -offset))
  # modulated carrier and its derivatives
  x <- d * sin(w * t)
  dx <- dd * sin(w * t) + w * d * cos(w * t)
  d2x <- (d2d - w^2 * d) * sin(w * t) + 2 * w * dd * cos(w * t)
  d3x <- (d3d - 3 * w^2 * dd) * sin(w * t) + (3 * w * d2d - w^3 * d) * cos(w * t)

  y <- sgolayfilt(x, sgolay(8, 41, 0, dt))
  expect_equal(norm(y - x, '2') / norm(x, '2'), 0, tolerance = 5e-6)

  y <- sgolayfilt(x, sgolay(8, 41, 1, dt))
  expect_equal(norm(y - dx, '2') / norm(dx, '2'), 0, tolerance = 5e-6)

  y <- sgolayfilt(x,sgolay(8, 41, 2, dt))
  expect_equal(norm(y - d2x, '2') / norm(d2x, '2'), 0, tolerance = 1e-5)

  y <- sgolayfilt(x, sgolay(8, 41, 3, dt))
  expect_equal(norm(y - d3x, '2') / norm(d3x, '2'), 0, tolerance = 1e-4)
})

# -----------------------------------------------------------------------
# fir2()

test_that("parameters to fir2() are correct", {
  expect_error(fir2())
  expect_error(fir2(-1))
  expect_error(fir2(0.5))
  expect_error(fir2(1))
  expect_error(fir2(1, -1))
  expect_error(fir2(1, 0.5))
  expect_error(fir2(1, 2))
  expect_error(fir2(1, c(0,1), c(1,0), 4, 5, hamming(-1)))
  expect_error(fir2(1, c(0,1), c(1,0), 4, 5, hamming(1), 7))
})

test_that("fir2() tests are correct", {

  # Test that the grid size is rounded up to the next power of 2
  f <- c(0, 0.6, 0.6, 1); m <- c(1, 1, 0, 0)
  b9  <- fir2 (30, f, m, 9)
  b16 <- fir2 (30, f, m, 16)
  b17 <- fir2 (30, f, m, 17)
  b32 <- fir2 (30, f, m, 32)
  expect_equal(b9,  b16, tolerance = tol)
  expect_equal(b17, b32, tolerance = tol)
  expect_false(isTRUE(all.equal(b16, b17, tolerance = tol)))
  
  # Test expected magnitudes of passbands, stopbands, and cutoff frequencies
  f <- c(0, 0.7, 0.7, 1); m <- c(0, 0, 1, 1)
  b <- fir2 (50, f, m)
  h <- abs(freqz (b, c(0, 0.7, 1), fs = 2)$h)
  expect_lte(h[1], 3e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1, tolerance = 2e-3)
  
  f <- c(0, 0.25, 0.25, 0.75, 0.75, 1); m <- c(0, 0, 1, 1, 0, 0)
  b <- fir2 (50, f, m)
  h <- abs (freqz (b, c(0, 0.25, 0.5, 0.75, 1), fs = 2)$h)
  expect_lte(h[1], 3e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1, tolerance = 2e-3)
  expect_lte(h[4], 1 / sqrt (2))
  expect_lte(h[5], 3e-3)
  
  f <- c(0, 0.45, 0.45, 0.55, 0.55, 1); m <- c(1, 1, 0, 0, 1, 1)
  b <- fir2 (50, f, m)
  h <- abs (freqz (b, c(0, 0.45, 0.5, 0.55, 1), fs = 2)$h)
  expect_equal(h[1], 1, tolerance = 2e-3)
  expect_lte(h[2], 1 / sqrt(2))
  expect_equal(h[3], 1e-1, tolerance = 2e-2)
  expect_lte(h[4], 1 / sqrt (2))
  expect_equal(h[5], 1, tolerance = 2e-3)
  
})

# -----------------------------------------------------------------------
# fir1()

test_that("parameters to fir1() are correct", {
  expect_error(fir1())
  expect_error(fir1(-1))
  expect_error(fir1(0.5))
  expect_error(fir1(1))
  expect_error(fir1(1, -1))
  expect_error(fir1(1, 2))
  expect_error(fir1(1, 0.5, 'invalid'))
  expect_error(fir1(1, 0.5, "low", 'invalid(1)'))
  expect_error(fir1(1, 0.5, "low", hamming(2), 'invalid'))
  expect_error(fir1(1, 0.5, "low", hamming(2), 'scale', 6))
})

test_that("fir1() tests are correct", {
  
  b <- fir1(30, 0.3)
  h <- abs (freqz (b, c(0, 0.3, 1), fs = 2)$h)
  expect_equal(h[1], 1, tolerance = 1e-2)
  expect_true(all(h[2:3] <= 1 / sqrt(2)))
  
  b <- fir1(30, 0.7, "high")
  h <- abs (freqz (b, c(0, 0.7, 1), fs = 2)$h)
  expect_equal(h[3], 1, tolerance = 1e-2)
  expect_true(all(h[1:2] <= c(3e-3, 1 / sqrt(2))))
  
  b <- fir1 (30, c(0.3, 0.7), 'pass')
  h <- abs (freqz (b, c(0, 0.3, 0.5, 0.7, 1), fs = 2)$h)
  expect_equal(h[3], 1, tolerance = 1e-3)
  expect_true(all(h[-3] <= c(3e-3, 1 / sqrt(2), 1 / sqrt(2), 3e-3)))
  
  b <- fir1(50, c(0.3, 0.7), "stop")
  h <- abs (freqz (b, c(0, 0.3, 0.5, 0.7, 1), fs = 2)$h)
  expect_equal(h[c(1, 5)], c(1, 1), tolerance = 1e-3)
  expect_true(all(h[2:4] <= c(1 / sqrt(2), 3e-3, 1 / sqrt(2))))
})

# -----------------------------------------------------------------------
# firls()

test_that("parameters to firls() are correct", {
  expect_error(firls())
  expect_error(firls(-1))
  expect_error(firls(0.5))
  expect_error(firls(1))
  expect_error(firls(1, 2))
  expect_error(firls(1, 2, 3))
  expect_error(firls(1, 0.5, c(1, 1)))
  expect_error(firls(1, c(0.1, 0.5), c(1, 0), c(3, 3)))
  expect_error(firls(1, 2, 3, 4, 5))
})

# -----------------------------------------------------------------------
# firpm()

cases <- list(
  type_I = list(n=30L, f=c(0,.3,.4,1), a=c(1,0),
                class="symmetric", sign=1),
  type_II = list(n=31L, f=c(0,.3,.4,1), a=c(1,0),
                 class="symmetric", sign=1),
  type_III = list(n=30L, f=c(.1,.9), a=c(1,1),
                  class="antisymmetric", sign=-1),
  type_IV = list(n=31L, f=c(.1,1), a=c(1,1),
                 class="antisymmetric", sign=-1),
  differentiator_even = list(n=30L, f=c(0,.9), a=c(0,.9),
                             class="differentiator", sign=-1),
  differentiator_odd = list(n=31L, f=c(0,.9), a=c(0,.9),
                            class="differentiator", sign=-1),
  hilbert = list(n=30L, f=c(.1,.9), a=c(1,1),
                 class="hilbert", sign=-1)
)

call_firpm <- function(spec, output="coefficients", control=c(0,0,0)) {
  gsignal::firpm(spec$n, spec$f, spec$a,
                 class=spec$class, control=control, output=output)
}

# Absolute tolerance appropriate for the previously verified test designs.
tol <- 1e-12

test_that("default output, full output, and linear phase agree", {
  for (nm in names(cases)) {
    s <- cases[[nm]]
    b <- call_firpm(s)
    full <- call_firpm(s, "full")
    expect_type(b, "double")
    expect_length(b, s$n+1L)
    expect_true(all(is.finite(b)))
    expect_equal(b, full$b, tolerance=0)
    expect_equal(b, s$sign * rev(b), tolerance=tol)
    expect_true(is.finite(full$minimax))
    expect_true(is.numeric(full$iterations))
    expect_true(is.numeric(full$result))
    expect_true(is.character(full$message))
  }
})

test_that("full diagnostics have consistent shapes and indexing", {
  for (nm in names(cases)) {
    z <- call_firpm(cases[[nm]], "full")
    r <- z$res
    expect_named(r, c("fgrid","des","wt","H","error","iextr","fextr"),
                 ignore.order=TRUE)
    m <- length(r$fgrid)
    expect_gt(m, 0L)
    expect_length(r$des, m)
    expect_length(r$wt, m)
    expect_length(r$H, m)
    expect_length(r$error, m)
    expect_true(is.complex(r$H))
    expect_true(all(is.finite(r$fgrid)))
    expect_true(all(diff(r$fgrid) >= 0))
    expect_true(all(r$fgrid >= 0 & r$fgrid <= 1))
    expect_length(r$fextr, length(r$iextr))
    expect_true(all(r$iextr >= 1 & r$iextr <= m))
    expect_equal(as.numeric(r$fextr), as.numeric(r$fgrid[r$iextr]),
                 tolerance=tol)
  }
})

test_that("amplitude and weight encodings are equivalent", {
  s <- cases$type_I
  a_edge <- c(1,1,0,0)
  w_band <- c(1,10)
  w_edge <- c(1,1,10,10)
  b1 <- gsignal::firpm(s$n,s$f,s$a,w_band,output="coefficients")
  b2 <- gsignal::firpm(s$n,s$f,a_edge,w_edge,output="coefficients")
  expect_equal(b1,b2,tolerance=tol)
})

test_that("control defaults and zero padding agree", {
  s <- cases$type_I
  expect_equal(call_firpm(s,control=0),
               call_firpm(s,control=c(0,0,0)),tolerance=0)
  expect_equal(call_firpm(s,control=c(0,0)),
               call_firpm(s,control=c(0,0,0)),tolerance=0)
})

test_that("invalid inputs are rejected", {
  base <- cases$type_I
  expect_error(gsignal::firpm(0,base$f,base$a), "n")
  expect_error(gsignal::firpm(2.5,base$f,base$a), "n")
  expect_error(gsignal::firpm(base$n,c(0,.4,.3,1),base$a), "f")
  expect_error(gsignal::firpm(base$n,c(-.1,.3,.4,1),base$a), "f")
  expect_error(gsignal::firpm(base$n,base$f,c(1,0,0)), "a")
  expect_error(gsignal::firpm(base$n,base$f,base$a,w=c(1,0)), "w")
  expect_error(gsignal::firpm(base$n,base$f,base$a,control=-1), "control")
  expect_error(gsignal::firpm(base$n,base$f,base$a,control=c(1,2,3,4)), "control")
  expect_error(gsignal::firpm(base$n,base$f,base$a,class="other"))
  expect_error(gsignal::firpm(base$n,base$f,base$a,output="other"))
})

test_that("callback errors are recoverable", {
  
  result <- local({
    
    tryCatch(
      firpm(
        30,
        c(0, 1),
        function(n, f, g, w) {
          stop("bad callback")
        }
      ),
      error = function(e) {
        expect_match(
          conditionMessage(e),
          "R response function raised an error"
        )
        
        "RECOVERED"
      }
    )
    
  })
  
  expect_identical(result, "RECOVERED")
  
})

test_that("invalid callback results are rejected", {
  
  # Incorrect amplitude length
  expect_error(
    firpm(
      30,
      c(0, 1),
      function(n, f, g, w) numeric(0)
    ),
    "amplitude"
  )
  
  # Negative weights
  expect_error(
    firpm(
      30,
      c(0, 1),
      function(n, f, g, w) {
        list(
          desired = rep(1, length(g)),
          weight = rep(-1, length(g))
        )
      }
    ),
    "weight"
  )
  
})



