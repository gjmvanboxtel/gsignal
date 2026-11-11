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
  firpm(spec$n, spec$f, spec$a,
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
  b1 <- firpm(s$n,s$f,s$a,w_band,output="coefficients")
  b2 <- firpm(s$n,s$f,a_edge,w_edge,output="coefficients")
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
  expect_error(firpm(0,base$f,base$a), "n")
  expect_error(firpm(2.5,base$f,base$a), "n")
  expect_error(firpm(base$n,c(0,.4,.3,1),base$a), "f")
  expect_error(firpm(base$n,c(-.1,.3,.4,1),base$a), "f")
  expect_error(firpm(base$n,base$f,c(1,0,0)), "a")
  expect_error(firpm(base$n,base$f,base$a,w=c(1,0)), "w")
  expect_error(firpm(base$n,base$f,base$a,control=-1), "control")
  expect_error(firpm(base$n,base$f,base$a,control=c(1,2,3,4)), "control")
  expect_error(firpm(base$n,base$f,base$a,class="other"))
  expect_error(firpm(base$n,base$f,base$a,output="other"))
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


# -----------------------------------------------------------------------
# remez()

# Return ordinary numeric coefficients (remez returns an Ma object).
remez_coeff <- function(...) as.numeric(remez(...))

# Evaluate the signed, phase-corrected amplitude of an antisymmetric FIR.
# The sign is calibrated against the desired positive response at midband.
odd_amplitude <- function(b, f) {
  k <- seq_along(b) - 1L
  delay <- (length(b) - 1) / 2
  vapply(f, function(fi) {
    omega <- pi * fi
    h <- sum(b * exp(-1i * omega * k))
    Im(h * exp(1i * omega * delay))
  }, numeric(1))
}

expect_remez_shape <- function(b, order, antisymmetric = FALSE) {
  expect_length(b, order + 1L)
  expect_true(all(is.finite(b)))
  target <- if (antisymmetric) -rev(b) else rev(b)
  expect_equal(b, target, tolerance = 1e-10)
  if (antisymmetric && order %% 2L == 0L) {
    expect_equal(b[(length(b) + 1L) %/% 2L], 0,
                           tolerance = 1e-12)
  }
}

test_that("Type I and II low-pass filters are symmetric", {
  b1 <- remez_coeff(30, c(0, .3, .4, 1), c(1, 1, 0, 0), density = 32)
  b2 <- remez_coeff(29, c(0, .3, .4, .95), c(1, 1, 0, 0), density = 32)
  expect_remez_shape(b1, 30)
  expect_remez_shape(b2, 29)
})


test_that("Type III and IV Hilbert transformers are antisymmetric", {
  b3 <- remez_coeff(30, c(.1, .9), c(1, 1),
                    ftype = "hilbert", density = 32)
  b4 <- remez_coeff(29, c(.1, .9), c(1, 1),
                    ftype = "hilbert", density = 32)
  expect_remez_shape(b3, 30, antisymmetric = TRUE)
  expect_remez_shape(b4, 29, antisymmetric = TRUE)
})

test_that("Type III and IV differentiators are antisymmetric", {
  for (order in c(30L, 29L)) {
    for (lo in c(0, .1)) {
      b <- remez_coeff(order, c(lo, .9), c(lo, .9),
                       ftype = "differentiator", density = 32)
      expect_remez_shape(b, order, antisymmetric = TRUE)
    }
  }
})

test_that("Type III differentiator has small weighted minimax error", {
  b <- remez_coeff(30, c(.1, .9), c(.1, .9),
                   ftype = "differentiator", density = 32)
  f <- seq(.1, .9, length.out = 4097L)
  amp <- odd_amplitude(b, f)
  sgn <- sign(amp[which.min(abs(f - .5))])
  expect_true(sgn != 0)
  err <- (sgn * amp - f) / f
  expect_lt(max(abs(err)), .0042)
  # Approximate balance of positive and negative ripple.
  expect_lt(abs(max(err) + min(err)), 5e-4)
})

test_that("Type IV differentiator approximates the desired slope", {
  b <- remez_coeff(29, c(.1, .9), c(.1, .9),
                   ftype = "differentiator", density = 32)
  f <- seq(.1, .9, length.out = 4097L)
  amp <- odd_amplitude(b, f)
  sgn <- sign(amp[which.min(abs(f - .5))])
  expect_true(sgn != 0)
  err <- (sgn * amp - f) / f
  # Broad correctness bound, not a hard-coded optimum.
  expect_lt(max(abs(err)), .02)
})

test_that("differentiator coefficients stabilize with grid density", {
  for (order in c(30L, 29L)) {
    b32 <- remez_coeff(order, c(.1, .9), c(.1, .9),
                       ftype = "differentiator", density = 32)
    b64 <- remez_coeff(order, c(.1, .9), c(.1, .9),
                       ftype = "differentiator", density = 64)
    expect_lt(max(abs(b32 - b64)), 1e-4)
  }
})

test_that("a multiband design produces finite symmetric taps", {
  b <- remez_coeff(30, c(0, .2, .3, .5, .6, 1),
                   c(0, 0, 1, 1, 0, 0), density = 32)
  expect_remez_shape(b, 30)
})

test_that("nonconstant desired response is evaluated consistently", {
  # Sloped bandpass-type design exercises frequency-based interpolation
  # without the differentiator-specific inverse-frequency weighting.
  b32 <- remez_coeff(30, c(0, .9), c(0, .9), density = 32)
  b64 <- remez_coeff(30, c(0, .9), c(0, .9), density = 64)
  expect_remez_shape(b32, 30)
  expect_lt(max(abs(b32 - b64)), 1e-3)
})

test_that("R wrapper rejects invalid inputs", {
  f <- c(0, .3, .4, 1)
  a <- c(1, 1, 0, 0)
  expect_error(remez(3, f, a), "n|length|order")
  expect_error(remez(30, c(0, .3, 1), c(1, 1, 0)), "even")
  expect_error(
    remez(30, f, c(1, 0)),
    "amplitudes matching f"
  )
  expect_error(remez(30, f, a, w = c(1)), "length")
  expect_error(remez(30, f, a, w = c(1, 0)), "positive")
  expect_error(remez(30, f, a, w = c(1, -1)), "positive")
  expect_error(remez(30, f, a, w = c(1, Inf)), "finite")
  expect_error(remez(30, f, a, density = 15), "density")
  expect_error(remez(30, c(0, .3, .4, 1.1), a), "range")
})

# Add the following tests after confirming the intended policy for
# zero-width bands, extremely narrow bands, and Type II at Nyquist.
# These cases should never trigger an out-of-bounds access or process crash.