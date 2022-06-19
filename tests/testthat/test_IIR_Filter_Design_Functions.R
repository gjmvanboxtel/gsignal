# gsignal IIR filter design functions
library(gsignal)
library(testthat)

tol <- 1e-6

# -----------------------------------------------------------------------
# cheb()

test_that("parameters to cheb() are correct", {
  expect_error(cheb())
  expect_error(cheb(0.5))
  expect_error(cheb(-1L))
  expect_error(cheb(array(1L, c(1, 4))))
})

test_that("cheb() tests are correct", {
  expect_equal(cheb(1, 1), 1)
  expect_equal(cheb(2, 1), 1)
  expect_equal(cheb(5, 2), 362)
  expect_equal(cheb(5, c(2,3)), c(362, 3363))
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
               tolerance = tol)
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
  expect_equal(zpg$z, c(0-1.414214i, 0+1.414214i), tolerance = tol)
  expect_equal(zpg$p, c(-0.5+1.322876i, -0.5-1.322876i), tolerance = tol)
  expect_equal(zpg$g, 1L)

  zpg <- besself(1, c(1, 2), 'pass')
  expect_equal(zpg$z, 0L)
  expect_equal(zpg$p, c(-0.5+1.322876i, -0.5-1.322876i), tolerance = tol)
  expect_equal(zpg$g, 1L)

  zpg <- besself(2, 1, 'low')
  expect_equal(zpg$z, complex(0))
  expect_equal(zpg$p, c(-0.8660254+0.5i, -0.8660254-0.5i), tolerance = tol)
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, 1, 'high')
  expect_equal(zpg$z, c(0L, 0L))
  expect_equal(zpg$p, c(-0.8660254-0.5i, -0.8660254+0.5i), tolerance = tol)
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, c(1, 2), 'stop')
  expect_equal(zpg$z, c(0-1.414214i, 0+1.414214i, 0-1.414214i, 0+1.414214i), tolerance = tol)
  expect_equal(zpg$p, c(-0.354087+1.121579i, -0.354087-1.121579i, -0.511939-1.621579i, -0.511939+1.621579i), tolerance = tol)
  expect_equal(zpg$g, 1L)
  
  zpg <- besself(2, c(1, 2), 'pass')
  expect_equal(zpg$z, c(0L, 0L))
  expect_equal(zpg$p, c( -0.354087-1.121579i, -0.354087+1.121579i, -0.511939+1.621579i, -0.511939-1.621579i), tolerance = tol)
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
  expect_equal(res$p, 1 / 3, tolerance = tol)
  expect_equal(res$g, 1/3)

  res <- sftrans(1, 3, 1, 1, FALSE)
  expect_equal(res$z, 1)
  expect_equal(res$p, 3)
  expect_equal(res$g, 1)
  
})

# -----------------------------------------------------------------------
# buttord()

test_that("parameters to buttord() are correct", {
  expect_error(buttord())
  expect_error(buttord(.1))
  expect_error(buttord(.1, .2))
  expect_error(buttord(c(.1, .1), c(.2, .2), 3, 4))
  expect_error(buttord(c(.1, .2), c(.5, .6), 3, 4))
  expect_error(buttord(c(.1, .5), c(.2, .6), 3, 4))
  expect_error(buttord(.1, .2, 3, 4, 5))
  expect_error(buttord(1, 2, 3, 4, 's', 6))
})

test_that("buttord() tests are correct", {
  
  # Analog band-pass
  res <- buttord(2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10436), 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), c(61903, 63775))
  expect_equal(round(res$Wc_s), c(61575, 64114))
    
  # Analog band-pass
  res <- buttord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9582, 11000), 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), c(61903, 63775))
  expect_equal(round(res$Wc_s), c(61575, 64115))

  # Analog band-pass
  res <- buttord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10437), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61850, 63830))
  expect_equal(round(res$Wc_s), c(61848, 63831))

  # Analog band-pass
  res <- buttord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9581, 11000), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61850, 63830))
  expect_equal(round(res$Wc_s), c(61847, 63832))

  # Analog high-pass
  res <- buttord (2 * pi * 13583, 2 * pi * 4000, 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), 72081)
  expect_equal(round(res$Wc_s), 53101)

  # Analog high-pass
  res <- buttord (2 * pi * 13584, 2 * pi * 4000, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 68140)
  expect_equal(round(res$Wc_s), 68138)

  # Analog low-pass
  res <- buttord (2 * pi * 4000, 2 * pi * 13583, 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), 29757)
  expect_equal(round(res$Wc_s), 40394)

  # Analog low-pass
  res <- buttord (2 * pi * 4000, 2 * pi * 13584, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 31481)
  expect_equal(round(res$Wc_s), 31482)

  # Analog notch (narrow band-stop)
  res <- buttord (2 * pi * c(9000, 10436), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), c(60607, 65138))
  expect_equal(round(res$Wc_s), c(61184, 64524))
  
  # Analog notch (narrow band-stop)
  res <- buttord (2 * pi * c(9582, 11000), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 4)
  expect_equal(round(res$Wc), c(60606, 65139))
  expect_equal(round(res$Wc_s), c(61184, 64524))

  # Analog notch (narrow band-stop)
  res <- buttord (2 * pi * c(9000, 10437), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(60722, 65015))
  expect_equal(round(res$Wc_s), c(60726, 65011))

  # Analog notch (narrow band-stop)
  res <- buttord (2 * pi * c(9581, 11000), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(60721, 65016))
  expect_equal(round(res$Wc_s), c(60726, 65011))

  # Digital band-pass
  fs <- 44100
  res <- buttord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10051), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), c(9477, 9773))
  expect_equal(round(Wc_s), c(9425, 9826))

  # Digital band-pass
  fs <- 44100
  res <- buttord (2 / fs * c(9500, 9750), 2 / fs * c(9204, 10700), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), c(9477, 9773))
  expect_equal(round(Wc_s), c(9425, 9826))

  # Digital band-pass
  fs <- 44100
  res <- buttord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10052), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9469, 9782))
  expect_equal(round(Wc_s), c(9468, 9782))

  # Digital band-pass
  fs <- 44100
  res <- buttord (2 / fs * c(9500, 9750), 2 / fs * c(9203, 10700), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9469, 9782))
  expect_equal(round(Wc_s), c(9468, 9782))
  
  # Digital high-pass
  fs <- 44100
  res <- buttord (2 / fs * 10987, 2 / fs * 4000, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), 9808)
  expect_equal(round(Wc_s), 7780)
  
  # Digital high-pass
  fs <- 44100
  res <- buttord (2 / fs * 10988, 2 / fs * 4000, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 9421)
  expect_equal(round(Wc_s), 9421)
  
  # Digital low-pass
  fs <- 44100
  res <- buttord (2 / fs * 4000, 2 / fs * 10987, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), 4686)
  expect_equal(round(Wc_s), 6176)

  # Digital low-pass
  fs <- 44100
  res <- buttord (2 / fs * 4000, 2 / fs * 10988, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 4936)
  expect_equal(round(Wc_s), 4936)

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- buttord (2 / fs * c(8500, 10833), 2 / fs * c(9875, 10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), c(9369, 10640))
  expect_equal(round(Wc_s), c(9605, 10400))

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- buttord (2 / fs * c(9183, 11000), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 4)
  expect_equal(round(Wc), c(9370, 10640))
  expect_equal(round(Wc_s), c(9605, 10400))  

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- buttord (2 / fs * c(8500, 10834), 2 / fs * c(9875, 10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9421, 10587))
  expect_equal(round(Wc_s), c(9422, 10587))  

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- buttord (2 / fs * c(9182, 11000), 2 / fs * c(9875, 10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9421, 10587))
  expect_equal(round(Wc_s), c(9422, 10587))  

})

# -----------------------------------------------------------------------
# butter()

test_that("parameters to butter() are correct", {
  expect_error(butter())
  expect_error(butter(1))
  expect_error(butter(1, 2, 3, 4, 5))
  expect_error(butter(.5, .2))
  expect_error(butter(3, .2, "invalid"))
  expect_error(butter(9, .6, "stop"))
  expect_error(butter(9, .6, "pass"))
  expect_error(butter(9, .6, "pass", "q"))
  expect_error(butter(9, .6, "pass", "z", "invalid"))
})

test_that("butter() tests are correct", {
  
  # shared sf, sf2, off_db
  off_db <- 0.5
  fs <- 6000; fs2 <- fs / 2
  sinetone <- function(f, r, s, a) a * sin(2 * pi * f * seq(0, s, length.out = r * s))
  data <- cbind(sinetone(5,fs,10,1), sinetone(10,fs,10,1), sinetone(50,fs,10,1), sinetone(200,fs,10,1), sinetone(400,fs,10,1))
  l <- nrow(data)
  
  # Test low pass order 1 with 3dB @ 50Hz
  bf <-  butter ( 1, 50 / fs2 )
  filtered <- NULL; for (i in 1:5) filtered <- cbind(filtered, filter(bf, data[, i]))
  damp_db <- NULL; for (i in 1:5) damp_db <- cbind(damp_db, 20 * log10(max(filtered[(l - fs):l, i])))
  expect_equal(c(damp_db[4] - damp_db[5], damp_db[1:3]), c(6, 0, 0, -3), tolerance = off_db)

  # Test low pass order 4 with 3dB @ 50Hz
  bf <- butter(4, 50 / fs2)
  filtered <- NULL; for (i in 1:5) filtered <- cbind(filtered, filter(bf, data[, i]))
  damp_db <- NULL; for (i in 1:5) damp_db <- cbind(damp_db, 20 * log10(max(filtered[(l - fs):l, i])))
  expect_equal(c(damp_db[4] - damp_db[5], damp_db[1:3]), c(24, 0, 0, -3), tolerance = off_db)

  # Test high pass order 1 with 3dB @ 50Hz
  bf <- butter(1, 50 / fs2, "high")
  filtered <- NULL; for (i in 1:5) filtered <- cbind(filtered, filter(bf, data[, i]))
  damp_db <- NULL; for (i in 1:5) damp_db <- cbind(damp_db, 20 * log10(max(filtered[(l - fs):l, i])))
  expect_equal(c(damp_db[2] - damp_db[1], damp_db[3:5]), c(6, -3, 0, 0), tolerance = off_db)
  
  # Test high pass order 4 with 3dB @ 50Hz
  bf <- butter(4, 50 / fs2, "high")
  filtered <- NULL; for (i in 1:5) filtered <- cbind(filtered, filter(bf, data[, i]))
  damp_db <- NULL; for (i in 1:5) damp_db <- cbind(damp_db, 20 * log10(max(filtered[(l - fs):l, i])))
  expect_equal(c(damp_db[2] - damp_db[1], damp_db[3:5]), c(24, -3, 0, 0), tolerance = off_db)

  # Test outut formats
  zpg <- butter(3, 0.05, output = "Zpg")
  expect_equal(as.Arma(zpg), butter(3, 0.05))
  sos <- butter(3, 0.05, output = "Sos")
  expect_equal(as.Arma(sos), butter(3, 0.05))
})

# -----------------------------------------------------------------------
# cheb1ord()

test_that("parameters to cheb1ord() are correct", {
  expect_error(cheb1ord())
  expect_error(cheb1ord(.1))
  expect_error(cheb1ord(.1, .2))
  expect_error(cheb1ord(c(.1, .1), c(.2, .2), 3, 4))
  expect_error(cheb1ord(c(.1, .2), c(.5, .6), 3, 4))
  expect_error(cheb1ord(c(.1, .5), c(.2, .6), 3, 4))
  expect_error(cheb1ord(.1, .2, 3, 4, 5))
  expect_error(cheb1ord(1, 2, 3, 4, 's', 6))
})

test_that("cheb1ord() tests are correct", {
  
  # Analog band-pass
  res <- cheb1ord(2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10437), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(62046, 63627))
  expect_equal(round(res$Wc_s), c(61652, 64035))
  
  # Analog band-pass
  res <- cheb1ord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9581, 12000), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(62046, 63627))
  expect_equal(round(res$Wc_s), c(61651, 64036))
  
  # Analog high-pass
  res <- cheb1ord (2 * pi * 13584, 2 * pi * 4000, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 85351)
  expect_equal(round(res$Wc_s), 56700)
  
  # Analog high-pass
  res <- cheb1ord (2 * pi * 13584, 2 * pi * 4000, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 85351)
  expect_equal(round(res$Wc_s), 56700)
  
  # Analog low-pass
  res <- cheb1ord (2 * pi * 4000, 2 * pi * 13584, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 25133)
  expect_equal(round(res$Wc_s), 37832)

  # Analog notch (narrow band-stop)
  res <- cheb1ord (2 * pi * c(9000, 10437), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(60201, 65578))
  expect_equal(round(res$Wc_s), c(61074, 64640))
  
  # Analog notch (narrow band-stop)
  res <- cheb1ord (2 * pi * c(9581, 12000), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(60199, 65580))
  expect_equal(round(res$Wc_s), c(61074, 64640))
  
  # Digital band-pass
  fs <- 44100
  res <- cheb1ord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10052), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9500, 9750))
  expect_equal(round(Wc_s), c(9437, 9814))
  
  # Digital band-pass
  fs <- 44100
  res <- cheb1ord (2 / fs * c(9500, 9750), 2 / fs * c(9182, 12000), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9500, 9750))
  expect_equal(round(Wc_s), c(9428, 9823))
  
  # Digital high-pass
  fs <- 44100
  res <- cheb1ord (2 / fs * 10988, 2 / fs * 4000, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 10988)
  expect_equal(round(Wc_s), 8197)
  
  # Digital low-pass
  fs <- 44100
  res <- cheb1ord (2 / fs * 4000, 2 / fs * 10988, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 4000)
  expect_equal(round(Wc_s), 5829)
  
  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- cheb1ord (2 / fs * c(8500, 10834), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9182, 10834))
  expect_equal(round(Wc_s), c(9475, 10532))
  
  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- cheb1ord (2 / fs * c(9182, 12000), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9182, 10834))
  expect_equal(round(Wc_s), c(9475, 10532))  

})

# -----------------------------------------------------------------------
# cheby1()

test_that("parameters to cheby1() are correct", {
  expect_error(cheby1())
  expect_error(cheby1(1))
  expect_error(cheby1(1, 2, 3, 4, 5))
  expect_error(cheby1(.5, .2))
  expect_error(cheby1(3, .2, 0.5, "invalid"))
  expect_error(cheby1(9, .6, 0.5, "stop"))
  expect_error(cheby1(9, .6, 0.5, "pass"))
  expect_error(cheby1(9, .6, 0.5, "pass", "q"))
})

# -----------------------------------------------------------------------
# cheb2ord()

test_that("parameters to cheb2ord() are correct", {
  expect_error(cheb2ord())
  expect_error(cheb2ord(.1))
  expect_error(cheb2ord(.1, .2))
  expect_error(cheb2ord(c(.1, .1), c(.2, .2), 3, 4))
  expect_error(cheb2ord(c(.1, .2), c(.5, .6), 3, 4))
  expect_error(cheb2ord(c(.1, .5), c(.2, .6), 3, 4))
  expect_error(cheb2ord(.1, .2, 3, 4, 5))
  expect_error(cheb2ord(1, 2, 3, 4, 's', 6))
})

test_that("cheb2ord() tests are correct", {
  
  # Analog band-pass
  res <- cheb2ord(2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10437), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61074, 64640))
  expect_equal(round(res$Wc_s), c(60201, 65578))
  
  # Analog band-pass
  res <- cheb2ord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9581, 12000), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61074, 64640))
  expect_equal(round(res$Wc_s), c(60199, 65580))
  
  # Analog high-pass
  res <- cheb2ord (2 * pi * 13584, 2 * pi * 4000, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 37832)
  expect_equal(round(res$Wc_s), 25133)
  
  # Analog low-pass
  res <- cheb2ord (2 * pi * 4000, 2 * pi * 13584, 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 56700)
  expect_equal(round(res$Wc_s), 85351)
  
  # Analog notch (narrow band-stop)
  res <- cheb2ord (2 * pi * c(9000, 10437), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61652, 64035))
  expect_equal(round(res$Wc_s), c(62046, 63627))
  
  # Analog notch (narrow band-stop)
  res <- cheb2ord (2 * pi * c(9581, 12000), 2 * pi * c(9875, 10126.5823), 1, 26, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(61651, 64036))
  expect_equal(round(res$Wc_s), c(62046, 63627))
  
  # Digital band-pass
  fs <- 44100
  res <- cheb2ord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10052), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9344, 9908))
  expect_equal(round(Wc_s), c(9203, 10052))
  
  # Digital band-pass
  fs <- 44100
  res <- cheb2ord (2 / fs * c(9500, 9750), 2 / fs * c(9182, 12000), 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9344, 9908))
  expect_equal(round(Wc_s), c(9182, 10073))
  
  # Digital high-pass
  fs <- 44100
  res <- cheb2ord (2 / fs * 10988, 2 / fs * 4000, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 5829)
  expect_equal(round(Wc_s), 4000)
  
  # Digital low-pass
  fs <- 44100
  res <- cheb2ord (2 / fs * 4000, 2 / fs * 10988, 1, 26)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 8197)
  expect_equal(round(Wc_s), 10988)
  
  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- cheb2ord (2 / fs * c(8500, 10834), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9804, 10198))
  expect_equal(round(Wc_s), c(9875, 10127))
  
  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- cheb2ord (2 / fs * c(9182, 12000), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  Wc_s <- res$Wc_s * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9804, 10198))
  expect_equal(round(Wc_s), c(9875, 10127))  
  
})

# -----------------------------------------------------------------------
# cheby2()

test_that("parameters to cheby2() are correct", {
  expect_error(cheby2())
  expect_error(cheby2(1))
  expect_error(cheby2(1, 2, 3, 4, 5))
  expect_error(cheby2(.5, .2))
  expect_error(cheby2(3, .2, 0.5, "invalid"))
  expect_error(cheby2(9, .6, 0.5, "stop"))
  expect_error(cheby2(9, .6, 0.5, "pass"))
  expect_error(cheby2(9, .6, 0.5, "pass", "q"))
})

# -----------------------------------------------------------------------
# ellipord()

test_that("parameters to ellipord() are correct", {
  expect_error(ellipord())
  expect_error(ellipord(.1))
  expect_error(ellipord(.1, .2))
  expect_error(ellipord(c(.1, .1), c(.2, .2), 3, 4))
  expect_error(ellipord(c(.1, .2), c(.5, .6), 3, 4))
  expect_error(ellipord(c(.1, .5), c(.2, .6), 3, 4))
  expect_error(ellipord(.1, .2, 3, 4, 5))
  expect_error(ellipord(1, 2, 3, 4, 's', 6))
})

test_that("ellipord() tests are correct", {
  
  # Analog band-pass
  res <- ellipord(2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10657), 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), c(62046, 63627))

  # Analog band-pass
  res <- ellipord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9384, 12000), 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), c(62046, 63627))

  # Analog band-pass
  res <- ellipord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9000, 10656), 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(62046, 63627))

  # Analog band-pass
  res <- ellipord (2 * pi * c(9875, 10126.5823), 2 * pi * c(9385, 12000), 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(62046, 63627))
  
  # Analog high-pass
  res <- ellipord (2 * pi * 20224, 2 * pi * 4000, 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), 127071)

  # Analog high-pass
  res <- ellipord (2 * pi * 20223, 2 * pi * 4000, 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 127065)
  
  # Analog low-pass
  res <- ellipord (2 * pi * 4000, 2 * pi * 20224, 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), 25133)

  # Analog low-pass
  res <- ellipord (2 * pi * 4000, 2 * pi * 20223, 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), 25133)
  
  # Analog notch (narrow band-stop)
  res <- ellipord (2 * pi * c(9000, 10657), 2 * pi * c(9875, 10126.5823), 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), c(58958, 66960))

  # Analog notch (narrow band-stop)
  res <- ellipord (2 * pi * c(9384, 12000), 2 * pi * c(9875, 10126.5823), 3, 40, "s")
  expect_equal(res$n, 2)
  expect_equal(round(res$Wc), c(58961 , 66956))

  # Analog notch (narrow band-stop)
  res <- ellipord (2 * pi * c(9000, 10656), 2 * pi * c(9875, 10126.5823), 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(58964, 66954))

  # Analog notch (narrow band-stop)
  res <- ellipord (2 * pi * c(9385, 12000), 2 * pi * c(9875, 10126.5823), 3, 40, "s")
  expect_equal(res$n, 3)
  expect_equal(round(res$Wc), c(58968, 66949))
  
  # Digital band-pass
  fs <- 44100
  res <- ellipord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10261), 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), c(9500, 9750))

  # Digital band-pass
  fs <- 44100
  res <- ellipord (2 / fs * c(9500, 9750), 2 / fs * c(9000, 10700), 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), c(9500, 9750))

  # Digital band-pass
  fs <- 44100
  res <- ellipord (2 / fs * c(9500, 9750), 2 / fs * c(8500, 10260), 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9500, 9750))

  # Digital band-pass
  fs <- 44100
  res <- ellipord (2 / fs * c(9500, 9750), 2 / fs * c(9001, 10700), 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(9500, 9750))
  
  # Digital high-pass
  fs <- 44100
  res <- ellipord (2 / fs * 13713, 2 / fs * 4000, 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), 13713)

  # Digital high-pass
  fs <- 44100
  res <- ellipord (2 / fs * 13712, 2 / fs * 4000, 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 13712)
  
  # Digital low-pass
  fs <- 44100
  res <- ellipord (2 / fs * 4000, 2 / fs * 13713, 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), 4000)

  # Digital low-pass
  fs <- 44100
  res <- ellipord (2 / fs * 4000, 2 / fs * 13712, 3, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), 4000)
  
  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- ellipord (2 / fs * c(8500, 11073), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), c(8952, 11073))

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- ellipord (2 / fs * c(8952, 12000), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 2)
  expect_equal(round(Wc), c(8952, 11073))

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- ellipord (2 / fs * c(8500, 11072), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(8953, 11072))

  # Digital notch (narrow band-stop)
  fs <- 44100
  res <- ellipord (2 / fs * c(8953, 12000), 2 / fs * c(9875,  10126.5823), 0.5, 40)
  Wc <- res$Wc * fs / 2
  expect_equal(res$n, 3)
  expect_equal(round(Wc), c(8953, 11072))
  
})

# -----------------------------------------------------------------------
# ellip()

test_that("parameters to ellip() are correct", {
  expect_error(ellip())
  expect_error(ellip(1))
  expect_error(ellip(1, 2))
  expect_error(ellip(1, 2, 3))
  expect_error(ellip(1, 2, 3, 4, 5, 6, 7))
  expect_error(ellip(0.5, 2, 40, .2))
  expect_error(ellip(3, 2, 40, .2, "invalid"))
  expect_error(ellip(3, 2, 40, .2, "low", "invalid"))
})

# -----------------------------------------------------------------------
# pei_tseng_notch()

test_that("parameters to pei_tseng_notch() are correct", {
  expect_error(pei_tseng_notch())
  expect_error(pei_tseng_notch(1))
  expect_error(pei_tseng_notch(1, 2, 3))
  expect_error(pei_tseng_notch(c(1, 2), 3))
  expect_error(pei_tseng_notch(-1, 1))
  expect_error(pei_tseng_notch(1, -1))
  expect_error(pei_tseng_notch(3, "invalid"))
  expect_error(pei_tseng_notch("invalid", 2))
  expect_error(pei_tseng_notch(matrix(0, 2, 2), 2))
  expect_error(pei_tseng_notch(1, matrix(0L, 2, 2)))
  
})

test_that("pei_tseng_notch() tests are correct", {
  
  sinetone <- function(f, r, s, a) a * sin(2 * pi * f * seq(0, s, length.out = r * s))
  
  ## 2Hz bandwidth
  
  fs <- 800; nyq <- fs / 2
  data <- cbind(sinetone(49, fs, 10, 1), sinetone(50, fs, 10, 1), sinetone(51, fs, 10, 1))
  l <- nrow(data)
  ba <-  pei_tseng_notch( 50 / nyq, 2 / nyq)
  filtered <- filter(ba, data)
  damp_db <- apply(filtered, 2, function(x) 20 * log10(max(x[(l - 1000):l])))
  expect_equal(as.vector(damp_db), c(-3.037382, -44.16588, -3.065681), tolerance = tol)

  ## 1Hz bandwidth
  data <- cbind(sinetone(49.5, fs, 10, 1), sinetone(50, fs, 10, 1), sinetone(50.5, fs, 10, 1))
  l <- nrow(data)
  ba <-  pei_tseng_notch( 50 / nyq, 1 / nyq)
  filtered <- filter(ba, data)
  damp_db <- apply(filtered, 2, function(x) 20 * log10(max(x[(l - 1000):l])))
  expect_equal(as.vector(damp_db), c(-3.064986, -38.10409, -2.997267), tolerance = tol)
})
  
# -----------------------------------------------------------------------
# cheb2ap()

test_that("parameters to cheb2ap() are correct", {
  expect_error(cheb2ap())
  expect_error(cheb2ap(1))
  expect_error(cheb2ap(-1, 3))
  expect_error(cheb2ap(3, -1))
})

