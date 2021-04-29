# czt.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2004 Daniel Gunyan
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Version history
# 20201011  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Chirp Z-transform
#'
#' Compute the Chirp Z-transform along a spiral contour on the z-plane.
#'
#' The chirp Z-transform (CZT) is a generalization of the discrete Fourier
#' transform (DFT). While the DFT samples the Z plane at uniformly-spaced points
#' along the unit circle, the chirp Z-transform samples along spiral arcs in the
#' Z-plane, corresponding to straight lines in the S plane. The DFT, real DFT,
#' and zoom DFT can be calculated as special cases of the CZT[1]. For the
#' specific case of the DFT, \code{a = 0}, \code{m = NCOL(x)}, and \code{w = 2 *
#' pi / m}[2, p. 656].
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param m transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#' @param w ratio between spiral contour points in each step (i.e., radius
#'   increases exponentially, and angle increases linearly), specified as a
#'   complex scalar. Default: \code{exp(0-1i * 2 * pi / m)}.
#' @param a initial spiral contour point, specified as a complex scalar.
#'   Default: 1.
#'
#' @return Chirp Z-transform, returned as a vector or matrix.
#'
#' @examples
#' fs <- 1000                                           # sampling frequency
#' secs <- 10                                           # number of seconds
#' t <- seq(0, secs, 1/fs)                              # time series
#' x <- sin(100 * 2 * pi * t) + runif(length(t))        # 100 Hz signal + noise
#' m <- 32                                              # n of points desired
#' f0 <- 75; f1 <- 175;                                 # desired freq range
#' w <- exp(-1i * 2 * pi * (f1 - f0) / ((m - 1) * fs))  # freq step of f1-f0/m
#' a <- exp(1i * 2 * pi * f0 / fs);                     # starting at freq f0
#' y <- czt(x, m, w, a)
#'
#' # compare DFT and FFT
#' fs <- 1000
#' h <- as.numeric(fir1(100, 125/(fs / 2), type = "low"))
#' m <- 1024
#' y <- stats::fft(postpad(h, m))
#'
#' f1 <- 75; f2 <- 175;
#' w <- exp(-1i * 2 * pi * (f2 - f1) / (m * fs))
#' a <- exp(1i * 2 * pi * f1 / fs)
#' z <- czt(h, m, w, a)
#'
#' fn <- seq(0, m - 1, 1) / m
#' fy <- fs * fn
#' fz = (f2 - f1) * fn + f1
#' plot(fy, 10 * log10(abs(y)), type = "l", xlim = c(50, 200),
#'   xlab = "Frequency", ylab = "Magnitude (dB")
#' lines(fz, 10 * log10(abs(z)), col = "red")
#' legend("topright", legend = c("FFT", "CZT"), col=1:2, lty = 1)
#'
#' @author Daniel Gunyan.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Chirp_Z-transform}\cr
#'   [2]Oppenheim, A.V., Schafer, R.W., and Buck, J.R. (1999). Discrete-Time
#'   Signal Processing, 2nd edition. Prentice-Hall.
#'
#' @export

czt <- function(x, m = NROW(x),
                w = exp(complex(real = 0, imaginary = -2 * pi / m)),
                a = 1) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !(is.numeric(x) || is.complex(x))) {
    stop("x must be a numeric or complex vector or matrix")
  }

  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
    vec <- TRUE
  } else {
    vec <- FALSE
  }
  n <- nrow(x)

  if (!isPosscal(m) || !isWhole(m)) {
    stop("m must be a positive integer")
  }

  if (!(is.numeric(w) || is.complex(w)) || length(w) != 1) {
    stop("w must be a single complex value")
  }
  w <- as.complex(w)

  if (!(is.numeric(a) || is.complex(a)) || length(a) != 1) {
    stop("a must be a single complex value")
  }
  a <- as.complex(a)

  # indexing to make the statements a little more compact
  N <- seq(0, n - 1, 1) + n
  NM <- seq(- (n - 1), (m - 1), 1) + n
  M <- seq(0, m - 1, 1) + n

  nfft <- nextpow2(n + m - 1) # fft pad
  W2 <- w ^ ((seq(- (n - 1), max(m - 1, n - 1), 1)^2) / 2) # chirp

  fg <- stats::mvfft(postpad(x * (a ^ - (N - n)) * W2[N], nfft))
  fw <- stats::fft(postpad(1 / W2[NM], nfft))
  gg <- imvfft(fg * fw)
  y <-  gg[M, ] * W2[M]

  if (vec) {
    y <- as.vector(y)
  }
  y
}
