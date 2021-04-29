# hilbert.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2007 Peter L. Soendergaard
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
# 20200709  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Hilbert transform
#'
#' Computes the extension of a real valued signal to an analytic signal.
#'
#' The function returns returns a complex helical sequence, sometimes called the
#' analytic signal, from a real data sequence. The analytic signal has a real
#' part, which is the original data, and an imaginary part, which contains the
#' Hilbert transform. The imaginary part is a version of the original real
#' sequence with a 90 degrees phase shift. Sines are therefore transformed to
#' cosines, and conversely, cosines are transformed to sines. The
#' Hilbert-transformed series has the same amplitude and frequency content as
#' the original sequence. The transform includes phase information that depends
#' on the phase of the original.
#'
#' @param x Input array, specified as a vector or a matrix. In case of a matrix,
#'   the Hilbert transform of all columns is computed.
#' @param n  use an n-point FFT to compute the Hilbert transform. The input data
#'   is zero-padded or truncated to length n, as appropriate.
#'
#' @return Analytic signal, of length \code{n}, returned as a complex vector or
#'   matrix, the real part of which contains the original signal, and the
#'   imaginary part of which contains the Hilbert transform of \code{x}.
#'
#' @examples
#' ## notice that the imaginary signal is phase-shifted 90 degrees
#' t <- seq(0, 10, length = 256)
#' z <- hilbert(sin(2 * pi * 0.5 * t))
#' plot(t, Re(z), type = "l", col="blue")
#' lines (t, Im(z), col = "red")
#' legend('topright', lty = 1, legend = c("Real", "Imag"),
#'        col = c("blue", "red"))
#'
#' ## the magnitude of the hilbert transform eliminates the carrier
#' t <- seq(0, 10, length = 1024)
#' x <- 5 * cos(0.2 * t) * sin(100 * t)
#' plot(t, x, type = "l", col = "green")
#' lines (t, abs(hilbert(x)), col = "blue")
#' legend('topright', lty = 1, legend = c("x", "|hilbert(x)|"),
#'         col = c("green", "blue"))
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'  Peter L. Soendergaard.\cr
#'  Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @references \url{https://en.wikipedia.org/wiki/Hilbert_transform},
#'   \url{https://en.wikipedia.org/wiki/Analytic_signal}
#'
#' @export

hilbert <- function(x, n = ifelse(is.vector(x), length(x), nrow(x))) {

  # check arguments
  if (!(is.vector(x) || is.matrix(x))) {
    stop("x must be a vector or a matrix")
  }
  if (is.character(x)) {
    stop("x must be a numeric vector or matrix")
  }
  if (!is.numeric(x)) {
    warning("imaginary parts discarded in coercion")
    x <- Re(x)
  }
  if (!isPosscal(n)) {
    stop("n must be a positive scalar")
  }

  # pad input to length n
  x <- postpad(x, n)

  # construct multiplication vector
  if (n %% 2 == 0) {
    v <- c(1, rep(2, n / 2 - 1), 1, rep(0, n / 2 - 1))
  } else {
    v <- c(1, rep(2, (n - 1) / 2), rep(0, (n - 1) / 2))
  }

  # compute the Hilbert transform
  if (is.vector(x)) {
    X <- stats::fft(x)
    Xv <- X * v
    y <- ifft(Xv)
  } else {
    X <- stats::mvfft(x)
    Xv <- X * v
    y <- imvfft(Xv)
  }
  y
}
