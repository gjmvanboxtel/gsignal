# fht.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2008 Muthiah Annamalai <muthiah.annamalai@uta.edu>
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
# 20201020  GvB       setup for gsignal v0.1.0
# 20201023  GvB       corrected padding
#------------------------------------------------------------------------------

#' Fast Hartley Transform
#'
#' Compute the (inverse) Hartley transform of a signal using FFT
#'
#' The Hartley transform is an integral transform closely related to the Fourier
#' transform, but which transforms real-valued functions to real-valued
#' functions. Compared to the Fourier transform, the Hartley transform has
#' the advantages of transforming real functions to real functions (as opposed
#' to requiring complex numbers) and of being its own inverse [1].
#'
#' This function implements the Hartley transform by calculating the difference
#' between the real- and imaginary-valued parts of the Fourier-transformed
#' signal [1]. The forward and inverse Hartley transforms are the same (except
#' for a scale factor of 1/N for the inverse Hartley transform), but implemented
#' using different functions.
#'
#' @param x input data, specified as a numeric vector or matrix. In case of a
#'   vector it represents a single signal; in case of a matrix each column is a
#'   signal.
#' @param n transform length, specified as a positive integer scalar. Default:
#'   \code{NROW(x)}.
#'
#' @return (inverse) Hartley transform, returned as a vector or matrix.
#'
#' @examples
#' # FHT of a 2.5 Hz signal with offset
#' fs <- 100
#' secs <- 10
#' freq <- 2.5
#' t <- seq(0, secs - 1 / fs, 1 / fs)
#' x <- 5 * t + 50 * cos(freq * 2 * pi * t)
#' X <- fht(x)
#' op <- par(mfrow = c(2, 1))
#' plot(t, x, type = "l", xlab = "", ylab = "", main = "Signal")
#' f <- seq(0, fs - (1 / fs), length.out = length(t))
#' to <- which(f >= 5)[1]
#' plot(f[1:to], X[1:to], type = "l", xlab = "", ylab = "",
#'      main = "Hartley Transform")
#' par(op)
#'
#' @author Muthiah Annamalai, \email{muthiah.annamalai@@uta.edu}.\
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Hartley_transform}
#'
#' @seealso \code{\link{fft}}
#'
#' @rdname fht
#' @export

 fht <- function(x, n = NROW(x)) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x)) {
    stop("x must be a numeric or vector or matrix")
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  if (n != nr) {
    x <- postpad(x, n)
  }

  Y <- stats::mvfft(x)
  y <- Re(Y) - Im(Y)

  if (vec) {
    y <- as.vector(y)
  }
  y
}

#' @rdname fht
#' @export

ifht <- function(x, n = NROW(x)) {

  # check parameters
  if (!(is.vector(x) || is.matrix(x)) || !is.numeric(x)) {
    stop("x must be a numeric or vector or matrix")
  }

  if (is.vector(x)) {
    vec <- TRUE
    x <- as.matrix(x, ncol = 1)
  } else {
    vec <- FALSE
  }
  nr <- nrow(x)

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }
  if (n != nr) {
    x <- postpad(x, n)
  }

  Y <- imvfft(x)
  y <- Re(Y) + Im(Y)

  if (vec) {
    y <- as.vector(y)
  }
  y
 }
