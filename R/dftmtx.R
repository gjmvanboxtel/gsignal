# dftmtx.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave code:
# Copyright (C) 2003 David Bateman <adb014@gmail.com>
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
# 20201016  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Discrete Fourier Transform Matrix
#'
#' Compute the discrete Fourier transform matrix
#'
#' A discrete Fourier transform matrix is a complex matrix whose matrix product
#' with a vector computes the discrete Fourier transform of the vector.
#' \code{dftmtx} takes the FFT of the identity matrix to generate the transform
#' matrix. For a column vector \code{x}, \code{y <- dftmtx(n) * x} is the same
#' as \code{y <- fft(x, postpad(x, n)}. The inverse discrete Fourier transform
#' matrix is \code{inv <- Conj(dftmtx(n)) / n}.
#'
#' In general this is less efficient than calling the \code{fft} and \code{ifft}
#' functions directly.
#'
#' @param n Size of Fourier transformation matrix, specified as a positive
#'   integer.
#'
#' @return Fourier transform matrix.
#'
#' @examples
#' x <- seq_len(256)
#' y1 <- stats::fft(x)
#' n <- length(x)
#' y2 <- drop(x %*% dftmtx(n))
#' mx <- max(abs(y1 - y2))
#'
#' @author David Bateman, \email{adb014@@gmail.com}.\cr Conversion to R by Geert
#'   van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @seealso \code{\link[stats]{fft}}, \code{\link{ifft}}
#'
#' @export

 dftmtx <- function(n) {

  if (!isPosscal(n) || !isWhole(n)) {
    stop("n must be a positive integer")
  }

  y <- stats::mvfft(diag(1, n))
  y
}
