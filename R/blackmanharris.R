# blackmanharris.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2007 Sylvain Pelissier <sylvain.pelissier@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; see the file COPYING. If not, see
# <https://www.gnu.org/licenses/>.
#
# 20191210 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Blackman-Harris window
#'
#' Return the filter coefficients of a minimum four-term Blackman-Harris window.
#'
#' The Blackman window is a member of the family of cosine sum windows. It is a
#' generalization of the Hamming family, produced by adding more shifted sinc
#' functions, meant to minimize side-lobe levels.
#'
#' @param n Window length, specified as a positive integer.
#' @param method Character string. Window sampling method, specified as:
#' \describe{
#'   \item{"symmetric" (Default)}{Use this option when using windows for filter
#'   design.}
#'   \item{"periodic"}{This option is useful for spectral analysis because it
#'   enables a windowed signal to have the perfect periodic extension implicit
#'   in the discrete Fourier transform. When "periodic" is specified, the
#'   function computes a window of length \code{n + 1} and returns the first
#'   \code{n} points.}
#' }
#'
#' @return Blackman-Harris window, returned as a vector.
#'
#' @examples
#'
#' b <- blackmanharris(64)
#' plot (b, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' bs = blackmanharris(64,'symmetric')
#' bp = blackmanharris(63,'periodic')
#' plot (bs, type = "l", xlab = "Samples", ylab =" Amplitude")
#' lines(bp, col="red")
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

blackmanharris <- function(n, method = c("symmetric", "periodic")) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0) {
    stop("n must be an integer strictly positive")
  }
  method <- match.arg(method)

  if (method == "periodic") {
    N <- n
  } else if (method == "symmetric") {
    N <- n - 1
  } else {
    stop("method must be either 'periodic' or 'symmetric'")
  }

  if (n == 1) {
    w <- 1
  } else {
    a0 <- 0.35875
    a1 <- 0.48829
    a2 <- 0.14128
    a3 <- 0.01168
    k <- 0:(n - 1)
    w <- a0 - a1 * cos(2 * pi * k / N) +
      a2 * cos(4 * pi * k / N) - a3 * cos(6 * pi * k / N)
  }
  w
}
