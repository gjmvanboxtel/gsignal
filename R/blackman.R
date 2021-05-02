# blackman.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995-2017 Andreas Weingessel
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

#' Blackman window
#'
#' Return the filter coefficients of a Blackman window.
#'
#' The Blackman window is a member of the family of cosine sum windows.
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
#' @return Blackman window, returned as a vector.
#'
#' @examples
#'
#' h <- blackman(64)
#' plot (h, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' bs = blackman(64,'symmetric')
#' bp = blackman(63,'periodic')
#' plot (bs, type = "l", xlab = "Samples", ylab =" Amplitude")
#' lines(bp, col="red")
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

blackman <- function(n, method = c("symmetric", "periodic")) {

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
    n <- n - 1
    k <- (0:n) / N
    w <- 0.42 - 0.5 * cos(2 * pi * k) + 0.08 * cos(4 * pi * k)
  }
  w
}
