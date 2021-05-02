# flattopwin.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Author: Paul Kienzle <pkienzle@users.sf.net> (2004)
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
# 20191211 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Flat top window
#'
#' Return the filter coefficients of a flat top window.
#'
#' The Flat Top window is defined by the function:
#' \deqn{f(w) = 1 - 1.93 cos(2 \pi w) + 1.29 cos(4 \pi w) - 0.388 cos(6 \pi w) +
#' 0.0322 cos(8 \pi w)}
#' where \code{w = i/(n-1)} for \code{i=0:n-1} for a symmetric window, or
#' \code{w = i/n} for \code{i=0:n-1} for a periodic window. The default is
#' symmetric. The returned window is normalized to a peak of 1 at w = 0.5.
#'
#' Flat top windows have very low passband ripple (< 0.01 dB) and are used
#' primarily for calibration purposes. Their bandwidth is approximately 2.5
#' times wider than a Hann window.
#'
#' @param n Window length, specified as a positive integer.
#' @param method Character string. Window sampling method, specified as:
#' \describe{
#'   \item{"symmetric"}{(Default). Use this option when using windows for filter
#'   design.}
#'   \item{"periodic"}{This option is useful for spectral analysis because it
#'   enables a windowed signal to have the perfect periodic extension implicit
#'   in the discrete Fourier transform. When 'periodic' is specified, the
#'   function computes a window of length \code{n + 1} and returns the first
#'   \code{n} points.}
#' }
#'
#' @return Flat top window, returned as a vector.
#'
#' @examples
#'
#' ft <- flattopwin(64)
#' plot (ft, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

flattopwin <- function(n, method = c("symmetric", "periodic")) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
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
    x <- 2 * pi * (0:(n - 1)) / N
    w <- (1 - 1.93 * cos(x) + 1.29 * cos(2 * x) -
            0.388 * cos(3 * x) + 0.0322 * cos(4 * x)) / 4.6402
  }
  w
}
