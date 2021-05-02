# hanning.R
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
# 20191209 GbB            First version for v0.1.0
# 20200413 GvB            corrected definition of hanning
# 20200606 GvB            exported hanning from namespace
#------------------------------------------------------------------------------

#' Hann window
#'
#' Return the filter coefficients of a Hann window of length \code{n}.
#'
#' The Hann window is a member of the family of cosine sum windows. It was named
#' after Julius von Hann, and is sometimes referred to as Hanning, presumably
#' due to its linguistic and formulaic similarities to Hamming window.
#'
#' @param n Window length, specified as a positive integer.
#' @param method Character string. Window sampling method, specified as:
#' \describe{
#'   \item{"symmetric"}{(Default). Use this option when using windows for filter
#'   design.}
#'   \item{"periodic"}{This option is useful for spectral analysis because it
#'   enables a windowed signal to have the perfect periodic extension implicit
#'   in the discrete Fourier transform. When \code{"periodic"} is specified, the
#'   function computes a window of length \code{n + 1} and returns the first
#'   \code{n} points.}
#' }
#'
#' @return Hann window, returned as a vector.
#'
#' @examples
#'
#' h <- hann(64)
#' plot (h, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' hs = hann(64,'symmetric')
#' hp = hann(63,'periodic')
#' plot (hs, type = "l", xlab = "Samples", ylab =" Amplitude")
#' lines(hp, col="red")
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @rdname hann
#' @export

hann <- function(n, method = c("symmetric", "periodic")) {

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
    n <- n - 1
    w <- 0.5 - 0.5 * cos(2 * pi * (0:n) / N)
  }
  w
}

#' @rdname hann
#' @export
hanning <- function(n, method = c("symmetric", "periodic")) hann(n, method)
