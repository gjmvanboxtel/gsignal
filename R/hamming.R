# hamming.R
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
# 20191209 Geert van Boxtel          First version for v0.1.0
# 20231116 Geert van Boxtel          use coefficients 0.54 - 0.46 as in
#                                    Matlab/Octave
#------------------------------------------------------------------------------

#' Hamming window
#'
#' Return the filter coefficients of a Hamming window of length \code{n}.
#'
#' The Hamming window is a member of the family of cosine sum windows.
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
#' @return Hamming window, returned as a vector. If you specify a one-point
#'   window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' h <- hamming(64)
#' plot (h, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' hs = hamming(64,'symmetric')
#' hp = hamming(63,'periodic')
#' plot (hs, type = "l", xlab = "Samples", ylab =" Amplitude")
#' lines(hp, col="red")
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

hamming <- function(n, method = c("symmetric", "periodic")) {

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
    w <- 0.54 - 0.46 * cos(2 * pi * (0:n) / N)
  }
  w
}
