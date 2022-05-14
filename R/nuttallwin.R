# nuttallwin.R
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
# 20191215 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Nuttall-defined minimum 4-term Blackman-Harris window
#'
#' Return the filter coefficients of a Blackman-Harris window defined by Nuttall
#' of length \code{n}.
#'
#' The window is minimum in the sense that its maximum sidelobes are minimized.
#' The coefficients for this window differ from the Blackman-Harris window
#' coefficients computed with \code{blackmanharris} and produce slightly lower
#' sidelobes.
#'
#' @param n Window length, specified as a positive integer.
#' @param method Character string. Window sampling method, specified as:
#' \describe{
#'   \item{"symmetric"}{(Default). Use this option when using windows for filter
#'   design.}
#'   \item{"periodic"}{This option is useful for spectral analysis because it
#'   enables a windowed signal to have the perfect periodic extension implicit
#'   in the discrete Fourier transform. When \code{periodic} is specified, the
#'   function computes a window of length \code{n + 1} and returns the first
#'   \code{n} points.}
#' }
#'
#' @return Nuttall-defined Blackman-Harris window, returned as a vector.
#'
#' @examples
#'
#' n <- nuttallwin(64)
#' plot (n, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @seealso \code{\link{blackman}}, \code{\link{blackmanharris}}

#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

nuttallwin <- function(n, method = c("symmetric", "periodic")) {

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
    a0 <- 0.355768
    a1 <- 0.487396
    a2 <- 0.144232
    a3 <- 0.012604
    k <- (-N / 2):((n - 1) / 2)
    w <- a0 + a1 * cos(2 * pi * k / N) + a2 *
      cos(4 * pi * k / N) + a3 * cos(6 * pi * k / N)
  }
  w
}
