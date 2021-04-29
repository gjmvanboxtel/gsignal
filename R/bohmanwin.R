# bohmanwin.R
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
# 20191209 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Bohman window
#'
#' Return the filter coefficients of a Bohman window.
#'
#' A Bohman window is the convolution of two half-duration cosine lobes. In the
#' time domain, it is the product of a triangular window and a single cycle of a
#' cosine with a term added to set the first derivative to zero at the boundary.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return Bohman window, returned as a vector. If you specify a one-point
#'   window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' b <- bohmanwin(64)
#' plot (b, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @seealso \code{\link{triang}}
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

bohmanwin <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0) {
    stop("n must be an integer strictly positive")
  }

  if (n == 1) {
    w <- 1
  } else {
    N <- n - 1
    k <- (-N / 2):(N / 2)
    w <- (1 - 2 * abs(k) / N) *
      cos(2 * pi * abs(k) / N) + (1 / pi) *
      sin(2 * pi * abs(k) / N)
    w[1] <- w[length(w)] <- 0
  }
  w
}
