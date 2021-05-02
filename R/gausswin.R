# gausswin.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
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

#' Gaussian window
#'
#' Return the filter coefficients of a Gaussian window of length \code{n}.
#'
#' The width of the window is inversely proportional to the parameter \code{a}.
#' Use larger \code{a} for a narrower window. Use larger \code{m} for a smoother
#' curve.
#' \deqn{w = e^{(-(a*x)^{2}/2 )}}
#' for \code{x <- seq(-(n - 1) / n, (n - 1) / n, by = n)}.
#'
#' The exact correspondence with the standard deviation of a Gaussian
#' probability density function is \eqn{\sigma = (n - 1) / (2a)}.
#'
#' @param n Window length, specified as a positive integer.
#' @param a Width factor, specified as a positive real scalar. \code{a} is
#'   inversely proportional to the width of the window. Default: 2.5.
#'
#' @return Gaussian convolution window, returned as a vector.
#'
#' @examples
#'
#' g1 <- gausswin(64)
#' g2 <- gausswin(64, 5)
#' plot (g1, type = "l", xlab = "Samples", ylab =" Amplitude", ylim = c(0, 1))
#' lines(g2, col = "red")
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

gausswin <- function(n, a = 2.5) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  if (!isScalar(a))
    stop("a must be a scalar")

  if (n == 1) {
    w <- 1
  } else {
    w <- exp(-0.5 * (a / (n - 1) * seq(- (n - 1), (n - 1), 2))^2)
  }
  w
}
