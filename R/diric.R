# diric.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
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
# 20191123 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Dirichlet function
#'
#' Compute the Dirichlet or periodic sinc function.
#'
#' \code{y <- diric(x, n)} returns the Dirichlet Function of degree \code{n}
#' evaluated at the elements of the input array \code{x}.
#'
#' The Dirichlet function, or periodic sinc function, has period \eqn{2 \pi} for
#' odd \eqn{N} and period \eqn{4 \pi} for even \eqn{N}. Its maximum value is 1
#' for all N, and its minimum value is -1 for even N. The magnitude of the
#' function is 1 / N times the magnitude of the discrete-time Fourier transform
#' of the N-point rectangular window.
#'
#' @param x Input array, specified as a real scalar, vector, matrix, or
#'   multidimensional array. When \code{x} is non-scalar, \code{diric} is an
#'   element-wise operation.
#' @param n Function degree, specified as a positive integer scalar.
#'
#' @return Output array, returned as a real-valued scalar, vector, matrix, or
#'   multidimensional array of the same size as x.
#'
#' @examples
#'
#' ## Compute and plot the Dirichlet function between -2pi and 2pi for N = 7
#' ## and N = 8. The function has a period of 2pi for odd N and 4pi for even N.
#' x <- seq(-2*pi, 2*pi, len = 301)
#' d7 <- diric(x, 7)
#' d8 <- diric(x, 8)
#' op <- par(mfrow = c(2,1))
#' plot(x/pi, d7, type="l", main = "Dirichlet function",
#'      xlab = "", ylab = "N = 7")
#' plot(x/pi, d8, type="l", ylab = "N = 8", xlab = expression(x / pi))
#' par(op)
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

diric <- function(x, n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  y <- sin(n * x / 2) / (n * sin(x / 2))
  y[x %% (2 * pi) == 0] <- (-1) ^ ((n - 1) * x[x %% (2 * pi) == 0] / (2 * pi))
  y
}
