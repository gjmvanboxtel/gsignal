# rectwin.R
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

#' Rectangular window
#'
#' Return the filter coefficients of a rectangular window of length \code{n}.
#'
#' The output of the rectwin function with input \code{n} can also be created
#' using the \code{rep} function: w <- rep(1L, n)
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return rectangular window, returned as a vector.
#'
#' @examples
#'
#' r <- rectwin(64)
#' plot (r, type = "l", xlab = "Samples", ylab =" Amplitude", ylim = c(0, 1))
#'
#' @seealso \code{\link{boxcar}}
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

rectwin <- function(n) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  w <- rep(1L, n)
  w
}
