# boxcar.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
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

#' Rectangular window
#'
#' Return the filter coefficients of a boxcar (rectangular) window.
#'
#' The rectangular window (sometimes known as the boxcar or Dirichlet window) is
#' the simplest window, equivalent to replacing all but \code{n} values of a
#' data sequence by zeros, making it appear as though the waveform suddenly
#' turns on and off. Other windows are designed to moderate these sudden
#' changes, which reduces scalloping loss and improves dynamic range.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return rectangular window, returned as a vector.
#'
#' @examples
#'
#' b <- boxcar(64)
#' plot (b, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @seealso \code{\link{triang}}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

boxcar <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0) {
    stop("n must be an integer strictly positive")
  }

  w <- rep(1L, n)
  w
}
