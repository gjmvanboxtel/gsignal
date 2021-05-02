# bartlett.R
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
#------------------------------------------------------------------------------

#' Bartlett window
#'
#' Return the filter coefficients of a Bartlett (triangular) window.
#'
#' The Bartlett window is very similar to a triangular window as returned by the
#' \code{\link{triang}} function. However, the Bartlett window always has zeros
#' at the first and last samples, while the triangular window is nonzero at
#' those points.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return Bartlett window, returned as a vector. If you specify a one-point
#'   window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' bw <- bartlett(64)
#' plot (bw, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @seealso \code{\link{triang}}
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

bartlett <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0) {
    stop("n must be an integer strictly positive")
  }

  if (n == 1) {
    w <- 1
  } else {
    n <- n - 1
    m <- trunc(n / 2)
    w <- c(2 * (0:m) / n, 2 - 2 * ((m + 1):n) / n)
  }
  w
}
