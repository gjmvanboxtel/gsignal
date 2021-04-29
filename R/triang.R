# triang.R
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

#' Triangular window
#'
#' Return the filter coefficients of a triangular window of length \code{n}.
#'
#' Unlike the Bartlett window, \code{triang} does not go to zero at the edges of
#' the window. For odd \code{n}, \code{triang(n)} is equal to \code{bartlett(m +
#' 2)} except for the zeros at the edges of the window.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return triangular window, returned as a vector. If you specify a one-point
#'   window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' t <- triang(64)
#' plot (t, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @seealso \code{\link{bartlett}}
#'
#' @author Andreas Weingessel, \email{Andreas.Weingessel@@ci.tuwien.ac.at}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

triang <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  w <- 1 - abs(seq(- (n - 1), (n - 1), by = 2) / (n + n %% 2))
  w
}
