# buttap.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2013 Carne Draug <carandraug+dev@gmail.com>
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
# 20200519 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Chebyshev Type II filter prototype
#'
#' Return the poles and gain of an analog Chebyshev Type II lowpass filter
#' prototype.
#'
#' This function exists for compatibility with 'Matlab' and 'Octave' only, and
#' is equivalent to \code{cheby2(n, Rp, 1, "low", "s")}.
#'
#' @param n Order of the filter.
#' @param Rs dB of stop-band ripple.
#'
#' @return list of class \code{\link{Zpg}} containing poles and gain of the
#'   filter
#'
#' @examples
#' ## 9th order Chebyshev type II low-pass analog filter
#' zp <- cheb2ap(9, 30)
#' w <- seq(0, 4, length.out = 128)
#' freqs(zp, w)
#'
#' @author Carne Draug, \email{carandraug+dev@@gmail.com}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cheb2ap <- function(n, Rs) {

  if (!isPosscal(n) || ! isWhole(n))
    stop("n must be an integer strictly positive")
  if (!isPosscal(Rs) || !is.numeric(Rs)) {
    stop("passband ripple Rp must a non-negative scalar")
  }

  cheby2(n, Rs, 1, "low", "s", "Zpg")

}
