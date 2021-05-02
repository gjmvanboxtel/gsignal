# cheb1ap.R
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

#' Chebyshev Type I filter prototype
#'
#' Return the poles and gain of an analog Chebyshev Type I lowpass filter
#' prototype.
#'
#' This function exists for compatibility with 'Matlab' and 'OCtave' only, and
#' is equivalent to \code{cheby1(n, Rp, 1, "low", "s")}.
#'
#' @param n Order of the filter.
#' @param Rp dB of pass-band ripple.
#'
#' @return List of class \code{\link{Zpg}} containing the poles and gain of the
#'   filter.
#'
#' @examples
#' ## 9th order Chebyshev type I low-pass analog filter
#' zp <- cheb1ap(9, .1)
#' w <- seq(0, 4, length.out = 128)
#' freqs(zp, w)
#'
#' @author Carne Draug, \email{carandraug+dev@@gmail.com}.\cr
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

cheb1ap <- function(n, Rp) {

  if (!isPosscal(n) || ! isWhole(n))
    stop("n must be an integer strictly positive")
  if (!isPosscal(Rp) || !is.numeric(Rp)) {
    stop("passband ripple Rp must a non-negative scalar")
  }

  cheby1(n, Rp, 1, "low", "s", "Zpg")

}
