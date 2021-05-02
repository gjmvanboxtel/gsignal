# ellipap.R
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
# 20200527 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Low-pass analog elliptic filter
#'
#' Return the zeros, poles and gain of an analog elliptic low-pass filter
#' prototype.
#'
#' This function exists for compatibility with 'Matlab' and 'OCtave' only, and
#' is equivalent to \code{ellip(n, Rp, Rs, 1, "low", "s")}.
#'
#' @param n Order of the filter.
#' @param Rp dB of passband ripple.
#' @param Rs dB of stopband ripple.
#'
#' @return list of class \code{\link{Zpg}} containing zeros, poles and gain of
#'   the filter.
#'
#' @examples
#' ## 9th order elliptic low-pass analog filter
#' zp <- ellipap(9, .1, 40)
#' w <- seq(0, 4, length.out = 128)
#' freqs(zp, w)
#'
#' @author Carne Draug, \email{carandraug+dev@@gmail.com}.
#'  Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

ellipap <- function(n, Rp, Rs) {

  if (!isPosscal(n) || ! isWhole(n))
    stop("n must be an integer strictly positive")
  if (!isPosscal(Rp) || !is.numeric(Rp)) {
    stop("passband ripple Rp must a non-negative scalar")
  }
  if (!isPosscal(Rs) || !is.numeric(Rs)) {
    stop("stopband ripple Rs must a non-negative scalar")
  }

  ellip(n, Rp, Rs, 1, "low", "s", "Zpg")

}
