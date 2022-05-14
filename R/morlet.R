# morlet.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
# Copyright (C) 2007 Sylvain Pelissier
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
# 20191126 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Morlet Wavelet
#'
#' Compute the Morlet wavelet on a regular grid.
#'
#' The code \code{m <- morlet(lb, ub, n)} returns values of the Morlet wavelet
#' on an \code{n}-point regular grid in the interval \code{c(lb, ub)}.
#'
#' The Morlet waveform is defined as
#' \deqn{\psi(x) = e^{-x^{2}/2} cos (5x)}
#'
#' @param lb,ub Lower and upper bounds of the interval to evaluate the wavelet
#'   on. Default: -4 to 4.
#' @param n Number of points on the grid between \code{lb} and \code{ub} (length
#'   of the wavelet). Default: 1000.
#'
#' @return A list containing 2 variables; \code{x}, the grid on which the Morlet
#'   wavelet was evaluated, and \code{psi} (\eqn{\Psi}), the evaluated wavelet
#'   on the grid \code{x}.
##'
#' @examples
#'
#' m <- morlet(-4, 4, 1000)
#' plot(m$x, m$psi, type="l", main = "Morlet Wavelet", xlab = "", ylab = "")
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#'   Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

morlet <- function(lb = -4, ub = 4, n = 1000) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  x <- seq(lb, ub, length.out =  n)
  psi <- cos(5 * x) * exp(-x^2 / 2)
  list(x = x, psi = psi)

}
