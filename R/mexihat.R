# mexihat.R
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
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#
# 20191126 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Mexicat Hat
#'
#' Generate a Mexican Hat (Ricker) wavelet sampled on a regular grid.
#'
#' The Mexican Hat or Ricker wavelet is the negative normalized second
#' derivative of a Gaussian function, i.e., up to scale and normalization, the
#' second Hermite function. It is a special case of the family of continuous
#' wavelets (wavelets used in a continuous wavelet transform) known as Hermitian
#' wavelets. The Ricker wavelet is frequently employed to model seismic data,
#' and as a broad spectrum source term in computational electrodynamics. It is
#' usually only referred to as the Mexican hat wavelet in the Americas, due to
#' taking the shape of a sombrero when used as a 2D image processing kernel. It
#' is also known as the Marr wavelet (source: Wikipedia)
#'
#' @param lb,ub Lower and upper bounds of the interval to evaluate the wavelet
#'   on. Default: -5 to 5.
#' @param n Number of points on the grid between \code{lb} and \code{ub} (length
#'   of the wavelet). Default: 1000.
#'
#' @return A list containing 2 variables; \code{x}, the grid on which the
#'   complex Mexican Hat wavelet was evaluated, and \code{psi} (\eqn{\Psi}), the
#'   evaluated wavelet on the grid \code{x}.
#'
#' @examples
#'
#' mh <- mexihat(-5, 5, 1000)
#' plot(mh$x, mh$psi, type="l", main = "Mexican Hat Wavelet",
#'      xlab = "", ylab = "")
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

mexihat <- function(lb = -5, ub = 5, n = 1000) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  x <- seq(lb, ub, length.out =  n)
  psi <- (1 - x^2) * (2 / (sqrt(3) * pi^0.25)) * exp(-x^2 / 2)
  list(x = x, psi = psi)

}
