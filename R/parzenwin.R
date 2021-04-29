# parzenwin.R
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

#' Parzen (de la Vallée Poussin) window
#'
#' Return the filter coefficients of a Parzen window of length \code{n}.
#'
#' Parzen windows are piecewise-cubic approximations of Gaussian windows.
#'
#' @param n Window length, specified as a positive integer.
#'
#' @return Parzen window, returned as a vector.
#'
#' @examples
#'
#' p <- parzenwin(64)
#' g <- gausswin(64)
#' plot (p, type = "l", xlab = "Samples", ylab =" Amplitude", ylim = c(0, 1))
#' lines(g, col = "red")
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

parzenwin <- function(n) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")

  N <- n - 1
  k <- (- (N / 2)):(N / 2)
  k1 <- k[which(abs(k) <= (N / 4))]
  k2 <- k[which(k > (N / 4))]
  k3 <- k[which(k < (-N / 4))]

  w1 <- 1 - 6 * (abs(k1) / (n / 2))^2 + 6 * (abs(k1) / (n / 2))^3
  w2 <- 2 * (1 - abs(k2) / (n / 2))^3
  w3 <- 2 * (1 - abs(k3) / (n / 2))^3
  w <- c(w3, w1, w2)
  w
}
