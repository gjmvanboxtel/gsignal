# meyeraux.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Matlab/Octave signal package:
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
# 20191123 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Meyer wavelet auxiliary function
#'
#' Compute the Meyer wavelet auxiliary function.
#'
#' The code \code{y = meyeraux(x)} returns values of the auxiliary function used
#' for Meyer wavelet generation evaluated at the elements of \code{x}. The input
#' \code{x} is a vector or matrix of real values. The function is \deqn{y =
#' 35x^{4} - 84x^{5} + 70x^{6} - 20x^{7}.} \code{x} and \code{y} have the same
#' dimensions. The range of \code{meyeraux} is the closed interval c(0, 1).
#'
#' @param x Input array, specified as a real scalar, vector, matrix, or
#'   multidimensional array.
#'
#' @return Output array, returned as a real-valued scalar, vector, matrix, or
#'   multidimensional array of the same size as x.
#'
#' @examples
#'
#' x <- seq(0, 1, length.out = 100)
#' y <- meyeraux(x)
#' plot(x, y, type="l", main = "Meyer wavelet auxiliary function",
#'      xlab = "", ylab = "")
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

meyeraux <- function(x) {

  y <- 35 * x^4 - 84 * x^5 + 70 * x^6 - 20 * x^7
  y
}
