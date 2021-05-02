# shanwavf.R
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
# 20191130 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Complex Shannon Wavelet
#'
#' Compute the Complex Shannon wavelet.
#'
#' The complex Shannon wavelet is defined by a bandwidth parameter \code{fb}, a
#' wavelet center frequency \code{fc}, and the expression
#' \deqn{\psi(x) = (fb^{0.5} * (sinc(fb * x) * e^{2 * 1i * pi * fc * x}))}
#' on an \code{n}-point regular grid in the interval of \code{lb} to \code{ub}.
#'
#' @param lb,ub Lower and upper bounds of the interval to evaluate the waveform
#'   on. Default: -8 to 8.
#' @param n Number of points on the grid between \code{lb} and \code{ub} (length
#'   of the wavelet). Default: 1000.
#' @param fb Time-decay parameter of the wavelet (bandwidth in the frequency
#'   domain). Must be a positive scalar. Default: 5.
#' @param fc Center frequency of the wavelet. Must be a positive scalar.
#'   Default: 1.
#'
#' @return A list containing 2 variables; \code{x}, the grid on which the
#'   complex Shannon wavelet was evaluated, and \code{psi} (\eqn{\Psi}), the
#'   evaluated wavelet on the grid \code{x}.
#'
#' @examples
#'
#' fb <- 1
#' fc <- 1.5
#' lb <- -20
#' ub <- 20
#' n <- 1000
#' sw <- shanwavf(lb, ub, n, fb, fc)
#' op <- par(mfrow = c(2,1))
#' plot(sw$x, Re(sw$psi), type="l", main = "Complex Shannon Wavelet",
#'      xlab = "real part", ylab = "")
#' plot(sw$x, Im(sw$psi), type="l", xlab = "imaginary part", ylab = "")
#' par(op)
#'
#' @author Sylvain Pelissier, \email{sylvain.pelissier@@gmail.com}.\cr
#' Conversion to R by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

shanwavf <- function(lb = -8, ub = 8, n = 1000, fb = 5, fc = 1) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  if (!isPosscal(fb) || fb <= 0)
    stop("fb must be a positive scalar > 0")
  if (!isPosscal(fc) || fc <= 0)
    stop("fc must be a positive scalar > 0")

  x <- seq(lb, ub, length.out = n)
  psi <- (fb^0.5) * (sinc(fb * x) * exp(2 * 1i * pi * fc * x))

  list(x = x, psi = psi)
}
