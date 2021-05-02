# kaiser.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1995, 1996, 1997 Kurt Hornik <Kurt.Hornik@ci.tuwien.ac.at>
# Copyright (C) 2000 Paul Kienzle <pkienzle@users.sf.net>
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
# 20191215 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Kaiser window
#'
#' Return the filter coefficients of a kaiser window of length \code{n}.
#'
#' The Kaiser, or Kaiser-Bessel, window is a simple approximation of the DPSS
#' window using Bessel functions, discovered by James Kaiser.
#' \if{latex}{
#'   \deqn{w(x) = \frac{besselI(0, \beta \cdot \sqrt{(1 -
#'   (2*x/m)^{2}))}}{besselI(0, \beta)}; -m/2 <= x <= m/2}
#' }
#' \if{html}{\preformatted{
#'         besselI(0, Beta * sqrt(1-(2*x/m)^2))
#' k(x) =  -------------------------------------,  -m/2 <= x <= m/2
#'         besselO(0, Beta)
#' }}
#' The variable parameter \eqn{\beta} determines the trade-off between main lobe
#' width and side lobe levels of the spectral leakage pattern. Increasing
#' \eqn{\beta} widens the main lobe and decreases the amplitude of the side
#' lobes (i.e., increases the attenuation).
#'
#' @param n Window length, specified as a positive integer.
#' @param beta Shape factor, specified as a positive real scalar. The parameter
#'   \code{beta} affects the side lobe attenuation of the Fourier transform of
#'   the window. Default: 0.5
#'
#' @return Kaiser window, returned as a vector.
#'
#' @examples
#'
#' k <- kaiser(200, 2.5)
#' plot (k, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#' @author Kurt Hornik, \email{Kurt.Hornik@@ci.tuwien.ac.at},\cr Paul Kienzle,
#'   \email{pkienzle@@users.sf.net}.\cr Conversion to R by Geert van Boxtel
#'   \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

kaiser <- function(n, beta = 0.5) {

  if (!isPosscal(n) || !isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  if (!isScalar(beta) || !is.double(beta))
    stop("beta must be a real scalar")

  if (n == 1) {
    w <- 1
  } else {
    N <- n - 1
    k <- (0:N)
    k <- 2 * beta / N * sqrt(k * (N - k))
    w <- besselI(k, 0) / besselI(beta, 0)
  }
  w
}
