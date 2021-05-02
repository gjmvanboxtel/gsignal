# chebwin.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2002 André Carezia <acarezia@uol.com.br>
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
# 20191211 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Chebyshev window
#'
#' Return the filter coefficients of a Dolph-Chebyshev window.
#'
#' The window is described in frequency domain by the expression:
#' \if{latex}{
#'   \deqn{W(k) = \frac{Cheb(m - 1, \beta \cdot cos(\pi \cdot k / m))}{Cheb(m -
#'   1, \beta)}}
#' }
#' \if{html}{\preformatted{
#'                 Cheb(m - 1, Beta * cos(\pi * k / m))
#'          W(k) = ------------------------------------
#'                        Cheb(m - 1, Beta)
#' }}
#' with
#' \if{latex}{
#'   \deqn{\beta = cosh(1 / (m - 1) \cdot acosh(10^{(at / 20)})}
#' }
#' \if{html}{\preformatted{
#'   Beta = cosh(1 / (m - 1) * acosh(10^(at / 20))
#' }}
#' and and \eqn{Cheb(m, x)} denoting the \eqn{m}-th order Chebyshev polynomial
#' calculated at the point \eqn{x}.
#'
#' Note that the denominator in W(k) above is not computed, and after the
#' inverse Fourier transform the window is scaled by making its maximum value
#' unitary.
#'
#' @param n Window length, specified as a positive integer.
#' @param at Stop-band attenuation in dB. Default: 100.
#'
#' @return Chebyshev window, returned as a vector. If you specify a one-point
#'   window \code{(n = 1)}, the value 1 is returned.
#'
#' @examples
#'
#' cw <- chebwin(64)
#' plot (cw, type = "l", xlab = "Samples", ylab =" Amplitude")
#'
#'
#' @author André Carezia, \email{acarezia@@uol.com.br}.\cr
#' Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

chebwin <- function(n, at = 100) {

  if (!isPosscal(n) || ! isWhole(n) || n <= 0)
    stop("n must be an integer strictly positive")
  if (!isScalar(at) || ! is.double(at) || n <= 0)
    stop("at must be a real scalar")

  if (n == 1) {
    w <- 1
  } else {
    ## beta calculation
    gamma <- 10 ^ (-at / 20)
    beta <- cosh(1 / (n - 1) * acosh(1 / gamma))

    ## freq. scale
    k <- 0:(n - 1)
    x <- beta * cos(pi * k / n)

    ## Chebyshev window (freq. domain)
    p <- cheb(n - 1, x)

    ## inverse Fourier transform
    if (n %% 2) {
      w <- Re(stats::fft(p))
      M <- (n + 1) / 2
      w <- w[1:M] / w[1]
      w <- c(w[M:2], w)
    } else {
      ## half-sample delay (even order)
      p <- p * exp(1i * pi / n * (0:(n - 1)))
      w <- Re(stats::fft(p))
      M <- n / 2 + 1
      w <- w / w[2]
      w <- c(w[M:2], w[2:M])
    }
  }

  w <- w / max(w)
  w
}
