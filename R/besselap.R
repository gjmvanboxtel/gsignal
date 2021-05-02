# besselap.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2009 Thomas Sailer
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
# 20200426 Geert van Boxtel          First version for v0.1.0
# 20205001 Geert van Boxtel          return Zpg$z = complex(0) instead of NULL
#------------------------------------------------------------------------------

#' Bessel analog low-pass filter prototype
#'
#' Return the poles and gain of a Bessel analog low-pass filter prototype.
#'
#' The transfer function is
#' \if{latex}{
#'   \deqn{H(s) = \frac{k}{(s-p(1))(s-p(2))...(s-p(n))}}
#' }
#' \if{html}{\preformatted{
#'                      k
#'  H(s) = -----------------------------
#'          (s-p(1))(s-p(2))...(s-p(n))
#'
#' }}
#' \code{besselap} normalizes the poles and gain so that at low frequency and
#' high frequency the Bessel prototype is asymptotically equivalent to the
#' Butterworth prototype of the same order. The magnitude of the filter is less
#' than \eqn{1/\sqrt{2}} at the unity cutoff frequency \eqn{\Omega_c = 1}.
#'
#' Analog Bessel filters are characterized by a group delay that is maximally
#' flat at zero frequency and almost constant throughout the passband. The group
#' delay at zero frequency is
#' \if{latex}{
#'   \deqn{\left( \frac{(2n)!}{2^{n}n!} \right) ^{1/n}}
#' }
#' \if{html}{\preformatted{
#'    /  (2n!) \ 2
#'    | ------ |
#'    \ 2^n n! /
#' }}
#'
#' @param n order of the filter; must be < 25.
#'
#' @return List of class \code{\link{Zpg}} containing poles and gain of the
#'   filter
#'
#' @examples
#' ## 6th order Bessel low-pass analog filter
#' zp <- besselap(6)
#' w <- seq(0, 4, length.out = 128)
#' freqs(zp, w)
#'
#' @references \url{https://en.wikipedia.org/wiki/Bessel_polynomials}
#'
#' @author Thomas Sailer, email{t.sailer@@alumni.ethz.ch}.\cr Conversion to R by
#'   Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

besselap <- function(n) {

  if (!isPosscal(n) || ! isWhole(n)) {
    stop("n must be an integer strictly positive")
  }

  if (n == 1) {
    p <- -1
  } else {
    p0 <- 1
    p1 <- rep(1L, 2)
    for (nn in 2:n) {
      px <- (2 * nn - 1) * p1
      py <- c(p0, 0, 0)
      px <- prepad(px, max(length(px), length(py)), 0)
      py <- prepad(py, length(px))
      p0 <- p1
      p1 <- px + py
    }
    ## p1 now contains the reverse bessel polynomial for n

    ## scale it by replacing s->s/w0 so that the gain becomes 1
    p1 <- p1 * p1[length(p1)] ^ (seq(length(p1) - 1, 0, -1) / (length(p1) - 1))
    p <- pracma::roots(p1)
  }
  Zpg(complex(0), p, 1)
}
