# ncauer.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2001 Paulo Neis <p_neis@yahoo.com.br>
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

#' ncauer analog filter design
#'
#' Compute the transfer function coefficients of a Cauer analog filter.
#'
#' Cauer filters have equal maximum ripple in the passband and the stopband. The
#' Cauer filter has a faster transition from the passband to the stopband than
#' any other class of network synthesis filter. The term Cauer filter can be
#' used interchangeably with elliptical filter, but the general case of
#' elliptical filters can have unequal ripples in the passband and stopband. An
#' elliptical filter in the limit of zero ripple in the passband is identical to
#' a Chebyshev Type 2 filter. An elliptical filter in the limit of zero ripple
#' in the stopband is identical to a Chebyshev Type 1 filter. An elliptical
#' filter in the limit of zero ripple in both passbands is identical to a
#' Butterworth filter. The filter is named after Wilhelm Cauer and the transfer
#' function is based on elliptic rational functions.Cauer-type filters use
#' generalized continued fractions.[1]
#'
#' @param Rp dB of passband ripple.
#' @param Rs dB of stopband ripple.
#' @param n filter order.
#'
#' @return A list of class Zpg with the following list elements:
#' \describe{
#'   \item{zero}{complex vector of the zeros of the model}
#'   \item{pole}{complex vector of the poles of the model}
#'   \item{gain}{gain of the model}
#' }
#'
#' @examples
#' zpg <- ncauer(1, 40, 5)
#' freqz(zpg)
#' zplane(zpg)
#'
#' @references [1]
#'  \url{https://en.wikipedia.org/wiki/Network_synthesis_filters#Cauer_filter}
#'
#' @seealso \code{\link{Zpg}}, \code{\link{filter}}, \code{\link{ellip}}
#'
#' @author Paulo Neis, \email{p_neis@@yahoo.com.br}.\cr
#' Conversion to R Tom Short,\cr
#'  adapted by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @export

ncauer <- function(Rp, Rs, n) {

  # check input arguments
  if (!isPosscal(n) || !isWhole(n)) {
    stop("filter order n must be a positive integer")
  }
  if (!isPosscal(Rp) || !is.numeric(Rp)) {
    stop("passband ripple Rp must a non-negative scalar")
  }
  if (!isPosscal(Rs) || !is.numeric(Rs)) {
    stop("stopband ripple Rp must a non-negative scalar")
  }

  ## Calculate the stop band edge for the Cauer filter.
  ellip_ws <- function(n, rp, rs) {
    ellip_ws_min <- function(kl) {
      int <- pracma::ellipke(c(kl, 1 - kl))$k
      ql <- int[1]
      q <- int[2]
      abs((ql / q) - x)
    }
    kl0 <- ((10 ^ (0.1 * rp) - 1) / (10 ^ (0.1 * rs) - 1))
    k0 <- 1 - kl0
    int <- pracma::ellipke(c(kl0, k0))$k
    ql0 <- int[1]
    q0 <- int[2]
    x <- n * ql0 / q0
    kl <- stats::optimize(ellip_ws_min,
                          interval = c(.Machine$double.eps,
                                       1 - .Machine$double.eps))$minimum
    ws <- sqrt(1 / kl)
    ws
  }

  ## Cutoff frequency = 1:
  wp <- 1

  ## Stop band edge ws:
  ws <- ellip_ws(n, Rp, Rs)

  k  <- wp / ws
  k1 <- sqrt(1 - k^2)
  q0 <- (1 / 2) * ((1 - sqrt(k1)) / (1 + sqrt(k1)))
  q <-  q0 + 2 * q0^5 + 15 * q0^9 + 150 * q0^13

  l <- (1 / (2 * n)) * log((10 ^ (0.05 * Rp) + 1) / (10 ^ (0.05 * Rp) - 1))
  sig01 <- 0
  sig02 <- 0
  for (m in 0:30) {
    sig01 <- sig01 + (-1)^m * q ^ (m * (m + 1)) * sinh((2 * m + 1) * l)
  }
  for (m in 1:30) {
    sig02 <- sig02 + (-1)^m * q ^ (m^2) * cosh(2 * m * l)
  }
  sig0 <- abs((2 * q ^ (1 / 4) * sig01) / (1 + 2 * sig02))

  w <- sqrt((1 + k * sig0^2) * (1 + sig0^2 / k))

  r <- (n - (n %% 2)) / 2

  wi <- matrix(0, 1, r)
  for (ii in 1:r) {
    mu <- ii - (1 - (n %% 2)) / 2
    soma1 <- 0
    for (m in 0:30) {
      soma1 <- soma1 + 2 * q ^ (1 / 4) * ((-1)^m * q ^ (m * (m + 1)) *
                                            sin(((2 * m + 1) * pi * mu) / n))
    }
    soma2 <- 0
    for (m in 1:30) {
      soma2 <- soma2 + 2 * ((-1)^m * q ^ (m^2) * cos((2 * m * pi * mu) / n))
    }
    wi[ii] <- soma1 / (1 + soma2)
  }

  Vi <- sqrt((1 - (k * (wi^2))) * (1 - (wi^2) / k))
  A0i <- 1 / wi^2
  sqrA0i <- 1 / wi
  B0i <- ((sig0 * Vi)^2 + (w * wi)^2) / ((1 + sig0^2 * wi^2)^2)
  ## not used: ## B1i <- (2 * sig0 * Vi) / (1 + sig0^2 * wi^2)

  ##Gain T0:
  if (n %% 2) { # odd
    T0 <- sig0 * prod(B0i / A0i) * sqrt(ws)
  } else {
    T0 <- 10 ^ (-0.05 * Rp) * prod(B0i / A0i)
  }
  ##zeros:
  zer <- c(1i * sqrA0i, -1i * sqrA0i)

  ##poles:
  pol <- c((-2 * sig0 * Vi + 2 * 1i * wi * w) / (2 * (1 + sig0^2 * wi^2)),
           (-2 * sig0 * Vi - 2 * 1i * wi * w) / (2 * (1 + sig0^2 * wi^2)))

  ##If n odd, there is a real pole  -sig0:
  if (n %% 2) {   # odd
    pol <- c(pol, -sig0)
  }

  pole <- sqrt(ws) * pol
  zero <- sqrt(ws) * zer

  Zpg(z = zero, p = pole, g = T0)
}
