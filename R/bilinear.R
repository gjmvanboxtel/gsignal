# bilinear.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
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
# 20200501 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Bilinear transformation
#'
#' Transform a s-plane (analog) filter specification into a z-plane (digital)
#' specification.
#'
#' Given a piecewise flat filter design, you can transform it from the s-plane
#' to the z-plane while maintaining the band edges by means of the bilinear
#' transform. This maps the left hand side of the s-plane into the interior of
#' the unit circle. The mapping is highly non-linear, so you must design your
#' filter with band edges in the s-plane positioned at \eqn{2/T tan(wT / 2)} so
#' that they will be positioned at \code{w} after the bilinear transform is
#' complete.
#'
#' The bilinear transform is:
#' \deqn{z = (1 + sT / 2) / (1 - sT / 2)}
#' \deqn{s = (T / 2) (z - 1) / (z + 1)}
#'
#' Please note that a pole and a zero at the same place exactly cancel. This is
#' significant since the bilinear transform creates numerous extra poles and
#' zeros, most of which cancel. Those which do not cancel have a “fill-in”
#' effect, extending the shorter of the sets to have the same number of as the
#' longer of the sets of poles and zeros (or at least split the difference in
#' the case of the band pass filter). There may be other opportunistic
#' cancellations, but it will not check for them.
#'
#' Also note that any pole on the unit circle or beyond will result in an
#' unstable filter. Because of cancellation, this will only happen if the number
#' of poles is smaller than the number of zeros. The analytic design methods all
#' yield more poles than zeros, so this will not be a problem.
#'
#' @param Sz In the generic case, a model to be transformed. In the default
#'   case, a vector containing the zeros in a pole-zero-gain model.
#' @param Sp a vector containing the poles in a pole-zero-gain model.
#' @param Sg a vector containing the gain in a pole-zero-gain model.
#' @param T the sampling frequency represented in the z plane. Default:
#'  \code{2 * tan(1 / 2)}.
#' @param ...	arguments passed to the generic function.
#'
#' @return For the default case or for bilinear.Zpg, an object of class
#'   \code{'Zpg'}, containing the list elements:
#' \describe{
#'   \item{z}{complex vector of the zeros of the transformed model}
#'   \item{p}{complex vector of the poles of the transformed model}
#'   \item{g}{gain of the transformed model}
#' }
#' For bilinear.Arma, an object of class \code{'Arma'}, containing the list
#' elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#'
#' @examples
#' ## 6th order Bessel low-pass analog filter
#' zp <- besselap(6)
#' w <- seq(0, 4, length.out = 128)
#' freqs(zp, w)
#' zzp <- bilinear(zp)
#' freqz(zzp)
#'
#' @references \url{https://en.wikipedia.org/wiki/Bilinear_transform}
#'
#' @author Paul Kienzle \email{pkienzle@@users.sf.net}. Conversion to R by Tom
#'   Short, adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname bilinear
#' @export

bilinear <- function(Sz, ...) UseMethod("bilinear")

#' @rdname bilinear
#' @export

bilinear.Zpg <- function(Sz, T = 2 * tan(1 / 2), ...)
  bilinear(Sz$z, Sz$p, Sz$g, T)

#' @rdname bilinear
#' @export

bilinear.Arma <- function(Sz, T = 2 * tan(1 / 2), ...)
  as.Arma(bilinear(as.Zpg(Sz), T))

#' @rdname bilinear
#' @export

bilinear.default <- function(Sz, Sp, Sg, T = 2 * tan(1 / 2), ...)  {
  p <- length(Sp)
  z <- length(Sz)
  if (z > p || p == 0)
    stop("must have at least as many poles as zeros in s-plane")
  ## ----------------  -------------------------  ------------------------
  ## Bilinear          zero: (2+xT)/(2-xT)        pole: (2+xT)/(2-xT)
  ##      2 z-1        pole: -1                   zero: -1
  ## S -> - ---        gain: (2-xT)/T             gain: (2-xT)/T
  ##      T z+1
  ## ----------------  -------------------------  ------------------------
  Zg <- Re(Sg * prod((2 - Sz * T) / T) / prod((2 - Sp * T) / T))
  Zp <- (2 + Sp * T) / (2 - Sp * T)
  if (is.null(Sz))
    Zz <- rep.int(-1, length(Zp))
  else {
    Zz <- (2 + Sz * T) / (2 - Sz * T)
    Zz <- c(Zz, rep.int(-1, p - z))
  }
  Zpg(z = Zz, p = Zp, g = Zg)
}
