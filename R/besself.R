# besself.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
# Copyright (C) 2009 Thomas Sailer <t.sailer@alumni.ethz.ch>
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
# 20200428 Geert van Boxtel          First version for v0.1.0
#------------------------------------------------------------------------------

#' Bessel analog filter design
#'
#' Compute the transfer function coefficients of an analog Bessel filter.
#'
#' Bessel filters are characterized by an almost constant group delay across the
#' entire passband, thus preserving the wave shape of filtered signals in the
#' passband.
#'
#' Lowpass Bessel filters have a monotonically decreasing magnitude response, as
#' do lowpass Butterworth filters. Compared to the Butterworth, Chebyshev, and
#' elliptic filters, the Bessel filter has the slowest rolloff and requires the
#' highest order to meet an attenuation specification.
#'
#' @note As the important characteristic of a Bessel filter is its
#'   maximally-flat group delay, and not the amplitude response, it is
#'   inappropriate to use the bilinear transform to convert the analog Bessel
#'   filter into a digital form (since this preserves the amplitude response but
#'   not the group delay) [1].
#'
#' @param n filter order.
#' @param w critical frequencies of the filter. \code{w} must be a scalar for
#'   low-pass and high-pass filters, and \code{w} must be a two-element vector
#'   c(low, high) specifying the lower and upper bands in radians/second.
#' @param type filter type, one of \code{"low"} (default), \code{"high"},
#'   \code{"stop"}, or \code{"pass"}.
#'
#' @return List of class \code{'\link{Zpg}'} containing poles and gain of the
#'   filter.
#'
#' @examples
#' w <- seq(0, 4, length.out = 128)
#'
#' ## 5th order Bessel low-pass analog filter
#' zp <- besself(5, 1.0)
#' freqs(zp, w)
#'
#' ## 5th order Bessel high-pass analog filter
#' zp <- besself(5, 1.0, 'high')
#' freqs(zp, w)
#'
#' ## 5th order Bessel band-pass analog filter
#' zp <- besself(5, c(1, 2), 'pass')
#' freqs(zp, w)
#'
#' ## 5th order Bessel band-stop analog filter
#' zp <- besself(5, c(1, 2), 'stop')
#' freqs(zp, w)
#'
#' @references [1] \url{https://en.wikipedia.org/wiki/Bessel_filter}
#'
#' @author Paul Kienzle, \email{pkienzle@@users.sf.net},\cr
#'   Doug Stewart, \email{dastew@@sympatico.ca},\cr
#'   Thomas Sailer, \email{t.sailer@@alumni.ethz.ch}.\cr
#'   Conversion to R by Geert van Boxtel, \email{G.J.M.vanBoxtel@@gmail.com}.
#
#' @export

besself <- function(n, w, type = c("low", "high", "stop", "pass")) {

  if (!isPosscal(n)) {
    stop("filter order n must be a positive integer")
  }

  type <- match.arg(type)
  stop <- type == "stop" || type == "high"

  if (length(w) != 1 && length(w) != 2) {
    stop("frequency must be given as w0 or c(w0, w1)")
  }
  if (!all(w >= 0)) {
    stop("critical frequencies must be in [0 Inf]")
  }
  if ((length(w) == 2) && (w[2] <= w[1])) {
    stop("w[1] must be less than w[2]")
  }

  # Generate splane poles for the prototype Bessel filter
  zpg <- besselap(n)

  # splane frequency transform
  zpg <- sftrans(zpg, w, stop)

  zpg
}
