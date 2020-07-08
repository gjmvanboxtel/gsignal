# ellip.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 2001 Paulo Neis <p_neis@yahoo.com.br>
# Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
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
# 20200708 GvB                       renamed IIRfspec to FilterSpecs
#---------------------------------------------------------------------------------------------------------------------------------

#' Elliptic filter design
#' 
#' Compute the transfer function coefficients of an elliptic filter
#' 
#' An elliptic filter is a filter with equalized ripple (equiripple) behavior in
#' both the passband and the stopband. The amount of ripple in each band is
#' independently adjustable, and no other filter of equal order can have a
#' faster transition in gain between the passband and the stopband, for the
#' given values of ripple.
#' 
#' As the ripple in the stopband approaches zero, the filter becomes a type I
#' Chebyshev filter. As the ripple in the passband approaches zero, the filter
#' becomes a type II Chebyshev filter and finally, as both ripple values
#' approach zero, the filter becomes a Butterworth filter.
#' 
#' Because \code{ellip} is generic, it can be extended to accept other inputs, using
#' \code{ellipord} to generate filter criteria for example.
#' 
#' @param n filter order.
#' @param Rp dB of passband ripple.
#' @param Rs dB of stopband ripple.
#' @param w critical frequencies of the filter. \code{w} must be a scalar for
#'   low-pass and high-pass filters, and \code{w} must be a two-element vector
#'   \code{c(low, high)} specifying the lower and upper bands in radians/second.
#'   For digital filters, w must be between 0 and 1 where 1 is the Nyquist
#'   frequency.
#' @param type filter type, one of "low", "high", "stop", or "pass".
#' @param plane "z" for a digital filter or "s" for an analog filter.
#' @param ... additional arguments passed to ellip, overriding those given by n
#'   of class \code{FilterSpecs}.
#' 
#' @return list of class \code{'\link{Arma}'} with list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#' 
#' @examples
#' # compare the frequency responses of 5th-order Butterworth and elliptic filters.
#' bf <- butter(5, 0.1)
#' ef <- ellip(5, 3, 40, 0.1)
#' bfr <- freqz(bf)
#' efr <- freqz(ef)
#' plot(bfr$w, 20 * log10(abs(bfr$h)), type = "l", ylim = c(-80, 0),
#'      xlab = "Frequency (Rad)", ylab = c("dB"))
#' lines(efr$w, 20 * log10(abs(efr$h)), col = "red")
#' 
#' @references \url{https://en.wikipedia.org/wiki/Elliptic_filter}
#' 
#' @seealso \code{\link{Arma}}, \code{\link{filter}}, \code{\link{butter}}, \code{\link{cheby1}}, \code{\link{ellipord}}
#' 
#' @author Original Octave code by Paulo Neis \email{p_neis@@yahoo.com.br},
#'   adapted by Doug Stewart \email{dastew@@sympatico.ca}. Port to R Tom Short,
#'   adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname ellip
#' @export

ellip <- function(n, ...) UseMethod("ellip")

#' @rdname ellip
#' @export

ellip.FilterSpecs <- function(n, Rp = n$Rp, Rs = n$Rs, w = n$Wc, type = n$type, plane = n$plane, ...)
  ellip(n$n, Rp, Rs, w, type, plane, ...)

#' @rdname ellip
#' @export

ellip.default <- function (n, Rp, Rs, w, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {
  
  # check input arguments
  type <- match.arg(type)
  plane <- match.arg(plane)
  if (!isPosscal(n) || !isWhole(n)) {
    stop("filter order n must be a positive integer")
  }
  if (!isPosscal(Rp) || !is.numeric(Rp)) {
    stop("passband ripple Rp must a non-negative scalar")
  }
  if (!isPosscal(Rs) || !is.numeric(Rs)) {
    stop("stopband ripple Rp must a non-negative scalar")
  }
  stop <- type == "stop" || type == "high"
  digital <- plane == "z"
  if (!is.vector(w) || (length(w) != 1 && length(w) != 2)) {
    stop("frequency w must be specified as a vector of length 1 or 2 (either w0 or c(w0, w1))")
  }
  if ((type == "stop" || type == "pass") &&  length(w) != 2) {
    stop("w must be two elements for stop and bandpass filters")
  }
  if (digital && !all(w >= 0 & w <= 1)) {
    stop("critical frequencies w must be in the range [0 1]")
  } else if (!digital && !all(w >= 0)) {
    stop("critical frequencies w must be in the range [0 Inf]")
  }

  ## Prewarp to the band edges to s plane
  if (digital) {
    T <- 2                    # sampling frequency of 2 Hz
    w <- 2 / T * tan (pi * w / T)
  }
  
  ## Generate splane poles, zeros and gain
  zpg <- ncauer(Rp, Rs, n)
  
  ## s-plane frequency transform
  zpg <- sftrans(zpg, w = w, stop = stop)
  
  ## Use bilinear transform to convert poles to the z plane
  if (digital) {
    zpg <- bilinear(zpg, T = T)
  }
  
  as.Arma(zpg)
  
}