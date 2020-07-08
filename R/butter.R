# butter.R
# Copyright (C) 2019 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Octave signal package:
# Copyright (C) 1999 Paul Kienzle <pkienzle@users.sf.net>
# Copyright (C) 2003 Doug Stewart <dastew@sympatico.ca>
# Copyright (C) 2011 Alexander Klein <alexander.klein@math.uni-giessen.de>
# Copyright (C) 2018 John W. Eaton
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
# 20200513 Geert van Boxtel           First version for v0.1.0
# 20200519 Geert van Boxtel           Added plane parameter to butter.IIRfspec
# 20200708 GvB                        renamed IIRfspec to FilterSpecs
#---------------------------------------------------------------------------------------------------------------------------------

#' Butterworth filter design
#' 
#' Compute the transfer function coefficients of a Butterworth filter
#' 
#' Butterworth filters have a magnitude response that is maximally flat in the
#' passband and monotonic overall. This smoothness comes at the price of
#' decreased rolloff steepness. Elliptic and Chebyshev filters generally provide
#' steeper rolloff for a given filter order.
#' 
#' Because butter is generic, it can be extended to accept other inputs, using
#' "buttord" to generate filter criteria for example.
#' 
#' @param n filter order.
#' @param w critical frequencies of the filter. \code{w} must be a scalar for
#'   low-pass and high-pass filters, and \code{w} must be a two-element vector
#'   c(low, high) specifying the lower and upper bands in radians/second. For
#'   digital filters, W must be between 0 and 1 where 1 is the Nyquist
#'   frequency.
#' @param type filter type, one of "low", "high", "stop", or "pass".
#' @param plane "z" for a digital filter or "s" for an analog filter.
#' @param ... additional arguments passed to butter, overriding those given by n
#'   of class \code{FilterSpecs}.
#' 
#' @return list of class \code{'\link{Arma}'} with list elements:
#' \describe{
#'   \item{b}{moving average (MA) polynomial coefficients}
#'   \item{a}{autoregressive (AR) polynomial coefficients}
#' }
#' 
#' @examples
#' ## 50 Hz notch filter
#' fs <- 256
#' bf <- butter(4, c(48, 52) / (fs / 2), "stop")
#' freqz(bf, fs = fs)
#' 
#' ## EEG alpha rhythm (8 - 12 Hz) bandpass filter 
#' fs <- 128
#' fpass <- c(8, 12)
#' wpass <- fpass / (fs / 2)
#' but <- butter(5, wpass, "pass")
#' freqz(but, fs = fs)
#' 
#' ## filter to remove vocals from songs, 25 dB attenuation in stop band
#' ## (not optimal with a Butterworth filter)
#' fs <- 44100
#' specs <- buttord(230/(fs/2), 450/(fs/2), 1, 25)
#' bf <- butter(specs)
#' freqz(bf, fs = fs)
#' zplane(bf)
#' 
#' @references \url{https://en.wikipedia.org/wiki/Butterworth_filter}
#' 
#' @seealso \code{\link{Arma}}, \code{\link{filter}}, \code{cheby1}, \code{ellip}, \code{\link{buttord}}
#' 
#' @author Original Octave code by Paul Kienzle \email{pkienzle@@users.sf.net},
#'   Doug Stewart \email{dastew@@sympatico.ca}, Alexander Klein
#'   \email{alexander.klein@@math.uni-giessen.de}, John W. Eaton. Port to R Tom
#'   Short, adapted by Geert van Boxtel \email{G.J.M.vanBoxtel@@gmail.com}.
#'
#' @rdname butter
#' @export

butter <- function(n, ...) UseMethod("butter")

#' @rdname butter
#' @export

butter.FilterSpecs <- function(n, ...)
  butter(n$n, n$Wc, n$type, n$plane, ...)

#' @rdname butter
#' @export

butter.default <- function (n, w, type = c("low", "high", "stop", "pass"), plane = c("z", "s"), ...) {
  
  # check input arguments
  type <- match.arg(type)
  plane <- match.arg(plane)
  if (!isPosscal(n) || !isWhole(n)) {
    stop("filter order n must be a positive integer")
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
    wc <- 2 / T * tan (pi * w / T)
  }
  
  ## Generate splane poles for the prototype Butterworth filter
  ## source: Kuc
  C <- 1                      # default cutoff frequency
  pole <- C * exp(1i * pi * (2 * 1:n + n - 1) / (2 * n))
  if (n %% 2 == 1){
    pole[(n + 1) / 2] <- -1   # pure real value at exp(i*pi)
  }
  zero <- numeric(0)
  gain <- C^n
  
  zpg <- Zpg(z = zero, p = pole, g = gain)
  
  ## splane frequency transform
  zpg <- sftrans(zpg, w = w, stop = stop)
  
  ## Use bilinear transform to convert poles to the z plane
  if (digital) {
    zpg <- bilinear(zpg, T = T)
  }
  
  as.Arma(zpg)
  
}
