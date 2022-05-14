# pyule.R
# Copyright (C) 2020 Geert van Boxtel <gjmvanboxtel@gmail.com>
# Original Octave function:
# Copyright (C) 2006 Peter V. Lanspeary <pvl@mecheng.adelaide.edu.au>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# Version history
# 20201106  GvB       setup for gsignal v0.1.0
#------------------------------------------------------------------------------

#' Autoregressive PSD estimate - Yule-Walker method
#'
#' Calculate Yule-Walker autoregressive power spectral density.
#'
#' @param x input data, specified as a numeric or complex vector or matrix. In
#'   case of a vector it represents a single signal; in case of a matrix each
#'   column is a signal.
#' @param p model order; number of poles in the AR model or limit to the number
#'   of poles if a valid criterion is provided. Must be < length(x) - 2.
#' @param freq vector of frequencies at which power spectral density is
#'   calculated, or a scalar indicating the number of uniformly distributed
#'   frequency values at which spectral density is calculated. Default: 256.
#' @param fs sampling frequency (Hz). Default: 1
#' @param range character string. one of:
#' \describe{
#'   \item{\code{"half"} or \code{"onesided"}}{frequency range of the spectrum
#'   is from zero up to but not including \code{fs / 2}. Power from negative
#'   frequencies is added to the positive side of the spectrum.}
#'   \item{\code{"whole"} or \code{"twosided"}}{frequency range of the spectrum
#'   is \code{-fs / 2} to \code{fs / 2}, with negative frequencies stored in
#'   "wrap around order" after the positive frequencies; e.g. frequencies for a
#'   10-point \code{'twosided'} spectrum are 0 0.1 0.2 0.3 0.4 0.5 -0.4 -0.3
#'   -0.2. -0.1.}
#'   \item{\code{"shift"} or \code{"centerdc"}}{same as \code{"whole"} but with
#'   the first half of the spectrum swapped with second half to put the
#'   zero-frequency value in the middle. If \code{freq} is vector,
#'   \code{"shift"} is ignored.}
#' }
#'   Default: If model coefficients \code{a} are real, the default range is
#'   \code{"half"}, otherwise the default range is \code{"whole"}.
#' @param method method used to calculate the power spectral density, either
#'   \code{"fft"} (use the Fast Fourier Transform) or \code{"poly"} (calculate
#'   the power spectrum as a polynomial). This argument is ignored if the
#'   \code{freq} argument is a vector. The default is \code{"poly"} unless the
#'   \code{freq} argument is an integer power of 2.
#'
#' @return An object of class "ar_psd" , which is a list containing two
#'   elements, \code{freq} and \code{psd} containing the frequency values and
#'   the estimates of power-spectral density, respectively.
#'
#' @note This function is a wrapper for \code{arburg} and \code{ar_psd}.
#'
#' @examples
#' A <- Arma(1, c(1, -2.7607, 3.8106, -2.6535, 0.9238))
#' y <- filter(A, 0.2 * rnorm(1024))
#' py <- pyulear(y, 4)
#'
#' @author Peter V. Lanspeary, \email{pvl@@mecheng.adelaide.edu.au}.\cr
#' Conversion to R by Geert van Boxtel, \email{gjmvanboxtel@@gmail.com}
#'
#' @seealso \code{\link{ar_psd}}, \code{\link{arburg}}
#'
#' @export

pyulear <- function(x, p, freq = 256, fs = 1, range = NULL,
                    method = if (length(freq) == 1 &&
                                 bitwAnd(freq, freq - 1) == 0)
                      "fft" else "poly")  {

  coefs <- aryule(x, p)
  if (is.null(range)) {
    range <- ifelse(is.numeric(coefs$a), "half", "whole")
  }
  rv <- ar_psd(coefs$a, coefs$e, freq, fs, range, method)
  rv
}
